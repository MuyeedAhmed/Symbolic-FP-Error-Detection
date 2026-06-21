import z3
import struct
import numpy as np
import itertools

from FPReduction import FPReduction


def _bits_to_f32(n: int) -> np.float32:
    return np.float32(struct.unpack("f", struct.pack("I", n & 0xFFFF_FFFF))[0])


def _extract_f32(model, fp_expr) -> np.float32:
    bv = model.eval(z3.fpToIEEEBV(fp_expr))
    return _bits_to_f32(bv.as_long())


def _fval(v: float, sort) -> z3.ExprRef:
    return z3.FPVal(float(np.float32(v)), sort)


'''
left_fold(v) != tree_fold(v)
'''
def prove_arithmetic_divergence(n: int = 4, hw_a: str = "left", hw_b: str = "tree", timeout_ms: int = 30_000,
                                min_val: float = 0.0, max_val: float | None = None) -> np.ndarray | None:
    FP32 = z3.Float32()
    rm   = z3.RNE()
    fold_a = FPReduction.z3_get(hw_a)
    fold_b = FPReduction.z3_get(hw_b)

    s = z3.Solver()
    s.set("timeout", timeout_ms)
    vs = [z3.FP(f"v{i}", FP32) for i in range(n)]
    for v in vs:
        s.add(z3.Not(z3.fpIsNaN(v)), z3.Not(z3.fpIsInf(v)))
        lo = max(0.0, float(min_val))
        hi = float(max_val) if max_val is not None else 1e6
        s.add(z3.fpGT(v, _fval(lo, FP32)), z3.fpLT(v, _fval(hi, FP32)))

    s.add(z3.Not(z3.fpEQ(fold_a(vs, rm), fold_b(vs, rm))))

    result = s.check()
    if result != z3.sat:
        print("Arith Div Unsat", result)
        return None

    m = s.model()
    vals = np.array([_extract_f32(m, v) for v in vs], dtype=np.float32)
    fa   = FPReduction.get(hw_a)(vals)
    fb   = FPReduction.get(hw_b)(vals)
    print(f"Arithmetic divergence: SAT  v={vals}")
    print(f"  {hw_a}_fold = {fa!r}")
    print(f"  {hw_b}_fold = {fb!r}")
    print(f"  delta = {abs(fa - fb):.3e}")
    return vals


'''
left_fold(d) < left_fold(e)   [x closer to c0 under CPU]
tree_fold(d) > tree_fold(e)   [x closer to c1 under GPU]
'''
def find_distance_flip(d: int = 4, hw_a: str = "left", hw_b: str = "tree", timeout_ms: int = 60_000,
                       extra_constraints: list | None = None, d_sq_base: np.ndarray | None = None) -> dict | None:
    FP32      = z3.Float32()
    rm        = z3.RNE()
    fold_a_z3 = FPReduction.z3_get(hw_a)
    fold_b_z3 = FPReduction.z3_get(hw_b)
    fold_a    = FPReduction.get(hw_a)
    fold_b    = FPReduction.get(hw_b)

    def _find_sqrt(target: np.float32, n_ulp: int = 2) -> list:
        if target <= 0:
            return [np.float32(0.0)] if target == 0 else []
        base = FPReduction.f32_to_bits(np.float32(np.sqrt(float(target))))
        out  = []
        for delta in range(-n_ulp, n_ulp + 1):
            a = FPReduction.bits_to_f32(base + delta)
            if np.isfinite(a) and np.float32(a * a) == target:
                out.append(a)
        return out

    def _reconstruct(d_sq, e_sq):
        ca = [_find_sqrt(d_sq[j]) for j in range(d)]
        cb = [_find_sqrt(e_sq[j]) for j in range(d)]
        if any(len(c) == 0 for c in ca + cb):
            return None
        for ac in itertools.product(*ca):
            for bc in itertools.product(*cb):
                a = np.array(ac, dtype=np.float32)
                b = np.array(bc, dtype=np.float32)
                if (np.all((a * a).astype(np.float32) == d_sq) and np.all((b * b).astype(np.float32) == e_sq)):
                    return a, b
        return None

    # ── Phase 1: d_sq from Claim 1 ──────────────────────────────────────────

    # Phase 1: anchor d_sq where folds diverge.
    # If the caller supplies d_sq_base (e.g. from the __main__ Claim-1 run),
    # use it directly — avoids re-running a potentially slow Z3 query and
    # ensures Phase 2 gets the same witness that the diagnostics saw.
    if d_sq_base is not None:
        d_sq = np.abs(d_sq_base).astype(np.float32)
    else:
        claim1 = prove_arithmetic_divergence(n=d, hw_a=hw_a, hw_b=hw_b,
                                             timeout_ms=timeout_ms)
        if claim1 is None:
            return None
        d_sq = np.abs(claim1).astype(np.float32)

    # d_sq as Z3 constants (used in Phase-2 solver and returned for margin API)
    d_vars = [z3.FPVal(float(d_sq[j]), FP32) for j in range(d)]

    # ── Phase 2: fix d_sq, search for e_sq ──────────────────────────────────

    # Variable names are stable across calls so extra_constraints (built
    # from e_vars returned by a previous call) remain valid constraints.
    e_extra = list(extra_constraints or [])

    for _retry in range(50):
        e_vars = [z3.FP(f"e{j}", FP32) for j in range(d)]

        s = z3.Solver()
        s.set("timeout", timeout_ms)
        for v in e_vars:
            s.add(z3.Not(z3.fpIsNaN(v)), z3.Not(z3.fpIsInf(v)))
            s.add(z3.fpGT(v, _fval(1e-7, FP32)))
            s.add(z3.fpLEQ(v, _fval(1e6, FP32)))

        # d_vars are constants → Z3 evaluates fold_a(d_vars) / fold_b(d_vars)
        # statically, leaving only a 4-variable formula for the flip condition.
        dist_A_c0 = fold_a_z3(d_vars,  rm)
        dist_A_c1 = fold_a_z3(e_vars,  rm)
        dist_B_c0 = fold_b_z3(d_vars,  rm)
        dist_B_c1 = fold_b_z3(e_vars,  rm)
        s.add(z3.Or(
            z3.And(z3.fpLT(dist_A_c0, dist_A_c1), z3.fpLT(dist_B_c1, dist_B_c0)),
            z3.And(z3.fpLT(dist_A_c1, dist_A_c0), z3.fpLT(dist_B_c0, dist_B_c1)),
        ))
        for c in e_extra:
            s.add(c)

        result = s.check()
        if result != z3.sat:
            return None   # timeout or UNSAT for the e_sq sub-problem

        m    = s.model()
        e_sq = np.array([_extract_f32(m, v) for v in e_vars], dtype=np.float32)

        # ── Reconstruct geometry ─────────────────────────────────────────────

        recon = _reconstruct(d_sq, e_sq)
        if recon is None:
            e_extra.append(z3.Or(*[z3.Not(z3.fpEQ(e_vars[j], _fval(e_sq[j], FP32)))
                                   for j in range(d)]))
            continue

        a_vals, b_vals = recon
        x   = np.zeros(d, dtype=np.float32)
        c0  = -a_vals
        c1  = -b_vals

        diffs0 = (a_vals ** 2).astype(np.float32)
        diffs1 = (b_vals ** 2).astype(np.float32)
        assert np.all(diffs0 == d_sq) and np.all(diffs1 == e_sq)

        d0_A, d1_A = fold_a(diffs0), fold_a(diffs1)
        d0_B, d1_B = fold_b(diffs0), fold_b(diffs1)
        label_A = 0 if d0_A < d1_A else 1
        label_B = 0 if d0_B < d1_B else 1

        if label_A == label_B:
            e_extra.append(z3.Or(*[z3.Not(z3.fpEQ(e_vars[j], _fval(e_sq[j], FP32))) for j in range(d)]))
            continue

        return {
            "d_vars": d_vars, "e_vars": e_vars,
            "d_sq": d_sq,     "e_sq": e_sq,
            "a": a_vals,      "b": b_vals,
            "x": x, "c0": c0, "c1": c1,
            "diffs0": diffs0, "diffs1": diffs1,
            "d0_A": d0_A, "d1_A": d1_A,
            "d0_B": d0_B, "d1_B": d1_B,
            "label_A": label_A, "label_B": label_B,
            "margin_A": float(abs(d0_A - d1_A)),
            "margin_B": float(abs(d0_B - d1_B)),
            "hw_a": hw_a, "hw_b": hw_b,
        }

    return None


def blocking_clause(witness: dict, FP32) -> z3.BoolRef:
    d_vars = witness["d_vars"]
    e_vars = witness["e_vars"]
    d_vals = witness["d_sq"]
    e_vals = witness["e_sq"]
    diffs  = [z3.Not(z3.fpEQ(v, z3.FPVal(float(val), FP32)))
              for v, val in zip(d_vars + e_vars, list(d_vals) + list(e_vals))]
    return z3.Or(*diffs)


def margin_constraint(min_margin: float, witness_template: dict, FP32) -> z3.BoolRef:
    rm   = z3.RNE()
    fabs = lambda x: z3.fpAbs(x)
    fsub = lambda x, y: z3.fpSub(rm, x, y)
    fval = lambda v: z3.FPVal(float(np.float32(v)), FP32)

    d_vars  = witness_template["d_vars"]
    e_vars  = witness_template["e_vars"]
    hw_a    = witness_template["hw_a"]
    hw_b    = witness_template["hw_b"]
    fold_az = FPReduction.z3_get(hw_a)
    fold_bz = FPReduction.z3_get(hw_b)

    # d_vars and e_vars are already squared distances (additive only)
    gap_A = fabs(fsub(fold_az(d_vars, rm), fold_az(e_vars, rm)))
    gap_B = fabs(fsub(fold_bz(d_vars, rm), fold_bz(e_vars, rm)))
    return z3.Or(
        z3.fpGEQ(gap_A, fval(min_margin)),
        z3.fpGEQ(gap_B, fval(min_margin)),
    )



if __name__ == "__main__":
    print("=" * 60)
    print("Z3 proofs: FP tiling sensitivity in KMeans (float32)")
    print("=" * 60)

    c1 = prove_arithmetic_divergence(n=4, hw_a="left", hw_b="tree")

    print("\nFinding distance-flip witness")
    w = find_distance_flip(d=4, hw_a="left", hw_b="tree", d_sq_base=c1)
    if w:
        print(f"  a   = {w['a']}")
        print(f"  b   = {w['b']}")
        print(f"  x   = {w['x']}   c0 = {w['c0']}   c1 = {w['c1']}")
        print(f"  dist(x,c0)  left={w['d0_A']:.10f}  tree={w['d0_B']:.10f}")
        print(f"  dist(x,c1)  left={w['d1_A']:.10f}  tree={w['d1_B']:.10f}")
        print(f"  HW_A (left) assigns x: cluster {w['label_A']}")
        print(f"  HW_B (tree) assigns x: cluster {w['label_B']}")
    else:
        print("  No witness found")
