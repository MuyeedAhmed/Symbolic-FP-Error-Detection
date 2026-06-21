import z3
import struct
import itertools
import numpy as np

from FPReduction import FPReduction


def _bits_to_f32(n):
    return np.float32(struct.unpack("f", struct.pack("I", n & 0xFFFF_FFFF))[0])


def _extract_f32(model, fp_expr):
    return _bits_to_f32(model.eval(z3.fpToIEEEBV(fp_expr)).as_long())


def _fval(v, sort):
    return z3.FPVal(float(np.float32(v)), sort)


def prove_arithmetic_divergence(n=4, hw_a="left", hw_b="tree", timeout_ms=30_000, min_val=0.0, max_val=None):
    FP32 = z3.Float32()
    rm   = z3.RNE()
    s    = z3.Solver()
    s.set("timeout", timeout_ms)

    vs = [z3.FP(f"v{i}", FP32) for i in range(n)]
    for v in vs:
        s.add(z3.Not(z3.fpIsNaN(v)), z3.Not(z3.fpIsInf(v)))
        s.add(z3.fpGT(v, _fval(max(0.0, float(min_val)), FP32)))
        s.add(z3.fpLT(v, _fval(float(max_val) if max_val else 1e6, FP32)))
    s.add(z3.Not(z3.fpEQ(FPReduction.z3_get(hw_a)(vs, rm),
                          FPReduction.z3_get(hw_b)(vs, rm))))

    if s.check() != z3.sat:
        return None

    m    = s.model()
    vals = np.array([_extract_f32(m, v) for v in vs], dtype=np.float32)
    fa   = FPReduction.get(hw_a)(vals)
    fb   = FPReduction.get(hw_b)(vals)
    print(f"Claim 1: v={vals}")
    print(f"  {hw_a}_fold = {fa!r}  {hw_b}_fold = {fb!r}  delta = {abs(fa - fb):.3e}")
    return vals


def find_distance_flip(d=4, hw_a="left", hw_b="tree", timeout_ms=60_000, extra_constraints=None, d_sq_base=None):
    FP32      = z3.Float32()
    rm        = z3.RNE()
    fold_a_z3 = FPReduction.z3_get(hw_a)
    fold_b_z3 = FPReduction.z3_get(hw_b)
    fold_a    = FPReduction.get(hw_a)
    fold_b    = FPReduction.get(hw_b)

    def _sqrt_cands(target, n_ulp=2):
        if target <= 0:
            return [np.float32(0.0)] if target == 0 else []
        base = FPReduction.f32_to_bits(np.float32(np.sqrt(float(target))))
        return [a for delta in range(-n_ulp, n_ulp + 1)
                if np.isfinite(a := FPReduction.bits_to_f32(base + delta))
                and np.float32(a * a) == target]

    def _reconstruct(d_sq, e_sq):
        ca = [_sqrt_cands(d_sq[j]) for j in range(d)]
        cb = [_sqrt_cands(e_sq[j]) for j in range(d)]
        if any(not c for c in ca + cb):
            return None
        for ac in itertools.product(*ca):
            for bc in itertools.product(*cb):
                a = np.array(ac, dtype=np.float32)
                b = np.array(bc, dtype=np.float32)
                if (np.all((a * a).astype(np.float32) == d_sq) and
                        np.all((b * b).astype(np.float32) == e_sq)):
                    return a, b
        return None

    if d_sq_base is not None:
        d_sq = np.abs(d_sq_base).astype(np.float32)
    else:
        c1 = prove_arithmetic_divergence(n=d, hw_a=hw_a, hw_b=hw_b, timeout_ms=timeout_ms)
        if c1 is None:
            return None
        d_sq = np.abs(c1).astype(np.float32)

    d_vars  = [z3.FPVal(float(d_sq[j]), FP32) for j in range(d)]
    e_extra = list(extra_constraints or [])

    for _ in range(50):
        e_vars = [z3.FP(f"e{j}", FP32) for j in range(d)]
        s = z3.Solver()
        s.set("timeout", timeout_ms)
        for v in e_vars:
            s.add(z3.Not(z3.fpIsNaN(v)), z3.Not(z3.fpIsInf(v)))
            s.add(z3.fpGT(v, _fval(1e-7, FP32)), z3.fpLEQ(v, _fval(1e6, FP32)))

        dist_A_c0 = fold_a_z3(d_vars, rm)
        dist_A_c1 = fold_a_z3(e_vars, rm)
        dist_B_c0 = fold_b_z3(d_vars, rm)
        dist_B_c1 = fold_b_z3(e_vars, rm)
        s.add(z3.Or(
            z3.And(z3.fpLT(dist_A_c0, dist_A_c1), z3.fpLT(dist_B_c1, dist_B_c0)),
            z3.And(z3.fpLT(dist_A_c1, dist_A_c0), z3.fpLT(dist_B_c0, dist_B_c1)),
        ))
        for c in e_extra:
            s.add(c)

        if s.check() != z3.sat:
            return None

        m    = s.model()
        e_sq = np.array([_extract_f32(m, v) for v in e_vars], dtype=np.float32)
        recon = _reconstruct(d_sq, e_sq)
        if recon is None:
            e_extra.append(z3.Or(*[z3.Not(z3.fpEQ(e_vars[j], _fval(e_sq[j], FP32)))
                                   for j in range(d)]))
            continue

        a_vals, b_vals = recon
        x  = np.zeros(d, dtype=np.float32)
        c0 = -a_vals
        c1 = -b_vals

        diffs0 = (a_vals ** 2).astype(np.float32)
        diffs1 = (b_vals ** 2).astype(np.float32)
        d0_A, d1_A = fold_a(diffs0), fold_a(diffs1)
        d0_B, d1_B = fold_b(diffs0), fold_b(diffs1)
        label_A = 0 if d0_A < d1_A else 1
        label_B = 0 if d0_B < d1_B else 1

        if label_A == label_B:
            e_extra.append(z3.Or(*[z3.Not(z3.fpEQ(e_vars[j], _fval(e_sq[j], FP32)))
                                   for j in range(d)]))
            continue

        return dict(
            d_vars=d_vars, e_vars=e_vars, d_sq=d_sq, e_sq=e_sq,
            a=a_vals, b=b_vals, x=x, c0=c0, c1=c1,
            diffs0=diffs0, diffs1=diffs1,
            d0_A=d0_A, d1_A=d1_A, d0_B=d0_B, d1_B=d1_B,
            label_A=label_A, label_B=label_B,
            margin_A=float(abs(d0_A - d1_A)), margin_B=float(abs(d0_B - d1_B)),
            hw_a=hw_a, hw_b=hw_b,
        )
    return None


def blocking_clause(witness, FP32):
    d_vars, e_vars = witness["d_vars"], witness["e_vars"]
    d_vals, e_vals = witness["d_sq"], witness["e_sq"]
    return z3.Or(*[z3.Not(z3.fpEQ(v, z3.FPVal(float(val), FP32)))
                   for v, val in zip(d_vars + e_vars, list(d_vals) + list(e_vals))])


def margin_constraint(min_margin, witness, FP32):
    rm       = z3.RNE()
    d_vars   = witness["d_vars"]
    e_vars   = witness["e_vars"]
    fold_az  = FPReduction.z3_get(witness["hw_a"])
    fold_bz  = FPReduction.z3_get(witness["hw_b"])
    fval     = lambda v: z3.FPVal(float(np.float32(v)), FP32)
    gap_A    = z3.fpAbs(z3.fpSub(rm, fold_az(d_vars, rm), fold_az(e_vars, rm)))
    gap_B    = z3.fpAbs(z3.fpSub(rm, fold_bz(d_vars, rm), fold_bz(e_vars, rm)))
    return z3.Or(z3.fpGEQ(gap_A, fval(min_margin)), z3.fpGEQ(gap_B, fval(min_margin)))
