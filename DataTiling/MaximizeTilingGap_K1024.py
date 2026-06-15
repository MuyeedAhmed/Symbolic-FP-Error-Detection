import z3
import numpy as np
from decimal import Decimal, getcontext
import pickle
import os
import time
import sys
import random

getcontext().prec = 100

def ReadAndPrintPickleSolutions(k=1024):
    pickle_path = f'Solutions/TilingGaps_K{k}.pkl'
    with open(pickle_path, 'rb') as f:
        witnesses = pickle.load(f)

    for i, w in enumerate(witnesses):
        print(f"\n{'='*20} Witness {i+1} (K={w.get('K', k)}) {'='*20}")
        print(f"alpha: {w['alpha']}")
        print(f"c_init: {w['c_init']}")
        print(f"Gap: {w['gap']}")
        
        kc1 = w.get('KC1', 'unknown')
        kc2 = w.get('KC2', 'unknown')
        print(f"Block Sizes: {kc1} vs {kc2}")

        a = np.array(w['a'])
        b = np.array(w['b'])

        print("\n--- A ---")
        non_zero_a = np.where(a != 0)[0]
        for idx in non_zero_a:
            print(f"Index {idx:4d}: {a[idx]}")

        print("\n--- B ---")
        non_zero_b = np.where(b != 0)[0]
        for idx in non_zero_b:
            print(f"Index {idx:4d}: {b[idx]}")
        # print("\nFull Vector a:", a.tolist())
        # print("\nFull Vector b:", b.tolist())


def fma_python(a, b, c):
    res = Decimal(float(a)) * Decimal(float(b)) + Decimal(float(c))
    return np.float32(float(res))

def fp32_to_python(fp_val_model):
    return float(z3.simplify(fp_val_model).as_decimal(20).replace('?', ''))

def get_bit_exact_float(m, z3_var):
    val = m.eval(z3_var)
    if val is None: return np.float32(0.0)
    s_val = str(val).replace('*(2**', '*2**').replace(')', '')
    try:
        return np.float32(eval(s_val))
    except:
        try:
            return np.float32(float(fp32_to_python(val)))
        except:
            return np.float32(0.0)

def compute_tiled_fma(a, b, c_init, alpha, kc):
    K = len(a)
    c_acc = np.float32(c_init)
    alpha_f = np.float32(alpha)
    for i in range(0, K, kc):
        block_acc = np.float32(0.0)
        for j in range(i, min(i + kc, K)):
            block_acc = fma_python(a[j], b[j], block_acc)
        c_acc = np.float32(np.float32(alpha_f * block_acc) + c_acc)
    return c_acc

def MaximizeGapK1024():
    K = 1024
    KC1 = 256
    KC2 = 512
    
    start_time = time.time()
    s = z3.Tactic('qffp').solver()
    s.set("timeout", 120000) # 120 seconds
    
    float_sort = z3.FPSort(8, 24)
    rm = z3.RoundNearestTiesToEven()
    
    sym_indices = [0, 256, 512, 768]
            
    a_sym = [z3.FP(f"a_{idx}", float_sort) for idx in sym_indices]
    b_sym = [z3.FP(f"b_{idx}", float_sort) for idx in sym_indices]
    c_init = z3.FP("c_init", float_sort)
    alpha = z3.FP("alpha", float_sort)
    
    a_fixed = np.zeros(K, dtype=np.float32)
    b_fixed = np.zeros(K, dtype=np.float32)
    
    def build_z3_block_acc(start, end):
        f_sum = np.float32(0.0)
        sym_pairs = []
        for idx in range(start, end):
            if idx in sym_indices:
                s_idx = sym_indices.index(idx)
                sym_pairs.append((a_sym[s_idx], b_sym[s_idx]))
            else:
                f_sum = fma_python(a_fixed[idx], b_fixed[idx], f_sum)
        
        acc = z3.FPVal(float(f_sum), float_sort)
        for sa, sb in sym_pairs:
            acc = z3.fpAdd(rm, z3.fpMul(rm, sa, sb), acc)
        return acc

    c_acc_a = c_init
    for i in range(0, K, KC1):
        block_acc = build_z3_block_acc(i, i + KC1)
        c_acc_a = z3.fpAdd(rm, z3.fpMul(rm, alpha, block_acc), c_acc_a)
    c_final_a = c_acc_a
    
    c_acc_b = c_init
    for i in range(0, K, KC2):
        block_acc = build_z3_block_acc(i, i + KC2)
        c_acc_b = z3.fpAdd(rm, z3.fpMul(rm, alpha, block_acc), c_acc_b)
    c_final_b = c_acc_b

    for x in a_sym + b_sym:
        s.add(z3.fpIsNormal(x))
        s.add(z3.fpAbs(x) <= z3.FPVal(1e5, float_sort))
    
    s.add(z3.fpIsNormal(c_init), z3.fpIsNormal(alpha))
    s.add(z3.fpAbs(c_init) >= z3.FPVal(1e3, float_sort))
    s.add(z3.fpAbs(c_init) <= z3.FPVal(1e6, float_sort))
    s.add(z3.fpAbs(alpha) >= z3.FPVal(0.1, float_sort))
    s.add(z3.fpAbs(alpha) <= z3.FPVal(10.0, float_sort))

    s.add(z3.Not(z3.fpIsInf(c_final_a)), z3.Not(z3.fpIsInf(c_final_b)))

    all_solutions = []
    target_gap = 1e-3
    print(f"Searching for K={K} divergence (KC1={KC1} vs KC2={KC2})...")
    
    for attempt in range(1, 6):
        s.push()
        seed = random.randint(0, 10000)
        s.set("smt.random_seed", seed)
        s.add(z3.fpGT(z3.fpAbs(z3.fpSub(rm, c_final_a, c_final_b)), z3.FPVal(target_gap, float_sort)))
        
        print(f"Attempt {attempt}: Target Gap > {target_gap:.2e}...")
        iter_start = time.time()
        res = s.check()
        iter_end = time.time()
        
        if res == z3.sat:
            m = s.model()
            alpha_v = get_bit_exact_float(m, alpha)
            c_init_v = get_bit_exact_float(m, c_init)
            a_v = a_fixed.copy()
            b_v = b_fixed.copy()
            for idx, sa, sb in zip(sym_indices, a_sym, b_sym):
                a_v[idx] = get_bit_exact_float(m, sa)
                b_v[idx] = get_bit_exact_float(m, sb)

            print(f"Full Vector a: {a_v.tolist()}")
            print(f"Full Vector b: {b_v.tolist()}")


            res_a = compute_tiled_fma(a_v, b_v, c_init_v, alpha_v, KC1)
            res_b = compute_tiled_fma(a_v, b_v, c_init_v, alpha_v, KC2)
            gap = abs(float(res_a) - float(res_b))
            
            all_solutions.append({
                'alpha': alpha_v, 'c_init': c_init_v, 
                'a': a_v, 'b': b_v, 'gap': gap,
                'K': K, 'KC1': KC1, 'KC2': KC2
            })
            print(f"SUCCESS! Verified Gap: {gap:.10e} (Time: {iter_end - iter_start:.2f}s)")
            s.pop()
            if gap > target_gap:
                target_gap = gap * 1.5
        else:
            s.pop()
            print(f"Attempt {attempt} failed (Result: {res}, Time: {iter_end - iter_start:.2f}s)")
            target_gap /= 10.0
            if target_gap < 1e-12: break
            
    os.makedirs('Solutions', exist_ok=True)
    with open(f'Solutions/TilingGaps_K{K}.pkl', 'wb') as f:
        pickle.dump(all_solutions, f)

    duration = time.time() - start_time
    max_gap = max([sol['gap'] for sol in all_solutions]) if all_solutions else 0
    print(f"\nSearch complete. Found {len(all_solutions)} solutions in {duration:.2f}s.")
    print(f"Maximum Gap found: {max_gap}")

if __name__ == "__main__":
    MaximizeGapK1024()
    ReadAndPrintPickleSolutions(1024)
