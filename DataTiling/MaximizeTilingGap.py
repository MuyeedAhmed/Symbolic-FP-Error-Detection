import z3
import numpy as np
from decimal import Decimal, getcontext
import pickle
import os
import time
import sys
import random

getcontext().prec = 100

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

def MaximizeGap(K=4):
    start_time = time.time()
    s = z3.Tactic('qffp').solver()
    s.set("timeout", 60000)
    
    float_sort = z3.FPSort(8, 24)
    rm = z3.RoundNearestTiesToEven()
    
    a = [z3.FP(f"a_{i}", float_sort) for i in range(K)]
    b = [z3.FP(f"b_{i}", float_sort) for i in range(K)]
    c_init = z3.FP("c_init", float_sort)
    alpha = z3.FP("alpha", float_sort)
    
    def get_z3_kernel_acc(a_sub, b_sub):
        acc = z3.FPVal(0.0, float_sort)
        for ai, bi in zip(a_sub, b_sub):
            # acc = z3.fpFMA(rm, ai, bi, acc)
            acc = z3.fpAdd(rm, z3.fpMul(rm, ai, bi), acc)
        return acc

    c_acc = c_init
    for i in range(0, K, 2):
        acc = get_z3_kernel_acc(a[i:i+2], b[i:i+2])
        c_acc = z3.fpAdd(rm, z3.fpMul(rm, alpha, acc), c_acc)
    c_final_a = c_acc
    
    acc_all = get_z3_kernel_acc(a[0:K], b[0:K])
    c_final_b = z3.fpAdd(rm, z3.fpMul(rm, alpha, acc_all), c_init)
    
    for x in a + b + [c_init, alpha]:
        s.add(z3.fpIsNormal(x))
        s.add(z3.fpAbs(x) <= z3.FPVal(1e6, float_sort))

    # s.add(z3.fpAbs(c_init) >= z3.FPVal(1e4, float_sort))
    # s.add(z3.fpAbs(alpha) <= z3.FPVal(1.0, float_sort))

    s.add(z3.Not(z3.fpIsInf(c_final_a)), z3.Not(z3.fpIsInf(c_final_b)))

    all_solutions = []
    target_gap = 1e-3
    found_any = False
    
    print(f"K={K}")
    
    for attempt in range(1, 11):
        s.push()
        seed = random.randint(0, 10000)
        s.set("smt.random_seed", seed)
        
        s.add(z3.fpGT(z3.fpAbs(z3.fpSub(rm, c_final_a, c_final_b)), z3.FPVal(target_gap, float_sort)))
        
        print(f"Attempt {attempt}: Checking for gap > {target_gap:.2e} (seed: {seed})...")
        iter_start = time.time()
        res = s.check()
        iter_end = time.time()
        
        if res == z3.sat:
            m = s.model()
            alpha_v = get_bit_exact_float(m, alpha)
            c_i_v = get_bit_exact_float(m, c_init)
            a_v = [get_bit_exact_float(m, x) for x in a]
            b_v = [get_bit_exact_float(m, x) for x in b]

            def py_kernel_acc(a_sub, b_sub):
                acc = np.float32(0.0)
                for ai, bi in zip(a_sub, b_sub):
                    acc = fma_python(ai, bi, acc)
                return acc

            py_c_acc = np.float32(c_i_v)
            for i in range(0, K, 2):
                acc_v = py_kernel_acc(a_v[i:i+2], b_v[i:i+2])
                term = np.float32(float(alpha_v) * float(acc_v))
                py_c_acc = np.float32(float(term) + float(py_c_acc))
            p_c_final_a = py_c_acc

            p_acc_all = py_kernel_acc(a_v, b_v)
            term_b = np.float32(float(alpha_v) * float(p_acc_all))
            p_c_final_b = np.float32(float(term_b) + float(c_i_v))
            
            gap = abs(float(p_c_final_a) - float(p_c_final_b))
            
            all_solutions.append({'alpha': alpha_v, 'c_init': c_i_v, 'a': a_v, 'b': b_v, 'gap': gap})
            print(f"SUCCESS! Gap: {gap:.10e}, Time: {iter_end - iter_start:.2f}s")
            
            s.pop()
            if gap > target_gap:
                target_gap = gap * 1.5
            found_any = True
        else:
            s.pop()
            print(f"Attempt {attempt} failed (Result: {res}, Time: {iter_end - iter_start:.2f}s)")
            if not found_any:
                target_gap /= 10.0
                if target_gap < 1e-12: break
            else:
                pass
            
    os.makedirs('Solutions', exist_ok=True)
    with open(f'Solutions/TilingGaps_K{K}.pkl', 'wb') as f:
        pickle.dump(all_solutions, f)

    duration = time.time() - start_time
    max_gap = max([sol['gap'] for sol in all_solutions]) if all_solutions else 0
    print(f"K={K}: Found {len(all_solutions)} solutions in {duration:.2f}s. Max Gap: {max_gap}")

if __name__ == "__main__":
    k_size = 4
    if len(sys.argv) > 1:
        k_size = int(sys.argv[1])
        
    MaximizeGap(K=k_size)
