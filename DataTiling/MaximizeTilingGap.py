import z3
import numpy as np
from decimal import Decimal, getcontext
import time

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

def MaximizeGap(K=4, num_iterations=20):
    
    s = z3.Solver()
    s.set("timeout", 300000)
    
    float_sort = z3.FPSort(8, 24)
    rm = z3.RoundNearestTiesToEven()
    
    a = [z3.FP(f"a_{i}", float_sort) for i in range(K)]
    b = [z3.FP(f"b_{i}", float_sort) for i in range(K)]
    c_init = z3.FP("c_init", float_sort)
    alpha = z3.FP("alpha", float_sort)
    
    def get_z3_kernel_acc(a_sub, b_sub):
        acc = z3.FPVal(0.0, float_sort)
        for ai, bi in zip(a_sub, b_sub):
            acc = z3.fpFMA(rm, ai, bi, acc)
        return acc

    acc1 = get_z3_kernel_acc(a[0:2], b[0:2])
    c_mid = z3.fpAdd(rm, z3.fpMul(rm, alpha, acc1), c_init)
    acc2 = get_z3_kernel_acc(a[2:4], b[2:4])
    c_final_a = z3.fpAdd(rm, z3.fpMul(rm, alpha, acc2), c_mid)
    
    acc_all = get_z3_kernel_acc(a[0:4], b[0:4])
    c_final_b = z3.fpAdd(rm, z3.fpMul(rm, alpha, acc_all), c_init)
    
    s.add(c_final_a != c_final_b)
    
    for x in a + b + [c_init, alpha]:
        s.add(z3.Not(z3.fpIsNaN(x)), z3.Not(z3.fpIsInf(x)), z3.Not(z3.fpIsZero(x)))

    max_gap = 0.0
    best_inputs = None
    
    iteration = 0

    while iteration < num_iterations:
        iteration += 1
        print(f"  Checking iteration {iteration}...")
        res = s.check()
        print(f"  Result: {res}")
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

            p_acc1 = py_kernel_acc(a_v[0:2], b_v[0:2])
            p_c_mid = np.float32(float(alpha_v) * float(p_acc1) + float(c_i_v))
            p_acc2 = py_kernel_acc(a_v[2:4], b_v[2:4])
            p_c_final_a = np.float32(float(alpha_v) * float(p_acc2) + float(p_c_mid))
            
            p_acc_all = py_kernel_acc(a_v, b_v)
            p_c_final_b = np.float32(float(alpha_v) * float(p_acc_all) + float(c_i_v))
            
            current_gap = abs(float(p_c_final_a) - float(p_c_final_b))
            print(f"  Counter-example {iteration}: gap = {current_gap}")
            
            if current_gap > max_gap:
                max_gap = current_gap
                best_inputs = (alpha_v, c_i_v, a_v, b_v, p_c_final_a, p_c_final_b)
            
            s.add(z3.Or([x != m[x] for x in a + b + [c_init, alpha]]))
        else:
            break
            
    if best_inputs:
        alpha_v, c_i_v, a_v, b_v, res_a, res_b = best_inputs
        print(f"\n--- Maximum Gap Found ---")
        print(f"  Inputs: alpha={alpha_v}, c_init={c_i_v}")
        print(f"  a: {a_v}")
        print(f"  b: {b_v}")
        print(f"  Result Config A (KC=2): {res_a}")
        print(f"  Result Config B (KC=4): {res_b}")
        print(f"  Absolute Execution Difference: {abs(res_a - res_b)}")
    else:
        print("Could not find any divergence.")

if __name__ == "__main__":
    MaximizeGap()
