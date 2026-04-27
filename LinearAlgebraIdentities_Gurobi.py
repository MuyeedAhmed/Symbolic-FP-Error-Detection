import gurobipy as gp
from gurobipy import GRB
import numpy as np
import time
import pickle
import os
from decimal import Decimal, getcontext

# getcontext().prec = 100

def mat_mul_gurobi(model, A, B, n):
    res = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            model.addConstr(res[i][j] == gp.quicksum(A[i][k] * B[k][j] for k in range(n)))
    return res

def transpose_gurobi(M):
    n = len(M)
    m = len(M[0])
    return [[M[j][i] for j in range(n)] for i in range(m)]

def add_gurobi(model, A, B, n):
    res = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            model.addConstr(res[i][j] == A[i][j] + B[i][j])
    return res

def det_gurobi(model, M, n):
    if n == 1: return M[0][0]
    if n == 2:
        res = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
        model.addConstr(res == M[0][0] * M[1][1] - M[0][1] * M[1][0])
        return res
    
    d = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
    terms = []
    for j in range(n):
        minor = [row[:j] + row[j+1:] for row in M[1:]]
        minor_det = det_gurobi(model, minor, n-1)
        term = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
        model.addConstr(term == M[0][j] * minor_det)
        if j % 2 == 1:
            terms.append(-term)
        else:
            terms.append(term)
    model.addConstr(d == gp.quicksum(terms))
    return d

def inv_gurobi(model, M, n):
    d = det_gurobi(model, M, n)
    
    is_pos = model.addVar(vtype=GRB.BINARY)
    is_neg = model.addVar(vtype=GRB.BINARY)
    model.addConstr(is_pos + is_neg == 1)
    model.addGenConstrIndicator(is_pos, True, d >= 1e-4)
    model.addGenConstrIndicator(is_neg, True, d <= -1e-4)

    if n == 1:
        inv_val = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
        model.addConstr(inv_val * d == 1)
        return [[inv_val]], d
    
    cofactors = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            minor = [row[:j] + row[j+1:] for row in (M[:i] + M[i+1:])]
            minor_det = det_gurobi(model, minor, n-1)
            if (i + j) % 2 == 1:
                model.addConstr(cofactors[i][j] == -minor_det)
            else:
                model.addConstr(cofactors[i][j] == minor_det)
    
    adj = [[cofactors[j][i] for j in range(n)] for i in range(n)]
    inv = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            model.addConstr(inv[i][j] * d == adj[i][j])
    return inv, d

def identity_gurobi(n):
    return [[(1.0 if i == j else 0.0) for j in range(n)] for i in range(n)]

def fresh_mat(model, name, n):
    mat = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, name=f"{name}_{i}_{j}") for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            v = mat[i][j]
            is_nonzero = model.addVar(vtype=GRB.BINARY, name=f"nz_{name}_{i}_{j}")
            is_pos = model.addVar(vtype=GRB.BINARY, name=f"pos_{name}_{i}_{j}")
            is_neg = model.addVar(vtype=GRB.BINARY, name=f"neg_{name}_{i}_{j}")
            model.addConstr(is_pos + is_neg == is_nonzero)
            model.addGenConstrIndicator(is_nonzero, False, v == 0)
            model.addGenConstrIndicator(is_pos, True, v >= 1e-4)
            model.addGenConstrIndicator(is_neg, True, v <= -1e-4)
    return mat

def get_vals(mat, n):
    if isinstance(mat, list):
        return [[(v.X if hasattr(v, 'X') else v) for v in row] for row in mat]
    return mat.X if hasattr(mat, 'X') else mat

def verify_identity(name, data, n):
    print(f"--- Verifying with NumPy (N={n}): {name} ---")
    LHS_val = np.array(data['LHS'])
    RHS_val = np.array(data['RHS'])
    
    diff = LHS_val - RHS_val
    max_diff = np.max(np.abs(diff))
    print(f"Max Difference: {max_diff}")
    if max_diff > 1e-7:
        print("Counter-example confirmed by NumPy.")
    else:
        print("NumPy agrees with the identity (within tolerance).")

def run_check(name, n, setup_fn):
    print(f"\n=== Checking Identity: {name} (N={n}) ===")
    m = gp.Model(name)
    m.setParam('NonConvex', 2)
    m.setParam('FeasibilityTol', 1e-9)
    m.setParam('NumericFocus', 3)
    m.setParam('ScaleFlag', 0)
    m.setParam('OutputFlag', 0)
    m.setParam('TimeLimit', 20) 
    inputs, LHS, RHS = setup_fn(m, n)
    
    diff = m.addVar(lb=-GRB.INFINITY, name="diff")
    if isinstance(LHS, list):
        m.addConstr(diff == LHS[0][0] - RHS[0][0])
    else:
        m.addConstr(diff == LHS - RHS)
        
    m.addConstr(diff * diff >= 1.0) 
    m.setObjective(diff * diff, GRB.MAXIMIZE)
    m.Params.SolutionLimit = 1
    
    start = time.time()
    m.optimize()
    duration = time.time() - start
    
    if m.SolCount > 0:
        print(f"Counter-example found in {duration:.2f}s!")
        res_data = {
            'name': name,
            'n': n,
            'inputs': {k: get_vals(v, n) for k, v in inputs.items()},
            'LHS': get_vals(LHS, n),
            'RHS': get_vals(RHS, n)
        }
        file_name = f'solutions/solution_{name.replace(" ", "_").replace("*", "x").replace("^", "").replace("(", "").replace(")", "").replace("/", "_")}.pkl'
        with open(file_name, 'wb') as f:
            pickle.dump(res_data, f)
        verify_identity(name, res_data, n)
    else:
        print(f"No counter-example found in {duration:.2f}s.")


def setup_transpose_sum(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    LHS = transpose_gurobi(add_gurobi(m, A, B, n))
    RHS = add_gurobi(m, transpose_gurobi(A), transpose_gurobi(B), n)
    return {"A": A, "B": B}, LHS, RHS

def setup_transpose_product(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    AB = mat_mul_gurobi(m, A, B, n)
    LHS = transpose_gurobi(AB)
    RHS = mat_mul_gurobi(m, transpose_gurobi(B), transpose_gurobi(A), n)
    return {"A": A, "B": B}, LHS, RHS

def setup_inverse_inverse(m, n):
    A = fresh_mat(m, "A", n)
    A_inv, _ = inv_gurobi(m, A, n)
    LHS, _ = inv_gurobi(m, A_inv, n)
    RHS = A
    return {"A": A}, LHS, RHS

def setup_inverse_product(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    AB = mat_mul_gurobi(m, A, B, n)
    LHS, _ = inv_gurobi(m, AB, n)
    A_inv, _ = inv_gurobi(m, A, n)
    B_inv, _ = inv_gurobi(m, B, n)
    RHS = mat_mul_gurobi(m, B_inv, A_inv, n)
    return {"A": A, "B": B}, LHS, RHS

def setup_transpose_inverse(m, n):
    A = fresh_mat(m, "A", n)
    AT = transpose_gurobi(A)
    LHS, _ = inv_gurobi(m, AT, n)
    A_inv, _ = inv_gurobi(m, A, n)
    RHS = transpose_gurobi(A_inv)
    return {"A": A}, LHS, RHS

def setup_inverse_identity(m, n):
    A = fresh_mat(m, "A", n)
    A_inv, _ = inv_gurobi(m, A, n)
    LHS = mat_mul_gurobi(m, A, A_inv, n)
    RHS = identity_gurobi(n)
    return {"A": A}, LHS, RHS

def setup_det_product(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    AB = mat_mul_gurobi(m, A, B, n)
    LHS = det_gurobi(m, AB, n)
    detA = det_gurobi(m, A, n)
    detB = det_gurobi(m, B, n)
    RHS = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
    m.addConstr(RHS == detA * detB)
    return {"A": A, "B": B}, LHS, RHS

def setup_det_inverse(m, n):
    A = fresh_mat(m, "A", n)
    A_inv, detA = inv_gurobi(m, A, n)
    LHS = det_gurobi(m, A_inv, n)
    RHS = m.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY)
    m.addConstr(RHS * detA == 1.0)
    return {"A": A}, LHS, RHS

def setup_associativity(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    C = fresh_mat(m, "C", n)
    BC = mat_mul_gurobi(m, B, C, n)
    LHS = mat_mul_gurobi(m, A, BC, n)
    AB = mat_mul_gurobi(m, A, B, n)
    RHS = mat_mul_gurobi(m, AB, C, n)
    return {"A": A, "B": B, "C": C}, LHS, RHS

def setup_distributivity(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    C = fresh_mat(m, "C", n)
    BplusC = add_gurobi(m, B, C, n)
    LHS = mat_mul_gurobi(m, A, BplusC, n)
    AB = mat_mul_gurobi(m, A, B, n)
    AC = mat_mul_gurobi(m, A, C, n)
    RHS = add_gurobi(m, AB, AC, n)
    return {"A": A, "B": B, "C": C}, LHS, RHS

def setup_inverse_triple_product(m, n):
    A = fresh_mat(m, "A", n)
    B = fresh_mat(m, "B", n)
    C = fresh_mat(m, "C", n)
    ABC = mat_mul_gurobi(m, A, mat_mul_gurobi(m, B, C, n), n)
    LHS, _ = inv_gurobi(m, ABC, n)
    A_inv, _ = inv_gurobi(m, A, n)
    B_inv, _ = inv_gurobi(m, B, n)
    C_inv, _ = inv_gurobi(m, C, n)
    RHS = mat_mul_gurobi(m, C_inv, mat_mul_gurobi(m, B_inv, A_inv, n), n)
    return {"A": A, "B": B, "C": C}, LHS, RHS

if __name__ == "__main__":
    run_check("(A+B)^T = A^T + B^T", 3, setup_transpose_sum)
    run_check("(A^-1)^-1 = A", 2, setup_inverse_inverse)
    run_check("(A*B)^-1 = B^-1 * A^-1", 2, setup_inverse_product)
    run_check("(A*B*C)^-1 = C^-1 * B^-1 * A^-1", 2, setup_inverse_triple_product)
    run_check("(A*B)^T = B^T * A^T", 3, setup_transpose_product)
    run_check("(A^T)^-1 = (A^-1)^T", 2, setup_transpose_inverse)
    run_check("A * A^-1 = I", 2, setup_inverse_identity)
    run_check("det(A*B) = det(A) * det(B)", 3, setup_det_product)
    run_check("det(A^-1) = 1/det(A)", 2, setup_det_inverse)
    run_check("A(BC) = (AB)C", 2, setup_associativity)
    run_check("A(B+C) = AB + AC", 5, setup_distributivity)
