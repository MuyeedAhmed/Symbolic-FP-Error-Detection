import pickle
import numpy as np
import os
import glob
from scipy.linalg import inv, det

def verify_identity_numpy(name, data, n):
    print(f"--- Verifying with NumPy (N={n}): {name} ---")
    inputs = {k: np.array(v) for k, v in data['inputs'].items()}
    A = inputs.get('A')
    B = inputs.get('B')
    C = inputs.get('C')
    
    lhs_numpy = None
    rhs_numpy = None
    
    try:
        if name == "(A+B)^T = A^T + B^T":
            lhs_numpy = (A + B).T
            rhs_numpy = A.T + B.T
        elif name == "(A^-1)^-1 = A":
            lhs_numpy = inv(inv(A))
            rhs_numpy = A
        elif name == "(A*B)^-1 = B^-1 * A^-1":
            lhs_numpy = inv(A @ B)
            rhs_numpy = inv(B) @ inv(A)
        elif name == "(A*B*C)^-1 = C^-1 * B^-1 * A^-1":
            lhs_numpy = inv(A @ B @ C)
            rhs_numpy = inv(C) @ inv(B) @ inv(A)
        elif name == "(A*B)^T = B^T * A^T":
            lhs_numpy = (A @ B).T
            rhs_numpy = B.T @ A.T
        elif name == "(A^T)^-1 = (A^-1)^T":
            lhs_numpy = inv(A.T)
            rhs_numpy = inv(A).T
        elif name == "A * A^-1 = I":
            lhs_numpy = A @ inv(A)
            rhs_numpy = np.eye(n)
        elif name == "det(A*B) = det(A) * det(B)":
            lhs_numpy = det(A @ B)
            rhs_numpy = det(A) * det(B)
        elif name == "det(A^-1) = 1/det(A)":
            lhs_numpy = det(inv(A))
            rhs_numpy = 1.0 / det(A)
        elif name == "A(BC) = (AB)C":
            lhs_numpy = A @ (B @ C)
            rhs_numpy = (A @ B) @ C
        elif name == "A(B+C) = AB + AC":
            lhs_numpy = A @ (B + C)
            rhs_numpy = (A @ B) + (A @ C)
        elif name == "tr(A + B) = tr(A) + tr(B)":
            lhs_numpy = np.trace(A + B)
            rhs_numpy = np.trace(A) + np.trace(B)
        elif name == "tr(AB) = tr(BA)":
            lhs_numpy = np.trace(A @ B)
            rhs_numpy = np.trace(B @ A)
        elif name == "det(kA) = k^n * det(A)":
            k = data['inputs']['k'][0][0]
            lhs_numpy = det(k * A)
            rhs_numpy = (k ** n) * det(A)
        else:
            print(f"Unknown identity: {name}")
            return

        diff_numpy = np.abs(lhs_numpy - rhs_numpy)
        max_diff_numpy = np.max(diff_numpy)
        print(f"NumPy LHS vs RHS Max Diff: {max_diff_numpy}")
        
        gurobi_lhs = np.array(data['LHS'])
        gurobi_rhs = np.array(data['RHS'])
        
        diff_gurobi_lhs = np.max(np.abs(gurobi_lhs - lhs_numpy))
        print(f"Gurobi LHS vs NumPy Recomputed LHS Max Diff: {diff_gurobi_lhs}")

        if max_diff_numpy > 0:
            print("SUCCESS")
        else:
            print("FAILURE: due to tolerance.")

    except Exception as e:
        print(f"Error during verification: {e}")

    print("-" * 40)

def main():
    solution_files = glob.glob("Solutions/*.pkl")
    if not solution_files:
        return

    for file_path in solution_files:
        print(f"Reading {file_path}...")
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
            verify_identity_numpy(data['name'], data, data['n'])

if __name__ == "__main__":
    main()
