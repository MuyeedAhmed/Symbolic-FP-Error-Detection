import pickle
import numpy as np
import os
import glob
from scipy.linalg import inv, det

def verify_identity_numpy(name, data, n):
    print(f"--- {name} (N={n}) ---")
    inputs = {k: np.array(v) for k, v in data['inputs'].items()}
    A = inputs.get('A')
    B = inputs.get('B')
    C = inputs.get('C')
    
    lhs_numpy = None
    rhs_numpy = None

    print(f"Inputs:")
    for k, v in inputs.items():
        print(f"  {k}: {v}")
    try:
        if name == "(A+B)^T = A^T + B^T":
            lhs_numpy = (A + B).T
            rhs_numpy = A.T + B.T
        elif name == "(A^-1)^-1 = A":
            lhs_numpy = inv(inv(A, check_finite=False), check_finite=False)
            rhs_numpy = A
        elif name == "(A*B)^-1 = B^-1 * A^-1":
            lhs_numpy = inv(A @ B, check_finite=False)
            rhs_numpy = inv(B, check_finite=False) @ inv(A, check_finite=False)
        elif name == "(A*B*C)^-1 = C^-1 * B^-1 * A^-1":
            lhs_numpy = inv(A @ B @ C, check_finite=False)
            rhs_numpy = inv(C, check_finite=False) @ inv(B, check_finite=False) @ inv(A, check_finite=False)
        elif name == "(A*B)^T = B^T * A^T":
            lhs_numpy = (A @ B).T
            rhs_numpy = B.T @ A.T
        elif name == "(A^T)^-1 = (A^-1)^T":
            lhs_numpy = inv(A.T, check_finite=False)
            rhs_numpy = inv(A, check_finite=False).T
        elif name == "A * A^-1 = I":
            lhs_numpy = A @ inv(A, check_finite=False)
            rhs_numpy = np.eye(n)
        elif name == "det(A*B) = det(A) * det(B)":
            lhs_numpy = det(A @ B, check_finite=False)
            rhs_numpy = det(A, check_finite=False) * det(B, check_finite=False)
        elif name == "det(A^-1) = 1/det(A)":
            lhs_numpy = det(inv(A, check_finite=False), check_finite=False)
            rhs_numpy = 1.0 / det(A, check_finite=False)
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
            k = data['inputs']['k']
            lhs_numpy = det(k * A, check_finite=False)
            rhs_numpy = (k ** n) * det(A, check_finite=False)
        else:
            print(f"Unknown identity: {name}")
            return
        # print(f"LHS: {lhs_numpy}")
        # print(f"\nRHS: {rhs_numpy}")
        diff_numpy = np.abs(lhs_numpy - rhs_numpy)
        max_diff_numpy = np.nanmax(diff_numpy) if isinstance(diff_numpy, np.ndarray) else np.abs(lhs_numpy - rhs_numpy)
        print(f"LHS vs RHS Max Diff: {max_diff_numpy}")
        
    except (ValueError, np.linalg.LinAlgError) as e:
        print(f"Error during verification: {e}")
    except Exception as e:
        print(f"Unexpected error: {e}")

    print("-" * 40)

def main():
    solution_files = glob.glob("Solutions/*.pkl")
    if not solution_files:
        return

    for file_path in solution_files:
        if "Gurobi_Solution" in file_path:
            print("Gurobi")
        else:
            print("Z3")
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
            verify_identity_numpy(data['name'], data, data['n'])

if __name__ == "__main__":
    main()
