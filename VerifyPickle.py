import pickle
import numpy as np
import os
import glob

def verify_identity(name, data, n):
    print(f"--- Verifying with NumPy (N={n}): {name} ---")
    LHS_val = np.array(data['LHS'])
    RHS_val = np.array(data['RHS'])
    
    diff = LHS_val - RHS_val
    if isinstance(diff, np.ndarray):
        max_diff = np.max(np.abs(diff))
    else:
        max_diff = abs(diff)
        
    print(f"Max Difference: {max_diff}")
    if max_diff > 1e-7:
        print("Counter-example confirmed by NumPy.")
    else:
        print("NumPy agrees with the identity (within tolerance).")
    
    for m_name, val in data['inputs'].items():
        print(f"Matrix {m_name}:\n{np.array(val)}")
    print(f"LHS:\n{LHS_val}")
    print(f"RHS:\n{RHS_val}")
    print("-" * 40)

def main():
    solution_files = glob.glob("solutions/*.pkl")
    if not solution_files:
        print("No solution files found in 'solutions/' directory.")
        return

    for file_path in solution_files:
        print(f"Reading {file_path}...")
        with open(file_path, 'rb') as f:
            data = pickle.load(f)
            verify_identity(data['name'], data, data['n'])

if __name__ == "__main__":
    main()
