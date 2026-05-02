# Symbolic FP Error Detection

This project explores floating-point error detection in linear algebra identities using symbolic solvers like Z3 and Gurobi.

## Matrix Identities


| Identity | Gurobi |  |  |  | Z3 |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | Solution? | Time (s) | Size of N | Max Diff | Solution? | Time (s) | Size of N | Max Diff |
| A * A^-1 = I | Yes | 7.48 | 2 | 0.00195 | Yes | 28.35 | 2 |  |
| (A^-1)^-1 = A | Yes | 0.06 | 2 | 1.36e-20 | Yes | 88.55 | 2 |  |
| (A*B)^-1 = B^-1 * A^-1 | Yes | 8.4 | 2 | 9.31e-10 | Yes | 255.39 | 2 |  |
| (A*B*C)^-1 = C^-1 * B^-1 * A^-1 |  |  |  |  |  |  |  |  |
| A(BC) = (AB)C | Yes | 0.33 | 2,3 | 1.82e-12 | Yes | 57.78 | 2 |  |
| det(A*B) = det(A) * det(B) | No |  |  |  |  |  |  |  |
| det(A^-1) = 1/det(A) | Yes | 0.74 | 2 | 3.64e-11 | Yes | 5.5 | 2 |  |
| A(B+C) = AB + AC |  |  |  |  |  | 9.25 | 2 |  |
| (A+B)^T = A^T + B^T | No |  |  |  | No |  |  |  |
| (A*B)^T = B^T * A^T | No |  |  |  |  |  |  |  |
| (A^T)^-1 = (A^-1)^T | No |  |  |  |  |  |  |  |
| tr(A + B) = tr(A) + tr(B) | No |  |  |  | Yes | 0.31 | 2 |  |
| tr(AB) = tr(BA) |  No |  |  |  | Yes | 6.98 | 2 |  |
| det(kA) = k^n * det(A) |  Yes |  | 2 | 1.11e-16 | Yes | 3.37 | 2 |  |



## How to run

1. Run Gurobi verification: `python LinearAlgebraIdentities_Gurobi.py`
2. Run Z3 verification: `python LinearAlgebraIdentities.py`
3. Verify stored counter-examples: `python VerifyPickle.py`
