# Symbolic FP Error Detection

This project explores floating-point error detection in linear algebra identities using symbolic solvers like Z3 and Gurobi.

## Matrix Identities


| Identity | Gurobi |  |  |  | Z3 |  |  |  |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
|  | Solution? | Time (s) | Size of N | Max Diff | Solution? | Time (s) | Size of N | Max Diff |
| A * A^-1 = I | Yes | 7.48 | 2 | 0.00195 | Yes | 28.35 | 2 | 1.11e-16 |
| (A^-1)^-1 = A | Yes | 0.06 | 2 | 1.36e-20 | Yes | 88.55 | 2 | 8.88e-16 |
| (A*B)^-1 = B^-1 * A^-1 | Yes | 8.4 | 2 | 9.31e-10 | Yes | 255.39 | 2 | 1.11e-16 |
| (A*B*C)^-1 = C^-1 * B^-1 * A^-1 |  |  |  |  |  |  |  |  |
| A(BC) = (AB)C | Yes | 0.33 | 2,3 | 1.82e-12 | Yes | 238.07 | 3 | 0.375 |
| det(A*B) = det(A) * det(B) | No |  |  |  | Yes | 7.3 | 2 | 8.88e-16 |
| det(A^-1) = 1/det(A) | Yes | 0.74 | 2 | 3.64e-11 | Yes | 5.5 | 2 | 1.11e-16 |
| A(B+C) = AB + AC |  |  |  |  | Yes | 346 | 3 | 4.66e-10 |
| (A+B)^T = A^T + B^T | No |  |  |  | No |  |  |  |
| (A*B)^T = B^T * A^T | No |  |  |  | No |  |  |  |
| (A^T)^-1 = (A^-1)^T | No |  |  |  |  |  |  |  |
| tr(A + B) = tr(A) + tr(B) | No |  |  |  | Yes | 0.31 | 2 | 0 |
| tr(AB) = tr(BA) |  No |  |  |  | Yes | 6.98 | 2 | 0 |
| det(kA) = k^n * det(A) |  Yes |  | 2 | 1.11e-16 | Yes | 154 | 3 | 1.36e-20 |



## How to run

1. Run Gurobi verification: `python LinearAlgebraIdentities_Gurobi.py`
2. Run Z3 verification: `python LinearAlgebraIdentities.py`
3. Verify stored counter-examples: `python VerifyPickle.py`
