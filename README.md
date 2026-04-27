# Symbolic FP Error Detection

This project explores floating-point error detection in linear algebra identities using symbolic solvers like Z3 and Gurobi.

## Matrix Identities - Gurobi

| Identity | Solution? | Time (s) | Size of N |
| :--- | :--- | :--- | :--- |
| `A * A^-1 = I` |  Yes | 7.48 | 2 |
| `(A^-1)^-1 = A` | Yes | 0.06 | 2 |
| `(A*B)^-1 = B^-1 * A^-1` | Yes | 8.4 | 2 |
| `(A*B*C)^-1 = C^-1 * B^-1 * A^-1` | | | |
| `A(BC) = (AB)C` | Yes | 0.33 | 2,3 |
| `det(A*B) = det(A) * det(B)` | Yes | 0.48 | 2,3 |
| `det(A^-1) = 1/det(A)` | Yes | 0.74 | 2 |
| `A(B+C) = AB + AC` | | | |
| `(A+B)^T = A^T + B^T` | No | | |
| `(A*B)^T = B^T * A^T` | | | |
| `(A^T)^-1 = (A^-1)^T` | No | | |

## Matrix Identities - Z3

| Identity | Solution? | Time (s) | Size of N |
| :--- | :--- | :--- | :--- |
| `A * A^-1 = I` | | | |
| `(A^-1)^-1 = A` | | | |
| `(A*B)^-1 = B^-1 * A^-1` | | | |
| `(A*B*C)^-1 = C^-1 * B^-1 * A^-1` | | | |
| `A(BC) = (AB)C` | | | |
| `det(A*B) = det(A) * det(B)` | | | |
| `det(A^-1) = 1/det(A)` | | | |
| `A(B+C) = AB + AC` | | | |
| `(A+B)^T = A^T + B^T` | | | |
| `(A*B)^T = B^T * A^T` | | | |
| `(A^T)^-1 = (A^-1)^T` | | | |

## How to run

1. Run Gurobi verification: `python LinearAlgebraIdentities_Gurobi.py`
2. Run Z3 verification: `python LinearAlgebraIdentities.py`
3. Verify stored counter-examples: `python VerifyPickle.py`
