# Symbolic FP Error Detection

This project explores floating-point error detection in linear algebra identities using symbolic solvers like Z3 and Gurobi.

## Matrix Identities - Z3 (Floating Point)

| Identity | Solution? | Time (s) | Size of N |
| :--- | :--- | :--- | :--- |
| `(A+B)^T = A^T + B^T` | | | |
| `(A^-1)^-1 = A` | Yes | 20 | 2 |
| `(A*B)^-1 = B^-1 * A^-1` | Yes | 20 | 2 |
| `(A*B*C)^-1 = C^-1 * B^-1 * A^-1` | | | |
| `(A*B)^T = B^T * A^T` | | | |
| `(A^T)^-1 = (A^-1)^T` | | | |
| `A * A^-1 = I` |  Yes | 20 | 2 |
| `det(A*B) = det(A) * det(B)` | Yes | 20 | 3 |
| `det(A^-1) = 1/det(A)` | Yes | 20 | 2 |
| `A(BC) = (AB)C` | | | |
| `A(B+C) = AB + AC` | | | |

## Matrix Identities - Gurobi

| Identity | Solution? | Time (s) | Size of N |
| :--- | :--- | :--- | :--- |
| `(A+B)^T = A^T + B^T` | | | |
| `(A^-1)^-1 = A` | | | |
| `(A*B)^-1 = B^-1 * A^-1` | | | |
| `(A*B*C)^-1 = C^-1 * B^-1 * A^-1` | | | |
| `(A*B)^T = B^T * A^T` | | | |
| `(A^T)^-1 = (A^-1)^T` | | | |
| `A * A^-1 = I` | | | |
| `det(A*B) = det(A) * det(B)` | | | |
| `det(A^-1) = 1/det(A)` | | | |
| `A(BC) = (AB)C` | | | |
| `A(B+C) = AB + AC` | | | |

## How to run

1. Run Gurobi verification: `python LinearAlgebraIdentities_Gurobi.py`
2. Run Z3 verification: `python LinearAlgebraIdentities.py`
3. Verify stored counter-examples: `python verify_solutions.py`
