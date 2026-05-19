# Extracted Scipy Functions
# Source: Usage/DynamicFanInResults/Scipy_Dynamic_FanIn_Merged_ExcludingEmpty.xlsx


################################################################################
# Method Path: scipy.sparse.linalg._eigen.lobpcg.lobpcg.lobpcg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_eigen/lobpcg/lobpcg.py
# Target:      dot
# Call Count:  261564
# Total Exec:  942
# Func Size:   941
################################################################################

def lobpcg(
    A,
    X,
    B=None,
    M=None,
    Y=None,
    tol=None,
    maxiter=None,
    largest=True,
    verbosityLevel=0,
    retLambdaHistory=False,
    retResidualNormsHistory=False,
    restartControl=20,
):
    """Locally Optimal Block Preconditioned Conjugate Gradient Method (LOBPCG).

    LOBPCG is a preconditioned eigensolver for large real symmetric and complex
    Hermitian definite generalized eigenproblems.

    Parameters
    ----------
    A : {sparse matrix, ndarray, LinearOperator, callable object}
        The Hermitian linear operator of the problem, usually given by a
        sparse matrix.  Often called the "stiffness matrix".
    X : ndarray, float32 or float64
        Initial approximation to the ``k`` eigenvectors (non-sparse).
        If `A` has ``shape=(n,n)`` then `X` must have ``shape=(n,k)``.
    B : {sparse matrix, ndarray, LinearOperator, callable object}
        Optional. By default ``B = None``, which is equivalent to identity.
        The right hand side operator in a generalized eigenproblem if present.
        Often called the "mass matrix". Must be Hermitian positive definite.
    M : {sparse matrix, ndarray, LinearOperator, callable object}
        Optional. By default ``M = None``, which is equivalent to identity.
        Preconditioner aiming to accelerate convergence.
    Y : ndarray, float32 or float64, default: None
        An ``n-by-sizeY`` ndarray of constraints with ``sizeY < n``.
        The iterations will be performed in the ``B``-orthogonal complement
        of the column-space of `Y`. `Y` must be full rank if present.
    tol : scalar, optional
        The default is ``tol=n*sqrt(eps)``.
        Solver tolerance for the stopping criterion.
    maxiter : int, default: 20
        Maximum number of iterations.
    largest : bool, default: True
        When True, solve for the largest eigenvalues, otherwise the smallest.
    verbosityLevel : int, optional
        By default ``verbosityLevel=0`` no output.
        Controls the solver standard/screen output.
    retLambdaHistory : bool, default: False
        Whether to return iterative eigenvalue history.
    retResidualNormsHistory : bool, default: False
        Whether to return iterative history of residual norms.
    restartControl : int, optional
        Iterations restart if the residuals jump ``2**restartControl`` times
        compared to the smallest recorded in ``retResidualNormsHistory``.
        The default is ``restartControl=20``, making the restarts rare for
        backward compatibility.

    Returns
    -------
    lambda : ndarray of the shape ``(k, )``.
        Array of ``k`` approximate eigenvalues.
    v : ndarray of the same shape as ``X.shape``.
        An array of ``k`` approximate eigenvectors.
    lambdaHistory : ndarray, optional.
        The eigenvalue history, if `retLambdaHistory` is ``True``.
    ResidualNormsHistory : ndarray, optional.
        The history of residual norms, if `retResidualNormsHistory`
        is ``True``.

    Notes
    -----
    The iterative loop runs ``maxit=maxiter`` (20 if ``maxit=None``)
    iterations at most and finishes earlier if the tolerance is met.
    Breaking backward compatibility with the previous version, LOBPCG
    now returns the block of iterative vectors with the best accuracy rather
    than the last one iterated, as a cure for possible divergence.

    If ``X.dtype == np.float32`` and user-provided operations/multiplications
    by `A`, `B`, and `M` all preserve the ``np.float32`` data type,
    all the calculations and the output are in ``np.float32``.

    The size of the iteration history output equals to the number of the best
    (limited by `maxit`) iterations plus 3: initial, final, and postprocessing.

    If both `retLambdaHistory` and `retResidualNormsHistory` are ``True``,
    the return tuple has the following format
    ``(lambda, V, lambda history, residual norms history)``.

    In the following ``n`` denotes the matrix size and ``k`` the number
    of required eigenvalues (smallest or largest).

    The LOBPCG code internally solves eigenproblems of the size ``3k`` on every
    iteration by calling the dense eigensolver `eigh`, so if ``k`` is not
    small enough compared to ``n``, it makes no sense to call the LOBPCG code.
    Moreover, if one calls the LOBPCG algorithm for ``5k > n``, it would likely
    break internally, so the code calls the standard function `eigh` instead.
    It is not that ``n`` should be large for the LOBPCG to work, but rather the
    ratio ``n / k`` should be large. It you call LOBPCG with ``k=1``
    and ``n=10``, it works though ``n`` is small. The method is intended
    for extremely large ``n / k``.

    The convergence speed depends basically on three factors:

    1. Quality of the initial approximations `X` to the seeking eigenvectors.
       Randomly distributed around the origin vectors work well if no better
       choice is known.

    2. Relative separation of the desired eigenvalues from the rest
       of the eigenvalues. One can vary ``k`` to improve the separation.

    3. Proper preconditioning to shrink the spectral spread.
       For example, a rod vibration test problem (under tests
       directory) is ill-conditioned for large ``n``, so convergence will be
       slow, unless efficient preconditioning is used. For this specific
       problem, a good simple preconditioner function would be a linear solve
       for `A`, which is easy to code since `A` is tridiagonal.

    References
    ----------
    .. [1] A. V. Knyazev (2001),
           Toward the Optimal Preconditioned Eigensolver: Locally Optimal
           Block Preconditioned Conjugate Gradient Method.
           SIAM Journal on Scientific Computing 23, no. 2,
           pp. 517-541. :doi:`10.1137/S1064827500366124`

    .. [2] A. V. Knyazev, I. Lashuk, M. E. Argentati, and E. Ovchinnikov
           (2007), Block Locally Optimal Preconditioned Eigenvalue Xolvers
           (BLOPEX) in hypre and PETSc. :arxiv:`0705.2626`

    .. [3] A. V. Knyazev's C and MATLAB implementations:
           https://github.com/lobpcg/blopex

    Examples
    --------
    Our first example is minimalistic - find the largest eigenvalue of
    a diagonal matrix by solving the non-generalized eigenvalue problem
    ``A x = lambda x`` without constraints or preconditioning.

    >>> import numpy as np
    >>> from scipy.sparse import diags_array
    >>> from scipy.sparse.linalg import LinearOperator, aslinearoperator
    >>> from scipy.sparse.linalg import lobpcg

    The square matrix size is

    >>> n = 100

    and its diagonal entries are 1, ..., 100 defined by

    >>> vals = np.arange(1, n + 1).astype(np.int16)

    The first mandatory input parameter in this test is
    the sparse diagonal matrix `A`
    of the eigenvalue problem ``A x = lambda x`` to solve.

    >>> A = diags_array(vals, offsets=0, shape=(n, n), dtype=None)
    >>> A.toarray()
    array([[  1,   0,   0, ...,   0,   0,   0],
           [  0,   2,   0, ...,   0,   0,   0],
           [  0,   0,   3, ...,   0,   0,   0],
           ...,
           [  0,   0,   0, ...,  98,   0,   0],
           [  0,   0,   0, ...,   0,  99,   0],
           [  0,   0,   0, ...,   0,   0, 100]], shape=(100, 100), dtype=int16)

    The second mandatory input parameter `X` is a 2D array with the
    row dimension determining the number of requested eigenvalues.
    `X` is an initial guess for targeted eigenvectors.
    `X` must have linearly independent columns.
    If no initial approximations available, randomly oriented vectors
    commonly work best, e.g., with components normally distributed
    around zero or uniformly distributed on the interval [-1 1].
    Setting the initial approximations to dtype ``np.float32``
    forces all iterative values to dtype ``np.float32`` speeding up
    the run while still allowing accurate eigenvalue computations.

    >>> k = 1
    >>> rng = np.random.default_rng()
    >>> X = rng.normal(size=(n, k))
    >>> X = X.astype(np.float32)

    >>> eigenvalues, _ = lobpcg(A, X, maxiter=60)
    >>> eigenvalues
    array([100.], dtype=float32)

    `lobpcg` needs only access the matrix product with `A` rather
    then the matrix itself. Since the matrix `A` is diagonal in
    this example, one can write a function of the matrix product
    ``A @ X`` using the diagonal values ``vals`` only, e.g., by
    element-wise multiplication with broadcasting in the lambda-function

    >>> A_lambda = lambda X: vals[:, np.newaxis] * X

    or the regular function

    >>> def A_matmat(X):
    ...     return vals[:, np.newaxis] * X

    and use the handle to one of these callables as an input

    >>> eigenvalues, _ = lobpcg(A_lambda, X, maxiter=60)
    >>> eigenvalues
    array([100.], dtype=float32)
    >>> eigenvalues, _ = lobpcg(A_matmat, X, maxiter=60)
    >>> eigenvalues
    array([100.], dtype=float32)

    The traditional callable `LinearOperator` is no longer
    necessary but still supported as the input to `lobpcg`.
    Specifying ``matmat=A_matmat`` explicitly improves performance. 

    >>> A_lo = LinearOperator((n, n), matvec=A_matmat, matmat=A_matmat, dtype=np.int16)
    >>> eigenvalues, _ = lobpcg(A_lo, X, maxiter=80)
    >>> eigenvalues
    array([100.], dtype=float32)

    The least efficient callable option is `aslinearoperator`:

    >>> eigenvalues, _ = lobpcg(aslinearoperator(A), X, maxiter=80)
    >>> eigenvalues
    array([100.], dtype=float32)

    We now switch to computing the three smallest eigenvalues specifying

    >>> k = 3
    >>> X = np.random.default_rng().normal(size=(n, k))

    and ``largest=False`` parameter

    >>> eigenvalues, _ = lobpcg(A, X, largest=False, maxiter=90)
    >>> print(eigenvalues)  
    [1. 2. 3.]

    The next example illustrates computing 3 smallest eigenvalues of
    the same matrix `A` given by the function handle ``A_matmat`` but
    with constraints and preconditioning.

    Constraints - an optional input parameter is a 2D array comprising
    of column vectors that the eigenvectors must be orthogonal to

    >>> Y = np.eye(n, 3)

    The preconditioner acts as the inverse of `A` in this example, but
    in the reduced precision ``np.float32`` even though the initial `X`
    and thus all iterates and the output are in full ``np.float64``.

    >>> inv_vals = 1./vals
    >>> inv_vals = inv_vals.astype(np.float32)
    >>> M = lambda X: inv_vals[:, np.newaxis] * X

    Let us now solve the eigenvalue problem for the matrix `A` first
    without preconditioning requesting 80 iterations

    >>> eigenvalues, _ = lobpcg(A_matmat, X, Y=Y, largest=False, maxiter=80)
    >>> eigenvalues
    array([4., 5., 6.])
    >>> eigenvalues.dtype
    dtype('float64')

    With preconditioning we need only 20 iterations from the same `X`

    >>> eigenvalues, _ = lobpcg(A_matmat, X, Y=Y, M=M, largest=False, maxiter=20)
    >>> eigenvalues
    array([4., 5., 6.])

    Note that the vectors passed in `Y` are the eigenvectors of the 3
    smallest eigenvalues. The results returned above are orthogonal to those.

    The primary matrix `A` may be indefinite, e.g., after shifting
    ``vals`` by 50 from 1, ..., 100 to -49, ..., 50, we still can compute
    the 3 smallest or largest eigenvalues.

    >>> vals = vals - 50
    >>> X = rng.normal(size=(n, k))
    >>> eigenvalues, _ = lobpcg(A_matmat, X, largest=False, maxiter=99)
    >>> eigenvalues
    array([-49., -48., -47.])
    >>> eigenvalues, _ = lobpcg(A_matmat, X, largest=True, maxiter=99)
    >>> eigenvalues
    array([50., 49., 48.])

    """
    blockVectorX = X
    bestblockVectorX = blockVectorX
    blockVectorY = Y
    residualTolerance = tol
    if maxiter is None:
        maxiter = 20

    bestIterationNumber = maxiter

    sizeY = 0
    if blockVectorY is not None:
        if len(blockVectorY.shape) != 2:
            warnings.warn(
                f"Expected rank-2 array for argument Y, instead got "
                f"{len(blockVectorY.shape)}, "
                f"so ignore it and use no constraints.",
                UserWarning, stacklevel=2
            )
            blockVectorY = None
        else:
            sizeY = blockVectorY.shape[1]

    # Block size.
    if blockVectorX is None:
        raise ValueError("The mandatory initial matrix X cannot be None")
    if len(blockVectorX.shape) != 2:
        raise ValueError("expected rank-2 array for argument X")

    n, sizeX = blockVectorX.shape

    # Data type of iterates, determined by X, must be inexact
    if not np.issubdtype(blockVectorX.dtype, np.inexact):
        warnings.warn(
            f"Data type for argument X is {blockVectorX.dtype}, "
            f"which is not inexact, so casted to np.float32.",
            UserWarning, stacklevel=2
        )
        blockVectorX = np.asarray(blockVectorX, dtype=np.float32)

    if retLambdaHistory:
        lambdaHistory = np.zeros((maxiter + 3, sizeX),
                                 dtype=blockVectorX.dtype)
    if retResidualNormsHistory:
        residualNormsHistory = np.zeros((maxiter + 3, sizeX),
                                        dtype=blockVectorX.dtype)

    if verbosityLevel:
        aux = "Solving "
        if B is None:
            aux += "standard"
        else:
            aux += "generalized"
        aux += " eigenvalue problem with"
        if M is None:
            aux += "out"
        aux += " preconditioning\n\n"
        aux += f"matrix size {n}\n"
        aux += f"block size {sizeX}\n\n"
        if blockVectorY is None:
            aux += "No constraints\n\n"
        else:
            if sizeY > 1:
                aux += f"{sizeY} constraints\n\n"
            else:
                aux += f"{sizeY} constraint\n\n"
        print(aux)

    if (n - sizeY) < (5 * sizeX):
        warnings.warn(
            f"The problem size {n} minus the constraints size {sizeY} "
            f"is too small relative to the block size {sizeX}. "
            f"Using a dense eigensolver instead of LOBPCG iterations."
            f"No output of the history of the iterations.",
            UserWarning, stacklevel=2
        )

        sizeX = min(sizeX, n)

        if blockVectorY is not None:
            raise NotImplementedError(
                "The dense eigensolver does not support constraints."
            )

        # Define the closed range of indices of eigenvalues to return.
        if largest:
            eigvals = (n - sizeX, n - 1)
        else:
            eigvals = (0, sizeX - 1)

        try:
            if isinstance(A, LinearOperator):
                A = A(np.eye(n, dtype=int))
            elif callable(A):
                A = A(np.eye(n, dtype=int))
                if A.shape != (n, n):
                    raise ValueError(
                        f"The shape {A.shape} of the primary matrix\n"
                        f"defined by a callable object is wrong.\n"
                        f"Expected {(n, n)}."
                    )
            elif issparse(A):
                A = A.toarray()
            else:
                A = np.asarray(A)
        except Exception as e:
            raise Exception(
                f"Primary MatMul call failed with error\n"
                f"{e}\n")

        if B is not None:
            try:
                if isinstance(B, LinearOperator):
                    B = B(np.eye(n, dtype=int))
                elif callable(B):
                    B = B(np.eye(n, dtype=int))
                    if B.shape != (n, n):
                        raise ValueError(
                            f"The shape {B.shape} of the secondary matrix\n"
                            f"defined by a callable object is wrong.\n"
                        )
                elif issparse(B):
                    B = B.toarray()
                else:
                    B = np.asarray(B)
            except Exception as e:
                raise Exception(
                    f"Secondary MatMul call failed with error\n"
                    f"{e}\n")

        try:
            vals, vecs = eigh(A,
                              B,
                              subset_by_index=eigvals,
                              check_finite=False)
            if largest:
                # Reverse order to be compatible with eigs() in 'LM' mode.
                vals = vals[::-1]
                vecs = vecs[:, ::-1]

            return vals, vecs
        except Exception as e:
            raise Exception(
                f"Dense eigensolver failed with error\n"
                f"{e}\n"
            )

    if (residualTolerance is None) or (residualTolerance <= 0.0):
        residualTolerance = np.sqrt(np.finfo(blockVectorX.dtype).eps) * n

    A = _makeMatMat(A)
    B = _makeMatMat(B)
    M = _makeMatMat(M)

    # Apply constraints to X.
    if blockVectorY is not None:

        if B is not None:
            blockVectorBY = B(blockVectorY)
            if blockVectorBY.shape != blockVectorY.shape:
                raise ValueError(
                    f"The shape {blockVectorY.shape} "
                    f"of the constraint not preserved\n"
                    f"and changed to {blockVectorBY.shape} "
                    f"after multiplying by the secondary matrix.\n"
                )
        else:
            blockVectorBY = blockVectorY

        # gramYBY is a dense array.
        gramYBY = blockVectorY.T.conj() @ blockVectorBY
        try:
            # gramYBY is a Cholesky factor from now on...
            gramYBY = cho_factor(gramYBY, overwrite_a=True)
        except LinAlgError as e:
            raise ValueError("Linearly dependent constraints") from e

        _applyConstraints(blockVectorX, gramYBY, blockVectorBY, blockVectorY)

    ##
    # B-orthonormalize X.
    blockVectorX, blockVectorBX, _ = _b_orthonormalize(
        B, blockVectorX, verbosityLevel=verbosityLevel)
    if blockVectorX is None:
        raise ValueError("Linearly dependent initial approximations")

    ##
    # Compute the initial Ritz vectors: solve the eigenproblem.
    blockVectorAX = A(blockVectorX)
    if blockVectorAX.shape != blockVectorX.shape:
        raise ValueError(
            f"The shape {blockVectorX.shape} "
            f"of the initial approximations not preserved\n"
            f"and changed to {blockVectorAX.shape} "
            f"after multiplying by the primary matrix.\n"
        )

    gramXAX = blockVectorX.T.conj() @ blockVectorAX

    _lambda, eigBlockVector = eigh(gramXAX, check_finite=False)
    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    if retLambdaHistory:
        lambdaHistory[0, :] = _lambda

    eigBlockVector = np.asarray(eigBlockVector[:, ii])
    blockVectorX = _matmul_inplace(
        blockVectorX, eigBlockVector,
        verbosityLevel=verbosityLevel
    )
    blockVectorAX = _matmul_inplace(
        blockVectorAX, eigBlockVector,
        verbosityLevel=verbosityLevel
    )
    if B is not None:
        blockVectorBX = _matmul_inplace(
            blockVectorBX, eigBlockVector,
            verbosityLevel=verbosityLevel
        )

    ##
    # Active index set.
    activeMask = np.ones((sizeX,), dtype=bool)

    ##
    # Main iteration loop.

    blockVectorP = None  # set during iteration
    blockVectorAP = None
    blockVectorBP = None

    smallestResidualNorm = np.abs(np.finfo(blockVectorX.dtype).max)

    iterationNumber = -1
    restart = True
    forcedRestart = False
    explicitGramFlag = False
    while iterationNumber < maxiter:
        iterationNumber += 1

        if B is not None:
            aux = blockVectorBX * _lambda[np.newaxis, :]
        else:
            aux = blockVectorX * _lambda[np.newaxis, :]

        blockVectorR = blockVectorAX - aux

        aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
        residualNorms = np.sqrt(np.abs(aux))
        if retResidualNormsHistory:
            residualNormsHistory[iterationNumber, :] = residualNorms
        residualNorm = np.sum(np.abs(residualNorms)) / sizeX

        if residualNorm < smallestResidualNorm:
            smallestResidualNorm = residualNorm
            bestIterationNumber = iterationNumber
            bestblockVectorX = blockVectorX
        elif residualNorm > 2**restartControl * smallestResidualNorm:
            forcedRestart = True
            blockVectorAX = A(blockVectorX)
            if blockVectorAX.shape != blockVectorX.shape:
                raise ValueError(
                    f"The shape {blockVectorX.shape} "
                    f"of the restarted iterate not preserved\n"
                    f"and changed to {blockVectorAX.shape} "
                    f"after multiplying by the primary matrix.\n"
                )
            if B is not None:
                blockVectorBX = B(blockVectorX)
                if blockVectorBX.shape != blockVectorX.shape:
                    raise ValueError(
                        f"The shape {blockVectorX.shape} "
                        f"of the restarted iterate not preserved\n"
                        f"and changed to {blockVectorBX.shape} "
                        f"after multiplying by the secondary matrix.\n"
                    )

        ii = np.where(residualNorms > residualTolerance, True, False)
        activeMask = activeMask & ii
        currentBlockSize = activeMask.sum()

        if verbosityLevel:
            print(f"iteration {iterationNumber}")
            print(f"current block size: {currentBlockSize}")
            print(f"eigenvalue(s):\n{_lambda}")
            print(f"residual norm(s):\n{residualNorms}")

        if currentBlockSize == 0:
            break

        activeBlockVectorR = _as2d(blockVectorR[:, activeMask])

        if iterationNumber > 0:
            activeBlockVectorP = _as2d(blockVectorP[:, activeMask])
            activeBlockVectorAP = _as2d(blockVectorAP[:, activeMask])
            if B is not None:
                activeBlockVectorBP = _as2d(blockVectorBP[:, activeMask])

        if M is not None:
            # Apply preconditioner T to the active residuals.
            activeBlockVectorR = M(activeBlockVectorR)

        ##
        # Apply constraints to the preconditioned residuals.
        if blockVectorY is not None:
            _applyConstraints(activeBlockVectorR,
                              gramYBY,
                              blockVectorBY,
                              blockVectorY)

        ##
        # B-orthogonalize the preconditioned residuals to X.
        if B is not None:
            activeBlockVectorR = activeBlockVectorR - (
                blockVectorX @
                (blockVectorBX.T.conj() @ activeBlockVectorR)
            )
        else:
            activeBlockVectorR = activeBlockVectorR - (
                blockVectorX @
                (blockVectorX.T.conj() @ activeBlockVectorR)
            )

        ##
        # B-orthonormalize the preconditioned residuals.
        aux = _b_orthonormalize(
            B, activeBlockVectorR, verbosityLevel=verbosityLevel)
        activeBlockVectorR, activeBlockVectorBR, _ = aux

        if activeBlockVectorR is None:
            warnings.warn(
                f"Failed at iteration {iterationNumber} with accuracies "
                f"{residualNorms}\n not reaching the requested "
                f"tolerance {residualTolerance}.",
                UserWarning, stacklevel=2
            )
            break
        activeBlockVectorAR = A(activeBlockVectorR)

        if iterationNumber > 0:
            if B is not None:
                aux = _b_orthonormalize(
                    B, activeBlockVectorP, activeBlockVectorBP,
                    verbosityLevel=verbosityLevel
                )
                activeBlockVectorP, activeBlockVectorBP, invR = aux
            else:
                aux = _b_orthonormalize(B, activeBlockVectorP,
                                        verbosityLevel=verbosityLevel)
                activeBlockVectorP, _, invR = aux
            # Function _b_orthonormalize returns None if Cholesky fails
            if activeBlockVectorP is not None:
                activeBlockVectorAP = _matmul_inplace(
                    activeBlockVectorAP, invR,
                    verbosityLevel=verbosityLevel
                )
                restart = forcedRestart
            else:
                restart = True

        ##
        # Perform the Rayleigh Ritz Procedure:
        # Compute symmetric Gram matrices:

        if activeBlockVectorAR.dtype == "float32":
            myeps = 1
        else:
            myeps = np.sqrt(np.finfo(activeBlockVectorR.dtype).eps)

        if residualNorms.max() > myeps and not explicitGramFlag:
            explicitGramFlag = False
        else:
            # Once explicitGramFlag, forever explicitGramFlag.
            explicitGramFlag = True

        # Shared memory assignments to simplify the code
        if B is None:
            blockVectorBX = blockVectorX
            activeBlockVectorBR = activeBlockVectorR
            if not restart:
                activeBlockVectorBP = activeBlockVectorP

        # Common submatrices:
        gramXAR = np.dot(blockVectorX.T.conj(), activeBlockVectorAR)
        gramRAR = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAR)

        gramDtype = activeBlockVectorAR.dtype
        if explicitGramFlag:
            gramRAR = (gramRAR + gramRAR.T.conj()) / 2
            gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)
            gramXAX = (gramXAX + gramXAX.T.conj()) / 2
            gramXBX = np.dot(blockVectorX.T.conj(), blockVectorBX)
            gramRBR = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorBR)
            gramXBR = np.dot(blockVectorX.T.conj(), activeBlockVectorBR)
        else:
            gramXAX = np.diag(_lambda).astype(gramDtype)
            gramXBX = np.eye(sizeX, dtype=gramDtype)
            gramRBR = np.eye(currentBlockSize, dtype=gramDtype)
            gramXBR = np.zeros((sizeX, currentBlockSize), dtype=gramDtype)

        if not restart:
            gramXAP = np.dot(blockVectorX.T.conj(), activeBlockVectorAP)
            gramRAP = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorAP)
            gramPAP = np.dot(activeBlockVectorP.T.conj(), activeBlockVectorAP)
            gramXBP = np.dot(blockVectorX.T.conj(), activeBlockVectorBP)
            gramRBP = np.dot(activeBlockVectorR.T.conj(), activeBlockVectorBP)
            if explicitGramFlag:
                gramPAP = (gramPAP + gramPAP.T.conj()) / 2
                gramPBP = np.dot(activeBlockVectorP.T.conj(),
                                 activeBlockVectorBP)
            else:
                gramPBP = np.eye(currentBlockSize, dtype=gramDtype)

            gramA = np.block(
                [
                    [gramXAX, gramXAR, gramXAP],
                    [gramXAR.T.conj(), gramRAR, gramRAP],
                    [gramXAP.T.conj(), gramRAP.T.conj(), gramPAP],
                ]
            )
            gramB = np.block(
                [
                    [gramXBX, gramXBR, gramXBP],
                    [gramXBR.T.conj(), gramRBR, gramRBP],
                    [gramXBP.T.conj(), gramRBP.T.conj(), gramPBP],
                ]
            )

            _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel)

            try:
                _lambda, eigBlockVector = eigh(gramA,
                                               gramB,
                                               check_finite=False)
            except LinAlgError as e:
                # raise ValueError("eigh failed in lobpcg iterations") from e
                if verbosityLevel:
                    warnings.warn(
                        f"eigh failed at iteration {iterationNumber} \n"
                        f"with error {e} causing a restart.\n",
                        UserWarning, stacklevel=2
                    )
                # try again after dropping the direction vectors P from RR
                restart = True

        if restart:
            gramA = np.block([[gramXAX, gramXAR], [gramXAR.T.conj(), gramRAR]])
            gramB = np.block([[gramXBX, gramXBR], [gramXBR.T.conj(), gramRBR]])

            _handle_gramA_gramB_verbosity(gramA, gramB, verbosityLevel)

            try:
                _lambda, eigBlockVector = eigh(gramA,
                                               gramB,
                                               check_finite=False)
            except LinAlgError as e:
                # raise ValueError("eigh failed in lobpcg iterations") from e
                warnings.warn(
                    f"eigh failed at iteration {iterationNumber} with error\n"
                    f"{e}\n",
                    UserWarning, stacklevel=2
                )
                break

        ii = _get_indx(_lambda, sizeX, largest)
        _lambda = _lambda[ii]
        eigBlockVector = eigBlockVector[:, ii]
        if retLambdaHistory:
            lambdaHistory[iterationNumber + 1, :] = _lambda

        # Compute Ritz vectors.
        if B is not None:
            if not restart:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:
                                                 sizeX + currentBlockSize]
                eigBlockVectorP = eigBlockVector[sizeX + currentBlockSize:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                pp += np.dot(activeBlockVectorP, eigBlockVectorP)

                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                app += np.dot(activeBlockVectorAP, eigBlockVectorP)

                bpp = np.dot(activeBlockVectorBR, eigBlockVectorR)
                bpp += np.dot(activeBlockVectorBP, eigBlockVectorP)
            else:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                bpp = np.dot(activeBlockVectorBR, eigBlockVectorR)

            blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
            blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app
            blockVectorBX = np.dot(blockVectorBX, eigBlockVectorX) + bpp

            blockVectorP, blockVectorAP, blockVectorBP = pp, app, bpp

        else:
            if not restart:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:
                                                 sizeX + currentBlockSize]
                eigBlockVectorP = eigBlockVector[sizeX + currentBlockSize:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                pp += np.dot(activeBlockVectorP, eigBlockVectorP)

                app = np.dot(activeBlockVectorAR, eigBlockVectorR)
                app += np.dot(activeBlockVectorAP, eigBlockVectorP)
            else:
                eigBlockVectorX = eigBlockVector[:sizeX]
                eigBlockVectorR = eigBlockVector[sizeX:]

                pp = np.dot(activeBlockVectorR, eigBlockVectorR)
                app = np.dot(activeBlockVectorAR, eigBlockVectorR)

            blockVectorX = np.dot(blockVectorX, eigBlockVectorX) + pp
            blockVectorAX = np.dot(blockVectorAX, eigBlockVectorX) + app

            blockVectorP, blockVectorAP = pp, app

    if B is not None:
        aux = blockVectorBX * _lambda[np.newaxis, :]
    else:
        aux = blockVectorX * _lambda[np.newaxis, :]

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))
    # Use old lambda in case of early loop exit.
    if retLambdaHistory:
        lambdaHistory[iterationNumber + 1, :] = _lambda
    if retResidualNormsHistory:
        residualNormsHistory[iterationNumber + 1, :] = residualNorms
    residualNorm = np.sum(np.abs(residualNorms)) / sizeX
    if residualNorm < smallestResidualNorm:
        smallestResidualNorm = residualNorm
        bestIterationNumber = iterationNumber + 1
        bestblockVectorX = blockVectorX

    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited at iteration {iterationNumber} with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.\n"
            f"Use iteration {bestIterationNumber} instead with accuracy \n"
            f"{smallestResidualNorm}.\n",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final iterative eigenvalue(s):\n{_lambda}")
        print(f"Final iterative residual norm(s):\n{residualNorms}")

    blockVectorX = bestblockVectorX
    # Making eigenvectors "exactly" satisfy the blockVectorY constrains
    if blockVectorY is not None:
        _applyConstraints(blockVectorX,
                          gramYBY,
                          blockVectorBY,
                          blockVectorY)

    # Making eigenvectors "exactly" othonormalized by final "exact" RR
    blockVectorAX = A(blockVectorX)
    if blockVectorAX.shape != blockVectorX.shape:
        raise ValueError(
            f"The shape {blockVectorX.shape} "
            f"of the postprocessing iterate not preserved\n"
            f"and changed to {blockVectorAX.shape} "
            f"after multiplying by the primary matrix.\n"
        )
    gramXAX = np.dot(blockVectorX.T.conj(), blockVectorAX)

    blockVectorBX = blockVectorX
    if B is not None:
        blockVectorBX = B(blockVectorX)
        if blockVectorBX.shape != blockVectorX.shape:
            raise ValueError(
                f"The shape {blockVectorX.shape} "
                f"of the postprocessing iterate not preserved\n"
                f"and changed to {blockVectorBX.shape} "
                f"after multiplying by the secondary matrix.\n"
            )

    gramXBX = np.dot(blockVectorX.T.conj(), blockVectorBX)
    _handle_gramA_gramB_verbosity(gramXAX, gramXBX, verbosityLevel)
    gramXAX = (gramXAX + gramXAX.T.conj()) / 2
    gramXBX = (gramXBX + gramXBX.T.conj()) / 2
    try:
        _lambda, eigBlockVector = eigh(gramXAX,
                                       gramXBX,
                                       check_finite=False)
    except LinAlgError as e:
        raise ValueError("eigh has failed in lobpcg postprocessing") from e

    ii = _get_indx(_lambda, sizeX, largest)
    _lambda = _lambda[ii]
    eigBlockVector = np.asarray(eigBlockVector[:, ii])

    blockVectorX = np.dot(blockVectorX, eigBlockVector)
    blockVectorAX = np.dot(blockVectorAX, eigBlockVector)

    if B is not None:
        blockVectorBX = np.dot(blockVectorBX, eigBlockVector)
        aux = blockVectorBX * _lambda[np.newaxis, :]
    else:
        aux = blockVectorX * _lambda[np.newaxis, :]

    blockVectorR = blockVectorAX - aux

    aux = np.sum(blockVectorR.conj() * blockVectorR, 0)
    residualNorms = np.sqrt(np.abs(aux))

    if retLambdaHistory:
        lambdaHistory[bestIterationNumber + 1, :] = _lambda
    if retResidualNormsHistory:
        residualNormsHistory[bestIterationNumber + 1, :] = residualNorms

    if retLambdaHistory:
        lambdaHistory = lambdaHistory[
            : bestIterationNumber + 2, :]
    if retResidualNormsHistory:
        residualNormsHistory = residualNormsHistory[
            : bestIterationNumber + 2, :]

    if np.max(np.abs(residualNorms)) > residualTolerance:
        warnings.warn(
            f"Exited postprocessing with accuracies \n"
            f"{residualNorms}\n"
            f"not reaching the requested tolerance {residualTolerance}.",
            UserWarning, stacklevel=2
        )

    if verbosityLevel:
        print(f"Final postprocessing eigenvalue(s):\n{_lambda}")
        print(f"Final residual norm(s):\n{residualNorms}")

    if retLambdaHistory:
        lambdaHistory = np.vsplit(lambdaHistory, np.shape(lambdaHistory)[0])
        lambdaHistory = [np.squeeze(i) for i in lambdaHistory]
    if retResidualNormsHistory:
        residualNormsHistory = np.vsplit(residualNormsHistory,
                                         np.shape(residualNormsHistory)[0])
        residualNormsHistory = [np.squeeze(i) for i in residualNormsHistory]

    if retLambdaHistory:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, lambdaHistory, residualNormsHistory
        else:
            return _lambda, blockVectorX, lambdaHistory
    else:
        if retResidualNormsHistory:
            return _lambda, blockVectorX, residualNormsHistory
        else:
            return _lambda, blockVectorX

################################################################################
# Method Path: scipy.sparse.linalg._isolve.iterative.gmres
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/iterative.py
# Target:      dot
# Call Count:  151367
# Total Exec:  1752
# Func Size:   261
################################################################################

def gmres(A, b, x0=None, *, rtol=1e-5, atol=0., restart=None, maxiter=None, M=None,
          callback=None, callback_type=None):
    """
    Solve ``Ax = b`` with the Generalized Minimal RESidual method.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution (a vector of zeros by default).
    rtol, atol : float
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    restart : int, optional
        Number of iterations between restarts. Larger values increase
        iteration cost, but may be necessary for convergence.
        If omitted, ``min(20, n)`` is used.
    maxiter : int, optional
        Maximum number of iterations (restart cycles).  Iteration will stop
        after maxiter steps even if the specified tolerance has not been
        achieved. See `callback_type`.
    M : {sparse array, ndarray, LinearOperator}
        Inverse of the preconditioner of `A`.  `M` should approximate the
        inverse of `A` and be easy to solve for (see Notes).  Effective
        preconditioning dramatically improves the rate of convergence,
        which implies that fewer iterations are needed to reach a given
        error tolerance.  By default, no preconditioner is used.
        In this implementation, left preconditioning is used,
        and the preconditioned residual is minimized. However, the final
        convergence is tested with respect to the ``b - A @ x`` residual.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as ``callback(args)``, where ``args`` are selected by `callback_type`.
    callback_type : {'x', 'pr_norm', 'legacy'}, optional
        Callback function argument requested:
          - ``x``: current iterate (ndarray), called on every restart
          - ``pr_norm``: relative (preconditioned) residual norm (float),
            called on every inner iteration
          - ``legacy`` (default): same as ``pr_norm``, but also changes the
            meaning of `maxiter` to count inner iterations instead of restart
            cycles.

        This keyword has no effect if `callback` is not set.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations

    See Also
    --------
    LinearOperator

    Notes
    -----
    A preconditioner, P, is chosen such that P is close to A but easy to solve
    for. The preconditioner parameter required by this routine is
    ``M = P^-1``. The inverse should preferably not be calculated
    explicitly.  Rather, use the following template to produce M::

      # Construct a linear operator that computes P^-1 @ x.
      import scipy.sparse.linalg as spla
      M_x = lambda x: spla.spsolve(P, x)
      M = spla.LinearOperator((n, n), M_x)

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import gmres
    >>> A = csc_array([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> x, exitCode = gmres(A, b, atol=1e-5)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    if callback is not None and callback_type is None:
        # Warn about 'callback_type' semantic changes.
        # Probably should be removed only in far future, Scipy 2.0 or so.
        msg = ("scipy.sparse.linalg.gmres called without specifying "
               "`callback_type`. The default value will be changed in"
               " a future release. For compatibility, specify a value "
               "for `callback_type` explicitly, e.g., "
               "``gmres(..., callback_type='pr_norm')``, or to retain the "
               "old behavior ``gmres(..., callback_type='legacy')``"
               )
        warnings.warn(msg, category=DeprecationWarning, stacklevel=3)

    if callback_type is None:
        callback_type = 'legacy'

    if callback_type not in ('x', 'pr_norm', 'legacy'):
        raise ValueError(f"Unknown callback_type: {callback_type!r}")

    if callback is None:
        callback_type = None

    A, M, x, b = make_system(A, M, x0, b)
    matvec = A.matvec
    psolve = M.matvec
    n = len(b)
    bnrm2 = np.linalg.norm(b)

    atol, _ = _get_atol_rtol('gmres', bnrm2, atol, rtol)

    if bnrm2 == 0:
        return b, 0

    eps = np.finfo(x.dtype.char).eps

    dotprod = np.vdot if np.iscomplexobj(x) else np.dot

    if maxiter is None:
        maxiter = n*10

    if restart is None:
        restart = 20
    restart = min(restart, n)

    Mb_nrm2 = np.linalg.norm(psolve(b))

    # ====================================================
    # =========== Tolerance control from gh-8400 =========
    # ====================================================
    # Tolerance passed to GMRESREVCOM applies to the inner
    # iteration and deals with the left-preconditioned
    # residual.
    ptol_max_factor = 1.
    ptol = Mb_nrm2 * min(ptol_max_factor, atol / bnrm2)
    presid = 0.
    # ====================================================
    lartg = get_lapack_funcs('lartg', dtype=x.dtype)

    # allocate internal variables
    v = np.empty([restart+1, n], dtype=x.dtype)
    h = np.zeros([restart, restart+1], dtype=x.dtype)
    givens = np.zeros([restart, 2], dtype=x.dtype)

    # legacy iteration count
    inner_iter = 0

    for iteration in range(maxiter):
        if iteration == 0:
            r = b - matvec(x) if x.any() else b.copy()
            if np.linalg.norm(r) < atol:  # Are we done?
                return x, 0

        v[0, :] = psolve(r)
        tmp = np.linalg.norm(v[0, :])
        v[0, :] *= (1 / tmp)
        # RHS of the Hessenberg problem
        S = np.zeros(restart+1, dtype=x.dtype)
        S[0] = tmp

        breakdown = False
        for col in range(restart):
            av = matvec(v[col, :])
            w = psolve(av)

            # Modified Gram-Schmidt
            h0 = np.linalg.norm(w)
            for k in range(col+1):
                tmp = dotprod(v[k, :], w)
                h[col, k] = tmp
                w -= tmp*v[k, :]

            h1 = np.linalg.norm(w)
            h[col, col + 1] = h1
            v[col + 1, :] = w[:]

            # Exact solution indicator
            if h1 <= eps*h0:
                h[col, col + 1] = 0
                breakdown = True
            else:
                v[col + 1, :] *= (1 / h1)

            # apply past Givens rotations to current h column
            for k in range(col):
                c, s = givens[k, 0], givens[k, 1]
                n0, n1 = h[col, [k, k+1]]
                h[col, [k, k + 1]] = [c*n0 + s*n1, -s.conj()*n0 + c*n1]

            # get and apply current rotation to h and S
            c, s, mag = lartg(h[col, col], h[col, col+1])
            givens[col, :] = [c, s]
            h[col, [col, col+1]] = mag, 0

            # S[col+1] component is always 0
            tmp = -np.conjugate(s)*S[col]
            S[[col, col + 1]] = [c*S[col], tmp]
            presid = np.abs(tmp)
            inner_iter += 1

            if callback_type in ('legacy', 'pr_norm'):
                callback(presid / bnrm2)
            # Legacy behavior
            if callback_type == 'legacy' and inner_iter == maxiter:
                break
            if presid <= ptol or breakdown:
                break

        # Solve h(col, col) upper triangular system and allow pseudo-solve
        # singular cases as in (but without the f2py copies):
        # y = trsv(h[:col+1, :col+1].T, S[:col+1])

        if h[col, col] == 0:
            S[col] = 0

        y = np.zeros([col+1], dtype=x.dtype)
        y[:] = S[:col+1]
        for k in range(col, 0, -1):
            if y[k] != 0:
                y[k] /= h[k, k]
                tmp = y[k]
                y[:k] -= tmp*h[k, :k]
        if y[0] != 0:
            y[0] /= h[0, 0]

        x += y @ v[:col+1, :]

        r = b - matvec(x)
        rnorm = np.linalg.norm(r)

        # Legacy exit
        if callback_type == 'legacy' and inner_iter == maxiter:
            return x, 0 if rnorm <= atol else maxiter

        if callback_type == 'x':
            callback(x)

        if rnorm <= atol:
            break
        elif breakdown:
            # Reached breakdown (= exact solution), but the external
            # tolerance check failed. Bail out with failure.
            break
        elif presid <= ptol:
            # Inner loop passed but outer didn't
            ptol_max_factor = max(eps, 0.25 * ptol_max_factor)
        else:
            ptol_max_factor = min(1.0, 1.5 * ptol_max_factor)

        ptol = presid * min(ptol_max_factor, atol / rnorm)

    info = 0 if (rnorm <= atol) else maxiter
    return x, info

################################################################################
# Method Path: scipy.sparse.linalg._funm_multiply_krylov._funm_multiply_krylov_arnoldi
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_funm_multiply_krylov.py
# Target:      dot
# Call Count:  139074
# Total Exec:  633
# Func Size:   50
################################################################################

def _funm_multiply_krylov_arnoldi(A, b, bnorm, V, H, m):
    """
    The Arnoldi iteration for constructing the basis V and the projection H = V * A V
    for the Krylov subspace Km(A, b) of order m.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose matrix function is of interest.
    b : ndarray
        The vector b to multiply the f(A) with.
    V : ndarray
        The n x (m + 1) matrix whose columns determines the basis for
        Krylov subspace Km(A, b).
    H : ndarray
        A (m + 1) x m upper Hessenberg matrix representing the projection of A
        onto Km(A, b).
    m : int
        The order of the Krylov subspace.

    Returns
    -------
    breakdown : bool
        Indicate if the Arnoldi broke down or not

    iter : int
        Returns the last valid iteration.

    """

    dotprod = np.vdot if np.iscomplexobj(b) else np.dot
    norm_tol = np.finfo(b.dtype.char).eps ** 2
    V[:, 0] = b / bnorm

    for k in range(0, m):
        V[:, k + 1] = A.dot(V[:, k])

        # Uses the modified Gram-Schmift process to orthogonalize V[:, k + 1]
        # against the previous basis vectors
        for i in range(0, k + 1):
            H[i, k] = dotprod(V[:, i], V[:, k + 1])
            V[:, k + 1] = V[:, k + 1] - H[i, k] * V[:, i]

        H[k + 1, k] = norm(V[:, k + 1])
        if H[k + 1, k] < norm_tol:
            return True, k

        V[:, k + 1] = V[:, k + 1] / H[k + 1, k]

    return False, m

################################################################################
# Method Path: scipy.sparse.linalg._expm_multiply._expm_multiply_interval_core_1
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_expm_multiply.py
# Target:      dot
# Call Count:  115852
# Total Exec:  349
# Func Size:   27
################################################################################

def _expm_multiply_interval_core_1(A, X, h, mu, m_star, s, q, tol):
    """
    A helper function, for the case q > s and q % s == 0.
    """
    d = q // s
    input_shape = X.shape[1:]
    K_shape = (m_star + 1, ) + input_shape
    K = np.empty(K_shape, dtype=X.dtype)
    for i in range(s):
        Z = X[i*d]
        K[0] = Z
        high_p = 0
        for k in range(1, d+1):
            F = K[0]
            c1 = _exact_inf_norm(F)
            for p in range(1, m_star+1):
                if p > high_p:
                    K[p] = h * A.dot(K[p-1]) / float(p)
                coeff = float(pow(k, p))
                F = F + coeff * K[p]
                inf_norm_K_p_1 = _exact_inf_norm(K[p])
                c2 = coeff * inf_norm_K_p_1
                if c1 + c2 <= tol * _exact_inf_norm(F):
                    break
                c1 = c2
            X[k + i*d] = np.exp(k*h*mu) * F
    return X, 1

################################################################################
# Method Path: scipy.sparse.linalg._onenormest.vectors_are_parallel
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_onenormest.py
# Target:      dot
# Call Count:  109657
# Total Exec:  109657
# Func Size:   8
################################################################################

def vectors_are_parallel(v, w):
    # Columns are considered parallel when they are equal or negative.
    # Entries are required to be in {-1, 1},
    # which guarantees that the magnitudes of the vectors are identical.
    if v.ndim != 1 or v.shape != w.shape:
        raise ValueError('expected conformant vectors with entries in {-1,1}')
    n = v.shape[0]
    return np.dot(v, w) == n

################################################################################
# Method Path: scipy.optimize._trustregion_constr.minimize_trustregion_constr.matvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/minimize_trustregion_constr.py
# Target:      dot
# Call Count:  100146
# Total Exec:  50073
# Func Size:   2
################################################################################

def matvec(p):
            return self.hessp(x, p, *args)

################################################################################
# Method Path: scipy.sparse.linalg._matfuncs._onenorm_matrix_power_nnm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_matfuncs.py
# Target:      dot
# Call Count:  99934
# Total Exec:  4146
# Func Size:   31
################################################################################

def _onenorm_matrix_power_nnm(A, p):
    """
    Compute the 1-norm of a non-negative integer power of a non-negative matrix.

    Parameters
    ----------
    A : a square ndarray or matrix or sparse arrays
        Input matrix with non-negative entries.
    p : non-negative integer
        The power to which the matrix is to be raised.

    Returns
    -------
    out : float
        The 1-norm of the matrix power p of A.

    """
    # Check input
    if int(p) != p or p < 0:
        raise ValueError('expected non-negative integer p')
    p = int(p)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected A to be like a square matrix')

    # Explicitly make a column vector so that this works when A is a
    # numpy matrix (in addition to ndarray and sparse arrays).
    v = np.ones((A.shape[0], 1), dtype=float)
    M = A.T
    for i in range(p):
        v = M.dot(v)
    return np.max(v)

################################################################################
# Method Path: scipy._external.pyprima.common.linalg.inprod
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/subprojects/pyprima/pyprima/pyprima/src/pyprima/common/linalg.py
# Target:      dot
# Call Count:  83572
# Total Exec:  83572
# Func Size:   7
################################################################################

def inprod(x, y):
    if not USE_NAIVE_MATH:
        return np.dot(x, y)
    result = 0
    for i in range(len(x)):
        result += x[i] * y[i]
    return result

################################################################################
# Method Path: scipy.optimize._differentialevolution.fun
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_differentialevolution.py
# Target:      dot
# Call Count:  68947
# Total Exec:  68947
# Func Size:   18
################################################################################

def fun(x):
                x = np.asarray(x)
                return np.atleast_1d(constraint.fun(x))

################################################################################
# Method Path: scipy.optimize._trustregion_constr.qp_subproblem.projected_cg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/qp_subproblem.py
# Target:      dot
# Call Count:  65876
# Total Exec:  6938
# Func Size:   227
################################################################################

def projected_cg(H, c, Z, Y, b, trust_radius=np.inf,
                 lb=None, ub=None, tol=None,
                 max_iter=None, max_infeasible_iter=None,
                 return_all=False):
    """Solve EQP problem with projected CG method.

    Solve equality-constrained quadratic programming problem
    ``min 1/2 x.T H x + x.t c``  subject to ``A x + b = 0`` and,
    possibly, to trust region constraints ``||x|| < trust_radius``
    and box constraints ``lb <= x <= ub``.

    Parameters
    ----------
    H : LinearOperator (or sparse array or ndarray), shape (n, n)
        Operator for computing ``H v``.
    c : array_like, shape (n,)
        Gradient of the quadratic objective function.
    Z : LinearOperator (or sparse array or ndarray), shape (n, n)
        Operator for projecting ``x`` into the null space of A.
    Y : LinearOperator,  sparse array, ndarray, shape (n, m)
        Operator that, for a given a vector ``b``, compute smallest
        norm solution of ``A x + b = 0``.
    b : array_like, shape (m,)
        Right-hand side of the constraint equation.
    trust_radius : float, optional
        Trust radius to be considered. By default, uses ``trust_radius=inf``,
        which means no trust radius at all.
    lb : array_like, shape (n,), optional
        Lower bounds to each one of the components of ``x``.
        If ``lb[i] = -Inf`` the lower bound for the i-th
        component is just ignored (default).
    ub : array_like, shape (n, ), optional
        Upper bounds to each one of the components of ``x``.
        If ``ub[i] = Inf`` the upper bound for the i-th
        component is just ignored (default).
    tol : float, optional
        Tolerance used to interrupt the algorithm.
    max_iter : int, optional
        Maximum algorithm iterations. Where ``max_inter <= n-m``.
        By default, uses ``max_iter = n-m``.
    max_infeasible_iter : int, optional
        Maximum infeasible (regarding box constraints) iterations the
        algorithm is allowed to take.
        By default, uses ``max_infeasible_iter = n-m``.
    return_all : bool, optional
        When ``true``, return the list of all vectors through the iterations.

    Returns
    -------
    x : array_like, shape (n,)
        Solution of the EQP problem.
    info : dict
        Dictionary containing the following:

            - niter : Number of iterations.
            - stop_cond : Reason for algorithm termination:
                1. Iteration limit was reached;
                2. Reached the trust-region boundary;
                3. Negative curvature detected;
                4. Tolerance was satisfied.
            - allvecs : List containing all intermediary vectors (optional).
            - hits_boundary : True if the proposed step is on the boundary
              of the trust region.

    Notes
    -----
    Implementation of Algorithm 6.2 on [1]_.

    In the absence of spherical and box constraints, for sufficient
    iterations, the method returns a truly optimal result.
    In the presence of those constraints, the value returned is only
    an inexpensive approximation of the optimal value.

    References
    ----------
    .. [1] Gould, Nicholas IM, Mary E. Hribar, and Jorge Nocedal.
           "On the solution of equality constrained quadratic
            programming problems arising in optimization."
            SIAM Journal on Scientific Computing 23.4 (2001): 1376-1395.
    """
    CLOSE_TO_ZERO = 1e-25

    n, = np.shape(c)  # Number of parameters
    m, = np.shape(b)  # Number of constraints

    # Initial Values
    x = Y.dot(-b)
    r = Z.dot(H.dot(x) + c)
    g = Z.dot(r)
    p = -g

    # Store ``x`` value
    if return_all:
        allvecs = [x]
    # Values for the first iteration
    H_p = H.dot(p)
    rt_g = norm(g)**2  # g.T g = r.T Z g = r.T g (ref [1]_ p.1389)

    # If x > trust-region the problem does not have a solution.
    tr_distance = trust_radius - norm(x)
    if tr_distance < 0:
        raise ValueError("Trust region problem does not have a solution.")
    # If x == trust_radius, then x is the solution
    # to the optimization problem, since x is the
    # minimum norm solution to Ax=b.
    elif tr_distance < CLOSE_TO_ZERO:
        info = {'niter': 0, 'stop_cond': 2, 'hits_boundary': True}
        if return_all:
            allvecs.append(x)
            info['allvecs'] = allvecs
        return x, info

    # Set default tolerance
    if tol is None:
        tol = max(min(0.01 * np.sqrt(rt_g), 0.1 * rt_g), CLOSE_TO_ZERO)
    # Set default lower and upper bounds
    if lb is None:
        lb = np.full(n, -np.inf)
    if ub is None:
        ub = np.full(n, np.inf)
    # Set maximum iterations
    if max_iter is None:
        max_iter = n-m
    max_iter = min(max_iter, n-m)
    # Set maximum infeasible iterations
    if max_infeasible_iter is None:
        max_infeasible_iter = n-m

    hits_boundary = False
    stop_cond = 1
    counter = 0
    last_feasible_x = np.zeros_like(x)
    k = 0
    for i in range(max_iter):
        # Stop criteria - Tolerance : r.T g < tol
        if rt_g < tol:
            stop_cond = 4
            break
        k += 1
        # Compute curvature
        pt_H_p = H_p.dot(p)
        # Stop criteria - Negative curvature
        if pt_H_p <= 0:
            if np.isinf(trust_radius):
                raise ValueError("Negative curvature not allowed "
                                 "for unrestricted problems.")
            else:
                # Find intersection with constraints
                _, alpha, intersect = box_sphere_intersections(
                    x, p, lb, ub, trust_radius, entire_line=True)
                # Update solution
                if intersect:
                    x = x + alpha*p
                # Reinforce variables are inside box constraints.
                # This is only necessary because of roundoff errors.
                x = reinforce_box_boundaries(x, lb, ub)
                # Attribute information
                stop_cond = 3
                hits_boundary = True
                break

        # Get next step
        alpha = rt_g / pt_H_p
        x_next = x + alpha*p

        # Stop criteria - Hits boundary
        if np.linalg.norm(x_next) >= trust_radius:
            # Find intersection with box constraints
            _, theta, intersect = box_sphere_intersections(x, alpha*p, lb, ub,
                                                           trust_radius)
            # Update solution
            if intersect:
                x = x + theta*alpha*p
            # Reinforce variables are inside box constraints.
            # This is only necessary because of roundoff errors.
            x = reinforce_box_boundaries(x, lb, ub)
            # Attribute information
            stop_cond = 2
            hits_boundary = True
            break

        # Check if ``x`` is inside the box and start counter if it is not.
        if inside_box_boundaries(x_next, lb, ub):
            counter = 0
        else:
            counter += 1
        # Whenever outside box constraints keep looking for intersections.
        if counter > 0:
            _, theta, intersect = box_sphere_intersections(x, alpha*p, lb, ub,
                                                           trust_radius)
            if intersect:
                last_feasible_x = x + theta*alpha*p
                # Reinforce variables are inside box constraints.
                # This is only necessary because of roundoff errors.
                last_feasible_x = reinforce_box_boundaries(last_feasible_x,
                                                           lb, ub)
                counter = 0
        # Stop after too many infeasible (regarding box constraints) iteration.
        if counter > max_infeasible_iter:
            break
        # Store ``x_next`` value
        if return_all:
            allvecs.append(x_next)

        # Update residual
        r_next = r + alpha*H_p
        # Project residual g+ = Z r+
        g_next = Z.dot(r_next)
        # Compute conjugate direction step d
        rt_g_next = norm(g_next)**2  # g.T g = r.T g (ref [1]_ p.1389)
        beta = rt_g_next / rt_g
        p = - g_next + beta*p
        # Prepare for next iteration
        x = x_next
        g = g_next
        r = g_next
        rt_g = norm(g)**2  # g.T g = r.T Z g = r.T g (ref [1]_ p.1389)
        H_p = H.dot(p)

    if not inside_box_boundaries(x, lb, ub):
        x = last_feasible_x
        hits_boundary = True
    info = {'niter': k, 'stop_cond': stop_cond,
            'hits_boundary': hits_boundary}
    if return_all:
        info['allvecs'] = allvecs
    return x, info

################################################################################
# Method Path: scipy.linalg._expm_frechet.expm_frechet_algo_64
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_expm_frechet.py
# Target:      dot
# Call Count:  63192
# Total Exec:  3303
# Func Size:   54
################################################################################

def expm_frechet_algo_64(A, E):
    n = A.shape[0]
    s = None
    ident = np.identity(n)
    A_norm_1 = scipy.linalg.norm(A, 1)
    m_pade_pairs = (
            (3, _diff_pade3),
            (5, _diff_pade5),
            (7, _diff_pade7),
            (9, _diff_pade9))
    for m, pade in m_pade_pairs:
        if A_norm_1 <= ell_table_61[m]:
            U, V, Lu, Lv = pade(A, E, ident)
            s = 0
            break
    if s is None:
        # scaling
        s = max(0, int(np.ceil(np.log2(A_norm_1 / ell_table_61[13]))))
        A = A * 2.0**-s
        E = E * 2.0**-s
        # pade order 13
        A2 = np.dot(A, A)
        M2 = np.dot(A, E) + np.dot(E, A)
        A4 = np.dot(A2, A2)
        M4 = np.dot(A2, M2) + np.dot(M2, A2)
        A6 = np.dot(A2, A4)
        M6 = np.dot(A4, M2) + np.dot(M4, A2)
        b = (64764752532480000., 32382376266240000., 7771770303897600.,
                1187353796428800., 129060195264000., 10559470521600.,
                670442572800., 33522128640., 1323241920., 40840800., 960960.,
                16380., 182., 1.)
        W1 = b[13]*A6 + b[11]*A4 + b[9]*A2
        W2 = b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident
        Z1 = b[12]*A6 + b[10]*A4 + b[8]*A2
        Z2 = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
        W = np.dot(A6, W1) + W2
        U = np.dot(A, W)
        V = np.dot(A6, Z1) + Z2
        Lw1 = b[13]*M6 + b[11]*M4 + b[9]*M2
        Lw2 = b[7]*M6 + b[5]*M4 + b[3]*M2
        Lz1 = b[12]*M6 + b[10]*M4 + b[8]*M2
        Lz2 = b[6]*M6 + b[4]*M4 + b[2]*M2
        Lw = np.dot(A6, Lw1) + np.dot(M6, W1) + Lw2
        Lu = np.dot(A, Lw) + np.dot(E, W)
        Lv = np.dot(A6, Lz1) + np.dot(M6, Z1) + Lz2
    # factor once and solve twice
    lu_piv = scipy.linalg.lu_factor(-U + V)
    R = scipy.linalg.lu_solve(lu_piv, U + V)
    L = scipy.linalg.lu_solve(lu_piv, Lu + Lv + np.dot((Lu - Lv), R))
    # squaring
    for k in range(s):
        L = np.dot(R, L) + np.dot(L, R)
        R = np.dot(R, R)
    return R, L

################################################################################
# Method Path: scipy.optimize._trustregion_constr.equality_constrained_sqp.equality_constrained_sqp
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/equality_constrained_sqp.py
# Target:      dot
# Call Count:  58240
# Total Exec:  1692
# Func Size:   215
################################################################################

def equality_constrained_sqp(fun_and_constr, grad_and_jac, lagr_hess,
                             x0, fun0, grad0, constr0,
                             jac0, stop_criteria,
                             state,
                             initial_penalty,
                             initial_trust_radius,
                             factorization_method,
                             trust_lb=None,
                             trust_ub=None,
                             scaling=default_scaling):
    """Solve nonlinear equality-constrained problem using trust-region SQP.

    Solve optimization problem:

        minimize fun(x)
        subject to: constr(x) = 0

    using Byrd-Omojokun Trust-Region SQP method described in [1]_. Several
    implementation details are based on [2]_ and [3]_, p. 549.

    References
    ----------
    .. [1] Lalee, Marucha, Jorge Nocedal, and Todd Plantenga. "On the
           implementation of an algorithm for large-scale equality
           constrained optimization." SIAM Journal on
           Optimization 8.3 (1998): 682-706.
    .. [2] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal.
           "An interior point algorithm for large-scale nonlinear
           programming." SIAM Journal on Optimization 9.4 (1999): 877-900.
    .. [3] Nocedal, Jorge, and Stephen J. Wright. "Numerical optimization"
           Second Edition (2006).
    """
    PENALTY_FACTOR = 0.3  # Rho from formula (3.51), reference [2]_, p.891.
    LARGE_REDUCTION_RATIO = 0.9
    INTERMEDIARY_REDUCTION_RATIO = 0.3
    SUFFICIENT_REDUCTION_RATIO = 1e-8  # Eta from reference [2]_, p.892.
    TRUST_ENLARGEMENT_FACTOR_L = 7.0
    TRUST_ENLARGEMENT_FACTOR_S = 2.0
    MAX_TRUST_REDUCTION = 0.5
    MIN_TRUST_REDUCTION = 0.1
    SOC_THRESHOLD = 0.1
    TR_FACTOR = 0.8  # Zeta from formula (3.21), reference [2]_, p.885.
    BOX_FACTOR = 0.5

    n, = np.shape(x0)  # Number of parameters

    # Set default lower and upper bounds.
    if trust_lb is None:
        trust_lb = np.full(n, -np.inf)
    if trust_ub is None:
        trust_ub = np.full(n, np.inf)

    # Initial values
    x = np.copy(x0)
    trust_radius = initial_trust_radius
    penalty = initial_penalty
    # Compute Values
    f = fun0
    c = grad0
    b = constr0
    A = jac0
    S = scaling(x)
    # Get projections
    try:
        Z, LS, Y = projections(A, factorization_method)
    except ValueError as e:
        if str(e) == "expected square matrix":
            # can be the case if there are more equality
            # constraints than independent variables
            raise ValueError(
                "The 'expected square matrix' error can occur if there are"
                " more equality constraints than independent variables."
                " Consider how your constraints are set up, or use"
                " factorization_method='SVDFactorization'."
            ) from e
        else:
            raise e

    # Compute least-square lagrange multipliers
    v = -LS.dot(c)
    # Compute Hessian
    H = lagr_hess(x, v)

    # Update state parameters
    optimality = norm(c + A.T.dot(v), np.inf)
    constr_violation = norm(b, np.inf) if len(b) > 0 else 0
    cg_info = {'niter': 0, 'stop_cond': 0,
               'hits_boundary': False}

    last_iteration_failed = False
    while not stop_criteria(state, x, last_iteration_failed,
                            optimality, constr_violation,
                            trust_radius, penalty, cg_info):
        # Normal Step - `dn`
        # minimize 1/2*||A dn + b||^2
        # subject to:
        # ||dn|| <= TR_FACTOR * trust_radius
        # BOX_FACTOR * lb <= dn <= BOX_FACTOR * ub.
        dn = modified_dogleg(A, Y, b,
                             TR_FACTOR*trust_radius,
                             BOX_FACTOR*trust_lb,
                             BOX_FACTOR*trust_ub)

        # Tangential Step - `dt`
        # Solve the QP problem:
        # minimize 1/2 dt.T H dt + dt.T (H dn + c)
        # subject to:
        # A dt = 0
        # ||dt|| <= sqrt(trust_radius**2 - ||dn||**2)
        # lb - dn <= dt <= ub - dn
        c_t = H.dot(dn) + c
        b_t = np.zeros_like(b)
        trust_radius_t = np.sqrt(trust_radius**2 - np.linalg.norm(dn)**2)
        lb_t = trust_lb - dn
        ub_t = trust_ub - dn
        dt, cg_info = projected_cg(H, c_t, Z, Y, b_t,
                                   trust_radius_t,
                                   lb_t, ub_t)

        # Compute update (normal + tangential steps).
        d = dn + dt

        # Compute second order model: 1/2 d H d + c.T d + f.
        quadratic_model = 1/2*(H.dot(d)).dot(d) + c.T.dot(d)
        # Compute linearized constraint: l = A d + b.
        linearized_constr = A.dot(d)+b
        # Compute new penalty parameter according to formula (3.52),
        # reference [2]_, p.891.
        vpred = norm(b) - norm(linearized_constr)
        # Guarantee `vpred` always positive,
        # regardless of roundoff errors.
        vpred = max(1e-16, vpred)
        previous_penalty = penalty
        if quadratic_model > 0:
            new_penalty = quadratic_model / ((1-PENALTY_FACTOR)*vpred)
            penalty = max(penalty, new_penalty)
        # Compute predicted reduction according to formula (3.52),
        # reference [2]_, p.891.
        predicted_reduction = -quadratic_model + penalty*vpred

        # Compute merit function at current point
        merit_function = f + penalty*norm(b)
        # Evaluate function and constraints at trial point
        x_next = x + S.dot(d)
        f_next, b_next = fun_and_constr(x_next)
        # Compute merit function at trial point
        merit_function_next = f_next + penalty*norm(b_next)
        # Compute actual reduction according to formula (3.54),
        # reference [2]_, p.892.
        actual_reduction = merit_function - merit_function_next
        # Compute reduction ratio
        reduction_ratio = actual_reduction / predicted_reduction

        # Second order correction (SOC), reference [2]_, p.892.
        if reduction_ratio < SUFFICIENT_REDUCTION_RATIO and \
           norm(dn) <= SOC_THRESHOLD * norm(dt):
            # Compute second order correction
            y = -Y.dot(b_next)
            # Make sure increment is inside box constraints
            _, t, intersect = box_intersections(d, y, trust_lb, trust_ub)
            # Compute tentative point
            x_soc = x + S.dot(d + t*y)
            f_soc, b_soc = fun_and_constr(x_soc)
            # Recompute actual reduction
            merit_function_soc = f_soc + penalty*norm(b_soc)
            actual_reduction_soc = merit_function - merit_function_soc
            # Recompute reduction ratio
            reduction_ratio_soc = actual_reduction_soc / predicted_reduction
            if intersect and reduction_ratio_soc >= SUFFICIENT_REDUCTION_RATIO:
                x_next = x_soc
                f_next = f_soc
                b_next = b_soc
                reduction_ratio = reduction_ratio_soc

        # Readjust trust region step, formula (3.55), reference [2]_, p.892.
        if reduction_ratio >= LARGE_REDUCTION_RATIO:
            trust_radius = max(TRUST_ENLARGEMENT_FACTOR_L * norm(d),
                               trust_radius)
        elif reduction_ratio >= INTERMEDIARY_REDUCTION_RATIO:
            trust_radius = max(TRUST_ENLARGEMENT_FACTOR_S * norm(d),
                               trust_radius)
        # Reduce trust region step, according to reference [3]_, p.696.
        elif reduction_ratio < SUFFICIENT_REDUCTION_RATIO:
            trust_reduction = ((1-SUFFICIENT_REDUCTION_RATIO) /
                               (1-reduction_ratio))
            new_trust_radius = trust_reduction * norm(d)
            if new_trust_radius >= MAX_TRUST_REDUCTION * trust_radius:
                trust_radius *= MAX_TRUST_REDUCTION
            elif new_trust_radius >= MIN_TRUST_REDUCTION * trust_radius:
                trust_radius = new_trust_radius
            else:
                trust_radius *= MIN_TRUST_REDUCTION

        # Update iteration
        if reduction_ratio >= SUFFICIENT_REDUCTION_RATIO:
            x = x_next
            f, b = f_next, b_next
            c, A = grad_and_jac(x)
            S = scaling(x)
            # Get projections
            Z, LS, Y = projections(A, factorization_method)
            # Compute least-square lagrange multipliers
            v = -LS.dot(c)
            # Compute Hessian
            H = lagr_hess(x, v)
            # Set Flag
            last_iteration_failed = False
            # Optimality values
            optimality = norm(c + A.T.dot(v), np.inf)
            constr_violation = norm(b, np.inf) if len(b) > 0 else 0
        else:
            penalty = previous_penalty
            last_iteration_failed = True

    return x, state

################################################################################
# Method Path: scipy.sparse.linalg._expm_multiply._expm_multiply_simple_core
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_expm_multiply.py
# Target:      dot
# Call Count:  56185
# Total Exec:  1249
# Func Size:   24
################################################################################

def _expm_multiply_simple_core(A, B, t, mu, m_star, s, tol=None, balance=False):
    """
    A helper function.
    """
    if balance:
        raise NotImplementedError
    if tol is None:
        u_d = 2 ** -53
        tol = u_d
    F = B
    eta = np.exp(t*mu / float(s))
    for i in range(s):
        c1 = _exact_inf_norm(B)
        for j in range(m_star):
            coeff = t / float(s*(j+1))
            B = coeff * A.dot(B)
            c2 = _exact_inf_norm(B)
            F = F + B
            if c1 + c2 <= tol * _exact_inf_norm(F):
                break
            c1 = c2
        F = eta * F
        B = F
    return F

################################################################################
# Method Path: scipy.linalg._matfuncs_inv_ssq._MatrixM1PowerOperator._matmat
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs_inv_ssq.py
# Target:      dot
# Call Count:  49217
# Total Exec:  14970
# Func Size:   4
################################################################################

def _matmat(self, X):
        for i in range(self._p):
            X = self._A.dot(X) - X
        return X

################################################################################
# Method Path: scipy.optimize._trustregion_constr.projections.null_space
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/projections.py
# Target:      dot
# Call Count:  46141
# Total Exec:  31041
# Func Size:   29
################################################################################

def null_space(x):
        v = factor(A.dot(x))
        z = x - A.T.dot(v)

        # Iterative refinement to improve roundoff
        # errors described in [2]_, algorithm 5.1.
        k = 0
        while orthogonality(A, z) > orth_tol:
            if k >= max_refin:
                break
            # z_next = z - A.T inv(A A.T) A z
            v = factor(A.dot(z))
            z = z - A.T.dot(v)
            k += 1

        return z

################################################################################
# Method Path: scipy.sparse.linalg._eigen.lobpcg.lobpcg._b_orthonormalize
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_eigen/lobpcg/lobpcg.py
# Target:      inv
# Call Count:  40257
# Total Exec:  40258
# Func Size:   47
################################################################################

def _b_orthonormalize(B, blockVectorV, blockVectorBV=None,
                      verbosityLevel=0):
    """in-place B-orthonormalize the given block vector using Cholesky."""
    if blockVectorBV is None:
        if B is None:
            blockVectorBV = blockVectorV
        else:
            try:
                blockVectorBV = B(blockVectorV)
            except Exception as e:
                if verbosityLevel:
                    warnings.warn(
                        f"Secondary MatMul call failed with error\n"
                        f"{e}\n",
                        UserWarning, stacklevel=3
                    )
                    return None, None, None
            if blockVectorBV.shape != blockVectorV.shape:
                raise ValueError(
                    f"The shape {blockVectorV.shape} "
                    f"of the orthogonalized matrix not preserved\n"
                    f"and changed to {blockVectorBV.shape} "
                    f"after multiplying by the secondary matrix.\n"
                )

    VBV = blockVectorV.T.conj() @ blockVectorBV
    try:
        # VBV is a Cholesky factor from now on...
        VBV = cholesky(VBV, overwrite_a=True)
        VBV = inv(VBV, overwrite_a=True)
        blockVectorV = _matmul_inplace(
            blockVectorV, VBV,
            verbosityLevel=verbosityLevel
        )
        if B is not None:
            blockVectorBV = _matmul_inplace(
                blockVectorBV, VBV,
                verbosityLevel=verbosityLevel
            )
        return blockVectorV, blockVectorBV, VBV
    except LinAlgError:
        if verbosityLevel:
            warnings.warn(
                "Cholesky has failed.",
                UserWarning, stacklevel=3
            )
        return None, None, None

################################################################################
# Method Path: scipy.optimize._trustregion_constr.canonical_constraint.matvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/canonical_constraint.py
# Target:      dot
# Call Count:  37752
# Total Exec:  17504
# Func Size:   5
################################################################################

def matvec(p):
                result = np.zeros_like(p, dtype=float)
                for h in hess_all:
                    result += h.dot(p)
                return result

################################################################################
# Method Path: scipy.optimize._linesearch.derphi
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linesearch.py
# Target:      dot
# Call Count:  31130
# Total Exec:  31130
# Func Size:   5
################################################################################

def derphi(s):
        gval[0] = fprime(xk + s*pk, *args)
        gc[0] += 1
        return np.dot(gval[0], pk)

################################################################################
# Method Path: scipy.optimize._trustregion_constr.tr_interior_point.BarrierSubproblem.matvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/tr_interior_point.py
# Target:      dot
# Call Count:  30605
# Total Exec:  30605
# Func Size:   7
################################################################################

def matvec(vec):
            return diag_elements*vec

################################################################################
# Method Path: scipy.sparse.linalg._matfuncs._smart_matrix_product
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_matfuncs.py
# Target:      dot
# Call Count:  27449
# Total Exec:  27601
# Func Size:   41
################################################################################

def _smart_matrix_product(A, B, alpha=None, structure=None):
    """
    A matrix product that knows about sparse and structured matrices.

    Parameters
    ----------
    A : 2d ndarray
        First matrix.
    B : 2d ndarray
        Second matrix.
    alpha : float
        The matrix product will be scaled by this constant.
    structure : str, optional
        A string describing the structure of both matrices `A` and `B`.
        Only `upper_triangular` is currently supported.

    Returns
    -------
    M : 2d ndarray
        Matrix product of A and B.

    """
    if len(A.shape) != 2:
        raise ValueError('expected A to be a rectangular matrix')
    if len(B.shape) != 2:
        raise ValueError('expected B to be a rectangular matrix')
    f = None
    if structure == UPPER_TRIANGULAR:
        if (not issparse(A) and not issparse(B)
                and not is_pydata_spmatrix(A) and not is_pydata_spmatrix(B)):
            f, = scipy.linalg.get_blas_funcs(('trmm',), (A, B))
    if f is not None:
        if alpha is None:
            alpha = 1.
        out = f(alpha, A, B)
    else:
        if alpha is None:
            out = A.dot(B)
        else:
            out = alpha * A.dot(B)
    return out

################################################################################
# Method Path: scipy.integrate._ivp.rk.rk_step
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  25803
# Total Exec:  4280
# Func Size:   58
################################################################################

def rk_step(fun, t, y, f, h, A, B, C, K):
    """Perform a single Runge-Kutta step.

    This function computes a prediction of an explicit Runge-Kutta method and
    also estimates the error of a less accurate method.

    Notation for Butcher tableau is as in [1]_.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t : float
        Current time.
    y : ndarray, shape (n,)
        Current state.
    f : ndarray, shape (n,)
        Current value of the derivative, i.e., ``fun(x, y)``.
    h : float
        Step to use.
    A : ndarray, shape (n_stages, n_stages)
        Coefficients for combining previous RK stages to compute the next
        stage. For explicit methods the coefficients at and above the main
        diagonal are zeros.
    B : ndarray, shape (n_stages,)
        Coefficients for combining RK stages for computing the final
        prediction.
    C : ndarray, shape (n_stages,)
        Coefficients for incrementing time for consecutive RK stages.
        The value for the first stage is always zero.
    K : ndarray, shape (n_stages + 1, n)
        Storage array for putting RK stages here. Stages are stored in rows.
        The last row is a linear combination of the previous rows with
        coefficients

    Returns
    -------
    y_new : ndarray, shape (n,)
        Solution at t + h computed with a higher accuracy.
    f_new : ndarray, shape (n,)
        Derivative ``fun(t + h, y_new)``.

    References
    ----------
    .. [1] E. Hairer, S. P. Norsett G. Wanner, "Solving Ordinary Differential
           Equations I: Nonstiff Problems", Sec. II.4.
    """
    K[0] = f
    for s, (a, c) in enumerate(zip(A[1:], C[1:]), start=1):
        dy = np.dot(K[:s].T, a[:s]) * h
        K[s] = fun(t + c * h, y + dy)

    y_new = y + h * np.dot(K[:-1].T, B)
    f_new = fun(t + h, y_new)

    K[-1] = f_new

    return y_new, f_new

################################################################################
# Method Path: scipy.optimize._lsq.common.evaluate_quadratic
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/common.py
# Target:      dot
# Call Count:  25099
# Total Exec:  7586
# Func Size:   37
################################################################################

def evaluate_quadratic(J, g, s, diag=None):
    """Compute values of a quadratic function arising in least squares.

    The function is 0.5 * s.T * (J.T * J + diag) * s + g.T * s.

    Parameters
    ----------
    J : ndarray, sparse array or LinearOperator, shape (m, n)
        Jacobian matrix, affects the quadratic term.
    g : ndarray, shape (n,)
        Gradient, defines the linear term.
    s : ndarray, shape (k, n) or (n,)
        Array containing steps as rows.
    diag : ndarray, shape (n,), optional
        Addition diagonal part, affects the quadratic term.
        If None, assumed to be 0.

    Returns
    -------
    values : ndarray with shape (k,) or float
        Values of the function. If `s` was 2-D, then ndarray is
        returned, otherwise, float is returned.
    """
    if s.ndim == 1:
        Js = J.dot(s)
        q = np.dot(Js, Js)
        if diag is not None:
            q += np.dot(s * diag, s)
    else:
        Js = J.dot(s.T)
        q = np.sum(Js**2, axis=0)
        if diag is not None:
            q += np.sum(diag * s**2, axis=1)

    l = np.dot(s, g)

    return 0.5 * q + l

################################################################################
# Method Path: scipy.optimize._trustregion_constr.projections.orthogonality
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/projections.py
# Target:      dot
# Call Count:  24087
# Total Exec:  35878
# Func Size:   33
################################################################################

def orthogonality(A, g):
    """Measure orthogonality between a vector and the null space of a matrix.

    Compute a measure of orthogonality between the null space
    of the (possibly sparse) matrix ``A`` and a given vector ``g``.

    The formula is a simplified (and cheaper) version of formula (3.13)
    from [1]_.
    ``orth =  norm(A g, ord=2)/(norm(A, ord='fro')*norm(g, ord=2))``.

    References
    ----------
    .. [1] Gould, Nicholas IM, Mary E. Hribar, and Jorge Nocedal.
           "On the solution of equality constrained quadratic
            programming problems arising in optimization."
            SIAM Journal on Scientific Computing 23.4 (2001): 1376-1395.
    """
    # Compute vector norms
    norm_g = np.linalg.norm(g)
    # Compute Froebnius norm of the matrix A
    if issparse(A):
        norm_A = scipy.sparse.linalg.norm(A, ord='fro')
    else:
        norm_A = np.linalg.norm(A, ord='fro')

    # Check if norms are zero
    if norm_g == 0 or norm_A == 0:
        return 0

    norm_A_g = np.linalg.norm(A.dot(g))
    # Orthogonality measure
    orth = norm_A_g / (norm_A*norm_g)
    return orth

################################################################################
# Method Path: scipy.sparse.linalg._funm_multiply_krylov._funm_multiply_krylov_lanczos
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_funm_multiply_krylov.py
# Target:      dot
# Call Count:  24000
# Total Exec:  600
# Func Size:   52
################################################################################

def _funm_multiply_krylov_lanczos(A, b, bnorm, V, H, m):
    """
    The Lanczos iteration for constructing the basis V and the projection H = V * A V
    for the Krylov subspace Km(A, b) of order m. A must be Hermitian.

    Parameters
    ----------
    A : transposable linear operator
        The operator whose matrix function is of interest.
    b : ndarray
        The vector b to multiply the f(A) with.
    V : ndarray
        The n x (m + 1) matrix whose columns determines the basis for
        Krylov subspace Km(A, b).
    H : ndarray
        A (m + 1) x m upper Hessenberg matrix representing the projection of A
        onto Km(A, b).
    m : int
        The order of the Krylov subspace.

    Returns
    -------
    breakdown : bool
        Indicate if the Arnoldi broke down or not

    iter : int
        Returns the last valid iteration.

    """
    dotprod = np.vdot if np.iscomplexobj(b) else np.dot
    norm_tol = np.finfo(b.dtype.char).eps ** 2
    V[:, 0] = b / bnorm

    for k in range(0, m):
        if k > 0:
            V[:, k + 1] = A.dot(V[:, k]) - H[k, k - 1] * V[:, k - 1]
        else:
            V[:, k + 1] = A.dot(V[:, k])

        H[k, k] = dotprod(V[:, k + 1], V[:, k])
        V[:, k + 1] = V[:, k + 1] - H[k, k] * V[:, k]

        H[k + 1, k] = norm(V[:, k + 1])

        if H[k + 1, k] < norm_tol:
            return True, k

        V[:, k + 1] = V[:, k + 1] / H[k + 1, k]
        if k < m - 1:
            H[k, k + 1] = H[k + 1, k]

    return False, m

################################################################################
# Method Path: scipy.optimize._lsq.common.build_quadratic_1d
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/common.py
# Target:      dot
# Call Count:  22596
# Total Exec:  5137
# Func Size:   49
################################################################################

def build_quadratic_1d(J, g, s, diag=None, s0=None):
    """Parameterize a multivariate quadratic function along a line.

    The resulting univariate quadratic function is given as follows::

        f(t) = 0.5 * (s0 + s*t).T * (J.T*J + diag) * (s0 + s*t) +
               g.T * (s0 + s*t)

    Parameters
    ----------
    J : ndarray, sparse array or LinearOperator shape (m, n)
        Jacobian matrix, affects the quadratic term.
    g : ndarray, shape (n,)
        Gradient, defines the linear term.
    s : ndarray, shape (n,)
        Direction vector of a line.
    diag : None or ndarray with shape (n,), optional
        Addition diagonal part, affects the quadratic term.
        If None, assumed to be 0.
    s0 : None or ndarray with shape (n,), optional
        Initial point. If None, assumed to be 0.

    Returns
    -------
    a : float
        Coefficient for t**2.
    b : float
        Coefficient for t.
    c : float
        Free term. Returned only if `s0` is provided.
    """
    v = J.dot(s)
    a = np.dot(v, v)
    if diag is not None:
        a += np.dot(s * diag, s)
    a *= 0.5

    b = np.dot(g, s)

    if s0 is not None:
        u = J.dot(s0)
        b += np.dot(u, v)
        c = 0.5 * np.dot(u, u) + np.dot(g, s0)
        if diag is not None:
            b += np.dot(s0 * diag, s)
            c += 0.5 * np.dot(s0 * diag, s0)
        return a, b, c
    else:
        return a, b

################################################################################
# Method Path: scipy.sparse.linalg._interface._CustomLinearOperator.__mul__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_interface.py
# Target:      dot
# Call Count:  20902
# Total Exec:  20902
# Func Size:   6
################################################################################

def __mul__(self, x):
        """Multiplication.

        Used by the ``*`` operator. Equivalent to `dot`.
        """
        return self.dot(x)

################################################################################
# Method Path: scipy.optimize._optimize._minimize_bfgs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  15253
# Total Exec:  609
# Func Size:   182
################################################################################

def _minimize_bfgs(fun, x0, args=(), jac=None, callback=None,
                   gtol=1e-5, norm=np.inf, eps=_epsilon, maxiter=None,
                   disp=False, return_all=False, finite_diff_rel_step=None,
                   xrtol=0, c1=1e-4, c2=0.9,
                   hess_inv0=None, workers=None, **unknown_options):
    """
    Minimization of scalar function of one or more variables using the
    BFGS algorithm.

    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    gtol : float
        Terminate successfully if gradient norm is less than `gtol`.
    norm : float
        Order of norm (Inf is max, -Inf is min).
    eps : float or ndarray
        If `jac is None` the absolute step size used for numerical
        approximation of the jacobian via forward differences.
    return_all : bool, optional
        Set to True to return a list of the best solution at each of the
        iterations.
    finite_diff_rel_step : None or array_like, optional
        If ``jac in ['2-point', '3-point', 'cs']`` the relative step size to
        use for numerical approximation of the jacobian. The absolute step
        size is computed as ``h = rel_step * sign(x) * max(1, abs(x))``,
        possibly adjusted to fit into the bounds. For ``jac='3-point'``
        the sign of `h` is ignored. If None (default) then step is selected
        automatically.
    xrtol : float, default: 0
        Relative tolerance for `x`. Terminate successfully if step size is
        less than ``xk * xrtol`` where ``xk`` is the current parameter vector.
    c1 : float, default: 1e-4
        Parameter for Armijo condition rule.
    c2 : float, default: 0.9
        Parameter for curvature condition rule.
    hess_inv0 : None or ndarray, optional
        Initial inverse hessian estimate, shape (n, n). If None (default) then
        the identity matrix is used.
    workers : int, map-like callable, optional
        A map-like callable, such as `multiprocessing.Pool.map` for evaluating
        any numerical differentiation in parallel.
        This evaluation is carried out as ``workers(fun, iterable)``.

        .. versionadded:: 1.16.0

    Notes
    -----
    Parameters `c1` and `c2` must satisfy ``0 < c1 < c2 < 1``.

    If minimization doesn't complete successfully, with an error message of
    ``Desired error not necessarily achieved due to precision loss``, then
    consider setting `gtol` to a higher value. This precision loss typically
    occurs when the (finite difference) numerical differentiation cannot provide
    sufficient precision to satisfy the `gtol` termination criterion.
    This can happen when working in single precision and a callable jac is not
    provided. For single precision problems a `gtol` of 1e-3 seems to work.
    """
    _check_unknown_options(unknown_options)
    _check_positive_definite(hess_inv0)
    retall = return_all

    x0 = asarray(x0).flatten()
    if x0.ndim == 0:
        x0 = x0.reshape((1,))
    if maxiter is None:
        maxiter = len(x0) * 200

    sf = _prepare_scalar_function(fun, x0, jac, args=args, epsilon=eps,
                                  finite_diff_rel_step=finite_diff_rel_step,
                                  workers=workers)

    f = sf.fun
    myfprime = sf.grad

    old_fval = f(x0)
    gfk = myfprime(x0)

    k = 0
    N = len(x0)
    I = np.eye(N, dtype=int)
    Hk = I if hess_inv0 is None else hess_inv0

    # Sets the initial step guess to dx ~ 1
    old_old_fval = old_fval + np.linalg.norm(gfk) / 2

    xk = x0
    if retall:
        allvecs = [x0]
    warnflag = 0
    gnorm = vecnorm(gfk, ord=norm)
    while (gnorm > gtol) and (k < maxiter):
        pk = -np.dot(Hk, gfk)
        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     _line_search_wolfe12(f, myfprime, xk, pk, gfk,
                                          old_fval, old_old_fval, amin=1e-100,
                                          amax=1e100, c1=c1, c2=c2)
        except _LineSearchError:
            # Line search failed to find a better solution.
            warnflag = 2
            break

        sk = alpha_k * pk
        xkp1 = xk + sk

        if retall:
            allvecs.append(xkp1)
        xk = xkp1
        if gfkp1 is None:
            gfkp1 = myfprime(xkp1)

        yk = gfkp1 - gfk
        gfk = gfkp1
        k += 1
        intermediate_result = OptimizeResult(x=xk, fun=old_fval)
        if _call_callback_maybe_halt(callback, intermediate_result):
            break
        gnorm = vecnorm(gfk, ord=norm)
        if (gnorm <= gtol):
            break

        #  See Chapter 5 in  P.E. Frandsen, K. Jonasson, H.B. Nielsen,
        #  O. Tingleff: "Unconstrained Optimization", IMM, DTU.  1999.
        #  These notes are available here:
        #  http://www2.imm.dtu.dk/documents/ftp/publlec.html
        if (alpha_k*vecnorm(pk) <= xrtol*(xrtol + vecnorm(xk))):
            break

        if not np.isfinite(old_fval):
            # We correctly found +-Inf as optimal value, or something went
            # wrong.
            warnflag = 2
            break

        rhok_inv = np.dot(yk, sk)
        # this was handled in numeric, let it remains for more safety
        # Cryptic comment above is preserved for posterity. Future reader:
        # consider change to condition below proposed in gh-1261/gh-17345.
        if rhok_inv == 0.:
            rhok = 1000.0
            if disp:
                msg = "Divide-by-zero encountered: rhok assumed large"
                _print_success_message_or_warn(True, msg)
        else:
            rhok = 1. / rhok_inv

        A1 = I - sk[:, np.newaxis] * yk[np.newaxis, :] * rhok
        A2 = I - yk[:, np.newaxis] * sk[np.newaxis, :] * rhok
        Hk = np.dot(A1, np.dot(Hk, A2)) + (rhok * sk[:, np.newaxis] *
                                                 sk[np.newaxis, :])

    fval = old_fval

    if warnflag == 2:
        msg = _status_message['pr_loss']
    elif k >= maxiter:
        warnflag = 1
        msg = _status_message['maxiter']
    elif np.isnan(gnorm) or np.isnan(fval) or np.isnan(xk).any():
        warnflag = 3
        msg = _status_message['nan']
    else:
        msg = _status_message['success']

    if disp:
        _print_success_message_or_warn(warnflag, msg)
        print(f"         Current function value: {fval:f}")
        print(f"         Iterations: {k:d}")
        print(f"         Function evaluations: {sf.nfev:d}")
        print(f"         Gradient evaluations: {sf.ngev:d}")

    result = OptimizeResult(fun=fval, jac=gfk, hess_inv=Hk, nfev=sf.nfev,
                            njev=sf.ngev, status=warnflag,
                            success=(warnflag == 0), message=msg, x=xk,
                            nit=k)
    if retall:
        result['allvecs'] = allvecs
    return result

################################################################################
# Method Path: scipy.optimize._hessian_update_strategy.BFGS._update_implementation
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_hessian_update_strategy.py
# Target:      dot
# Call Count:  14883
# Total Exec:  7397
# Func Size:   44
################################################################################

def _update_implementation(self, delta_x, delta_grad):
        # Auxiliary variables w and z
        if self.approx_type == 'hess':
            w = delta_x
            z = delta_grad
        else:
            w = delta_grad
            z = delta_x
        # Do some common operations
        wz = np.dot(w, z)
        Mw = self @ w
        wMw = Mw.dot(w)
        # Guarantee that wMw > 0 by reinitializing matrix.
        # While this is always true in exact arithmetic,
        # indefinite matrix may appear due to roundoff errors.
        if wMw <= 0.0:
            scale = self._auto_scale(delta_x, delta_grad)
            # Reinitialize matrix
            if self.approx_type == 'hess':
                self.B = scale * np.eye(self.n, dtype=float)
            else:
                self.H = scale * np.eye(self.n, dtype=float)
            # Do common operations for new matrix
            Mw = self @ w
            wMw = Mw.dot(w)
        # Check if curvature condition is violated
        if wz <= self.min_curvature * wMw:
            # If the option 'skip_update' is set
            # we just skip the update when the condition
            # is violated.
            if self.exception_strategy == 'skip_update':
                return
            # If the option 'damp_update' is set we
            # interpolate between the actual BFGS
            # result and the unmodified matrix.
            elif self.exception_strategy == 'damp_update':
                update_factor = (1-self.min_curvature) / (1 - wz/wMw)
                z = update_factor*z + (1-update_factor)*Mw
                wz = np.dot(w, z)
        # Update matrix
        if self.approx_type == 'hess':
            self._update_hessian(wz, Mw, wMw, z)
        else:
            self._update_inverse_hessian(wz, Mw, wMw, z)

################################################################################
# Method Path: scipy.sparse.linalg._expm_multiply._expm_multiply_interval_core_2
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_expm_multiply.py
# Target:      dot
# Call Count:  14306
# Total Exec:  85
# Func Size:   34
################################################################################

def _expm_multiply_interval_core_2(A, X, h, mu, m_star, s, q, tol):
    """
    A helper function, for the case q > s and q % s > 0.
    """
    d = q // s
    j = q // d
    r = q - d * j
    input_shape = X.shape[1:]
    K_shape = (m_star + 1, ) + input_shape
    K = np.empty(K_shape, dtype=X.dtype)
    for i in range(j + 1):
        Z = X[i*d]
        K[0] = Z
        high_p = 0
        if i < j:
            effective_d = d
        else:
            effective_d = r
        for k in range(1, effective_d+1):
            F = K[0]
            c1 = _exact_inf_norm(F)
            for p in range(1, m_star+1):
                if p == high_p + 1:
                    K[p] = h * A.dot(K[p-1]) / float(p)
                    high_p = p
                coeff = float(pow(k, p))
                F = F + coeff * K[p]
                inf_norm_K_p_1 = _exact_inf_norm(K[p])
                c2 = coeff * inf_norm_K_p_1
                if c1 + c2 <= tol * _exact_inf_norm(F):
                    break
                c1 = c2
            X[k + i*d] = np.exp(k*h*mu) * F
    return X, 2

################################################################################
# Method Path: scipy.optimize._trustregion_constr.projections.row_space
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/projections.py
# Target:      dot
# Call Count:  14239
# Total Exec:  13577
# Func Size:   8
################################################################################

def row_space(x):
        return A.T.dot(factor(x))

################################################################################
# Method Path: scipy.optimize._linprog_ip._get_delta
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  13933
# Total Exec:  959
# Func Size:   199
################################################################################

def _get_delta(A, b, c, x, y, z, tau, kappa, gamma, eta, sparse=False,
               lstsq=False, sym_pos=True, cholesky=True, pc=True, ip=False,
               permc_spec='MMD_AT_PLUS_A'):
    """
    Given standard form problem defined by ``A``, ``b``, and ``c``;
    current variable estimates ``x``, ``y``, ``z``, ``tau``, and ``kappa``;
    algorithmic parameters ``gamma and ``eta;
    and options ``sparse``, ``lstsq``, ``sym_pos``, ``cholesky``, ``pc``
    (predictor-corrector), and ``ip`` (initial point improvement),
    get the search direction for increments to the variable estimates.

    Parameters
    ----------
    As defined in [4], except:
    sparse : bool
        True if the system to be solved is sparse. This is typically set
        True when the original ``A_ub`` and ``A_eq`` arrays are sparse.
    lstsq : bool
        True if the system is ill-conditioned and/or (nearly) singular and
        thus a more robust least-squares solver is desired. This is sometimes
        needed as the solution is approached.
    sym_pos : bool
        True if the system matrix is symmetric positive definite
        Sometimes this needs to be set false as the solution is approached,
        even when the system should be symmetric positive definite, due to
        numerical difficulties.
    cholesky : bool
        True if the system is to be solved by Cholesky, rather than LU,
        decomposition. This is typically faster unless the problem is very
        small or prone to numerical difficulties.
    pc : bool
        True if the predictor-corrector method of Mehrota is to be used. This
        is almost always (if not always) beneficial. Even though it requires
        the solution of an additional linear system, the factorization
        is typically (implicitly) reused so solution is efficient, and the
        number of algorithm iterations is typically reduced.
    ip : bool
        True if the improved initial point suggestion due to [4] section 4.3
        is desired. It's unclear whether this is beneficial.
    permc_spec : str (default = 'MMD_AT_PLUS_A')
        (Has effect only with ``sparse = True``, ``lstsq = False``, ``sym_pos =
        True``.) A matrix is factorized in each iteration of the algorithm.
        This option specifies how to permute the columns of the matrix for
        sparsity preservation. Acceptable values are:

        - ``NATURAL``: natural ordering.
        - ``MMD_ATA``: minimum degree ordering on the structure of A^T A.
        - ``MMD_AT_PLUS_A``: minimum degree ordering on the structure of A^T+A.
        - ``COLAMD``: approximate minimum degree column ordering.

        This option can impact the convergence of the
        interior point algorithm; test different values to determine which
        performs best for your problem. For more information, refer to
        ``scipy.sparse.linalg.splu``.

    Returns
    -------
    Search directions as defined in [4]

    References
    ----------
    .. [4] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    if A.shape[0] == 0:
        # If there are no constraints, some solvers fail (understandably)
        # rather than returning empty solution. This gets the job done.
        lstsq, sym_pos, cholesky = False, True, False
    n_x = len(x)

    # [4] Equation 8.8
    r_P = b * tau - A.dot(x)
    r_D = c * tau - A.T.dot(y) - z
    r_G = c.dot(x) - b.transpose().dot(y) + kappa
    mu = (x.dot(z) + tau * kappa) / (n_x + 1)

    #  Assemble M from [4] Equation 8.31
    Dinv = x / z

    if sparse:
        M = A.dot(sps.diags_array(Dinv, format="csc").dot(A.T))
    else:
        M = A.dot(Dinv.reshape(-1, 1) * A.T)
    solve = _get_solver(M, sparse, lstsq, sym_pos, cholesky, permc_spec)

    # pc: "predictor-corrector" [4] Section 4.1
    # In development this option could be turned off
    # but it always seems to improve performance substantially
    n_corrections = 1 if pc else 0

    i = 0
    alpha, d_x, d_z, d_tau, d_kappa = 0, 0, 0, 0, 0
    while i <= n_corrections:
        # Reference [4] Eq. 8.6
        rhatp = eta(gamma) * r_P
        rhatd = eta(gamma) * r_D
        rhatg = eta(gamma) * r_G

        # Reference [4] Eq. 8.7
        rhatxs = gamma * mu - x * z
        rhattk = gamma * mu - tau * kappa

        if i == 1:
            if ip:  # if the correction is to get "initial point"
                # Reference [4] Eq. 8.23
                rhatxs = ((1 - alpha) * gamma * mu -
                          x * z - alpha**2 * d_x * d_z)
                rhattk = ((1 - alpha) * gamma * mu -
                    tau * kappa -
                    alpha**2 * d_tau * d_kappa)
            else:  # if the correction is for "predictor-corrector"
                # Reference [4] Eq. 8.13
                rhatxs -= d_x * d_z
                rhattk -= d_tau * d_kappa

        # sometimes numerical difficulties arise as the solution is approached
        # this loop tries to solve the equations using a sequence of functions
        # for solve. For dense systems, the order is:
        # 1. scipy.linalg.cho_factor/scipy.linalg.cho_solve,
        # 2. scipy.linalg.solve w/ sym_pos = True,
        # 3. scipy.linalg.solve w/ sym_pos = False, and if all else fails
        # 4. scipy.linalg.lstsq
        # For sparse systems, the order is:
        # 1. sksparse.cholmod.cholesky (if available)
        # 2. scipy.sparse.linalg.factorized (if umfpack available)
        # 3. scipy.sparse.linalg.splu
        # 4. scipy.sparse.linalg.lsqr
        solved = False
        while not solved:
            try:
                # [4] Equation 8.28
                p, q = _sym_solve(Dinv, A, c, b, solve)
                # [4] Equation 8.29
                u, v = _sym_solve(Dinv, A, rhatd -
                                  (1 / x) * rhatxs, rhatp, solve)
                if np.any(np.isnan(p)) or np.any(np.isnan(q)):
                    raise LinAlgError
                solved = True
            except (LinAlgError, ValueError, TypeError) as e:
                # Usually this doesn't happen. If it does, it happens when
                # there are redundant constraints or when approaching the
                # solution. If so, change solver.
                if cholesky:
                    cholesky = False
                    warn(
                        "Solving system with option 'cholesky':True "
                        "failed. It is normal for this to happen "
                        "occasionally, especially as the solution is "
                        "approached. However, if you see this frequently, "
                        "consider setting option 'cholesky' to False.",
                        OptimizeWarning, stacklevel=5)
                elif sym_pos:
                    sym_pos = False
                    warn(
                        "Solving system with option 'sym_pos':True "
                        "failed. It is normal for this to happen "
                        "occasionally, especially as the solution is "
                        "approached. However, if you see this frequently, "
                        "consider setting option 'sym_pos' to False.",
                        OptimizeWarning, stacklevel=5)
                elif not lstsq:
                    lstsq = True
                    warn(
                        "Solving system with option 'sym_pos':False "
                        "failed. This may happen occasionally, "
                        "especially as the solution is "
                        "approached. However, if you see this frequently, "
                        "your problem may be numerically challenging. "
                        "If you cannot improve the formulation, consider "
                        "setting 'lstsq' to True. Consider also setting "
                        "`presolve` to True, if it is not already.",
                        OptimizeWarning, stacklevel=5)
                else:
                    raise e
                solve = _get_solver(M, sparse, lstsq, sym_pos,
                                    cholesky, permc_spec)
        # [4] Results after 8.29
        d_tau = ((rhatg + 1 / tau * rhattk - (-c.dot(u) + b.dot(v))) /
                 (1 / tau * kappa + (-c.dot(p) + b.dot(q))))
        d_x = u + p * d_tau
        d_y = v + q * d_tau

        # [4] Relations between  after 8.25 and 8.26
        d_z = (1 / x) * (rhatxs - z * d_x)
        d_kappa = 1 / tau * (rhattk - kappa * d_tau)

        # [4] 8.12 and "Let alpha be the maximal possible step..." before 8.23
        alpha = _get_step(x, d_x, z, d_z, tau, d_tau, kappa, d_kappa, 1)
        if ip:  # initial point - see [4] 4.4
            gamma = 10
        else:  # predictor-corrector, [4] definition after 8.12
            beta1 = 0.1  # [4] pg. 220 (Table 8.1)
            gamma = (1 - alpha)**2 * min(beta1, (1 - alpha))
        i += 1

    return d_x, d_y, d_z, d_tau, d_kappa

################################################################################
# Method Path: scipy.optimize._differentiable_functions.VectorFunction._update_hess
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_differentiable_functions.py
# Target:      dot
# Call Count:  12718
# Total Exec:  8612
# Func Size:   18
################################################################################

def _update_hess(self):
        if not self.H_updated:
            if callable(self._orig_hess):
                self.H = self.hess_wrapped(xp_copy(self.x), self.v)
                self._nhev += 1
            elif self._orig_hess in FD_METHODS:
                self._update_jac()
                self.H = self.hess_wrapped(xp_copy(self.x), self.v, J0=self.J)
            elif isinstance(self._orig_hess, HessianUpdateStrategy):
                self._update_jac()
                # When v is updated before x was updated, then x_prev and
                # J_prev are None and we need this check.
                if self.x_prev is not None and self.J_prev is not None:
                    delta_x = self.x - self.x_prev
                    delta_g = self.J.T.dot(self.v) - self.J_prev.T.dot(self.v)
                    self.H.update(delta_x, delta_g)

            self.H_updated = True

################################################################################
# Method Path: scipy.integrate._ivp.radau.solve_collocation_system
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/radau.py
# Target:      dot
# Call Count:  12515
# Total Exec:  1719
# Func Size:   89
################################################################################

def solve_collocation_system(fun, t, y, h, Z0, scale, tol,
                             LU_real, LU_complex, solve_lu):
    """Solve the collocation system.

    Parameters
    ----------
    fun : callable
        Right-hand side of the system.
    t : float
        Current time.
    y : ndarray, shape (n,)
        Current state.
    h : float
        Step to try.
    Z0 : ndarray, shape (3, n)
        Initial guess for the solution. It determines new values of `y` at
        ``t + h * C`` as ``y + Z0``, where ``C`` is the Radau method constants.
    scale : ndarray, shape (n)
        Problem tolerance scale, i.e. ``rtol * abs(y) + atol``.
    tol : float
        Tolerance to which solve the system. This value is compared with
        the normalized by `scale` error.
    LU_real, LU_complex
        LU decompositions of the system Jacobians.
    solve_lu : callable
        Callable which solves a linear system given a LU decomposition. The
        signature is ``solve_lu(LU, b)``.

    Returns
    -------
    converged : bool
        Whether iterations converged.
    n_iter : int
        Number of completed iterations.
    Z : ndarray, shape (3, n)
        Found solution.
    rate : float
        The rate of convergence.
    """
    n = y.shape[0]
    M_real = MU_REAL / h
    M_complex = MU_COMPLEX / h

    W = TI.dot(Z0)
    Z = Z0

    F = np.empty((3, n))
    ch = h * C

    dW_norm_old = None
    dW = np.empty_like(W)
    converged = False
    rate = None
    for k in range(NEWTON_MAXITER):
        for i in range(3):
            F[i] = fun(t + ch[i], y + Z[i])

        if not np.all(np.isfinite(F)):
            break

        f_real = F.T.dot(TI_REAL) - M_real * W[0]
        f_complex = F.T.dot(TI_COMPLEX) - M_complex * (W[1] + 1j * W[2])

        dW_real = solve_lu(LU_real, f_real)
        dW_complex = solve_lu(LU_complex, f_complex)

        dW[0] = dW_real
        dW[1] = dW_complex.real
        dW[2] = dW_complex.imag

        dW_norm = norm(dW / scale)
        if dW_norm_old is not None:
            rate = dW_norm / dW_norm_old

        if (rate is not None and (rate >= 1 or
                rate ** (NEWTON_MAXITER - k) / (1 - rate) * dW_norm > tol)):
            break

        W += dW
        Z = T.dot(W)

        if (dW_norm == 0 or
                rate is not None and rate / (1 - rate) * dW_norm < tol):
            converged = True
            break

        dW_norm_old = dW_norm

    return converged, k + 1, Z, rate

################################################################################
# Method Path: scipy.sparse.linalg._isolve.iterative.bicgstab
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/iterative.py
# Target:      dot
# Call Count:  11046
# Total Exec:  451
# Func Size:   148
################################################################################

def bicgstab(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=None, M=None,
             callback=None):
    """
    Solve ``Ax = b`` with the BIConjugate Gradient STABilized method.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` and ``A^T x`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    maxiter : int
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse array, ndarray, LinearOperator}
        Preconditioner for `A`. It should approximate the
        inverse of `A` (see Notes). Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as ``callback(xk)``, where ``xk`` is the current solution vector.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : parameter breakdown

    Notes
    -----
    The preconditioner `M` should be a matrix such that ``M @ A`` has a smaller
    condition number than `A`, see [1]_ .

    References
    ----------
    .. [1] "Preconditioner", Wikipedia, 
           https://en.wikipedia.org/wiki/Preconditioner
    .. [2] "Biconjugate gradient stabilized method", 
           Wikipedia, https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import bicgstab
    >>> R = np.array([[4, 2, 0, 1],
    ...               [3, 0, 0, 2],
    ...               [0, 1, 1, 1],
    ...               [0, 2, 1, 0]])
    >>> A = csc_array(R)
    >>> b = np.array([-1, -0.5, -1, 2])
    >>> x, exit_code = bicgstab(A, b, atol=1e-5)
    >>> print(exit_code)  # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A, M, x, b = make_system(A, M, x0, b)
    bnrm2 = np.linalg.norm(b)

    atol, _ = _get_atol_rtol('bicgstab', bnrm2, atol, rtol)

    if bnrm2 == 0:
        return b, 0

    n = len(b)

    dotprod = np.vdot if np.iscomplexobj(x) else np.dot

    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec

    # These values make no sense but coming from original Fortran code
    # sqrt might have been meant instead.
    rhotol = np.finfo(x.dtype.char).eps**2
    omegatol = rhotol

    # Dummy values to initialize vars, silence linter warnings
    rho_prev, omega, alpha, p, v = None, None, None, None, None

    r = b - matvec(x) if x.any() else b.copy()
    rtilde = r.copy()

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return x, 0

        rho = dotprod(rtilde, r)
        if np.abs(rho) < rhotol:  # rho breakdown
            return x, -10

        if iteration > 0:
            if np.abs(omega) < omegatol:  # omega breakdown
                return x, -11

            beta = (rho / rho_prev) * (alpha / omega)
            p -= omega*v
            p *= beta
            p += r
        else:  # First spin
            s = np.empty_like(r)
            p = r.copy()

        phat = psolve(p)
        v = matvec(phat)
        rv = dotprod(rtilde, v)
        if rv == 0:
            return x, -11
        alpha = rho / rv
        r -= alpha*v
        s[:] = r[:]

        if np.linalg.norm(s) < atol:
            x += alpha*phat
            return x, 0

        shat = psolve(s)
        t = matvec(shat)
        omega = dotprod(t, s) / dotprod(t, t)
        x += alpha*phat
        x += omega*shat
        r -= omega*t
        rho_prev = rho

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return x, maxiter

################################################################################
# Method Path: scipy.optimize._trustregion_constr.qp_subproblem.modified_dogleg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/qp_subproblem.py
# Target:      dot
# Call Count:  10942
# Total Exec:  6934
# Func Size:   94
################################################################################

def modified_dogleg(A, Y, b, trust_radius, lb, ub):
    """Approximately  minimize ``1/2*|| A x + b ||^2`` inside trust-region.

    Approximately solve the problem of minimizing ``1/2*|| A x + b ||^2``
    subject to ``||x|| < Delta`` and ``lb <= x <= ub`` using a modification
    of the classical dogleg approach.

    Parameters
    ----------
    A : LinearOperator (or sparse array or ndarray), shape (m, n)
        Matrix ``A`` in the minimization problem. It should have
        dimension ``(m, n)`` such that ``m < n``.
    Y : LinearOperator (or sparse array or ndarray), shape (n, m)
        LinearOperator that apply the projection matrix
        ``Q = A.T inv(A A.T)`` to the vector. The obtained vector
        ``y = Q x`` being the minimum norm solution of ``A y = x``.
    b : array_like, shape (m,)
        Vector ``b``in the minimization problem.
    trust_radius: float
        Trust radius to be considered. Delimits a sphere boundary
        to the problem.
    lb : array_like, shape (n,)
        Lower bounds to each one of the components of ``x``.
        It is expected that ``lb <= 0``, otherwise the algorithm
        may fail. If ``lb[i] = -Inf``, the lower
        bound for the ith component is just ignored.
    ub : array_like, shape (n, )
        Upper bounds to each one of the components of ``x``.
        It is expected that ``ub >= 0``, otherwise the algorithm
        may fail. If ``ub[i] = Inf``, the upper bound for the ith
        component is just ignored.

    Returns
    -------
    x : array_like, shape (n,)
        Solution to the problem.

    Notes
    -----
    Based on implementations described in pp. 885-886 from [1]_.

    References
    ----------
    .. [1] Byrd, Richard H., Mary E. Hribar, and Jorge Nocedal.
           "An interior point algorithm for large-scale nonlinear
           programming." SIAM Journal on Optimization 9.4 (1999): 877-900.
    """
    # Compute minimum norm minimizer of 1/2*|| A x + b ||^2.
    newton_point = -Y.dot(b)
    # Check for interior point
    if inside_box_boundaries(newton_point, lb, ub)  \
       and norm(newton_point) <= trust_radius:
        x = newton_point
        return x

    # Compute gradient vector ``g = A.T b``
    g = A.T.dot(b)
    # Compute Cauchy point
    # `cauchy_point = g.T g / (g.T A.T A g)``.
    A_g = A.dot(g)
    cauchy_point = -np.dot(g, g) / np.dot(A_g, A_g) * g
    # Origin
    origin_point = np.zeros_like(cauchy_point)

    # Check the segment between cauchy_point and newton_point
    # for a possible solution.
    z = cauchy_point
    p = newton_point - cauchy_point
    _, alpha, intersect = box_sphere_intersections(z, p, lb, ub,
                                                   trust_radius)
    if intersect:
        x1 = z + alpha*p
    else:
        # Check the segment between the origin and cauchy_point
        # for a possible solution.
        z = origin_point
        p = cauchy_point
        _, alpha, _ = box_sphere_intersections(z, p, lb, ub,
                                               trust_radius)
        x1 = z + alpha*p

    # Check the segment between origin and newton_point
    # for a possible solution.
    z = origin_point
    p = newton_point
    _, alpha, _ = box_sphere_intersections(z, p, lb, ub,
                                           trust_radius)
    x2 = z + alpha*p

    # Return the best solution among x1 and x2.
    if norm(A.dot(x1) + b) < norm(A.dot(x2) + b):
        return x1
    else:
        return x2

################################################################################
# Method Path: scipy.optimize._hessian_update_strategy.BFGS.__matmul__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_hessian_update_strategy.py
# Target:      dot
# Call Count:  10215
# Total Exec:  10215
# Func Size:   2
################################################################################

def __matmul__(self, p):
        return self.dot(p)

################################################################################
# Method Path: scipy.optimize._trustregion_constr.qp_subproblem.sphere_intersections
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/qp_subproblem.py
# Target:      dot
# Call Count:  9786
# Total Exec:  3268
# Func Size:   84
################################################################################

def sphere_intersections(z, d, trust_radius,
                         entire_line=False):
    """Find the intersection between segment (or line) and spherical constraints.

    Find the intersection between the segment (or line) defined by the
    parametric  equation ``x(t) = z + t*d`` and the ball
    ``||x|| <= trust_radius``.

    Parameters
    ----------
    z : array_like, shape (n,)
        Initial point.
    d : array_like, shape (n,)
        Direction.
    trust_radius : float
        Ball radius.
    entire_line : bool, optional
        When ``True``, the function returns the intersection between the line
        ``x(t) = z + t*d`` (``t`` can assume any value) and the ball
        ``||x|| <= trust_radius``. When ``False``, the function returns the intersection
        between the segment ``x(t) = z + t*d``, ``0 <= t <= 1``, and the ball.

    Returns
    -------
    ta, tb : float
        The line/segment ``x(t) = z + t*d`` is inside the ball for
        for ``ta <= t <= tb``.
    intersect : bool
        When ``True``, there is an intersection between the line/segment
        and the sphere. On the other hand, when ``False``, there is no
        intersection.
    """
    # Special case when d=0
    if norm(d) == 0:
        return 0, 0, False
    # Check for inf trust_radius
    if np.isinf(trust_radius):
        if entire_line:
            ta = -np.inf
            tb = np.inf
        else:
            ta = 0
            tb = 1
        intersect = True
        return ta, tb, intersect

    a = np.dot(d, d)
    b = 2 * np.dot(z, d)
    c = np.dot(z, z) - trust_radius**2
    discriminant = b*b - 4*a*c
    if discriminant < 0:
        intersect = False
        return 0, 0, intersect
    sqrt_discriminant = np.sqrt(discriminant)

    # The following calculation is mathematically
    # equivalent to:
    # ta = (-b - sqrt_discriminant) / (2*a)
    # tb = (-b + sqrt_discriminant) / (2*a)
    # but produce smaller round off errors.
    # Look at Matrix Computation p.97
    # for a better justification.
    aux = b + copysign(sqrt_discriminant, b)
    ta = -aux / (2*a)
    tb = -2*c / aux
    ta, tb = sorted([ta, tb])

    if entire_line:
        intersect = True
    else:
        # Checks to see if intersection happens
        # within vectors length.
        if tb < 0 or ta > 1:
            intersect = False
            ta = 0
            tb = 0
        else:
            intersect = True
            # Restrict intersection interval
            # between 0 and 1.
            ta = max(0, ta)
            tb = min(1, tb)

    return ta, tb, intersect

################################################################################
# Method Path: scipy.stats._multivariate.multivariate_t_gen._logpdf
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  9636
# Total Exec:  9639
# Func Size:   42
################################################################################

def _logpdf(self, x, loc, prec_U, log_pdet, df, dim, rank, cov_object=None):
        """Utility method `pdf`, `logpdf` for parameters.

        Parameters
        ----------
        x : ndarray
            Points at which to evaluate the log of the probability density
            function.
        loc : ndarray
            Location of the distribution.
        prec_U : ndarray
            A decomposition such that `np.dot(prec_U, prec_U.T)` is the inverse
            of the shape matrix.
        log_pdet : float
            Logarithm of the determinant of the shape matrix.
        df : float
            Degrees of freedom of the distribution.
        dim : int
            Dimension of the quantiles x.
        rank : int
            Rank of the shape matrix.

        Notes
        -----
        As this function does no argument checking, it should not be called
        directly; use 'logpdf' instead.

        """
        if df == np.inf:
            return multivariate_normal._logpdf(x, loc, cov_object)

        dev = x - loc
        maha = np.square(np.dot(dev, prec_U)).sum(axis=-1)

        t = 0.5 * (df + dim)
        A = gammaln(t)
        B = gammaln(0.5 * df)
        C = dim/2. * np.log(df * np.pi)
        D = 0.5 * log_pdet
        E = -t * np.log(1 + (1./df) * maha)

        return _squeeze_output(A - B - C - D + E)

################################################################################
# Method Path: scipy.optimize._trustregion_constr.minimize_trustregion_constr.update_state_sqp
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/minimize_trustregion_constr.py
# Target:      dot
# Call Count:  8959
# Total Exec:  8617
# Func Size:   41
################################################################################

def update_state_sqp(state, x, last_iteration_failed, objective, prepared_constraints,
                     start_time, tr_radius, constr_penalty, cg_info):
    state.nit += 1
    state.nfev = objective.nfev
    state.njev = objective.ngev
    state.nhev = objective.nhev
    state.constr_nfev = [c.fun.nfev if isinstance(c.fun, VectorFunction) else 0
                         for c in prepared_constraints]
    state.constr_njev = [c.fun.njev if isinstance(c.fun, VectorFunction) else 0
                         for c in prepared_constraints]
    state.constr_nhev = [c.fun.nhev if isinstance(c.fun, VectorFunction) else 0
                         for c in prepared_constraints]

    if not last_iteration_failed:
        state.x = x
        state.fun = objective.f
        state.grad = objective.g
        state.v = [c.fun.v for c in prepared_constraints]
        state.constr = [c.fun.f for c in prepared_constraints]
        state.jac = [c.fun.J for c in prepared_constraints]
        # Compute Lagrangian Gradient
        state.lagrangian_grad = np.copy(state.grad)
        for c in prepared_constraints:
            state.lagrangian_grad += c.fun.J.T.dot(c.fun.v)
        state.optimality = np.linalg.norm(state.lagrangian_grad, np.inf)
        # Compute maximum constraint violation
        state.constr_violation = 0
        for i in range(len(prepared_constraints)):
            lb, ub = prepared_constraints[i].bounds
            c = state.constr[i]
            state.constr_violation = np.max([state.constr_violation,
                                             np.max(lb - c),
                                             np.max(c - ub)])

    state.execution_time = time.time() - start_time
    state.tr_radius = tr_radius
    state.constr_penalty = constr_penalty
    state.cg_niter += cg_info["niter"]
    state.cg_stop_cond = cg_info["stop_cond"]

    return state

################################################################################
# Method Path: scipy.optimize._trustregion_constr.projections.least_squares
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_constr/projections.py
# Target:      dot
# Call Count:  8875
# Total Exec:  8813
# Func Size:   7
################################################################################

def least_squares(x):
        return factor(A.dot(x))

################################################################################
# Method Path: scipy.optimize._linesearch.line_search_wolfe1
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linesearch.py
# Target:      dot
# Call Count:  8655
# Total Exec:  8655
# Func Size:   69
################################################################################

def line_search_wolfe1(f, fprime, xk, pk, gfk=None,
                       old_fval=None, old_old_fval=None,
                       args=(), c1=1e-4, c2=0.9, amax=50, amin=1e-8,
                       xtol=1e-14):
    """
    As `scalar_search_wolfe1` but do a line search to direction `pk`.

    This implements the DCSRCH function of MINPACK.

    Parameters
    ----------
    f : callable
        Function `f(x)`
    fprime : callable
        Gradient of `f`
    xk : array_like
        Current point
    pk : array_like
        Search direction
    gfk : array_like, optional
        Gradient of `f` at point `xk`
    old_fval : float, optional
        Value of `f` at point `xk`
    old_old_fval : float, optional
        Value of `f` at point preceding `xk`

    The rest of the parameters are the same as for `scalar_search_wolfe1`.

    Returns
    -------
    stp, f_count, g_count, fval, old_fval
        As in `line_search_wolfe1`
    gval : array
        Gradient of `f` at the final point

    Notes
    -----
    Parameters `c1` and `c2` must satisfy ``0 < c1 < c2 < 1``.

    References
    ----------
    .. [1] Jorge J. Moré and David J. Thuente. 1994. "Line search algorithms with
       guaranteed sufficient decrease". ACM Trans. Math. Softw. 20, 3 (Sept. 1994),
       286-307. doi:`10.1145/192115.192132`

    """
    if gfk is None:
        gfk = fprime(xk, *args)

    gval = [gfk]
    gc = [0]
    fc = [0]

    def phi(s):
        fc[0] += 1
        return f(xk + s*pk, *args)

    def derphi(s):
        gval[0] = fprime(xk + s*pk, *args)
        gc[0] += 1
        return np.dot(gval[0], pk)

    derphi0 = np.dot(gfk, pk)

    stp, fval, old_fval = scalar_search_wolfe1(
            phi, derphi, old_fval, old_old_fval, derphi0,
            c1=c1, c2=c2, amax=amax, amin=amin, xtol=xtol)

    return stp, fc[0], gc[0], fval, old_fval, gval[0]

################################################################################
# Method Path: scipy.optimize._lsq.trf.trf_bounds
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/trf.py
# Target:      dot
# Call Count:  8227
# Total Exec:  226
# Func Size:   207
################################################################################

def trf_bounds(fun, jac, x0, f0, J0, lb, ub, ftol, xtol, gtol, max_nfev,
               x_scale, loss_function, tr_solver, tr_options, verbose, 
               callback=None):
    x = x0.copy()

    f = f0
    f_true = f.copy()
    nfev = 1

    J = J0
    njev = 1
    m, n = J.shape

    if loss_function is not None:
        rho = loss_function(f)
        cost = 0.5 * np.sum(rho[0])
        J, f = scale_for_robust_loss_function(J, f, rho)
    else:
        cost = 0.5 * np.dot(f, f)

    g = compute_grad(J, f)

    jac_scale = isinstance(x_scale, str) and x_scale == 'jac'
    if jac_scale:
        scale, scale_inv = compute_jac_scale(J)
    else:
        scale, scale_inv = x_scale, 1 / x_scale

    v, dv = CL_scaling_vector(x, g, lb, ub)
    v[dv != 0] *= scale_inv[dv != 0]
    Delta = norm(x0 * scale_inv / v**0.5)
    if Delta == 0:
        Delta = 1.0

    g_norm = norm(g * v, ord=np.inf)

    f_augmented = np.zeros(m + n)
    if tr_solver == 'exact':
        J_augmented = np.empty((m + n, n))
    elif tr_solver == 'lsmr':
        reg_term = 0.0
        regularize = tr_options.pop('regularize', True)

    if max_nfev is None:
        max_nfev = x0.size * 100

    alpha = 0.0  # "Levenberg-Marquardt" parameter

    termination_status = None
    iteration = 0
    step_norm = None
    actual_reduction = None

    if verbose == 2:
        print_header_nonlinear()

    while True:
        v, dv = CL_scaling_vector(x, g, lb, ub)

        g_norm = norm(g * v, ord=np.inf)
        if g_norm < gtol:
            termination_status = 1

        if verbose == 2:
            print_iteration_nonlinear(iteration, nfev, cost, actual_reduction,
                                      step_norm, g_norm)

        if termination_status is not None or nfev == max_nfev:
            break

        # Now compute variables in "hat" space. Here, we also account for
        # scaling introduced by `x_scale` parameter. This part is a bit tricky,
        # you have to write down the formulas and see how the trust-region
        # problem is formulated when the two types of scaling are applied.
        # The idea is that first we apply `x_scale` and then apply Coleman-Li
        # approach in the new variables.

        # v is recomputed in the variables after applying `x_scale`, note that
        # components which were identically 1 not affected.
        v[dv != 0] *= scale_inv[dv != 0]

        # Here, we apply two types of scaling.
        d = v**0.5 * scale

        # C = diag(g * scale) Jv
        diag_h = g * dv * scale

        # After all this has been done, we continue normally.

        # "hat" gradient.
        g_h = d * g

        f_augmented[:m] = f
        if tr_solver == 'exact':
            J_augmented[:m] = J * d
            J_h = J_augmented[:m]  # Memory view.
            J_augmented[m:] = np.diag(diag_h**0.5)
            U, s, V = svd(J_augmented, full_matrices=False)
            V = V.T
            uf = U.T.dot(f_augmented)
        elif tr_solver == 'lsmr':
            J_h = right_multiplied_operator(J, d)

            if regularize:
                a, b = build_quadratic_1d(J_h, g_h, -g_h, diag=diag_h)
                to_tr = Delta / norm(g_h)
                ag_value = minimize_quadratic_1d(a, b, 0, to_tr)[1]
                reg_term = -ag_value / Delta**2

            lsmr_op = regularized_lsq_operator(J_h, (diag_h + reg_term)**0.5)
            gn_h = lsmr(lsmr_op, f_augmented, **tr_options)[0]
            S = np.vstack((g_h, gn_h)).T
            S, _ = qr(S, mode='economic')
            JS = J_h.dot(S)  # LinearOperator does dot too.
            B_S = np.dot(JS.T, JS) + np.dot(S.T * diag_h, S)
            g_S = S.T.dot(g_h)

        # theta controls step back step ratio from the bounds.
        theta = max(0.995, 1 - g_norm)

        actual_reduction = -1
        while actual_reduction <= 0 and nfev < max_nfev:
            if tr_solver == 'exact':
                p_h, alpha, n_iter = solve_lsq_trust_region(
                    n, m, uf, s, V, Delta, initial_alpha=alpha)
            elif tr_solver == 'lsmr':
                p_S, _ = solve_trust_region_2d(B_S, g_S, Delta)
                p_h = S.dot(p_S)

            p = d * p_h  # Trust-region solution in the original space.
            step, step_h, predicted_reduction = select_step(
                x, J_h, diag_h, g_h, p, p_h, d, Delta, lb, ub, theta)

            x_new = make_strictly_feasible(x + step, lb, ub, rstep=0)
            f_new = fun(x_new)
            nfev += 1

            step_h_norm = norm(step_h)

            if not np.all(np.isfinite(f_new)):
                Delta = 0.25 * step_h_norm
                continue

            # Usual trust-region step quality estimation.
            if loss_function is not None:
                cost_new = loss_function(f_new, cost_only=True)
            else:
                cost_new = 0.5 * np.dot(f_new, f_new)
            actual_reduction = cost - cost_new
            Delta_new, ratio = update_tr_radius(
                Delta, actual_reduction, predicted_reduction,
                step_h_norm, step_h_norm > 0.95 * Delta)

            step_norm = norm(step)
            termination_status = check_termination(
                actual_reduction, cost, step_norm, norm(x), ratio, ftol, xtol)
            if termination_status is not None:
                break

            alpha *= Delta / Delta_new
            Delta = Delta_new

        if actual_reduction > 0:
            x = x_new

            f = f_new
            f_true = f.copy()

            cost = cost_new

            J = jac(x)
            njev += 1

            if loss_function is not None:
                rho = loss_function(f)
                J, f = scale_for_robust_loss_function(J, f, rho)

            g = compute_grad(J, f)

            if jac_scale:
                scale, scale_inv = compute_jac_scale(J, scale_inv)
        else:
            step_norm = 0
            actual_reduction = 0
            
        iteration += 1
            
        # Call callback function and possibly stop optimization
        if callback is not None:
            intermediate_result = OptimizeResult(
                x=x, fun=f_true, nit=iteration, nfev=nfev)
            intermediate_result["cost"] = cost

            if _call_callback_maybe_halt(
                callback, intermediate_result
            ):
                termination_status = -2
                break

    if termination_status is None:
        termination_status = 0

    active_mask = find_active_constraints(x, lb, ub, rtol=xtol)
    return OptimizeResult(
        x=x, cost=cost, fun=f_true, jac=J, grad=g, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev,
        status=termination_status)

################################################################################
# Method Path: scipy.optimize._linprog_ip._sym_solve
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  7876
# Total Exec:  4040
# Func Size:   18
################################################################################

def _sym_solve(Dinv, A, r1, r2, solve):
    """
    An implementation of [4] equation 8.31 and 8.32

    References
    ----------
    .. [4] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """
    # [4] 8.31
    r = r2 + A.dot(Dinv * r1)
    v = solve(r)
    # [4] 8.32
    u = Dinv * (A.T.dot(v) - r1)
    return u, v

################################################################################
# Method Path: scipy.optimize._optimize._minimize_newtoncg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  7758
# Total Exec:  175
# Func Size:   192
################################################################################

def _minimize_newtoncg(fun, x0, args=(), jac=None, hess=None, hessp=None,
                       callback=None, xtol=1e-5, eps=_epsilon, maxiter=None,
                       disp=False, return_all=False, c1=1e-4, c2=0.9, workers=None,
                       **unknown_options):
    """
    Minimization of scalar function of one or more variables using the
    Newton-CG algorithm.

    Note that the `jac` parameter (Jacobian) is required.

    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    xtol : float
        Average relative error in solution `xopt` acceptable for
        convergence.
    maxiter : int
        Maximum number of iterations to perform.
    eps : float or ndarray
        If `hessp` is approximated, use this value for the step size.
    return_all : bool, optional
        Set to True to return a list of the best solution at each of the
        iterations.
    c1 : float, default: 1e-4
        Parameter for Armijo condition rule.
    c2 : float, default: 0.9
        Parameter for curvature condition rule.
    workers : int, map-like callable, optional
        A map-like callable, such as `multiprocessing.Pool.map` for evaluating
        any numerical differentiation in parallel.
        This evaluation is carried out as ``workers(fun, iterable)``.

        .. versionadded:: 1.16.0

    Notes
    -----
    Parameters `c1` and `c2` must satisfy ``0 < c1 < c2 < 1``.
    """
    _check_unknown_options(unknown_options)
    if jac is None:
        raise ValueError('Jacobian is required for Newton-CG method')
    fhess_p = hessp
    fhess = hess
    avextol = xtol
    epsilon = eps
    retall = return_all

    x0 = asarray(x0).flatten()
    # TODO: add hessp (callable or FD) to ScalarFunction?
    sf = _prepare_scalar_function(
        fun, x0, jac, args=args, epsilon=eps, hess=hess, workers=workers
    )
    f = sf.fun
    fprime = sf.grad
    _h = sf.hess(x0)

    # Logic for hess/hessp
    # - If a callable(hess) is provided, then use that
    # - If hess is a FD_METHOD, or the output from hess(x) is a LinearOperator
    #   then create a hessp function using those.
    # - If hess is None but you have callable(hessp) then use the hessp.
    # - If hess and hessp are None then approximate hessp using the grad/jac.

    if (hess in FD_METHODS or isinstance(_h, LinearOperator)):
        fhess = None

        def _hessp(x, p, *args):
            return sf.hess(x).dot(p)

        fhess_p = _hessp

    def terminate(warnflag, msg):
        if disp:
            _print_success_message_or_warn(warnflag, msg)
            print(f"         Current function value: {old_fval:f}")
            print(f"         Iterations: {k:d}")
            print(f"         Function evaluations: {sf.nfev:d}")
            print(f"         Gradient evaluations: {sf.ngev:d}")
            print(f"         Hessian evaluations: {hcalls:d}")
        fval = old_fval
        result = OptimizeResult(fun=fval, jac=gfk, nfev=sf.nfev,
                                njev=sf.ngev, nhev=hcalls, status=warnflag,
                                success=(warnflag == 0), message=msg, x=xk,
                                nit=k)
        if retall:
            result['allvecs'] = allvecs
        return result

    hcalls = 0
    if maxiter is None:
        maxiter = len(x0)*200
    cg_maxiter = 20*len(x0)

    xtol = len(x0) * avextol
    # Make sure we enter the while loop.
    update_l1norm = np.finfo(float).max
    xk = np.copy(x0)
    if retall:
        allvecs = [xk]
    k = 0
    gfk = None
    old_fval = f(x0)
    old_old_fval = None
    float64eps = np.finfo(np.float64).eps
    while update_l1norm > xtol:
        if k >= maxiter:
            msg = "Warning: " + _status_message['maxiter']
            return terminate(1, msg)
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - grad f(xk) starting from 0.
        b = -fprime(xk)
        maggrad = np.linalg.norm(b, ord=1)
        eta = min(0.5, math.sqrt(maggrad))
        termcond = eta * maggrad
        xsupi = zeros(len(x0), dtype=x0.dtype)
        ri = -b
        psupi = -ri
        i = 0
        dri0 = np.dot(ri, ri)

        if fhess is not None:             # you want to compute hessian once.
            A = sf.hess(xk)
            hcalls += 1

        for k2 in range(cg_maxiter):
            if np.add.reduce(np.abs(ri)) <= termcond:
                break
            if fhess is None:
                if fhess_p is None:
                    Ap = approx_fhess_p(xk, psupi, fprime, epsilon)
                else:
                    Ap = fhess_p(xk, psupi, *args)
                    hcalls += 1
            else:
                # hess was supplied as a callable or hessian update strategy, so
                # A is a dense numpy array or sparse array
                Ap = A.dot(psupi)
            # check curvature
            Ap = asarray(Ap).squeeze()  # get rid of matrices...
            curv = np.dot(psupi, Ap)
            if 0 <= curv <= 3 * float64eps:
                break
            elif curv < 0:
                if (i > 0):
                    break
                else:
                    # fall back to steepest descent direction
                    xsupi = dri0 / (-curv) * b
                    break
            alphai = dri0 / curv
            xsupi += alphai * psupi
            ri += alphai * Ap
            dri1 = np.dot(ri, ri)
            betai = dri1 / dri0
            psupi = -ri + betai * psupi
            i += 1
            dri0 = dri1          # update np.dot(ri,ri) for next time.
        else:
            # curvature keeps increasing, bail out
            msg = ("Warning: CG iterations didn't converge. The Hessian is not "
                   "positive definite.")
            return terminate(3, msg)

        pk = xsupi  # search direction is solution to system.
        gfk = -b    # gradient at xk

        try:
            alphak, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     _line_search_wolfe12(f, fprime, xk, pk, gfk,
                                          old_fval, old_old_fval, c1=c1, c2=c2)
        except _LineSearchError:
            # Line search failed to find a better solution.
            msg = "Warning: " + _status_message['pr_loss']
            return terminate(2, msg)

        update = alphak * pk
        xk += update        # upcast if necessary
        if retall:
            allvecs.append(xk)
        k += 1
        intermediate_result = OptimizeResult(x=xk, fun=old_fval)
        if _call_callback_maybe_halt(callback, intermediate_result):
            return terminate(5, "")
        update_l1norm = np.linalg.norm(update, ord=1)

    else:
        if np.isnan(old_fval) or np.isnan(update_l1norm):
            return terminate(3, _status_message['nan'])

        msg = _status_message['success']
        return terminate(0, msg)

################################################################################
# Method Path: scipy.optimize._lsq.common.compute_grad
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/common.py
# Target:      dot
# Call Count:  6931
# Total Exec:  6942
# Func Size:   6
################################################################################

def compute_grad(J, f):
    """Compute gradient of the least-squares cost function."""
    if isinstance(J, LinearOperator):
        return J.rmatvec(f)
    else:
        return J.T.dot(f)

################################################################################
# Method Path: scipy.sparse.linalg._isolve.iterative.cgs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/iterative.py
# Target:      dot
# Call Count:  6332
# Total Exec:  426
# Func Size:   157
################################################################################

def cgs(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=None, M=None, callback=None):
    """
    Solve ``Ax = b`` with the Conjugate Gradient Squared method.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real-valued N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    maxiter : int
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse array, ndarray, LinearOperator}
        Preconditioner for ``A``. It should approximate the
        inverse of `A` (see Notes). Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as ``callback(xk)``, where ``xk`` is the current solution vector.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : parameter breakdown

    Notes
    -----
    The preconditioner `M` should be a matrix such that ``M @ A`` has a smaller
    condition number than `A`, see [1]_.

    References
    ----------
    .. [1] "Preconditioner", Wikipedia, 
           https://en.wikipedia.org/wiki/Preconditioner
    .. [2] "Conjugate gradient squared", Wikipedia,
           https://en.wikipedia.org/wiki/Conjugate_gradient_squared_method

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import cgs
    >>> R = np.array([[4, 2, 0, 1],
    ...               [3, 0, 0, 2],
    ...               [0, 1, 1, 1],
    ...               [0, 2, 1, 0]])
    >>> A = csc_array(R)
    >>> b = np.array([-1, -0.5, -1, 2])
    >>> x, exit_code = cgs(A, b)
    >>> print(exit_code)  # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A, M, x, b = make_system(A, M, x0, b)
    bnrm2 = np.linalg.norm(b)

    atol, _ = _get_atol_rtol('cgs', bnrm2, atol, rtol)

    if bnrm2 == 0:
        return b, 0

    n = len(b)

    dotprod = np.vdot if np.iscomplexobj(x) else np.dot

    if maxiter is None:
        maxiter = n*10

    matvec = A.matvec
    psolve = M.matvec

    rhotol = np.finfo(x.dtype.char).eps**2

    r = b - matvec(x) if x.any() else b.copy()

    rtilde = r.copy()
    bnorm = np.linalg.norm(b)
    if bnorm == 0:
        bnorm = 1

    # Dummy values to initialize vars, silence linter warnings
    rho_prev, p, u, q = None, None, None, None

    for iteration in range(maxiter):
        rnorm = np.linalg.norm(r)
        if rnorm < atol:  # Are we done?
            return x, 0

        rho_cur = dotprod(rtilde, r)
        if np.abs(rho_cur) < rhotol:  # Breakdown case
            return x, -10

        if iteration > 0:
            beta = rho_cur / rho_prev

            # u = r + beta * q
            # p = u + beta * (q + beta * p);
            u[:] = r[:]
            u += beta*q

            p *= beta
            p += q
            p *= beta
            p += u

        else:  # First spin
            p = r.copy()
            u = r.copy()
            q = np.empty_like(r)

        phat = psolve(p)
        vhat = matvec(phat)
        rv = dotprod(rtilde, vhat)

        if rv == 0:  # Dot product breakdown
            return x, -11

        alpha = rho_cur / rv
        q[:] = u[:]
        q -= alpha*vhat
        uhat = psolve(u + q)
        x += alpha*uhat

        # Due to numerical error build-up the actual residual is computed
        # instead of the following two lines that were in the original
        # FORTRAN templates, still using a single matvec.

        # qhat = matvec(uhat)
        # r -= alpha*qhat
        r = b - matvec(x)

        rho_prev = rho_cur

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return x, maxiter

################################################################################
# Method Path: scipy.stats._multivariate.wishart_gen._rvs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  6118
# Total Exec:  16
# Func Size:   44
################################################################################

def _rvs(self, n, shape, dim, df, C, random_state):
        """Draw random samples from a Wishart distribution.

        Parameters
        ----------
        n : int
            Number of variates to generate
        shape : iterable
            Shape of the variates to generate
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        C : ndarray
            Cholesky factorization of the scale matrix, lower triangular.
        %(_doc_random_state)s

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'rvs' instead.

        """
        random_state = self._get_random_state(random_state)
        # Calculate the matrices A, which are actually lower triangular
        # Cholesky factorizations of a matrix B such that B ~ W(df, I)
        A = self._standard_rvs(n, shape, dim, df, random_state)

        # Calculate SA = C A A' C', where SA ~ W(df, scale)
        # Note: this is the product of a (lower) (lower) (lower)' (lower)'
        #       or, denoting B = AA', it is C B C' where C is the lower
        #       triangular Cholesky factorization of the scale matrix.
        #       this appears to conflict with the instructions in [1]_, which
        #       suggest that it should be D' B D where D is the lower
        #       triangular factorization of the scale matrix. However, it is
        #       meant to refer to the Bartlett (1933) representation of a
        #       Wishart random variate as L A A' L' where L is lower triangular
        #       so it appears that understanding D' to be upper triangular
        #       is either a typo in or misreading of [1]_.
        for index in np.ndindex(shape):
            CA = np.dot(C, A[index])
            A[index] = np.dot(CA, CA.T)

        return A

################################################################################
# Method Path: scipy.optimize._optimize.descent_condition
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  6042
# Total Exec:  3535
# Func Size:   16
################################################################################

def descent_condition(alpha, xkp1, fp1, gfkp1):
            # Polak-Ribiere+ needs an explicit check of a sufficient
            # descent condition, which is not guaranteed by strong Wolfe.
            #
            # See Gilbert & Nocedal, "Global convergence properties of
            # conjugate gradient methods for optimization",
            # SIAM J. Optimization 2, 21 (1992).
            cached_step[:] = polak_ribiere_powell_step(alpha, gfkp1)
            alpha, xk, pk, gfk, gnorm = cached_step

            # Accept step if it leads to convergence.
            if gnorm <= gtol:
                return True

            # Accept step if sufficient descent condition applies.
            return np.dot(pk, gfk) <= -sigma_3 * np.dot(gfk, gfk)

################################################################################
# Method Path: scipy.sparse.linalg._isolve.iterative.qmr
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/iterative.py
# Target:      dot
# Call Count:  5738
# Total Exec:  393
# Func Size:   204
################################################################################

def qmr(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=None, M1=None, M2=None,
        callback=None):
    """
    Solve ``Ax = b`` with the Quasi-Minimal Residual method.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real-valued N-by-N matrix of the linear system.
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` and ``A^T x`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    maxiter : int
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M1 : {sparse array, ndarray, LinearOperator}
        Left preconditioner for A.
    M2 : {sparse array, ndarray, LinearOperator}
        Right preconditioner for A. Used together with the left
        preconditioner M1.  The matrix M1@A@M2 should have better
        conditioned than A alone.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as callback(xk), where xk is the current solution vector.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : parameter breakdown

    See Also
    --------
    LinearOperator

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import qmr
    >>> A = csc_array([[3., 2., 0.], [1., -1., 0.], [0., 5., 1.]])
    >>> b = np.array([2., 4., -1.])
    >>> x, exitCode = qmr(A, b, atol=1e-5)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A_ = A
    A, M, x, b = make_system(A, None, x0, b)
    bnrm2 = np.linalg.norm(b)

    atol, _ = _get_atol_rtol('qmr', bnrm2, atol, rtol)

    if bnrm2 == 0:
        return b, 0

    if M1 is None and M2 is None:
        if hasattr(A_, 'psolve'):
            def left_psolve(b):
                return A_.psolve(b, 'left')

            def right_psolve(b):
                return A_.psolve(b, 'right')

            def left_rpsolve(b):
                return A_.rpsolve(b, 'left')

            def right_rpsolve(b):
                return A_.rpsolve(b, 'right')
            M1 = LinearOperator(A.shape,
                                matvec=left_psolve,
                                rmatvec=left_rpsolve)
            M2 = LinearOperator(A.shape,
                                matvec=right_psolve,
                                rmatvec=right_rpsolve)
        else:
            def id(b):
                return b
            M1 = LinearOperator(A.shape, matvec=id, rmatvec=id)
            M2 = LinearOperator(A.shape, matvec=id, rmatvec=id)

    n = len(b)
    if maxiter is None:
        maxiter = n*10

    dotprod = np.vdot if np.iscomplexobj(x) else np.dot

    rhotol = np.finfo(x.dtype.char).eps
    betatol = rhotol
    gammatol = rhotol
    deltatol = rhotol
    epsilontol = rhotol
    xitol = rhotol

    r = b - A.matvec(x) if x.any() else b.copy()

    vtilde = r.copy()
    y = M1.matvec(vtilde)
    rho = np.linalg.norm(y)
    wtilde = r.copy()
    z = M2.rmatvec(wtilde)
    xi = np.linalg.norm(z)
    gamma, eta, theta = 1, -1, 0
    v = np.empty_like(vtilde)
    w = np.empty_like(wtilde)

    # Dummy values to initialize vars, silence linter warnings
    epsilon, q, d, p, s = None, None, None, None, None

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return x, 0
        if np.abs(rho) < rhotol:  # rho breakdown
            return x, -10
        if np.abs(xi) < xitol:  # xi breakdown
            return x, -15

        v[:] = vtilde[:]
        v *= (1 / rho)
        y *= (1 / rho)
        w[:] = wtilde[:]
        w *= (1 / xi)
        z *= (1 / xi)
        delta = dotprod(z, y)

        if np.abs(delta) < deltatol:  # delta breakdown
            return x, -13

        ytilde = M2.matvec(y)
        ztilde = M1.rmatvec(z)

        if iteration > 0:
            ytilde -= (xi * delta / epsilon) * p
            p[:] = ytilde[:]
            ztilde -= (rho * (delta / epsilon).conj()) * q
            q[:] = ztilde[:]
        else:  # First spin
            p = ytilde.copy()
            q = ztilde.copy()

        ptilde = A.matvec(p)
        epsilon = dotprod(q, ptilde)
        if np.abs(epsilon) < epsilontol:  # epsilon breakdown
            return x, -14

        beta = epsilon / delta
        if np.abs(beta) < betatol:  # beta breakdown
            return x, -11

        vtilde[:] = ptilde[:]
        vtilde -= beta*v
        y = M1.matvec(vtilde)

        rho_prev = rho
        rho = np.linalg.norm(y)
        wtilde[:] = w[:]
        wtilde *= - beta.conj()
        wtilde += A.rmatvec(q)
        z = M2.rmatvec(wtilde)
        xi = np.linalg.norm(z)
        gamma_prev = gamma
        theta_prev = theta
        theta = rho / (gamma_prev * np.abs(beta))
        gamma = 1 / np.sqrt(1 + theta**2)

        if np.abs(gamma) < gammatol:  # gamma breakdown
            return x, -12

        eta *= -(rho_prev / beta) * (gamma / gamma_prev)**2

        if iteration > 0:
            d *= (theta_prev * gamma) ** 2
            d += eta*p
            s *= (theta_prev * gamma) ** 2
            s += eta*ptilde
        else:
            d = p.copy()
            d *= eta
            s = ptilde.copy()
            s *= eta

        x += d
        r -= s

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return x, maxiter

################################################################################
# Method Path: scipy.sparse.linalg._isolve.iterative.bicg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/iterative.py
# Target:      dot
# Call Count:  5722
# Total Exec:  393
# Func Size:   133
################################################################################

def bicg(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=None, M=None, callback=None):
    """
    Solve ``Ax = b`` with the BIConjugate Gradient method.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` and ``A^T x`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    maxiter : int
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse array, ndarray, LinearOperator}
        Preconditioner for `A`. It should approximate the
        inverse of `A` (see Notes). Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as ``callback(xk)``, where ``xk`` is the current solution vector.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations
            <0 : parameter breakdown

    Notes
    -----
    The preconditioner `M` should be a matrix such that ``M @ A`` has a smaller
    condition number than `A`, see [1]_ .

    References
    ----------
    .. [1] "Preconditioner", Wikipedia, 
           https://en.wikipedia.org/wiki/Preconditioner
    .. [2] "Biconjugate gradient method", Wikipedia, 
           https://en.wikipedia.org/wiki/Biconjugate_gradient_method

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import bicg
    >>> A = csc_array([[3, 2, 0], [1, -1, 0], [0, 5, 1.]])
    >>> b = np.array([2., 4., -1.])
    >>> x, exitCode = bicg(A, b, atol=1e-5)
    >>> print(exitCode)  # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A, M, x, b = make_system(A, M, x0, b)
    bnrm2 = np.linalg.norm(b)

    atol, _ = _get_atol_rtol('bicg', bnrm2, atol, rtol)

    if bnrm2 == 0:
        return b, 0

    n = len(b)
    dotprod = np.vdot if np.iscomplexobj(x) else np.dot

    if maxiter is None:
        maxiter = n*10

    matvec, rmatvec = A.matvec, A.rmatvec
    psolve, rpsolve = M.matvec, M.rmatvec

    rhotol = np.finfo(x.dtype.char).eps**2

    # Dummy values to initialize vars, silence linter warnings
    rho_prev, p, ptilde = None, None, None

    r = b - matvec(x) if x.any() else b.copy()
    rtilde = r.copy()

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return x, 0

        z = psolve(r)
        ztilde = rpsolve(rtilde)
        # order matters in this dot product
        rho_cur = dotprod(rtilde, z)

        if np.abs(rho_cur) < rhotol:  # Breakdown case
            return x, -10

        if iteration > 0:
            beta = rho_cur / rho_prev
            p *= beta
            p += z
            ptilde *= beta.conj()
            ptilde += ztilde
        else:  # First spin
            p = z.copy()
            ptilde = ztilde.copy()

        q = matvec(p)
        qtilde = rmatvec(ptilde)
        rv = dotprod(ptilde, q)

        if rv == 0:
            return x, -11

        alpha = rho_cur / rv
        x += alpha*p
        r -= alpha*q
        rtilde -= alpha.conj()*qtilde
        rho_prev = rho_cur

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return x, maxiter

################################################################################
# Method Path: scipy.sparse.linalg._matfuncs._expm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_matfuncs.py
# Target:      dot
# Call Count:  5480
# Total Exec:  3866
# Func Size:   88
################################################################################

def _expm(A, use_exact_onenorm):
    # Core of expm, separated to allow testing exact and approximate
    # algorithms.

    # Avoid indiscriminate asarray() to allow sparse or other strange arrays.
    if isinstance(A, list | tuple | np.matrix):
        A = np.asarray(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected a square matrix')

    # gracefully handle size-0 input,
    # carefully handling sparse scenario
    if A.shape == (0, 0):
        out = np.zeros([0, 0], dtype=A.dtype)
        if issparse(A) or is_pydata_spmatrix(A):
            return A.__class__(out)
        return out

    # Trivial case
    if A.shape == (1, 1):
        out = [[np.exp(A[0, 0])]]

        # Avoid indiscriminate casting to ndarray to
        # allow for sparse or other strange arrays
        if issparse(A) or is_pydata_spmatrix(A):
            return A.__class__(out)

        return np.array(out)

    # Ensure input is of float type, to avoid integer overflows etc.
    if ((isinstance(A, np.ndarray) or issparse(A) or is_pydata_spmatrix(A))
            and not np.issubdtype(A.dtype, np.inexact)):
        A = A.astype(float)

    # Detect upper triangularity.
    structure = UPPER_TRIANGULAR if _is_upper_triangular(A) else None

    if use_exact_onenorm == "auto":
        # Hardcode a matrix order threshold for exact vs. estimated one-norms.
        use_exact_onenorm = A.shape[0] < 200

    # Track functions of A to help compute the matrix exponential.
    h = _ExpmPadeHelper(
            A, structure=structure, use_exact_onenorm=use_exact_onenorm)

    # Try Pade order 3.
    eta_1 = max(h.d4_loose, h.d6_loose)
    if eta_1 < 1.495585217958292e-002 and _ell(h.A, 3) == 0:
        U, V = h.pade3()
        return _solve_P_Q(U, V, structure=structure)

    # Try Pade order 5.
    eta_2 = max(h.d4_tight, h.d6_loose)
    if eta_2 < 2.539398330063230e-001 and _ell(h.A, 5) == 0:
        U, V = h.pade5()
        return _solve_P_Q(U, V, structure=structure)

    # Try Pade orders 7 and 9.
    eta_3 = max(h.d6_tight, h.d8_loose)
    if eta_3 < 9.504178996162932e-001 and _ell(h.A, 7) == 0:
        U, V = h.pade7()
        return _solve_P_Q(U, V, structure=structure)
    if eta_3 < 2.097847961257068e+000 and _ell(h.A, 9) == 0:
        U, V = h.pade9()
        return _solve_P_Q(U, V, structure=structure)

    # Use Pade order 13.
    eta_4 = max(h.d8_loose, h.d10_loose)
    eta_5 = min(eta_3, eta_4)
    theta_13 = 4.25

    # Choose smallest s>=0 such that 2**(-s) eta_5 <= theta_13
    if eta_5 == 0:
        # Nilpotent special case
        s = 0
    else:
        s = max(int(np.ceil(np.log2(eta_5 / theta_13))), 0)
    s = s + _ell(2**-s * h.A, 13)
    U, V = h.pade13_scaled(s)
    X = _solve_P_Q(U, V, structure=structure)
    if structure == UPPER_TRIANGULAR:
        # Invoke Code Fragment 2.1.
        X = _fragment_2_1(X, h.A, s)
    else:
        # X = r_13(A)^(2^s) by repeated squaring.
        for i in range(s):
            X = X.dot(X)
    return X

################################################################################
# Method Path: scipy.optimize._lsq.trf.trf_no_bounds
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/trf.py
# Target:      dot
# Call Count:  5396
# Total Exec:  229
# Func Size:   173
################################################################################

def trf_no_bounds(fun, jac, x0, f0, J0, ftol, xtol, gtol, max_nfev,
                  x_scale, loss_function, tr_solver, tr_options, verbose, 
                  callback=None):
    x = x0.copy()

    f = f0
    f_true = f.copy()
    nfev = 1

    J = J0
    njev = 1
    m, n = J.shape

    if loss_function is not None:
        rho = loss_function(f)
        cost = 0.5 * np.sum(rho[0])
        J, f = scale_for_robust_loss_function(J, f, rho)
    else:
        cost = 0.5 * np.dot(f, f)

    g = compute_grad(J, f)

    jac_scale = isinstance(x_scale, str) and x_scale == 'jac'
    if jac_scale:
        scale, scale_inv = compute_jac_scale(J)
    else:
        scale, scale_inv = x_scale, 1 / x_scale

    Delta = norm(x0 * scale_inv)
    if Delta == 0:
        Delta = 1.0

    if tr_solver == 'lsmr':
        reg_term = 0
        damp = tr_options.pop('damp', 0.0)
        regularize = tr_options.pop('regularize', True)

    if max_nfev is None:
        max_nfev = x0.size * 100

    alpha = 0.0  # "Levenberg-Marquardt" parameter

    termination_status = None
    iteration = 0
    step_norm = None
    actual_reduction = None

    if verbose == 2:
        print_header_nonlinear()

    while True:
        g_norm = norm(g, ord=np.inf)
        if g_norm < gtol:
            termination_status = 1

        if verbose == 2:
            print_iteration_nonlinear(iteration, nfev, cost, actual_reduction,
                                      step_norm, g_norm)

        if termination_status is not None or nfev == max_nfev:
            break

        d = scale
        g_h = d * g

        if tr_solver == 'exact':
            J_h = J * d
            U, s, V = svd(J_h, full_matrices=False)
            V = V.T
            uf = U.T.dot(f)
        elif tr_solver == 'lsmr':
            J_h = right_multiplied_operator(J, d)

            if regularize:
                a, b = build_quadratic_1d(J_h, g_h, -g_h)
                to_tr = Delta / norm(g_h)
                ag_value = minimize_quadratic_1d(a, b, 0, to_tr)[1]
                reg_term = -ag_value / Delta**2

            damp_full = (damp**2 + reg_term)**0.5
            gn_h = lsmr(J_h, f, damp=damp_full, **tr_options)[0]
            S = np.vstack((g_h, gn_h)).T
            S, _ = qr(S, mode='economic')
            JS = J_h.dot(S)
            B_S = np.dot(JS.T, JS)
            g_S = S.T.dot(g_h)

        actual_reduction = -1
        while actual_reduction <= 0 and nfev < max_nfev:
            if tr_solver == 'exact':
                step_h, alpha, n_iter = solve_lsq_trust_region(
                    n, m, uf, s, V, Delta, initial_alpha=alpha)
            elif tr_solver == 'lsmr':
                p_S, _ = solve_trust_region_2d(B_S, g_S, Delta)
                step_h = S.dot(p_S)

            predicted_reduction = -evaluate_quadratic(J_h, g_h, step_h)
            step = d * step_h
            x_new = x + step
            f_new = fun(x_new)
            nfev += 1

            step_h_norm = norm(step_h)

            if not np.all(np.isfinite(f_new)):
                Delta = 0.25 * step_h_norm
                continue

            # Usual trust-region step quality estimation.
            if loss_function is not None:
                cost_new = loss_function(f_new, cost_only=True)
            else:
                cost_new = 0.5 * np.dot(f_new, f_new)
            actual_reduction = cost - cost_new

            Delta_new, ratio = update_tr_radius(
                Delta, actual_reduction, predicted_reduction,
                step_h_norm, step_h_norm > 0.95 * Delta)

            step_norm = norm(step)
            termination_status = check_termination(
                actual_reduction, cost, step_norm, norm(x), ratio, ftol, xtol)
            if termination_status is not None:
                break

            alpha *= Delta / Delta_new
            Delta = Delta_new

        if actual_reduction > 0:
            x = x_new

            f = f_new
            f_true = f.copy()

            cost = cost_new

            J = jac(x)
            njev += 1

            if loss_function is not None:
                rho = loss_function(f)
                J, f = scale_for_robust_loss_function(J, f, rho)

            g = compute_grad(J, f)

            if jac_scale:
                scale, scale_inv = compute_jac_scale(J, scale_inv)
        else:
            step_norm = 0
            actual_reduction = 0

        iteration += 1
        
        # Call callback function and possibly stop optimization
        if callback is not None:
            intermediate_result = OptimizeResult(
                x=x, fun=f_true, nit=iteration, nfev=nfev)
            intermediate_result["cost"] = cost

            if _call_callback_maybe_halt(
                callback, intermediate_result
            ):
                termination_status = -2
                break

    if termination_status is None:
        termination_status = 0

    active_mask = np.zeros_like(x)
    return OptimizeResult(
        x=x, cost=cost, fun=f_true, jac=J, grad=g, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev,
        status=termination_status)

################################################################################
# Method Path: scipy.optimize._differentiable_functions.IdentityVectorFunction.fun
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_differentiable_functions.py
# Target:      dot
# Call Count:  5183
# Total Exec:  7766
# Func Size:   6
################################################################################

def fun(self, x):
        if not np.array_equal(x, self.x):
            self._update_x(x)
        self._update_fun()
        return self.f

################################################################################
# Method Path: scipy.sparse.linalg._isolve.iterative.cg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/iterative.py
# Target:      dot
# Call Count:  5012
# Total Exec:  331
# Func Size:   121
################################################################################

def cg(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=None, M=None, callback=None):
    """
    Solve ``Ax = b`` with the Conjugate Gradient method, for a symmetric,
    positive-definite `A`.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        `A` must represent a hermitian, positive definite matrix.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-5``.
    maxiter : int
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse array, ndarray, LinearOperator}
        Preconditioner for `A`. `M` must represent a hermitian, positive definite
        matrix. It should approximate the inverse of `A` (see Notes).
        Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as ``callback(xk)``, where ``xk`` is the current solution vector.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:
            0  : successful exit
            >0 : convergence to tolerance not achieved, number of iterations

    Notes
    -----
    The preconditioner `M` should be a matrix such that ``M @ A`` has a smaller
    condition number than `A`, see [2]_.

    References
    ----------
    .. [1] "Conjugate Gradient Method, Wikipedia, 
           https://en.wikipedia.org/wiki/Conjugate_gradient_method
    .. [2] "Preconditioner", 
           Wikipedia, https://en.wikipedia.org/wiki/Preconditioner

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import cg
    >>> P = np.array([[4, 0, 1, 0],
    ...               [0, 5, 0, 0],
    ...               [1, 0, 3, 2],
    ...               [0, 0, 2, 4]])
    >>> A = csc_array(P)
    >>> b = np.array([-1, -0.5, -1, 2])
    >>> x, exit_code = cg(A, b, atol=1e-5)
    >>> print(exit_code)    # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A, M, x, b = make_system(A, M, x0, b)
    bnrm2 = np.linalg.norm(b)

    atol, _ = _get_atol_rtol('cg', bnrm2, atol, rtol)

    if bnrm2 == 0:
        return b, 0

    n = len(b)

    if maxiter is None:
        maxiter = n*10

    dotprod = np.vdot if np.iscomplexobj(x) else np.dot

    matvec = A.matvec
    psolve = M.matvec
    r = b - matvec(x) if x.any() else b.copy()

    # Dummy value to initialize var, silences warnings
    rho_prev, p = None, None

    for iteration in range(maxiter):
        if np.linalg.norm(r) < atol:  # Are we done?
            return x, 0

        z = psolve(r)
        rho_cur = dotprod(r, z)
        if iteration > 0:
            beta = rho_cur / rho_prev
            p *= beta
            p += z
        else:  # First spin
            p = np.empty_like(r)
            p[:] = z[:]

        q = matvec(p)
        alpha = rho_cur / dotprod(p, q)
        x += alpha*p
        r -= alpha*q
        rho_prev = rho_cur

        if callback:
            callback(x)

    else:  # for loop exhausted
        # Return incomplete progress
        return x, maxiter

################################################################################
# Method Path: scipy.sparse.linalg._isolve._gcrotmk.gcrotmk
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/_gcrotmk.py
# Target:      dot
# Call Count:  4863
# Total Exec:  610
# Func Size:   321
################################################################################

def gcrotmk(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=1000, M=None, callback=None,
            m=20, k=None, CU=None, discard_C=False, truncate='oldest'):
    """
    Solve ``Ax = b`` with the flexible GCROT(m,k) algorithm.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        `LinearOperator`.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``rtol=1e-5`` and ``atol=0.0``.
    maxiter : int, optional
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved. The
        default is ``1000``.
    M : {sparse array, ndarray, LinearOperator}, optional
        Preconditioner for `A`.  The preconditioner should approximate the
        inverse of `A`. gcrotmk is a 'flexible' algorithm and the preconditioner
        can vary from iteration to iteration. Effective preconditioning
        dramatically improves the rate of convergence, which implies that
        fewer iterations are needed to reach a given error tolerance.
    callback : function, optional
        User-supplied function to call after each iteration.  It is called
        as ``callback(xk)``, where ``xk`` is the current solution vector.
    m : int, optional
        Number of inner FGMRES iterations per each outer iteration.
        Default: 20
    k : int, optional
        Number of vectors to carry between inner FGMRES iterations.
        According to [2]_, good values are around `m`.
        Default: `m`
    CU : list of tuples, optional
        List of tuples ``(c, u)`` which contain the columns of the matrices
        C and U in the GCROT(m,k) algorithm. For details, see [2]_.
        The list given and vectors contained in it are modified in-place.
        If not given, start from empty matrices. The ``c`` elements in the
        tuples can be ``None``, in which case the vectors are recomputed
        via ``c = A u`` on start and orthogonalized as described in [3]_.
    discard_C : bool, optional
        Discard the C-vectors at the end. Useful if recycling Krylov subspaces
        for different linear systems.
    truncate : {'oldest', 'smallest'}, optional
        Truncation scheme to use. Drop: oldest vectors, or vectors with
        smallest singular values using the scheme discussed in [1,2].
        See [2]_ for detailed comparison.
        Default: 'oldest'

    Returns
    -------
    x : ndarray
        The solution found.
    info : int
        Provides convergence information:

        * 0  : successful exit
        * >0 : convergence to tolerance not achieved, number of iterations

    References
    ----------
    .. [1] E. de Sturler, ''Truncation strategies for optimal Krylov subspace
           methods'', SIAM J. Numer. Anal. 36, 864 (1999).
    .. [2] J.E. Hicken and D.W. Zingg, ''A simplified and flexible variant
           of GCROT for solving nonsymmetric linear systems'',
           SIAM J. Sci. Comput. 32, 172 (2010).
    .. [3] M.L. Parks, E. de Sturler, G. Mackey, D.D. Johnson, S. Maiti,
           ''Recycling Krylov subspaces for sequences of linear systems'',
           SIAM J. Sci. Comput. 28, 1651 (2006).

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import gcrotmk
    >>> R = np.random.randn(5, 5)
    >>> A = csc_array(R)
    >>> b = np.random.randn(5)
    >>> x, exit_code = gcrotmk(A, b, atol=1e-5)
    >>> print(exit_code)
    0
    >>> np.allclose(A.dot(x), b)
    True

    """
    A,M,x,b = make_system(A,M,x0,b)

    if not np.isfinite(b).all():
        raise ValueError("RHS must contain only finite numbers")

    if truncate not in ('oldest', 'smallest'):
        raise ValueError(f"Invalid value for 'truncate': {truncate!r}")

    matvec = A.matvec
    psolve = M.matvec

    if CU is None:
        CU = []

    if k is None:
        k = m

    axpy, dot, scal = None, None, None

    if x0 is None:
        r = b.copy()
    else:
        r = b - matvec(x)

    axpy, dot, scal, nrm2 = get_blas_funcs(['axpy', 'dot', 'scal', 'nrm2'], (x, r))

    b_norm = nrm2(b)

    # we call this to get the right atol/rtol and raise errors as necessary
    atol, rtol = _get_atol_rtol('gcrotmk', b_norm, atol, rtol)

    if b_norm == 0:
        x = b
        return (x, 0)

    if discard_C:
        CU[:] = [(None, u) for c, u in CU]

    # Reorthogonalize old vectors
    if CU:
        # Sort already existing vectors to the front
        CU.sort(key=lambda cu: cu[0] is not None)

        # Fill-in missing ones
        C = np.empty((A.shape[0], len(CU)), dtype=r.dtype, order='F')
        us = []
        j = 0
        while CU:
            # More memory-efficient: throw away old vectors as we go
            c, u = CU.pop(0)
            if c is None:
                c = matvec(u)
            C[:,j] = c
            j += 1
            us.append(u)

        # Orthogonalize
        Q, R, P = qr(C, overwrite_a=True, mode='economic', pivoting=True)
        del C

        # C := Q
        cs = list(Q.T)

        # U := U P R^-1,  back-substitution
        new_us = []
        for j in range(len(cs)):
            u = us[P[j]]
            for i in range(j):
                u = axpy(us[P[i]], u, u.shape[0], -R[i,j])
            if abs(R[j,j]) < 1e-12 * abs(R[0,0]):
                # discard rest of the vectors
                break
            u = scal(1.0/R[j,j], u)
            new_us.append(u)

        # Form the new CU lists
        CU[:] = list(zip(cs, new_us))[::-1]

    if CU:
        axpy, dot = get_blas_funcs(['axpy', 'dot'], (r,))

        # Solve first the projection operation with respect to the CU
        # vectors. This corresponds to modifying the initial guess to
        # be
        #
        #     x' = x + U y
        #     y = argmin_y || b - A (x + U y) ||^2
        #
        # The solution is y = C^H (b - A x)
        for c, u in CU:
            yc = dot(c, r)
            x = axpy(u, x, x.shape[0], yc)
            r = axpy(c, r, r.shape[0], -yc)

    # GCROT main iteration
    for j_outer in range(maxiter):
        # -- callback
        if callback is not None:
            callback(x)

        beta = nrm2(r)

        # -- check stopping condition
        beta_tol = max(atol, rtol * b_norm)

        if beta <= beta_tol and (j_outer > 0 or CU):
            # recompute residual to avoid rounding error
            r = b - matvec(x)
            beta = nrm2(r)

        if beta <= beta_tol:
            j_outer = -1
            break

        ml = m + max(k - len(CU), 0)

        cs = [c for c, u in CU]

        try:
            Q, R, B, vs, zs, y, pres = _fgmres(matvec,
                                               r/beta,
                                               ml,
                                               rpsolve=psolve,
                                               atol=max(atol, rtol*b_norm)/beta,
                                               cs=cs)
            y *= beta
        except LinAlgError:
            # Floating point over/underflow, non-finite result from
            # matmul etc. -- report failure.
            break

        #
        # At this point,
        #
        #     [A U, A Z] = [C, V] G;   G =  [ I  B ]
        #                                   [ 0  H ]
        #
        # where [C, V] has orthonormal columns, and r = beta v_0. Moreover,
        #
        #     || b - A (x + Z y + U q) ||_2 = || r - C B y - V H y - C q ||_2 = min!
        #
        # from which y = argmin_y || beta e_1 - H y ||_2, and q = -B y
        #

        #
        # GCROT(m,k) update
        #

        # Define new outer vectors

        # ux := (Z - U B) y
        ux = zs[0]*y[0]
        for z, yc in zip(zs[1:], y[1:]):
            ux = axpy(z, ux, ux.shape[0], yc)  # ux += z*yc
        by = B.dot(y)
        for cu, byc in zip(CU, by):
            c, u = cu
            ux = axpy(u, ux, ux.shape[0], -byc)  # ux -= u*byc

        # cx := V H y
        with np.errstate(invalid="ignore"):
            hy = Q.dot(R.dot(y))
        cx = vs[0] * hy[0]
        for v, hyc in zip(vs[1:], hy[1:]):
            cx = axpy(v, cx, cx.shape[0], hyc)  # cx += v*hyc

        # Normalize cx, maintaining cx = A ux
        # This new cx is orthogonal to the previous C, by construction
        try:
            alpha = 1/nrm2(cx)
            if not np.isfinite(alpha):
                raise FloatingPointError()
        except (FloatingPointError, ZeroDivisionError):
            # Cannot update, so skip it
            continue

        cx = scal(alpha, cx)
        ux = scal(alpha, ux)

        # Update residual and solution
        gamma = dot(cx, r)
        r = axpy(cx, r, r.shape[0], -gamma)  # r -= gamma*cx
        x = axpy(ux, x, x.shape[0], gamma)  # x += gamma*ux

        # Truncate CU
        if truncate == 'oldest':
            while len(CU) >= k and CU:
                del CU[0]
        elif truncate == 'smallest':
            if len(CU) >= k and CU:
                # cf. [1,2]
                D = solve(R[:-1,:].T, B.T).T
                W, sigma, V = svd(D)

                # C := C W[:,:k-1],  U := U W[:,:k-1]
                new_CU = []
                for j, w in enumerate(W[:,:k-1].T):
                    c, u = CU[0]
                    c = c * w[0]
                    u = u * w[0]
                    for cup, wp in zip(CU[1:], w[1:]):
                        cp, up = cup
                        c = axpy(cp, c, c.shape[0], wp)
                        u = axpy(up, u, u.shape[0], wp)

                    # Reorthogonalize at the same time; not necessary
                    # in exact arithmetic, but floating point error
                    # tends to accumulate here
                    for cp, up in new_CU:
                        alpha = dot(cp, c)
                        c = axpy(cp, c, c.shape[0], -alpha)
                        u = axpy(up, u, u.shape[0], -alpha)
                    alpha = nrm2(c)
                    c = scal(1.0/alpha, c)
                    u = scal(1.0/alpha, u)

                    new_CU.append((c, u))
                CU[:] = new_CU

        # Add new vector to CU
        CU.append((cx, ux))

    # Include the solution vector to the span
    CU.append((None, x.copy()))
    if discard_C:
        CU[:] = [(None, uz) for cz, uz in CU]

    return x, j_outer + 1

################################################################################
# Method Path: scipy.optimize._linprog_ip.mu
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  4456
# Total Exec:  2228
# Func Size:   2
################################################################################

def mu(x, tau, z, kappa):
        return (x.dot(z) + np.dot(tau, kappa)) / (len(x) + 1)

################################################################################
# Method Path: scipy.optimize._linprog_ip.r_g
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  4456
# Total Exec:  2228
# Func Size:   2
################################################################################

def r_g(x, y, kappa):
        return kappa + c.dot(x) - b.dot(y)

################################################################################
# Method Path: scipy.optimize._linprog_ip._indicators
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  4456
# Total Exec:  1114
# Func Size:   47
################################################################################

def _indicators(A, b, c, c0, x, y, z, tau, kappa):
    """
    Implementation of several equations from [4] used as indicators of
    the status of optimization.

    References
    ----------
    .. [4] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.

    """

    # residuals for termination are relative to initial values
    x0, y0, z0, tau0, kappa0 = _get_blind_start(A.shape)

    # See [4], Section 4 - The Homogeneous Algorithm, Equation 8.8
    def r_p(x, tau):
        return b * tau - A.dot(x)

    def r_d(y, z, tau):
        return c * tau - A.T.dot(y) - z

    def r_g(x, y, kappa):
        return kappa + c.dot(x) - b.dot(y)

    # np.dot unpacks if they are arrays of size one
    def mu(x, tau, z, kappa):
        return (x.dot(z) + np.dot(tau, kappa)) / (len(x) + 1)

    obj = c.dot(x / tau) + c0

    def norm(a):
        return np.linalg.norm(a)

    # See [4], Section 4.5 - The Stopping Criteria
    r_p0 = r_p(x0, tau0)
    r_d0 = r_d(y0, z0, tau0)
    r_g0 = r_g(x0, y0, kappa0)
    mu_0 = mu(x0, tau0, z0, kappa0)
    rho_A = norm(c.T.dot(x) - b.T.dot(y)) / (tau + norm(b.T.dot(y)))
    rho_p = norm(r_p(x, tau)) / max(1, norm(r_p0))
    rho_d = norm(r_d(y, z, tau)) / max(1, norm(r_d0))
    rho_g = norm(r_g(x, y, kappa)) / max(1, norm(r_g0))
    rho_mu = mu(x, tau, z, kappa) / mu_0
    return rho_p, rho_d, rho_A, rho_g, rho_mu, obj

################################################################################
# Method Path: scipy.optimize._trustregion.DoglegSubproblem.get_boundaries_intersections
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion.py
# Target:      dot
# Call Count:  4359
# Total Exec:  1453
# Func Size:   22
################################################################################

def get_boundaries_intersections(self, z, d, trust_radius):
        """
        Solve the scalar quadratic equation ``||z + t d|| == trust_radius``.
        This is like a line-sphere intersection.
        Return the two values of t, sorted from low to high.
        """
        a = np.dot(d, d)
        b = 2 * np.dot(z, d)
        c = np.dot(z, z) - trust_radius**2
        sqrt_discriminant = math.sqrt(b*b - 4*a*c)

        # The following calculation is mathematically
        # equivalent to:
        # ta = (-b - sqrt_discriminant) / (2*a)
        # tb = (-b + sqrt_discriminant) / (2*a)
        # but produce smaller round off errors.
        # Look at Matrix Computation p.97
        # for a better justification.
        aux = b + math.copysign(sqrt_discriminant, b)
        ta = -aux / (2*a)
        tb = -2*c / aux
        return sorted([ta, tb])

################################################################################
# Method Path: scipy.stats._mannwhitneyu._MWU.build_u_freqs_array
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_mannwhitneyu.py
# Target:      dot
# Call Count:  4270
# Total Exec:  920
# Func Size:   36
################################################################################

def build_u_freqs_array(self, maxu):
        """
        Build all the array of frequencies for u from 0 to maxu.
        Assumptions:
          n1 <= n2
          maxu <= n1 * n2 / 2
        """
        n1, n2 = self.n1, self.n2
        total = special.binom(n1 + n2, n1)

        if maxu + 1 <= self.configurations.size:
            return self.configurations[:maxu + 1] / total

        s_array = self.build_sigma_array(maxu)

        # Start working with ints, for maximum precision and efficiency:
        configurations = np.zeros(maxu + 1, dtype=np.uint64)
        configurations_is_uint = True
        uint_max = np.iinfo(np.uint64).max
        # How many ways to have U=0? 1
        configurations[0] = 1

        for u in np.arange(1, maxu + 1):
            coeffs = s_array[u - 1::-1]
            new_val = np.dot(configurations[:u], coeffs) / u
            if new_val > uint_max and configurations_is_uint:
                # OK, we got into numbers too big for uint64.
                # So now we start working with floats.
                # By doing this since the beginning, we would have lost precision.
                # (And working on python long ints would be unbearably slow)
                configurations = configurations.astype(float)
                configurations_is_uint = False
            configurations[u] = new_val

        self.configurations = configurations
        return configurations / total

################################################################################
# Method Path: scipy.integrate._ivp.rk.RK45._estimate_error
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  4172
# Total Exec:  4172
# Func Size:   2
################################################################################

def _estimate_error(self, K, h):
        return np.dot(K.T, self.E) * h

################################################################################
# Method Path: scipy.optimize._lsq.common.solve_lsq_trust_region
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/common.py
# Target:      dot
# Call Count:  4011
# Total Exec:  2925
# Func Size:   112
################################################################################

def solve_lsq_trust_region(n, m, uf, s, V, Delta, initial_alpha=None,
                           rtol=0.01, max_iter=10):
    """Solve a trust-region problem arising in least-squares minimization.

    This function implements a method described by J. J. More [1]_ and used
    in MINPACK, but it relies on a single SVD of Jacobian instead of series
    of Cholesky decompositions. Before running this function, compute:
    ``U, s, VT = svd(J, full_matrices=False)``.

    Parameters
    ----------
    n : int
        Number of variables.
    m : int
        Number of residuals.
    uf : ndarray
        Computed as U.T.dot(f).
    s : ndarray
        Singular values of J.
    V : ndarray
        Transpose of VT.
    Delta : float
        Radius of a trust region.
    initial_alpha : float, optional
        Initial guess for alpha, which might be available from a previous
        iteration. If None, determined automatically.
    rtol : float, optional
        Stopping tolerance for the root-finding procedure. Namely, the
        solution ``p`` will satisfy ``abs(norm(p) - Delta) < rtol * Delta``.
    max_iter : int, optional
        Maximum allowed number of iterations for the root-finding procedure.

    Returns
    -------
    p : ndarray, shape (n,)
        Found solution of a trust-region problem.
    alpha : float
        Positive value such that (J.T*J + alpha*I)*p = -J.T*f.
        Sometimes called Levenberg-Marquardt parameter.
    n_iter : int
        Number of iterations made by root-finding procedure. Zero means
        that Gauss-Newton step was selected as the solution.

    References
    ----------
    .. [1] More, J. J., "The Levenberg-Marquardt Algorithm: Implementation
           and Theory," Numerical Analysis, ed. G. A. Watson, Lecture Notes
           in Mathematics 630, Springer Verlag, pp. 105-116, 1977.
    """
    def phi_and_derivative(alpha, suf, s, Delta):
        """Function of which to find zero.

        It is defined as "norm of regularized (by alpha) least-squares
        solution minus `Delta`". Refer to [1]_.
        """
        denom = s**2 + alpha
        p_norm = norm(suf / denom)
        phi = p_norm - Delta
        phi_prime = -np.sum(suf ** 2 / denom**3) / p_norm
        return phi, phi_prime

    suf = s * uf

    # Check if J has full rank and try Gauss-Newton step.
    if m >= n:
        threshold = EPS * m * s[0]
        full_rank = s[-1] > threshold
    else:
        full_rank = False

    if full_rank:
        p = -V.dot(uf / s)
        if norm(p) <= Delta:
            return p, 0.0, 0

    alpha_upper = norm(suf) / Delta

    if full_rank:
        phi, phi_prime = phi_and_derivative(0.0, suf, s, Delta)
        alpha_lower = -phi / phi_prime
    else:
        alpha_lower = 0.0

    if initial_alpha is None or not full_rank and initial_alpha == 0:
        alpha = max(0.001 * alpha_upper, (alpha_lower * alpha_upper)**0.5)
    else:
        alpha = initial_alpha

    for it in range(max_iter):
        if alpha < alpha_lower or alpha > alpha_upper:
            alpha = max(0.001 * alpha_upper, (alpha_lower * alpha_upper)**0.5)

        phi, phi_prime = phi_and_derivative(alpha, suf, s, Delta)

        if phi < 0:
            alpha_upper = alpha

        ratio = phi / phi_prime
        alpha_lower = max(alpha_lower, alpha - ratio)
        alpha -= (phi + Delta) * ratio / Delta

        if np.abs(phi) < rtol * Delta:
            break

    p = -V.dot(suf / (s**2 + alpha))

    # Make the norm of p equal to Delta, p is changed only slightly during
    # this. It is done to prevent p lie outside the trust region (which can
    # cause problems later).
    p *= Delta / norm(p)

    return p, alpha, it + 1

################################################################################
# Method Path: scipy.signal._ltisys._YT_real
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  3811
# Total Exec:  491
# Func Size:   63
################################################################################

def _YT_real(ker_pole, Q, transfer_matrix, i, j):
    """
    Applies algorithm from YT section 6.1 page 19 related to real pairs
    """
    # step 1 page 19
    u = Q[:, -2, np.newaxis]
    v = Q[:, -1, np.newaxis]

    # step 2 page 19
    m = np.dot(np.dot(ker_pole[i].T, np.dot(u, v.T) -
        np.dot(v, u.T)), ker_pole[j])

    # step 3 page 19
    um, sm, vm = np.linalg.svd(m)
    # mu1, mu2 two first columns of U => 2 first lines of U.T
    mu1, mu2 = um.T[:2, :, np.newaxis]
    # VM is V.T with numpy we want the first two lines of V.T
    nu1, nu2 = vm[:2, :, np.newaxis]

    # what follows is a rough python translation of the formulas
    # in section 6.2 page 20 (step 4)
    transfer_matrix_j_mo_transfer_matrix_j = np.vstack((
            transfer_matrix[:, i, np.newaxis],
            transfer_matrix[:, j, np.newaxis]))

    if not np.allclose(sm[0], sm[1]):
        ker_pole_imo_mu1 = np.dot(ker_pole[i], mu1)
        ker_pole_i_nu1 = np.dot(ker_pole[j], nu1)
        ker_pole_mu_nu = np.vstack((ker_pole_imo_mu1, ker_pole_i_nu1))
    else:
        ker_pole_ij = np.vstack((
                                np.hstack((ker_pole[i],
                                           np.zeros(ker_pole[i].shape))),
                                np.hstack((np.zeros(ker_pole[j].shape),
                                                    ker_pole[j]))
                                ))
        mu_nu_matrix = np.vstack(
            (np.hstack((mu1, mu2)), np.hstack((nu1, nu2)))
            )
        ker_pole_mu_nu = np.dot(ker_pole_ij, mu_nu_matrix)
    transfer_matrix_ij = np.dot(np.dot(ker_pole_mu_nu, ker_pole_mu_nu.T),
                             transfer_matrix_j_mo_transfer_matrix_j)
    if not np.allclose(transfer_matrix_ij, 0):
        transfer_matrix_ij = (np.sqrt(2)*transfer_matrix_ij /
                              np.linalg.norm(transfer_matrix_ij))
        transfer_matrix[:, i] = transfer_matrix_ij[
            :transfer_matrix[:, i].shape[0], 0
            ]
        transfer_matrix[:, j] = transfer_matrix_ij[
            transfer_matrix[:, i].shape[0]:, 0
            ]
    else:
        # As in knv0 if transfer_matrix_j_mo_transfer_matrix_j is orthogonal to
        # Vect{ker_pole_mu_nu} assign transfer_matrixi/transfer_matrix_j to
        # ker_pole_mu_nu and iterate. As we are looking for a vector in
        # Vect{Matker_pole_MU_NU} (see section 6.1 page 19) this might help
        # (that's a guess, not a claim !)
        transfer_matrix[:, i] = ker_pole_mu_nu[
            :transfer_matrix[:, i].shape[0], 0
            ]
        transfer_matrix[:, j] = ker_pole_mu_nu[
            transfer_matrix[:, i].shape[0]:, 0
            ]

################################################################################
# Method Path: scipy.optimize._optimize.polak_ribiere_powell_step
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  3617
# Total Exec:  3617
# Func Size:   9
################################################################################

def polak_ribiere_powell_step(alpha, gfkp1=None):
            xkp1 = xk + alpha * pk
            if gfkp1 is None:
                gfkp1 = myfprime(xkp1)
            yk = gfkp1 - gfk
            beta_k = max(0, np.dot(yk, gfkp1) / deltak)
            pkp1 = -gfkp1 + beta_k * pk
            gnorm = vecnorm(gfkp1, ord=norm)
            return (alpha, xkp1, pkp1, gfkp1, gnorm)

################################################################################
# Method Path: scipy.optimize._lsq.common.solve_trust_region_2d
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/common.py
# Target:      dot
# Call Count:  3264
# Total Exec:  1683
# Func Size:   49
################################################################################

def solve_trust_region_2d(B, g, Delta):
    """Solve a general trust-region problem in 2 dimensions.

    The problem is reformulated as a 4th order algebraic equation,
    the solution of which is found by numpy.roots.

    Parameters
    ----------
    B : ndarray, shape (2, 2)
        Symmetric matrix, defines a quadratic term of the function.
    g : ndarray, shape (2,)
        Defines a linear term of the function.
    Delta : float
        Radius of a trust region.

    Returns
    -------
    p : ndarray, shape (2,)
        Found solution.
    newton_step : bool
        Whether the returned solution is the Newton step which lies within
        the trust region.
    """
    try:
        R, lower = cho_factor(B)
        p = -cho_solve((R, lower), g)
        if np.dot(p, p) <= Delta**2:
            return p, True
    except LinAlgError:
        pass

    a = B[0, 0] * Delta**2
    b = B[0, 1] * Delta**2
    c = B[1, 1] * Delta**2

    d = g[0] * Delta
    f = g[1] * Delta

    coeffs = np.array(
        [-b + d, 2 * (a - c + f), 6 * b, 2 * (-a + c + f), -b - d])
    t = np.roots(coeffs)  # Can handle leading zeros.
    t = np.real(t[np.isreal(t)])

    p = Delta * np.vstack((2 * t / (1 + t**2), (1 - t**2) / (1 + t**2)))
    value = 0.5 * np.sum(p * B.dot(p), axis=0) + np.dot(g, p)
    i = np.argmin(value)
    p = p[:, i]

    return p, False

################################################################################
# Method Path: scipy.optimize._lsq.dogbox.dogbox
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/dogbox.py
# Target:      dot
# Call Count:  2834
# Total Exec:  393
# Func Size:   195
################################################################################

def dogbox(fun, jac, x0, f0, J0, lb, ub, ftol, xtol, gtol, max_nfev, x_scale,
           loss_function, tr_solver, tr_options, verbose, callback=None):
    f = f0
    f_true = f.copy()
    nfev = 1

    J = J0
    njev = 1

    if loss_function is not None:
        rho = loss_function(f)
        cost = 0.5 * np.sum(rho[0])
        J, f = scale_for_robust_loss_function(J, f, rho)
    else:
        cost = 0.5 * np.dot(f, f)

    g = compute_grad(J, f)

    jac_scale = isinstance(x_scale, str) and x_scale == 'jac'
    if jac_scale:
        scale, scale_inv = compute_jac_scale(J)
    else:
        scale, scale_inv = x_scale, 1 / x_scale

    Delta = norm(x0 * scale_inv, ord=np.inf)
    if Delta == 0:
        Delta = 1.0

    on_bound = np.zeros_like(x0, dtype=int)
    on_bound[np.equal(x0, lb)] = -1
    on_bound[np.equal(x0, ub)] = 1

    x = x0
    step = np.empty_like(x0)

    if max_nfev is None:
        max_nfev = x0.size * 100

    termination_status = None
    iteration = 0
    step_norm = None
    actual_reduction = None

    if verbose == 2:
        print_header_nonlinear()

    while True:
        active_set = on_bound * g < 0
        free_set = ~active_set

        g_free = g[free_set]
        g_full = g.copy()
        g[active_set] = 0

        g_norm = norm(g, ord=np.inf)
        if g_norm < gtol:
            termination_status = 1

        if verbose == 2:
            print_iteration_nonlinear(iteration, nfev, cost, actual_reduction,
                                      step_norm, g_norm)

        if termination_status is not None or nfev == max_nfev:
            break

        x_free = x[free_set]
        lb_free = lb[free_set]
        ub_free = ub[free_set]
        scale_free = scale[free_set]

        # Compute (Gauss-)Newton and build quadratic model for Cauchy step.
        if tr_solver == 'exact':
            J_free = J[:, free_set]
            newton_step = lstsq(J_free, -f, rcond=-1)[0]

            # Coefficients for the quadratic model along the anti-gradient.
            a, b = build_quadratic_1d(J_free, g_free, -g_free)
        elif tr_solver == 'lsmr':
            Jop = aslinearoperator(J)

            # We compute lsmr step in scaled variables and then
            # transform back to normal variables, if lsmr would give exact lsq
            # solution, this would be equivalent to not doing any
            # transformations, but from experience it's better this way.

            # We pass active_set to make computations as if we selected
            # the free subset of J columns, but without actually doing any
            # slicing, which is expensive for sparse matrices and impossible
            # for LinearOperator.

            lsmr_op = lsmr_operator(Jop, scale, active_set)
            newton_step = -lsmr(lsmr_op, f, **tr_options)[0][free_set]
            newton_step *= scale_free

            # Components of g for active variables were zeroed, so this call
            # is correct and equivalent to using J_free and g_free.
            a, b = build_quadratic_1d(Jop, g, -g)

        actual_reduction = -1.0
        while actual_reduction <= 0 and nfev < max_nfev:
            tr_bounds = Delta * scale_free

            step_free, on_bound_free, tr_hit = dogleg_step(
                x_free, newton_step, g_free, a, b, tr_bounds, lb_free, ub_free)

            step.fill(0.0)
            step[free_set] = step_free

            if tr_solver == 'exact':
                predicted_reduction = -evaluate_quadratic(J_free, g_free,
                                                          step_free)
            elif tr_solver == 'lsmr':
                predicted_reduction = -evaluate_quadratic(Jop, g, step)

            # gh11403 ensure that solution is fully within bounds.
            x_new = np.clip(x + step, lb, ub)

            f_new = fun(x_new)
            nfev += 1

            step_h_norm = norm(step * scale_inv, ord=np.inf)

            if not np.all(np.isfinite(f_new)):
                Delta = 0.25 * step_h_norm
                continue

            # Usual trust-region step quality estimation.
            if loss_function is not None:
                cost_new = loss_function(f_new, cost_only=True)
            else:
                cost_new = 0.5 * np.dot(f_new, f_new)
            actual_reduction = cost - cost_new

            Delta, ratio = update_tr_radius(
                Delta, actual_reduction, predicted_reduction,
                step_h_norm, tr_hit
            )

            step_norm = norm(step)
            termination_status = check_termination(
                actual_reduction, cost, step_norm, norm(x), ratio, ftol, xtol)

            if termination_status is not None:
                break

        if actual_reduction > 0:
            on_bound[free_set] = on_bound_free

            x = x_new
            # Set variables exactly at the boundary.
            mask = on_bound == -1
            x[mask] = lb[mask]
            mask = on_bound == 1
            x[mask] = ub[mask]

            f = f_new
            f_true = f.copy()

            cost = cost_new

            J = jac(x)
            njev += 1

            if loss_function is not None:
                rho = loss_function(f)
                J, f = scale_for_robust_loss_function(J, f, rho)

            g = compute_grad(J, f)

            if jac_scale:
                scale, scale_inv = compute_jac_scale(J, scale_inv)
        else:
            step_norm = 0
            actual_reduction = 0

        iteration += 1
        
        # Call callback function and possibly stop optimization
        if callback is not None:
            intermediate_result = OptimizeResult(
                x=x, fun=f, nit=iteration, nfev=nfev)
            intermediate_result["cost"] = cost_new

            if _call_callback_maybe_halt(
                callback, intermediate_result
            ):
                termination_status = -2
                break

    if termination_status is None:
        termination_status = 0

    return OptimizeResult(
        x=x, cost=cost, fun=f_true, jac=J, grad=g_full, optimality=g_norm,
        active_mask=on_bound, nfev=nfev, njev=njev, status=termination_status)

################################################################################
# Method Path: scipy.optimize._hessian_update_strategy.SR1._update_implementation
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_hessian_update_strategy.py
# Target:      dot
# Call Count:  2814
# Total Exec:  2814
# Func Size:   21
################################################################################

def _update_implementation(self, delta_x, delta_grad):
        # Auxiliary variables w and z
        if self.approx_type == 'hess':
            w = delta_x
            z = delta_grad
        else:
            w = delta_grad
            z = delta_x
        # Do some common operations
        Mw = self @ w
        z_minus_Mw = z - Mw
        denominator = np.dot(w, z_minus_Mw)
        # If the denominator is too small
        # we just skip the update.
        if np.abs(denominator) <= self.min_denominator*norm(w)*norm(z_minus_Mw):
            return
        # Update matrix
        if self.approx_type == 'hess':
            self.B = self._syr(1/denominator, z_minus_Mw, a=self.B)
        else:
            self.H = self._syr(1/denominator, z_minus_Mw, a=self.H)

################################################################################
# Method Path: scipy.optimize._optimize._minimize_cg
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  2776
# Total Exec:  546
# Func Size:   160
################################################################################

def _minimize_cg(fun, x0, args=(), jac=None, callback=None,
                 gtol=1e-5, norm=np.inf, eps=_epsilon, maxiter=None,
                 disp=False, return_all=False, finite_diff_rel_step=None,
                 c1=1e-4, c2=0.4, workers=None,
                 **unknown_options):
    """
    Minimization of scalar function of one or more variables using the
    conjugate gradient algorithm.

    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter : int
        Maximum number of iterations to perform.
    gtol : float
        Gradient norm must be less than `gtol` before successful
        termination.
    norm : float
        Order of norm (Inf is max, -Inf is min).
    eps : float or ndarray
        If `jac is None` the absolute step size used for numerical
        approximation of the jacobian via forward differences.
    return_all : bool, optional
        Set to True to return a list of the best solution at each of the
        iterations.
    finite_diff_rel_step : None or array_like, optional
        If ``jac in ['2-point', '3-point', 'cs']`` the relative step size to
        use for numerical approximation of the jacobian. The absolute step
        size is computed as ``h = rel_step * sign(x) * max(1, abs(x))``,
        possibly adjusted to fit into the bounds. For ``jac='3-point'``
        the sign of `h` is ignored. If None (default) then step is selected
        automatically.
    c1 : float, default: 1e-4
        Parameter for Armijo condition rule.
    c2 : float, default: 0.4
        Parameter for curvature condition rule.
    workers : int, map-like callable, optional
        A map-like callable, such as `multiprocessing.Pool.map` for evaluating
        any numerical differentiation in parallel.
        This evaluation is carried out as ``workers(fun, iterable)``.

        .. versionadded:: 1.16.0

    Notes
    -----
    Parameters `c1` and `c2` must satisfy ``0 < c1 < c2 < 1``.
    """
    _check_unknown_options(unknown_options)

    retall = return_all

    x0 = asarray(x0).flatten()
    if maxiter is None:
        maxiter = len(x0) * 200

    sf = _prepare_scalar_function(fun, x0, jac=jac, args=args, epsilon=eps,
                                  finite_diff_rel_step=finite_diff_rel_step,
                                  workers=workers)

    f = sf.fun
    myfprime = sf.grad

    old_fval = f(x0)
    gfk = myfprime(x0)

    k = 0
    xk = x0
    # Sets the initial step guess to dx ~ 1
    old_old_fval = old_fval + np.linalg.norm(gfk) / 2

    if retall:
        allvecs = [xk]
    warnflag = 0
    pk = -gfk
    gnorm = vecnorm(gfk, ord=norm)

    sigma_3 = 0.01

    while (gnorm > gtol) and (k < maxiter):
        deltak = np.dot(gfk, gfk)

        cached_step = [None]

        def polak_ribiere_powell_step(alpha, gfkp1=None):
            xkp1 = xk + alpha * pk
            if gfkp1 is None:
                gfkp1 = myfprime(xkp1)
            yk = gfkp1 - gfk
            beta_k = max(0, np.dot(yk, gfkp1) / deltak)
            pkp1 = -gfkp1 + beta_k * pk
            gnorm = vecnorm(gfkp1, ord=norm)
            return (alpha, xkp1, pkp1, gfkp1, gnorm)

        def descent_condition(alpha, xkp1, fp1, gfkp1):
            # Polak-Ribiere+ needs an explicit check of a sufficient
            # descent condition, which is not guaranteed by strong Wolfe.
            #
            # See Gilbert & Nocedal, "Global convergence properties of
            # conjugate gradient methods for optimization",
            # SIAM J. Optimization 2, 21 (1992).
            cached_step[:] = polak_ribiere_powell_step(alpha, gfkp1)
            alpha, xk, pk, gfk, gnorm = cached_step

            # Accept step if it leads to convergence.
            if gnorm <= gtol:
                return True

            # Accept step if sufficient descent condition applies.
            return np.dot(pk, gfk) <= -sigma_3 * np.dot(gfk, gfk)

        try:
            alpha_k, fc, gc, old_fval, old_old_fval, gfkp1 = \
                     _line_search_wolfe12(f, myfprime, xk, pk, gfk, old_fval,
                                          old_old_fval, c1=c1, c2=c2, amin=1e-100,
                                          amax=1e100, extra_condition=descent_condition)
        except _LineSearchError:
            # Line search failed to find a better solution.
            warnflag = 2
            break

        # Reuse already computed results if possible
        if alpha_k == cached_step[0]:
            alpha_k, xk, pk, gfk, gnorm = cached_step
        else:
            alpha_k, xk, pk, gfk, gnorm = polak_ribiere_powell_step(alpha_k, gfkp1)

        if retall:
            allvecs.append(xk)
        k += 1
        intermediate_result = OptimizeResult(x=xk, fun=old_fval)
        if _call_callback_maybe_halt(callback, intermediate_result):
            break

    fval = old_fval
    if warnflag == 2:
        msg = _status_message['pr_loss']
    elif k >= maxiter:
        warnflag = 1
        msg = _status_message['maxiter']
    elif np.isnan(gnorm) or np.isnan(fval) or np.isnan(xk).any():
        warnflag = 3
        msg = _status_message['nan']
    else:
        msg = _status_message['success']

    if disp:
        _print_success_message_or_warn(warnflag, msg)
        print(f"         Current function value: {fval:f}")
        print(f"         Iterations: {k:d}")
        print(f"         Function evaluations: {sf.nfev:d}")
        print(f"         Gradient evaluations: {sf.ngev:d}")

    result = OptimizeResult(fun=fval, jac=gfk, nfev=sf.nfev,
                            njev=sf.ngev, status=warnflag,
                            success=(warnflag == 0), message=msg, x=xk,
                            nit=k)
    if retall:
        result['allvecs'] = allvecs
    return result

################################################################################
# Method Path: scipy.optimize._differentiable_functions._VectorHessWrapper.jac_dot_v
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_differentiable_functions.py
# Target:      dot
# Call Count:  2716
# Total Exec:  2716
# Func Size:   3
################################################################################

def jac_dot_v(self, x, v):
        self.njev += 1
        return self.jac(x).T.dot(v)

################################################################################
# Method Path: scipy.optimize._trustregion_exact.IterativeSubproblem.solve
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_exact.py
# Target:      dot
# Call Count:  2504
# Total Exec:  1031
# Func Size:   155
################################################################################

def solve(self, tr_radius):
        """Solve quadratic subproblem"""

        lambda_current, lambda_lb, lambda_ub = self._initial_values(tr_radius)
        n = self.dimension
        hits_boundary = True
        already_factorized = False
        self.niter = 0

        while self.niter < self.maxiter:

            # Compute Cholesky factorization
            if already_factorized:
                already_factorized = False
            else:
                H = self.hess+lambda_current*np.eye(n)
                U, info = self.cholesky(H, lower=False,
                                        overwrite_a=False,
                                        clean=True)

            self.niter += 1

            # Check if factorization succeeded
            if info == 0 and self.jac_mag > self.CLOSE_TO_ZERO:
                # Successful factorization

                # Solve `U.T U p = s`
                p = cho_solve((U, False), -self.jac)

                p_norm = norm(p)

                # Check for interior convergence
                if p_norm <= tr_radius and lambda_current == 0:
                    hits_boundary = False
                    break

                # Solve `U.T w = p`
                w = solve_triangular(U, p, trans='T')

                w_norm = norm(w)

                # Compute Newton step accordingly to
                # formula (4.44) p.87 from ref [2]_.
                delta_lambda = (p_norm/w_norm)**2 * (p_norm-tr_radius)/tr_radius
                lambda_new = lambda_current + delta_lambda

                if p_norm < tr_radius:  # Inside boundary
                    s_min, z_min = estimate_smallest_singular_value(U)

                    ta, tb = self.get_boundaries_intersections(p, z_min,
                                                               tr_radius)

                    # Choose `step_len` with the smallest magnitude.
                    # The reason for this choice is explained at
                    # ref [3]_, p. 6 (Immediately before the formula
                    # for `tau`).
                    step_len = min([ta, tb], key=abs)

                    # Compute the quadratic term  (p.T*H*p)
                    quadratic_term = np.dot(p, np.dot(H, p))

                    # Check stop criteria
                    relative_error = ((step_len**2 * s_min**2)
                                      / (quadratic_term + lambda_current*tr_radius**2))
                    if relative_error <= self.k_hard:
                        p += step_len * z_min
                        break

                    # Update uncertainty bounds
                    lambda_ub = lambda_current
                    lambda_lb = max(lambda_lb, lambda_current - s_min**2)

                    # Compute Cholesky factorization
                    H = self.hess + lambda_new*np.eye(n)
                    c, info = self.cholesky(H, lower=False,
                                            overwrite_a=False,
                                            clean=True)

                    # Check if the factorization have succeeded
                    #
                    if info == 0:  # Successful factorization
                        # Update damping factor
                        lambda_current = lambda_new
                        already_factorized = True
                        U = c
                    else:  # Unsuccessful factorization
                        # Update uncertainty bounds
                        lambda_lb = max(lambda_lb, lambda_new)

                        # Update damping factor
                        lambda_current = max(
                            np.sqrt(np.abs(lambda_lb * lambda_ub)),
                            lambda_lb + self.UPDATE_COEFF*(lambda_ub-lambda_lb)
                        )

                else:  # Outside boundary
                    # Check stop criteria
                    relative_error = abs(p_norm - tr_radius) / tr_radius
                    if relative_error <= self.k_easy:
                        break

                    # Update uncertainty bounds
                    lambda_lb = lambda_current

                    # Update damping factor
                    lambda_current = lambda_new

            elif info == 0 and self.jac_mag <= self.CLOSE_TO_ZERO:
                # jac_mag very close to zero

                # Check for interior convergence
                if lambda_current == 0:
                    p = np.zeros(n)
                    hits_boundary = False
                    break

                s_min, z_min = estimate_smallest_singular_value(U)
                step_len = tr_radius

                p = step_len * z_min
                # Check stop criteria
                if (step_len**2 * s_min**2
                    <= self.k_hard * lambda_current * tr_radius**2):
                    break

                # Update uncertainty bounds
                lambda_ub = lambda_current
                lambda_lb = max(lambda_lb, lambda_current - s_min**2)

                # Update damping factor
                lambda_current = max(
                    np.sqrt(lambda_lb * lambda_ub),
                    lambda_lb + self.UPDATE_COEFF*(lambda_ub-lambda_lb)
                )

            else:  # Unsuccessful factorization

                # Compute auxiliary terms
                delta, v = singular_leading_submatrix(H, U, info)
                v_norm = norm(v)

                # Update uncertainty interval
                lambda_lb = max(lambda_lb, lambda_current + delta/v_norm**2)

                # Update damping factor
                lambda_current = max(
                    np.sqrt(np.abs(lambda_lb * lambda_ub)),
                    lambda_lb + self.UPDATE_COEFF*(lambda_ub-lambda_lb)
                )

        self.lambda_lb = lambda_lb
        self.lambda_current = lambda_current
        self.previous_tr_radius = tr_radius

        return p, hits_boundary

################################################################################
# Method Path: scipy.optimize._linprog_util._postsolve
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_util.py
# Target:      dot
# Call Count:  2445
# Total Exec:  815
# Func Size:   102
################################################################################

def _postsolve(x, postsolve_args, complete=False):
    """
    Given solution x to presolved, standard form linear program x, add
    fixed variables back into the problem and undo the variable substitutions
    to get solution to original linear program. Also, calculate the objective
    function value, slack in original upper bound constraints, and residuals
    in original equality constraints.

    Parameters
    ----------
    x : 1-D array
        Solution vector to the standard-form problem.
    postsolve_args : tuple
        Data needed by _postsolve to convert the solution to the standard-form
        problem into the solution to the original problem, including:

    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, lower bounds in the 1st column, upper
            bounds in the 2nd column. The bounds are possibly tightened
            by the presolve procedure.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    revstack: list of functions
        the functions in the list reverse the operations of _presolve()
        the function signature is x_org = f(x_mod), where x_mod is the result
        of a presolve step and x_org the value at the start of the step
    complete : bool
        Whether the solution is was determined in presolve (``True`` if so)

    Returns
    -------
    x : 1-D array
        Solution vector to original linear programming problem
    fun: float
        optimal objective value for original problem
    slack : 1-D array
        The (non-negative) slack in the upper bound constraints, that is,
        ``b_ub - A_ub @ x``
    con : 1-D array
        The (nominally zero) residuals of the equality constraints, that is,
        ``b - A_eq @ x``
    """
    # note that all the inputs are the ORIGINAL, unmodified versions
    # no rows, columns have been removed

    c, A_ub, b_ub, A_eq, b_eq, bounds, x0, integrality = postsolve_args[0]
    revstack, C, b_scale = postsolve_args[1:]

    x = _unscale(x, C, b_scale)

    # Undo variable substitutions of _get_Abc()
    # if "complete", problem was solved in presolve; don't do anything here
    n_x = bounds.shape[0]
    if not complete and bounds is not None:  # bounds are never none, probably
        n_unbounded = 0
        for i, bi in enumerate(bounds):
            lbi = bi[0]
            ubi = bi[1]
            if lbi == -np.inf and ubi == np.inf:
                n_unbounded += 1
                x[i] = x[i] - x[n_x + n_unbounded - 1]
            else:
                if lbi == -np.inf:
                    x[i] = ubi - x[i]
                else:
                    x[i] += lbi
    # all the rest of the variables were artificial
    x = x[:n_x]

    # If there were variables removed from the problem, add them back into the
    # solution vector
    # Apply the functions in revstack (reverse direction)
    for rev in reversed(revstack):
        x = rev(x)

    fun = x.dot(c)
    with np.errstate(invalid="ignore"):
        slack = b_ub - A_ub.dot(x)  # report slack for ORIGINAL UB constraints
        # report residuals of ORIGINAL EQ constraints
        con = b_eq - A_eq.dot(x)

    return x, fun, slack, con

################################################################################
# Method Path: scipy.optimize._lbfgsb_py.LbfgsInvHessProduct._matvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lbfgsb_py.py
# Target:      dot
# Call Count:  2364
# Total Exec:  48
# Func Size:   33
################################################################################

def _matvec(self, x):
        """Efficient matrix-vector multiply with the BFGS matrices.

        This calculation is described in Section (4) of [1].

        Parameters
        ----------
        x : ndarray
            An array with shape (n,) or (n,1).

        Returns
        -------
        y : ndarray
            The matrix-vector product

        """
        s, y, n_corrs, rho = self.sk, self.yk, self.n_corrs, self.rho
        q = np.array(x, dtype=self.dtype, copy=True)
        if q.ndim == 2 and q.shape[1] == 1:
            q = q.reshape(-1)

        alpha = np.empty(n_corrs)

        for i in range(n_corrs-1, -1, -1):
            alpha[i] = rho[i] * np.dot(s[i], q)
            q = q - alpha[i]*y[i]

        r = q
        for i in range(n_corrs):
            beta = rho[i] * np.dot(y[i], r)
            r = r + s[i] * (alpha[i] - beta)

        return r

################################################################################
# Method Path: scipy.optimize._lbfgsb_py.LbfgsInvHessProduct._matmat
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lbfgsb_py.py
# Target:      dot
# Call Count:  2364
# Total Exec:  48
# Func Size:   36
################################################################################

def _matmat(self, X):
        """Efficient matrix-matrix multiply with the BFGS matrices.

        This calculation is described in Section (4) of [1].

        Parameters
        ----------
        X : ndarray
            An array with shape (n,m)

        Returns
        -------
        Y : ndarray
            The matrix-matrix product

        Notes
        -----
        This implementation is written starting from _matvec and broadcasting
        all expressions along the second axis of X.

        """
        s, y, n_corrs, rho = self.sk, self.yk, self.n_corrs, self.rho
        Q = np.array(X, dtype=self.dtype, copy=True)

        alpha = np.empty((n_corrs, Q.shape[1]))

        for i in range(n_corrs-1, -1, -1):
            alpha[i] = rho[i] * np.dot(s[i], Q)
            Q -= alpha[i]*y[i][:, np.newaxis]

        R = Q
        for i in range(n_corrs):
            beta = rho[i] * np.dot(y[i], R)
            R += s[i][:, np.newaxis] * (alpha[i] - beta)

        return R

################################################################################
# Method Path: scipy.integrate._ivp.radau.RadauDenseOutput._call_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/radau.py
# Target:      dot
# Call Count:  2298
# Total Exec:  2298
# Func Size:   16
################################################################################

def _call_impl(self, t):
        x = (t - self.t_old) / self.h
        if t.ndim == 0:
            p = np.tile(x, self.order + 1)
            p = np.cumprod(p)
        else:
            p = np.tile(x, (self.order + 1, 1))
            p = np.cumprod(p, axis=0)
        # Here we don't multiply by h, not a mistake.
        y = np.dot(self.Q, p)
        if y.ndim == 2:
            y += self.y_old[:, None]
        else:
            y += self.y_old

        return y

################################################################################
# Method Path: scipy.optimize._trustregion.DoglegSubproblem.__call__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion.py
# Target:      dot
# Call Count:  2246
# Total Exec:  1123
# Func Size:   2
################################################################################

def __call__(self, p):
        return self.fun + np.dot(self.jac, p) + 0.5 * np.dot(p, self.hessp(p))

################################################################################
# Method Path: scipy.optimize._lsq.common.intersect_trust_region
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/common.py
# Target:      dot
# Call Count:  2242
# Total Exec:  748
# Func Size:   37
################################################################################

def intersect_trust_region(x, s, Delta):
    """Find the intersection of a line with the boundary of a trust region.

    This function solves the quadratic equation with respect to t
    ||(x + s*t)||**2 = Delta**2.

    Returns
    -------
    t_neg, t_pos : tuple of float
        Negative and positive roots.

    Raises
    ------
    ValueError
        If `s` is zero or `x` is not within the trust region.
    """
    a = np.dot(s, s)
    if a == 0:
        raise ValueError("`s` is zero.")

    b = np.dot(x, s)

    c = np.dot(x, x) - Delta**2
    if c > 0:
        raise ValueError("`x` is not within the trust region.")

    d = np.sqrt(b*b - a*c)  # Root from one fourth of the discriminant.

    # Computations below avoid loss of significance, see "Numerical Recipes".
    q = -(b + copysign(d, b))
    t1 = q / a
    t2 = c / q

    if t1 < t2:
        return t1, t2
    else:
        return t2, t1

################################################################################
# Method Path: scipy.optimize._linprog_ip.r_d
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  2228
# Total Exec:  2228
# Func Size:   2
################################################################################

def r_d(y, z, tau):
        return c * tau - A.T.dot(y) - z

################################################################################
# Method Path: scipy.optimize._linprog_ip.r_p
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  2228
# Total Exec:  2228
# Func Size:   2
################################################################################

def r_p(x, tau):
        return b * tau - A.dot(x)

################################################################################
# Method Path: scipy.optimize._trustregion._minimize_trust_region
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion.py
# Target:      dot
# Call Count:  2097
# Total Exec:  94
# Func Size:   212
################################################################################

def _minimize_trust_region(fun, x0, args=(), jac=None, hess=None, hessp=None,
                           subproblem=None, initial_trust_radius=1.0,
                           max_trust_radius=1000.0, eta=0.15, gtol=1e-4,
                           maxiter=None, disp=False, return_all=False,
                           callback=None, inexact=True, workers=None,
                           subproblem_maxiter=None, **unknown_options):
    """
    Minimization of scalar function of one or more variables using a
    trust-region algorithm.

    Options for the trust-region algorithm are:
        initial_trust_radius : float
            Initial trust radius.
        max_trust_radius : float
            Never propose steps that are longer than this value.
        eta : float
            Trust region related acceptance stringency for proposed steps.
        gtol : float
            Gradient norm must be less than `gtol`
            before successful termination.
        maxiter : int
            Maximum number of iterations to perform.
        disp : bool
            If True, print convergence message.
        inexact : bool
            Accuracy to solve subproblems. If True requires less nonlinear
            iterations, but more vector products. Only effective for method
            trust-krylov.
        workers : int, map-like callable, optional
            A map-like callable, such as `multiprocessing.Pool.map` for evaluating
            any numerical differentiation in parallel.
            This evaluation is carried out as ``workers(fun, iterable)``.
            Only for 'trust-krylov', 'trust-ncg'.

            .. versionadded:: 1.16.0
        subproblem_maxiter : int, optional
            Maximum number of iterations to perform per subproblem. Only affects
            trust-exact. Default is 25.

            .. versionadded:: 1.17.0


    This function is called by the `minimize` function.
    It is not supposed to be called directly.
    """
    _check_unknown_options(unknown_options)

    if jac is None:
        raise ValueError('Jacobian is currently required for trust-region '
                         'methods')
    if hess is None and hessp is None:
        raise ValueError('Either the Hessian or the Hessian-vector product '
                         'is currently required for trust-region methods')
    if subproblem is None:
        raise ValueError('A subproblem solving strategy is required for '
                         'trust-region methods')
    if not (0 <= eta < 0.25):
        raise Exception('invalid acceptance stringency')
    if max_trust_radius <= 0:
        raise Exception('the max trust radius must be positive')
    if initial_trust_radius <= 0:
        raise ValueError('the initial trust radius must be positive')
    if initial_trust_radius >= max_trust_radius:
        raise ValueError('the initial trust radius must be less than the '
                         'max trust radius')

    # force the initial guess into a nice format
    x0 = np.asarray(x0).flatten()

    # A ScalarFunction representing the problem. This caches calls to fun, jac,
    # hess.
    # the workers kwd only has an effect for trust-ncg, trust-krylov when
    # estimating the Hessian with finite-differences. It's never used
    # during calculation of jacobian, because callables are required for all
    # methods.
    sf = _prepare_scalar_function(
        fun, x0, jac=jac, hess=hess, args=args, workers=workers
    )
    fun = sf.fun
    jac = sf.grad
    if callable(hess):
        hess = sf.hess
    elif callable(hessp):
        # this elif statement must come before examining whether hess
        # is estimated by FD methods or a HessianUpdateStrategy
        pass
    elif (hess in FD_METHODS or isinstance(hess, HessianUpdateStrategy)):
        # If the Hessian is being estimated by finite differences or a
        # Hessian update strategy then ScalarFunction.hess returns a
        # LinearOperator or a HessianUpdateStrategy. This enables the
        # calculation/creation of a hessp. BUT you only want to do this
        # if the user *hasn't* provided a callable(hessp) function.
        hess = None

        def hessp(x, p, *args):
            return sf.hess(x).dot(p)
    else:
        raise ValueError('Either the Hessian or the Hessian-vector product '
                         'is currently required for trust-region methods')

    # ScalarFunction doesn't represent hessp
    nhessp, hessp = _wrap_function(hessp, args)

    # limit the number of iterations
    if maxiter is None:
        maxiter = len(x0)*200

    # init the search status
    warnflag = 0

    # initialize the search
    trust_radius = initial_trust_radius
    x = x0
    if return_all:
        allvecs = [x]

    subproblem_init_kw = {}
    if hasattr(subproblem, 'MAXITER_DEFAULT'):
        subproblem_init_kw['maxiter'] = subproblem_maxiter

    m = subproblem(x, fun, jac, hess, hessp, **subproblem_init_kw)
    k = 0

    # search for the function min
    # do not even start if the gradient is small enough
    while m.jac_mag >= gtol:

        # Solve the sub-problem.
        # This gives us the proposed step relative to the current position
        # and it tells us whether the proposed step
        # has reached the trust region boundary or not.
        try:
            p, hits_boundary = m.solve(trust_radius)
        except np.linalg.LinAlgError:
            warnflag = 3
            break

        # calculate the predicted value at the proposed point
        predicted_value = m(p)

        # define the local approximation at the proposed point
        x_proposed = x + p
        m_proposed = subproblem(x_proposed, fun, jac, hess, hessp, **subproblem_init_kw)

        # evaluate the ratio defined in equation (4.4)
        actual_reduction = m.fun - m_proposed.fun
        predicted_reduction = m.fun - predicted_value
        if predicted_reduction <= 0:
            warnflag = 2
            break
        rho = actual_reduction / predicted_reduction

        # update the trust radius according to the actual/predicted ratio
        if rho < 0.25:
            trust_radius *= 0.25
        elif rho > 0.75 and hits_boundary:
            trust_radius = min(2*trust_radius, max_trust_radius)

        # if the ratio is high enough then accept the proposed step
        if rho > eta:
            x = x_proposed
            m = m_proposed

        # append the best guess, call back, increment the iteration count
        if return_all:
            allvecs.append(np.copy(x))
        k += 1

        intermediate_result = OptimizeResult(x=x, fun=m.fun)
        if _call_callback_maybe_halt(callback, intermediate_result):
            break

        # check if the gradient is small enough to stop
        if m.jac_mag < gtol:
            warnflag = 0
            break

        # check if we have looked at enough iterations
        if k >= maxiter:
            warnflag = 1
            break

    # print some stuff if requested
    status_messages = (
            _status_message['success'],
            _status_message['maxiter'],
            'A bad approximation caused failure to predict improvement.',
            'A linalg error occurred, such as a non-psd Hessian.',
            )
    if disp:
        if warnflag == 0:
            print(status_messages[warnflag])
        else:
            warnings.warn(status_messages[warnflag], RuntimeWarning, stacklevel=3)
        print(f"         Current function value: {m.fun:f}")
        print(f"         Iterations: {k:d}")
        print(f"         Function evaluations: {sf.nfev:d}")
        print(f"         Gradient evaluations: {sf.ngev:d}")
        print(f"         Hessian evaluations: {sf.nhev + nhessp[0]:d}")

    result = OptimizeResult(x=x, success=(warnflag == 0), status=warnflag,
                            fun=m.fun, jac=m.jac, nfev=sf.nfev, njev=sf.ngev,
                            nhev=sf.nhev + nhessp[0], nit=k,
                            message=status_messages[warnflag])

    if hess is not None:
        result['hess'] = m.hess

    if return_all:
        result['allvecs'] = allvecs

    return result

################################################################################
# Method Path: scipy.linalg._matfuncs_inv_ssq._remainder_matrix_power_triu
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs_inv_ssq.py
# Target:      dot
# Call Count:  2080
# Total Exec:  535
# Func Size:   77
################################################################################

def _remainder_matrix_power_triu(T, t):
    """
    Compute a fractional power of an upper triangular matrix.

    The fractional power is restricted to fractions -1 < t < 1.
    This uses algorithm (3.1) of [1]_.
    The Pade approximation itself uses algorithm (4.1) of [2]_.

    Parameters
    ----------
    T : (N, N) array_like
        Upper triangular matrix whose fractional power to evaluate.
    t : float
        Fractional power between -1 and 1 exclusive.

    Returns
    -------
    X : (N, N) array_like
        The fractional power of the matrix.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing Lin (2013)
           "An Improved Schur-Pade Algorithm for Fractional Powers
           of a Matrix and their Frechet Derivatives."

    .. [2] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    m_to_theta = {
            1: 1.51e-5,
            2: 2.24e-3,
            3: 1.88e-2,
            4: 6.04e-2,
            5: 1.24e-1,
            6: 2.00e-1,
            7: 2.79e-1,
            }
    n, n = T.shape
    T0 = T
    T0_diag = np.diag(T0)
    if np.array_equal(T0, np.diag(T0_diag)):
        U = np.diag(T0_diag ** t)
    else:
        R, s, m = _inverse_squaring_helper(T0, m_to_theta)

        # Evaluate the Pade approximation.
        # Note that this function expects the negative of the matrix
        # returned by the inverse squaring helper.
        U = _fractional_power_pade(-R, t, m)

        # Undo the inverse scaling and squaring.
        # Be less clever about this
        # if the principal branch does not exist at T0;
        # this happens when a diagonal entry of T0
        # is negative with imaginary part 0.
        eivals = np.diag(T0)
        has_principal_branch = all(x.real > 0 or x.imag != 0 for x in eivals)
        for i in range(s, -1, -1):
            if i < s:
                U = U.dot(U)
            else:
                if has_principal_branch:
                    p = t * np.exp2(-i)
                    U[np.diag_indices(n)] = T0_diag ** p
                    for j in range(n-1):
                        l1 = T0[j, j]
                        l2 = T0[j+1, j+1]
                        t12 = T0[j, j+1]
                        f12 = _fractional_power_superdiag_entry(l1, l2, t12, p)
                        U[j, j+1] = f12
    if not np.array_equal(U, np.triu(U)):
        raise Exception('U is not upper triangular')
    return U

################################################################################
# Method Path: scipy.optimize._linprog_rs._phase_two
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_rs.py
# Target:      dot
# Call Count:  2057
# Total Exec:  164
# Func Size:   81
################################################################################

def _phase_two(c, A, x, b, callback, postsolve_args, maxiter, tol, disp,
               maxupdate, mast, pivot, iteration=0, phase_one_n=None):
    """
    The heart of the simplex method. Beginning with a basic feasible solution,
    moves to adjacent basic feasible solutions successively lower reduced cost.
    Terminates when there are no basic feasible solutions with lower reduced
    cost or if the problem is determined to be unbounded.

    This implementation follows the revised simplex method based on LU
    decomposition. Rather than maintaining a tableau or an inverse of the
    basis matrix, we keep a factorization of the basis matrix that allows
    efficient solution of linear systems while avoiding stability issues
    associated with inverted matrices.
    """
    m, n = A.shape
    status = 0
    a = np.arange(n)                    # indices of columns of A
    ab = np.arange(m)                   # indices of columns of B
    if maxupdate:
        # basis matrix factorization object; similar to B = A[:, b]
        B = BGLU(A, b, maxupdate, mast)
    else:
        B = LU(A, b)

    for iteration in range(iteration, maxiter):

        if disp or callback is not None:
            _display_and_callback(phase_one_n, x, postsolve_args, status,
                                  iteration, disp, callback)

        bl = np.zeros(len(a), dtype=bool)
        bl[b] = 1

        xb = x[b]       # basic variables
        cb = c[b]       # basic costs

        try:
            v = B.solve(cb, transposed=True)    # similar to v = solve(B.T, cb)
        except LinAlgError:
            status = 4
            break

        # TODO: cythonize?
        c_hat = c - v.dot(A)    # reduced cost
        c_hat = c_hat[~bl]
        # Above is much faster than:
        # N = A[:, ~bl]                 # slow!
        # c_hat = c[~bl] - v.T.dot(N)
        # Can we perform the multiplication only on the nonbasic columns?

        if np.all(c_hat >= -tol):  # all reduced costs positive -> terminate
            break

        j = _select_enter_pivot(c_hat, bl, a, rule=pivot, tol=tol)
        u = B.solve(A[:, j])        # similar to u = solve(B, A[:, j])

        i = u > tol                 # if none of the u are positive, unbounded
        if not np.any(i):
            status = 3
            break

        th = xb[i]/u[i]
        l = np.argmin(th)           # implicitly selects smallest subscript
        th_star = th[l]             # step size

        x[b] = x[b] - th_star*u     # take step
        x[j] = th_star
        B.update(ab[i][l], j)       # modify basis
        b = B.b                     # similar to b[ab[i][l]] =

    else:
        # If the end of the for loop is reached (without a break statement),
        # then another step has been taken, so the iteration counter should
        # increment, info should be displayed, and callback should be called.
        iteration += 1
        status = 1
        if disp or callback is not None:
            _display_and_callback(phase_one_n, x, postsolve_args, status,
                                  iteration, disp, callback)

    return x, b, status, iteration

################################################################################
# Method Path: scipy.linalg._solvers.solve_sylvester
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_solvers.py
# Target:      dot
# Call Count:  1836
# Total Exec:  837
# Func Size:   83
################################################################################

def solve_sylvester(a, b, q):
    """
    Computes a solution (X) to the Sylvester equation :math:`AX + XB = Q`.

    Parameters
    ----------
    a : (M, M) array_like
        Leading matrix of the Sylvester equation
    b : (N, N) array_like
        Trailing matrix of the Sylvester equation
    q : (M, N) array_like
        Right-hand side

    Returns
    -------
    x : (M, N) ndarray
        The solution to the Sylvester equation.

    Raises
    ------
    LinAlgError
        If solution was not found

    Notes
    -----
    Computes a solution to the Sylvester matrix equation via the Bartels-
    Stewart algorithm. The A and B matrices first undergo Schur
    decompositions. The resulting matrices are used to construct an
    alternative Sylvester equation (``RY + YS^T = F``) where the R and S
    matrices are in quasi-triangular form (or, when R, S or F are complex,
    triangular form). The simplified equation is then solved using
    ``*TRSYL`` from LAPACK directly.

    .. versionadded:: 0.11.0

    Examples
    --------
    Given `a`, `b`, and `q` solve for `x`:

    >>> import numpy as np
    >>> from scipy import linalg
    >>> a = np.array([[-3, -2, 0], [-1, -1, 3], [3, -5, -1]])
    >>> b = np.array([[1]])
    >>> q = np.array([[1],[2],[3]])
    >>> x = linalg.solve_sylvester(a, b, q)
    >>> x
    array([[ 0.0625],
           [-0.5625],
           [ 0.6875]])
    >>> np.allclose(a.dot(x) + x.dot(b), q)
    True

    """
    # Accommodate empty a
    if a.size == 0 or b.size == 0:
        tdict = {'s': np.float32, 'd': np.float64,
                 'c': np.complex64, 'z': np.complex128}
        func, = get_lapack_funcs(('trsyl',), arrays=(a, b, q))
        return np.empty(q.shape, dtype=tdict[func.typecode])

    # Compute the Schur decomposition form of a
    r, u = schur(a, output='real')

    # Compute the Schur decomposition of b
    s, v = schur(b.conj().transpose(), output='real')

    # Construct f = u'*q*v
    f = np.dot(np.dot(u.conj().transpose(), q), v)

    # Call the Sylvester equation solver
    trsyl, = get_lapack_funcs(('trsyl',), (r, s, f))
    if trsyl is None:
        raise RuntimeError('LAPACK implementation does not contain a proper '
                           'Sylvester equation solver (TRSYL)')
    y, scale, info = trsyl(r, s, f, tranb='C')

    y = scale*y

    if info < 0:
        raise LinAlgError(f"Illegal value encountered in the {-info} term")

    return np.dot(np.dot(u, y), v.conj().transpose())

################################################################################
# Method Path: scipy.optimize._trustregion_ncg.CGSteihaugSubproblem.solve
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_ncg.py
# Target:      dot
# Call Count:  1792
# Total Exec:  411
# Func Size:   83
################################################################################

def solve(self, trust_radius):
        """
        Solve the subproblem using a conjugate gradient method.

        Parameters
        ----------
        trust_radius : float
            We are allowed to wander only this far away from the origin.

        Returns
        -------
        p : ndarray
            The proposed step.
        hits_boundary : bool
            True if the proposed step is on the boundary of the trust region.

        Notes
        -----
        This is algorithm (7.2) of Nocedal and Wright 2nd edition.
        Only the function that computes the Hessian-vector product is required.
        The Hessian itself is not required, and the Hessian does
        not need to be positive semidefinite.
        """

        # get the norm of jacobian and define the origin
        p_origin = np.zeros_like(self.jac)

        # define a default tolerance
        tolerance = min(0.5, math.sqrt(self.jac_mag)) * self.jac_mag

        # Stop the method if the search direction
        # is a direction of nonpositive curvature.
        if self.jac_mag < tolerance:
            hits_boundary = False
            return p_origin, hits_boundary

        # init the state for the first iteration
        z = p_origin
        r = self.jac
        d = -r

        # Search for the min of the approximation of the objective function.
        while True:

            # do an iteration
            Bd = self.hessp(d)
            dBd = np.dot(d, Bd)
            if dBd <= 0:
                # Look at the two boundary points.
                # Find both values of t to get the boundary points such that
                # ||z + t d|| == trust_radius
                # and then choose the one with the predicted min value.
                ta, tb = self.get_boundaries_intersections(z, d, trust_radius)
                pa = z + ta * d
                pb = z + tb * d
                if self(pa) < self(pb):
                    p_boundary = pa
                else:
                    p_boundary = pb
                hits_boundary = True
                return p_boundary, hits_boundary
            r_squared = np.dot(r, r)
            alpha = r_squared / dBd
            z_next = z + alpha * d
            if scipy.linalg.norm(z_next) >= trust_radius:
                # Find t >= 0 to get the boundary point such that
                # ||z + t d|| == trust_radius
                ta, tb = self.get_boundaries_intersections(z, d, trust_radius)
                p_boundary = z + tb * d
                hits_boundary = True
                return p_boundary, hits_boundary
            r_next = r + alpha * Bd
            r_next_squared = np.dot(r_next, r_next)
            if math.sqrt(r_next_squared) < tolerance:
                hits_boundary = False
                return z_next, hits_boundary
            beta_next = r_next_squared / r_squared
            d_next = -r_next + beta_next * d

            # update the state for the next iteration
            z = z_next
            r = r_next
            d = d_next

################################################################################
# Method Path: scipy.integrate._ivp.radau.Radau._step_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/radau.py
# Target:      dot
# Call Count:  1709
# Total Exec:  1697
# Func Size:   141
################################################################################

def _step_impl(self):
        t = self.t
        y = self.y
        f = self.f

        max_step = self.max_step
        atol = self.atol
        rtol = self.rtol

        min_step = 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t)
        if self.h_abs > max_step:
            h_abs = max_step
            h_abs_old = None
            error_norm_old = None
        elif self.h_abs < min_step:
            h_abs = min_step
            h_abs_old = None
            error_norm_old = None
        else:
            h_abs = self.h_abs
            h_abs_old = self.h_abs_old
            error_norm_old = self.error_norm_old

        J = self.J
        LU_real = self.LU_real
        LU_complex = self.LU_complex

        current_jac = self.current_jac
        jac = self.jac

        rejected = False
        step_accepted = False
        message = None
        while not step_accepted:
            if h_abs < min_step:
                return False, self.TOO_SMALL_STEP

            h = h_abs * self.direction
            t_new = t + h

            if self.direction * (t_new - self.t_bound) > 0:
                t_new = self.t_bound

            h = t_new - t
            h_abs = np.abs(h)

            if self.sol is None:
                Z0 = np.zeros((3, y.shape[0]))
            else:
                Z0 = self.sol(t + h * C).T - y

            scale = atol + np.abs(y) * rtol

            converged = False
            while not converged:
                if LU_real is None or LU_complex is None:
                    LU_real = self.lu(MU_REAL / h * self.I - J)
                    LU_complex = self.lu(MU_COMPLEX / h * self.I - J)

                converged, n_iter, Z, rate = solve_collocation_system(
                    self.fun, t, y, h, Z0, scale, self.newton_tol,
                    LU_real, LU_complex, self.solve_lu)

                if not converged:
                    if current_jac:
                        break

                    J = self.jac(t, y, f)
                    current_jac = True
                    LU_real = None
                    LU_complex = None

            if not converged:
                h_abs *= 0.5
                LU_real = None
                LU_complex = None
                continue

            y_new = y + Z[-1]
            ZE = Z.T.dot(E) / h
            error = self.solve_lu(LU_real, f + ZE)
            scale = atol + np.maximum(np.abs(y), np.abs(y_new)) * rtol
            error_norm = norm(error / scale)
            safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (2 * NEWTON_MAXITER
                                                       + n_iter)

            if rejected and error_norm > 1:
                error = self.solve_lu(LU_real, self.fun(t, y + error) + ZE)
                error_norm = norm(error / scale)

            if error_norm > 1:
                factor = predict_factor(h_abs, h_abs_old,
                                        error_norm, error_norm_old)
                h_abs *= max(MIN_FACTOR, safety * factor)

                LU_real = None
                LU_complex = None
                rejected = True
            else:
                step_accepted = True

        recompute_jac = jac is not None and n_iter > 2 and rate > 1e-3

        factor = predict_factor(h_abs, h_abs_old, error_norm, error_norm_old)
        factor = min(MAX_FACTOR, safety * factor)

        if not recompute_jac and factor < 1.2:
            factor = 1
        else:
            LU_real = None
            LU_complex = None

        f_new = self.fun(t_new, y_new)
        if recompute_jac:
            J = jac(t_new, y_new, f_new)
            current_jac = True
        elif jac is not None:
            current_jac = False

        self.h_abs_old = self.h_abs
        self.error_norm_old = error_norm

        self.h_abs = h_abs * factor

        self.y_old = y

        self.t = t_new
        self.y = y_new
        self.f = f_new

        self.Z = Z

        self.LU_real = LU_real
        self.LU_complex = LU_complex
        self.current_jac = current_jac
        self.J = J

        self.t_old = t
        self.sol = self._compute_dense_output()

        return step_accepted, message

################################################################################
# Method Path: scipy.optimize._nonlin.LowRankMatrix.matvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_nonlin.py
# Target:      dot
# Call Count:  1708
# Total Exec:  2238
# Func Size:   5
################################################################################

def matvec(self, v):
        """Evaluate w = M v"""
        if self.collapsed is not None:
            return np.dot(self.collapsed, v)
        return LowRankMatrix._matvec(v, self.alpha, self.cs, self.ds)

################################################################################
# Method Path: scipy.integrate._ivp.radau.Radau._compute_dense_output
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/radau.py
# Target:      dot
# Call Count:  1695
# Total Exec:  1695
# Func Size:   3
################################################################################

def _compute_dense_output(self):
        Q = np.dot(self.Z.T, P)
        return RadauDenseOutput(self.t_old, self.t, self.y_old, Q)

################################################################################
# Method Path: scipy.signal._ltisys._YT_complex
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  1540
# Total Exec:  220
# Func Size:   49
################################################################################

def _YT_complex(ker_pole, Q, transfer_matrix, i, j):
    """
    Applies algorithm from YT section 6.2 page 20 related to complex pairs
    """
    # step 1 page 20
    ur = np.sqrt(2)*Q[:, -2, np.newaxis]
    ui = np.sqrt(2)*Q[:, -1, np.newaxis]
    u = ur + 1j*ui

    # step 2 page 20
    ker_pole_ij = ker_pole[i]
    m = np.dot(np.dot(np.conj(ker_pole_ij.T), np.dot(u, np.conj(u).T) -
               np.dot(np.conj(u), u.T)), ker_pole_ij)

    # step 3 page 20
    e_val, e_vec = np.linalg.eig(m)
    # sort eigenvalues according to their module
    e_val_idx = np.argsort(np.abs(e_val))
    mu1 = e_vec[:, e_val_idx[-1], np.newaxis]
    mu2 = e_vec[:, e_val_idx[-2], np.newaxis]

    # what follows is a rough python translation of the formulas
    # in section 6.2 page 20 (step 4)

    # remember transfer_matrix_i has been split as
    # transfer_matrix[i]=real(transfer_matrix_i) and
    # transfer_matrix[j]=imag(transfer_matrix_i)
    transfer_matrix_j_mo_transfer_matrix_j = (
        transfer_matrix[:, i, np.newaxis] +
        1j*transfer_matrix[:, j, np.newaxis]
        )
    if not np.allclose(np.abs(e_val[e_val_idx[-1]]),
                              np.abs(e_val[e_val_idx[-2]])):
        ker_pole_mu = np.dot(ker_pole_ij, mu1)
    else:
        mu1_mu2_matrix = np.hstack((mu1, mu2))
        ker_pole_mu = np.dot(ker_pole_ij, mu1_mu2_matrix)
    transfer_matrix_i_j = np.dot(np.dot(ker_pole_mu, np.conj(ker_pole_mu.T)),
                              transfer_matrix_j_mo_transfer_matrix_j)

    if not np.allclose(transfer_matrix_i_j, 0):
        transfer_matrix_i_j = (transfer_matrix_i_j /
            np.linalg.norm(transfer_matrix_i_j))
        transfer_matrix[:, i] = np.real(transfer_matrix_i_j[:, 0])
        transfer_matrix[:, j] = np.imag(transfer_matrix_i_j[:, 0])
    else:
        # same idea as in YT_real
        transfer_matrix[:, i] = np.real(ker_pole_mu[:, 0])
        transfer_matrix[:, j] = np.imag(ker_pole_mu[:, 0])

################################################################################
# Method Path: scipy.stats._mstats_extras._hd_1D
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_mstats_extras.py
# Target:      dot
# Call Count:  1374
# Total Exec:  227
# Func Size:   29
################################################################################

def _hd_1D(data,prob,var):
        "Computes the HD quantiles for a 1D array. Returns nan for invalid data."
        xsorted = np.squeeze(np.sort(data.compressed().view(ndarray)))
        # Don't use length here, in case we have a numpy scalar
        n = xsorted.size

        hd = np.empty((2,len(prob)), float64)
        if n < 2:
            hd.flat = np.nan
            if var:
                return hd
            return hd[0]

        v = np.arange(n+1) / float(n)
        betacdf = beta.cdf
        for (i,p) in enumerate(prob):
            _w = betacdf(v, (n+1)*p, (n+1)*(1-p))
            w = _w[1:] - _w[:-1]
            hd_mean = np.dot(w, xsorted)
            hd[0,i] = hd_mean
            #
            hd[1,i] = np.dot(w, (xsorted-hd_mean)**2)
            #
        hd[0, prob == 0] = xsorted[0]
        hd[0, prob == 1] = xsorted[-1]
        if var:
            hd[1, prob == 0] = hd[1, prob == 1] = np.nan
            return hd
        return hd[0]

################################################################################
# Method Path: scipy.signal._ltisys.dlsim
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  1326
# Total Exec:  36
# Func Size:   118
################################################################################

def dlsim(system, u, t=None, x0=None):
    r"""Simulate output of a discrete-time linear system.

    Parameters
    ----------
    system : dlti | tuple
        An instance of the LTI class `dlti` or a tuple describing the system.
        The number of elements in the tuple determine the interpretation. I.e.:

        * ``system``: Instance of LTI class `dlti`. Note that derived instances, such
          as instances of `TransferFunction`, `ZerosPolesGain`, or `StateSpace`, are
          allowed as well.
        * ``(num, den, dt)``: Rational polynomial as described in `TransferFunction`.
          The coefficients of the polynomials should be specified in descending
          exponent order,  e.g., z² + 3z + 5 would be represented as ``[1, 3, 5]``.
        * ``(zeros, poles, gain, dt)``:  Zeros, poles, gain form as described
          in `ZerosPolesGain`.
        * ``(A, B, C, D, dt)``: State-space form as described in `StateSpace`.


    u : array_like
        An input array describing the input at each time `t` (interpolation is
        assumed between given times).  If there are multiple inputs, then each
        column of the rank-2 array represents an input.
    t : array_like, optional
        The time steps at which the input is defined.  If `t` is given, it
        must be the same length as `u`, and the final value in `t` determines
        the number of steps returned in the output.
    x0 : array_like, optional
        The initial conditions on the state vector (zero by default).

    Returns
    -------
    tout : ndarray
        Time values for the output, as a 1-D array.
    yout : ndarray
        System response, as a 1-D array.
    xout : ndarray, optional
        Time-evolution of the state-vector.  Only generated if the input is a
        `StateSpace` system.

    See Also
    --------
    lsim, dstep, dimpulse, cont2discrete

    Examples
    --------
    A simple integrator transfer function with a discrete time step of 1.0
    could be implemented as:

    >>> import numpy as np
    >>> from scipy import signal
    >>> tf = ([1.0,], [1.0, -1.0], 1.0)
    >>> t_in = [0.0, 1.0, 2.0, 3.0]
    >>> u = np.asarray([0.0, 0.0, 1.0, 1.0])
    >>> t_out, y = signal.dlsim(tf, u, t=t_in)
    >>> y.T
    array([[ 0.,  0.,  0.,  1.]])

    """
    # Convert system to dlti-StateSpace
    if isinstance(system, lti):
        raise AttributeError('dlsim can only be used with discrete-time dlti '
                             'systems.')
    elif not isinstance(system, dlti):
        system = dlti(*system[:-1], dt=system[-1])

    # Condition needed to ensure output remains compatible
    is_ss_input = isinstance(system, StateSpace)
    system = system._as_ss()

    u = np.atleast_1d(u)

    if u.ndim == 1:
        u = np.atleast_2d(u).T

    if t is None:
        out_samples = len(u)
        stoptime = (out_samples - 1) * system.dt
    else:
        stoptime = t[-1]
        out_samples = int(np.floor(stoptime / system.dt)) + 1

    # Pre-build output arrays
    xout = np.zeros((out_samples, system.A.shape[0]))
    yout = np.zeros((out_samples, system.C.shape[0]))
    tout = np.linspace(0.0, stoptime, num=out_samples)

    # Check initial condition
    if x0 is None:
        xout[0, :] = np.zeros((system.A.shape[1],))
    else:
        xout[0, :] = np.asarray(x0)

    # Pre-interpolate inputs into the desired time steps
    if t is None:
        u_dt = u
    else:
        if len(u.shape) == 1:
            u = u[:, np.newaxis]

        u_dt = make_interp_spline(t, u, k=1)(tout)

    # Simulate the system
    for i in range(0, out_samples - 1):
        xout[i+1, :] = (np.dot(system.A, xout[i, :]) +
                        np.dot(system.B, u_dt[i, :]))
        yout[i, :] = (np.dot(system.C, xout[i, :]) +
                      np.dot(system.D, u_dt[i, :]))

    # Last point
    yout[out_samples-1, :] = (np.dot(system.C, xout[out_samples-1, :]) +
                              np.dot(system.D, u_dt[out_samples-1, :]))

    if is_ss_input:
        return tout, yout, xout
    else:
        return tout, yout

################################################################################
# Method Path: scipy.linalg._solvers._solve_discrete_are
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_solvers.py
# Target:      dot
# Call Count:  1311
# Total Exec:  291
# Func Size:   88
################################################################################

def _solve_discrete_are(a, b, q, r, e, s, balanced):
    # Validate input arguments
    a, b, q, r, e, s, m, n, r_or_c, gen_are = _are_validate_args(
                                                     a, b, q, r, e, s, 'dare')

    # Form the matrix pencil
    H = np.zeros((2*m+n, 2*m+n), dtype=r_or_c)
    H[:m, :m] = a
    H[:m, 2*m:] = b
    H[m:2*m, :m] = -q
    H[m:2*m, m:2*m] = np.eye(m) if e is None else e.conj().T
    H[m:2*m, 2*m:] = 0. if s is None else -s
    H[2*m:, :m] = 0. if s is None else s.conj().T
    H[2*m:, 2*m:] = r

    J = np.zeros_like(H, dtype=r_or_c)
    J[:m, :m] = np.eye(m) if e is None else e
    J[m:2*m, m:2*m] = a.conj().T
    J[2*m:, m:2*m] = -b.conj().T

    if balanced:
        # xGEBAL does not remove the diagonals before scaling. Also
        # to avoid destroying the Symplectic structure, we follow Ref.3
        M = np.abs(H) + np.abs(J)
        np.fill_diagonal(M, 0.)
        _, (sca, _) = matrix_balance(M, separate=1, permute=0)
        # do we need to bother?
        if not np.allclose(sca, np.ones_like(sca)):
            # Now impose diag(D,inv(D)) from Benner where D is
            # square root of s_i/s_(n+i) for i=0,....
            sca = np.log2(sca)
            # NOTE: Py3 uses "Bankers Rounding: round to the nearest even" !!
            s = np.round((sca[m:2*m] - sca[:m])/2)
            sca = 2 ** np.r_[s, -s, sca[2*m:]]
            # Elementwise multiplication via broadcasting.
            elwisescale = sca[:, None] * np.reciprocal(sca)
            H *= elwisescale
            J *= elwisescale

    # Deflate the pencil by the R column ala Ref.1
    q_of_qr, _ = qr(H[:, -n:])
    H = q_of_qr[:, n:].conj().T.dot(H[:, :2*m])
    J = q_of_qr[:, n:].conj().T.dot(J[:, :2*m])

    # Decide on which output type is needed for QZ
    out_str = 'real' if r_or_c is float else 'complex'

    _, _, _, _, _, u = ordqz(H, J, sort='iuc',
                             overwrite_a=True,
                             overwrite_b=True,
                             check_finite=False,
                             output=out_str)

    # Get the relevant parts of the stable subspace basis
    if e is not None:
        u, _ = qr(np.vstack((e.dot(u[:m, :m]), u[m:, :m])))
    u00 = u[:m, :m]
    u10 = u[m:, :m]

    # Solve via back-substituion after checking the condition of u00
    up, ul, uu = lu(u00)

    if 1/cond(uu) < np.spacing(1.):
        raise LinAlgError('Failed to find a finite solution.')

    # Exploit the triangular structure
    x = solve_triangular(ul.conj().T,
                         solve_triangular(uu.conj().T,
                                          u10.conj().T,
                                          lower=True),
                         unit_diagonal=True,
                         ).conj().T.dot(up.conj().T)
    if balanced:
        x *= sca[:m, None] * sca[:m]

    # Check the deviation from symmetry for lack of success
    # See proof of Thm.5 item 3 in [2]
    u_sym = u00.conj().T.dot(u10)
    n_u_sym = norm(u_sym, 1)
    u_sym = u_sym - u_sym.conj().T
    sym_threshold = np.max([np.spacing(1000.), 0.1*n_u_sym])

    if norm(u_sym, 1) > sym_threshold:
        raise LinAlgError('The associated symplectic pencil has eigenvalues '
                          'too close to the unit circle')

    return (x + x.conj().T)/2

################################################################################
# Method Path: scipy.linalg._solvers._solve_continuous_are
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_solvers.py
# Target:      dot
# Call Count:  1283
# Total Exec:  284
# Func Size:   86
################################################################################

def _solve_continuous_are(a, b, q, r, e, s, balanced):
    # Validate input arguments
    a, b, q, r, e, s, m, n, r_or_c, gen_are = _are_validate_args(
                                                     a, b, q, r, e, s, 'care')

    H = np.empty((2*m+n, 2*m+n), dtype=r_or_c)
    H[:m, :m] = a
    H[:m, m:2*m] = 0.
    H[:m, 2*m:] = b
    H[m:2*m, :m] = -q
    H[m:2*m, m:2*m] = -a.conj().T
    H[m:2*m, 2*m:] = 0. if s is None else -s
    H[2*m:, :m] = 0. if s is None else s.conj().T
    H[2*m:, m:2*m] = b.conj().T
    H[2*m:, 2*m:] = r

    if gen_are and e is not None:
        J = block_diag(e, e.conj().T, np.zeros_like(r, dtype=r_or_c))
    else:
        J = block_diag(np.eye(2*m), np.zeros_like(r, dtype=r_or_c))

    if balanced:
        # xGEBAL does not remove the diagonals before scaling. Also
        # to avoid destroying the Symplectic structure, we follow Ref.3
        M = np.abs(H) + np.abs(J)
        np.fill_diagonal(M, 0.)
        _, (sca, _) = matrix_balance(M, separate=1, permute=0)
        # do we need to bother?
        if not np.allclose(sca, np.ones_like(sca)):
            # Now impose diag(D,inv(D)) from Benner where D is
            # square root of s_i/s_(n+i) for i=0,....
            sca = np.log2(sca)
            # NOTE: Py3 uses "Bankers Rounding: round to the nearest even" !!
            s = np.round((sca[m:2*m] - sca[:m])/2)
            sca = 2 ** np.r_[s, -s, sca[2*m:]]
            # Elementwise multiplication via broadcasting.
            elwisescale = sca[:, None] * np.reciprocal(sca)
            H *= elwisescale
            J *= elwisescale

    # Deflate the pencil to 2m x 2m ala Ref.1, eq.(55)
    q, r = qr(H[:, -n:])
    H = q[:, n:].conj().T.dot(H[:, :2*m])
    J = q[:2*m, n:].conj().T.dot(J[:2*m, :2*m])

    # Decide on which output type is needed for QZ
    out_str = 'real' if r_or_c is float else 'complex'

    _, _, _, _, _, u = ordqz(H, J, sort='lhp', overwrite_a=True,
                             overwrite_b=True, check_finite=False,
                             output=out_str)

    # Get the relevant parts of the stable subspace basis
    if e is not None:
        u, _ = qr(np.vstack((e.dot(u[:m, :m]), u[m:, :m])))
    u00 = u[:m, :m]
    u10 = u[m:, :m]

    # Solve via back-substituion after checking the condition of u00
    up, ul, uu = lu(u00)
    if 1/cond(uu) < np.spacing(1.):
        raise LinAlgError('Failed to find a finite solution.')

    # Exploit the triangular structure
    x = solve_triangular(ul.conj().T,
                         solve_triangular(uu.conj().T,
                                          u10.conj().T,
                                          lower=True),
                         unit_diagonal=True,
                         ).conj().T.dot(up.conj().T)
    if balanced:
        x *= sca[:m, None] * sca[:m]

    # Check the deviation from symmetry for lack of success
    # See proof of Thm.5 item 3 in [2]
    u_sym = u00.conj().T.dot(u10)
    n_u_sym = norm(u_sym, 1)
    u_sym = u_sym - u_sym.conj().T
    sym_threshold = np.max([np.spacing(1000.), 0.1*n_u_sym])

    if norm(u_sym, 1) > sym_threshold:
        raise LinAlgError('The associated Hamiltonian pencil has eigenvalues '
                          'too close to the imaginary axis')

    return (x + x.conj().T)/2

################################################################################
# Method Path: scipy.integrate._ivp.bdf.BDF._step_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/bdf.py
# Target:      dot
# Call Count:  1238
# Total Exec:  1224
# Func Size:   142
################################################################################

def _step_impl(self):
        t = self.t
        D = self.D

        max_step = self.max_step
        min_step = 10 * np.abs(np.nextafter(t, self.direction * np.inf) - t)
        if self.h_abs > max_step:
            h_abs = max_step
            change_D(D, self.order, max_step / self.h_abs)
            self.n_equal_steps = 0
        elif self.h_abs < min_step:
            h_abs = min_step
            change_D(D, self.order, min_step / self.h_abs)
            self.n_equal_steps = 0
        else:
            h_abs = self.h_abs

        atol = self.atol
        rtol = self.rtol
        order = self.order

        alpha = self.alpha
        gamma = self.gamma
        error_const = self.error_const

        J = self.J
        LU = self.LU
        current_jac = self.jac is None

        step_accepted = False
        while not step_accepted:
            if h_abs < min_step:
                return False, self.TOO_SMALL_STEP

            h = h_abs * self.direction
            t_new = t + h

            if self.direction * (t_new - self.t_bound) > 0:
                t_new = self.t_bound
                change_D(D, order, np.abs(t_new - t) / h_abs)
                self.n_equal_steps = 0
                LU = None

            h = t_new - t
            h_abs = np.abs(h)

            y_predict = np.sum(D[:order + 1], axis=0)

            scale = atol + rtol * np.abs(y_predict)
            psi = np.dot(D[1: order + 1].T, gamma[1: order + 1]) / alpha[order]

            converged = False
            c = h / alpha[order]
            while not converged:
                if LU is None:
                    LU = self.lu(self.I - c * J)

                converged, n_iter, y_new, d = solve_bdf_system(
                    self.fun, t_new, y_predict, c, psi, LU, self.solve_lu,
                    scale, self.newton_tol)

                if not converged:
                    if current_jac:
                        break
                    J = self.jac(t_new, y_predict)
                    LU = None
                    current_jac = True

            if not converged:
                factor = 0.5
                h_abs *= factor
                change_D(D, order, factor)
                self.n_equal_steps = 0
                LU = None
                continue

            safety = 0.9 * (2 * NEWTON_MAXITER + 1) / (2 * NEWTON_MAXITER
                                                       + n_iter)

            scale = atol + rtol * np.abs(y_new)
            error = error_const[order] * d
            error_norm = norm(error / scale)

            if error_norm > 1:
                factor = max(MIN_FACTOR,
                             safety * error_norm ** (-1 / (order + 1)))
                h_abs *= factor
                change_D(D, order, factor)
                self.n_equal_steps = 0
                # As we didn't have problems with convergence, we don't
                # reset LU here.
            else:
                step_accepted = True

        self.n_equal_steps += 1

        self.t = t_new
        self.y = y_new

        self.h_abs = h_abs
        self.J = J
        self.LU = LU

        # Update differences. The principal relation here is
        # D^{j + 1} y_n = D^{j} y_n - D^{j} y_{n - 1}. Keep in mind that D
        # contained difference for previous interpolating polynomial and
        # d = D^{k + 1} y_n. Thus this elegant code follows.
        D[order + 2] = d - D[order + 1]
        D[order + 1] = d
        for i in reversed(range(order + 1)):
            D[i] += D[i + 1]

        if self.n_equal_steps < order + 1:
            return True, None

        if order > 1:
            error_m = error_const[order - 1] * D[order]
            error_m_norm = norm(error_m / scale)
        else:
            error_m_norm = np.inf

        if order < MAX_ORDER:
            error_p = error_const[order + 1] * D[order + 2]
            error_p_norm = norm(error_p / scale)
        else:
            error_p_norm = np.inf

        error_norms = np.array([error_m_norm, error_norm, error_p_norm])
        with np.errstate(divide='ignore'):
            factors = error_norms ** (-1 / np.arange(order, order + 3))

        delta_order = np.argmax(factors) - 1
        order += delta_order
        self.order = order

        factor = min(MAX_FACTOR, safety * np.max(factors))
        self.h_abs *= factor
        change_D(D, order, factor)
        self.n_equal_steps = 0
        self.LU = None

        return True, None

################################################################################
# Method Path: scipy.sparse.linalg._funm_multiply_krylov.funm_multiply_krylov
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_funm_multiply_krylov.py
# Target:      dot
# Call Count:  1233
# Total Exec:  502
# Func Size:   192
################################################################################

def funm_multiply_krylov(f, A, b, *, assume_a = "general", t = 1.0, atol = 0.0,
                         rtol = 1e-6, restart_every_m = None, max_restarts = 20):
    """
    A restarted Krylov method for evaluating ``y = f(tA) b`` from [1]_ [2]_.

    Parameters
    ----------
    f : callable
        Callable object that computes the matrix function ``F = f(X)``.
    A : {sparse array, ndarray, LinearOperator}
        A real or complex N-by-N matrix.
        Alternatively, `A` can be a linear operator which can
        produce ``Ax`` using, e.g., ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        A vector to multiply the ``f(tA)`` with.
    assume_a : str, optional
        Indicate the structure of ``A``. The algorithm will use this information
        to select the appropriated code path. The available options are
        'hermitian'/'her' and 'general'/'gen'. If ommited, then it is assumed
        that ``A`` has a 'general' structure.
    t : float, optional
        The value to scale the matrix ``A`` with. The default is ``t = 1.0``
    atol, rtol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(||y_k - y_k-1||) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``atol=0.`` and ``rtol=1e-6``.
    restart_every_m : int
        If the iteration number reaches this value a restart is triggered.
        Larger values increase iteration cost but may be necessary for convergence.
        If omitted, ``min(20, n)`` is used.
    max_restarts : int, optional
        Maximum number of restart cycles. The algorithm will stop
        after max_restarts cycles even if the specified tolerance has not been
        achieved. The default is ``max_restarts=20``

    Returns
    -------
    y : ndarray
        The result of ``f(tA) b``.

    Notes
    -----
    The convergence of the Krylov method heavily depends on the spectrum
    of ``A`` and the function ``f``. With restarting, there are only formal
    proofs for functions of order 1 (e.g., ``exp``, ``sin``, ``cos``) and
    Stieltjes functions [2]_ [3]_, while the general case remains an open problem.

    References
    ----------
    .. [1] M. Afanasjew, M. Eiermann, O. G. Ernst, and S. Güttel,
           "Implementation of a restarted Krylov subspace method for the
           evaluation of matrix functions," Linear Algebra and its Applications,
           vol. 429, no. 10, pp. 2293-2314, Nov. 2008, :doi:`10.1016/j.laa.2008.06.029`.

    .. [2] M. Eiermann and O. G. Ernst, "A Restarted Krylov Subspace Method
           for the Evaluation of Matrix Functions," SIAM J. Numer. Anal., vol. 44,
           no. 6, pp. 2481-2504, Jan. 2006, :doi:`10.1137/050633846`.

    .. [3] A. Frommer, S. Güttel, and M. Schweitzer, "Convergence of Restarted
           Krylov Subspace Methods for Stieltjes Functions of Matrices," SIAM J.
           Matrix Anal. Appl., vol. 35, no. 4, pp. 1602-1624,
           Jan. 2014, :doi:`10.1137/140973463`.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csr_array
    >>> from scipy.sparse.linalg import funm_multiply_krylov
    >>> from scipy.linalg import expm, solve
    >>> A = csr_array([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> t = 0.1

    Compute ``y = exp(tA) b``.

    >>> y = funm_multiply_krylov(expm, A, b, t = t)
    >>> y
    array([3.6164913 , 3.88421511 , 0.96073457])

    >>> ref = expm(t * A.todense()) @ b
    >>> err = y - ref
    >>> err
    array([4.44089210e-16 , 0.00000000e+00 , 2.22044605e-16])

    Compute :math:`y = (A^3 - A) b`.

    >>> poly = lambda X : X @ X @ X - X
    >>> y = funm_multiply_krylov(poly, A, b)
    >>> y
    array([132. , 24. , 70.])

    >>> ref = poly(A.todense()) @ b
    >>> err = y - ref
    >>> err
    array([ 0.00000000e+00 , 7.10542736e-15 , -2.84217094e-14])

    Compute :math:`y = f(tA) b`, where  :math:`f(X) = X^{-1}(e^{X} - I)`. This is
    known as the "phi function" from the exponential integrator literature.

    >>> phim_1 = lambda X : solve(X, expm(X) - np.eye(X.shape[0]))
    >>> y = funm_multiply_krylov(phim_1, A, b, t = t)
    >>> y
    array([ 2.76984306 , 3.92769192 , -0.03111392])

    >>> ref = phim_1(t * A.todense()) @ b
    >>> err = y - ref
    >>> err
    array([ 0.00000000e+00 , 8.88178420e-16 , -4.60742555e-15])
    """

    if assume_a not in {'hermitian', 'general', 'her', 'gen'}:
        raise ValueError(f'scipy.sparse.linalg.funm_multiply_krylov: {assume_a} '
                         'is not a recognized matrix structure')
    is_hermitian = (assume_a == 'her') or (assume_a == 'hermitian')

    if len(b.shape) != 1:
        raise ValueError("scipy.sparse.linalg.funm_multiply_krylov: "
                         "argument 'b' must be a 1D array.")
    n = b.shape[0]

    if restart_every_m is None:
        restart_every_m = min(20, n)

    restart_every_m = int(restart_every_m)
    max_restarts = int(max_restarts)

    if restart_every_m <= 0:
        raise ValueError("scipy.sparse.linalg.funm_multiply_krylov: "
                         "argument 'restart_every_m' must be positive.")

    if max_restarts <= 0:
            raise ValueError("scipy.sparse.linalg.funm_multiply_krylov: "
                             "argument 'max_restarts' must be positive.")

    m = restart_every_m
    max_restarts = min(max_restarts, int(n / m) + 1)
    mmax = m * max_restarts

    bnorm = norm(b)
    atol, _ = _get_atol_rtol("funm_multiply_krylov", bnorm, atol, rtol)

    if bnorm == 0:
        y = np.array(b)
        return y

    # Preallocate the maximum memory space.
    # Using the column major order here since we work with
    # each individual column separately.
    internal_type = np.common_type(A, b)
    V = np.zeros((n, m + 1), dtype = internal_type, order = 'F')
    H = np.zeros((mmax + 1, mmax), dtype = internal_type, order = 'F')

    restart = 1

    if is_hermitian:
        breakdown, j = _funm_multiply_krylov_lanczos(A, b, bnorm, V,
                                                        H[:m + 1, :m], m)
    else:
        breakdown, j = _funm_multiply_krylov_arnoldi(A, b, bnorm, V,
                                                        H[:m + 1, :m], m)

    fH = f(t * H[:j, :j])
    y = bnorm * V[:, :j].dot(fH[:, 0])

    if breakdown:
        return y

    update_norm = norm(bnorm * fH[:, 0])

    while restart < max_restarts and update_norm > atol:
        begin = restart * m
        end = (restart + 1) * m

        if is_hermitian:
            breakdown, j = _funm_multiply_krylov_lanczos(A, V[:, m], 1, V,
                                                         H[begin:end + 1, begin:end], m)
        else:
            breakdown, j = _funm_multiply_krylov_arnoldi(A, V[:, m], 1, V,
                                                         H[begin:end + 1, begin:end], m)

        if breakdown:
            end = begin + j
            fH = f(t * H[:end, :end])
            y[:end] = y[:end] + bnorm * V[:, :m].dot(fH[begin:end, 0])
            return y

        fH = f(t * H[:end, :end])
        y = y + bnorm * V[:, :m].dot(fH[begin:end, 0])
        update_norm = norm(bnorm * fH[begin:end, 0])
        restart += 1

    return y

################################################################################
# Method Path: scipy.integrate._ivp.rk.RkDenseOutput._call_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  1195
# Total Exec:  1195
# Func Size:   15
################################################################################

def _call_impl(self, t):
        x = (t - self.t_old) / self.h
        if t.ndim == 0:
            p = np.tile(x, self.order + 1)
            p = np.cumprod(p)
        else:
            p = np.tile(x, (self.order + 1, 1))
            p = np.cumprod(p, axis=0)
        y = self.h * np.dot(self.Q, p)
        if y.ndim == 2:
            y += self.y_old[:, None]
        else:
            y += self.y_old

        return y

################################################################################
# Method Path: scipy.linalg._decomp_schur.rsf2csf
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_decomp_schur.py
# Target:      dot
# Call Count:  1125
# Total Exec:  599
# Func Size:   85
################################################################################

def rsf2csf(T, Z, check_finite=True):
    """
    Convert real Schur form to complex Schur form.

    Convert a quasi-diagonal real-valued Schur form to the upper-triangular
    complex-valued Schur form.

    Parameters
    ----------
    T : (M, M) array_like
        Real Schur form of the original array
    Z : (M, M) array_like
        Schur transformation matrix
    check_finite : bool, optional
        Whether to check that the input arrays contain only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.

    Returns
    -------
    T : (M, M) ndarray
        Complex Schur form of the original array
    Z : (M, M) ndarray
        Schur transformation matrix corresponding to the complex form

    See Also
    --------
    schur : Schur decomposition of an array

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import schur, rsf2csf
    >>> A = np.array([[0, 2, 2], [0, 1, 2], [1, 0, 1]])
    >>> T, Z = schur(A)
    >>> T
    array([[ 2.65896708,  1.42440458, -1.92933439],
           [ 0.        , -0.32948354, -0.49063704],
           [ 0.        ,  1.31178921, -0.32948354]])
    >>> Z
    array([[0.72711591, -0.60156188, 0.33079564],
           [0.52839428, 0.79801892, 0.28976765],
           [0.43829436, 0.03590414, -0.89811411]])
    >>> T2 , Z2 = rsf2csf(T, Z)
    >>> T2
    array([[2.65896708+0.j, -1.64592781+0.743164187j, -1.21516887+1.00660462j],
           [0.+0.j , -0.32948354+8.02254558e-01j, -0.82115218-2.77555756e-17j],
           [0.+0.j , 0.+0.j, -0.32948354-0.802254558j]])
    >>> Z2
    array([[0.72711591+0.j,  0.28220393-0.31385693j,  0.51319638-0.17258824j],
           [0.52839428+0.j,  0.24720268+0.41635578j, -0.68079517-0.15118243j],
           [0.43829436+0.j, -0.76618703+0.01873251j, -0.03063006+0.46857912j]])

    """
    if check_finite:
        Z, T = map(asarray_chkfinite, (Z, T))
    else:
        Z, T = map(asarray, (Z, T))

    for ind, X in enumerate([Z, T]):
        if X.ndim != 2 or X.shape[0] != X.shape[1]:
            raise ValueError(f"Input '{'ZT'[ind]}' must be square.")

    if T.shape[0] != Z.shape[0]:
        message = f"Input array shapes must match: Z: {Z.shape} vs. T: {T.shape}"
        raise ValueError(message)
    N = T.shape[0]
    t = _commonType(Z, T, array([3.0], 'F'))
    Z, T = _castCopy(t, Z, T)

    for m in range(N-1, 0, -1):
        if abs(T[m, m-1]) > eps*(abs(T[m-1, m-1]) + abs(T[m, m])):
            mu = eigvals(T[m-1:m+1, m-1:m+1]) - T[m, m]
            r = norm([mu[0], T[m, m-1]])
            c = mu[0] / r
            s = T[m, m-1] / r
            G = array([[c.conj(), s], [-s, c]], dtype=t)

            T[m-1:m+1, m-1:] = G.dot(T[m-1:m+1, m-1:])
            T[:m+1, m-1:m+1] = T[:m+1, m-1:m+1].dot(G.conj().T)
            Z[:, m-1:m+1] = Z[:, m-1:m+1].dot(G.conj().T)

        T[m, m-1] = 0.0
    return T, Z

################################################################################
# Method Path: scipy.optimize._trustregion.DoglegSubproblem.hessp
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion.py
# Target:      dot
# Call Count:  1055
# Total Exec:  2299
# Func Size:   5
################################################################################

def hessp(self, p):
        if self._hessp is not None:
            return self._hessp(self._x, p)
        else:
            return np.dot(self.hess, p)

################################################################################
# Method Path: scipy.optimize._linesearch._cubicmin
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linesearch.py
# Target:      dot
# Call Count:  1005
# Total Exec:  1005
# Func Size:   32
################################################################################

def _cubicmin(a, fa, fpa, b, fb, c, fc):
    """
    Finds the minimizer for a cubic polynomial that goes through the
    points (a,fa), (b,fb), and (c,fc) with derivative at a of fpa.

    If no minimizer can be found, return None.

    """
    # f(x) = A *(x-a)^3 + B*(x-a)^2 + C*(x-a) + D

    with np.errstate(divide='raise', over='raise', invalid='raise'):
        try:
            C = fpa
            db = b - a
            dc = c - a
            denom = (db * dc) ** 2 * (db - dc)
            d1 = np.empty((2, 2))
            d1[0, 0] = dc ** 2
            d1[0, 1] = -db ** 2
            d1[1, 0] = -dc ** 3
            d1[1, 1] = db ** 3
            [A, B] = np.dot(d1, np.asarray([fb - fa - C * db,
                                            fc - fa - C * dc]).flatten())
            A /= denom
            B /= denom
            radical = B * B - 3 * A * C
            xmin = a + (-B + np.sqrt(radical)) / (3 * A)
        except ArithmeticError:
            return None
    if not np.isfinite(xmin):
        return None
    return xmin

################################################################################
# Method Path: scipy.optimize._hessian_update_strategy.BFGS._auto_scale
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_hessian_update_strategy.py
# Target:      dot
# Call Count:  975
# Total Exec:  325
# Func Size:   13
################################################################################

def _auto_scale(self, delta_x, delta_grad):
        # Heuristic to scale matrix at first iteration.
        # Described in Nocedal and Wright "Numerical Optimization"
        # p.143 formula (6.20).
        s_norm2 = np.dot(delta_x, delta_x)
        y_norm2 = np.dot(delta_grad, delta_grad)
        ys = np.abs(np.dot(delta_grad, delta_x))
        if ys == 0.0 or y_norm2 == 0 or s_norm2 == 0:
            return 1
        if self.approx_type == 'hess':
            return y_norm2 / ys
        else:
            return ys / y_norm2

################################################################################
# Method Path: scipy.optimize._trustregion.hessp
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion.py
# Target:      dot
# Call Count:  912
# Total Exec:  912
# Func Size:   2
################################################################################

def hessp(self, p):
        if self._hessp is not None:
            return self._hessp(self._x, p)
        else:
            return np.dot(self.hess, p)

################################################################################
# Method Path: scipy.linalg._matfuncs_inv_ssq._logm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs_inv_ssq.py
# Target:      dot
# Call Count:  908
# Total Exec:  518
# Func Size:   47
################################################################################

def _logm(A):
    """
    Compute the matrix logarithm.

    See the logm docstring in matfuncs.py for more info.

    Notes
    -----
    In this function we look at triangular matrices that are similar
    to the input matrix. If any diagonal entry of such a triangular matrix
    is exactly zero then the original matrix is singular.
    The matrix logarithm does not exist for such matrices,
    but in such cases we will pretend that the diagonal entries that are zero
    are actually slightly positive by an ad-hoc amount, in the interest
    of returning something more useful than NaN. This will cause a warning.

    """
    A = np.asarray(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected a square matrix')

    # If the input matrix dtype is integer then copy to a float dtype matrix.
    if issubclass(A.dtype.type, np.integer):
        A = np.asarray(A, dtype=float)

    keep_it_real = np.isrealobj(A)
    try:
        if np.array_equal(A, np.triu(A)):
            A = _logm_force_nonsingular_triangular_matrix(A)
            if np.min(np.diag(A), initial=0.) < 0:
                A = A.astype(complex)
            return _logm_triu(A)
        else:
            if keep_it_real:
                T, Z = schur(A)
                if not np.array_equal(T, np.triu(T)):
                    T, Z = rsf2csf(T, Z)
            else:
                T, Z = schur(A, output='complex')
            T = _logm_force_nonsingular_triangular_matrix(T, inplace=True)
            U = _logm_triu(T)
            ZH = np.conjugate(Z).T
            return Z.dot(U).dot(ZH)
    except (SqrtmError, LogmError):
        X = np.empty_like(A)
        X.fill(np.nan)
        return X

################################################################################
# Method Path: scipy.sparse.linalg._isolve.lgmres.lgmres
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_isolve/lgmres.py
# Target:      dot
# Call Count:  892
# Total Exec:  437
# Func Size:   217
################################################################################

def lgmres(A, b, x0=None, *, rtol=1e-5, atol=0., maxiter=1000, M=None, callback=None,
           inner_m=30, outer_k=3, outer_v=None, store_outer_Av=True,
           prepend_outer_v=False):
    """
    Solve ``Ax = b`` with the LGMRES algorithm.

    The LGMRES algorithm [1]_ [2]_ is designed to avoid some problems
    in the convergence in restarted GMRES, and often converges in fewer
    iterations.

    Parameters
    ----------
    A : {sparse array, ndarray, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
        Alternatively, ``A`` can be a linear operator which can
        produce ``Ax`` using, e.g.,
        ``scipy.sparse.linalg.LinearOperator``.
    b : ndarray
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0 : ndarray
        Starting guess for the solution.
    rtol, atol : float, optional
        Parameters for the convergence test. For convergence,
        ``norm(b - A @ x) <= max(rtol*norm(b), atol)`` should be satisfied.
        The default is ``rtol=1e-5``, the default for ``atol`` is ``0.0``.
    maxiter : int, optional
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse array, ndarray, LinearOperator}, optional
        Preconditioner for A.  The preconditioner should approximate the
        inverse of A.  Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function, optional
        User-supplied function to call after each iteration.  It is called
        as callback(xk), where xk is the current solution vector.
    inner_m : int, optional
        Number of inner GMRES iterations per each outer iteration.
    outer_k : int, optional
        Number of vectors to carry between inner GMRES iterations.
        According to [1]_, good values are in the range of 1...3.
        However, note that if you want to use the additional vectors to
        accelerate solving multiple similar problems, larger values may
        be beneficial.
    outer_v : list of tuples, optional
        List containing tuples ``(v, Av)`` of vectors and corresponding
        matrix-vector products, used to augment the Krylov subspace, and
        carried between inner GMRES iterations. The element ``Av`` can
        be `None` if the matrix-vector product should be re-evaluated.
        This parameter is modified in-place by `lgmres`, and can be used
        to pass "guess" vectors in and out of the algorithm when solving
        similar problems.
    store_outer_Av : bool, optional
        Whether LGMRES should store also A@v in addition to vectors `v`
        in the `outer_v` list. Default is True.
    prepend_outer_v : bool, optional
        Whether to put outer_v augmentation vectors before Krylov iterates.
        In standard LGMRES, prepend_outer_v=False.

    Returns
    -------
    x : ndarray
        The converged solution.
    info : int
        Provides convergence information:

            - 0  : successful exit
            - >0 : convergence to tolerance not achieved, number of iterations
            - <0 : illegal input or breakdown

    Notes
    -----
    The LGMRES algorithm [1]_ [2]_ is designed to avoid the
    slowing of convergence in restarted GMRES, due to alternating
    residual vectors. Typically, it often outperforms GMRES(m) of
    comparable memory requirements by some measure, or at least is not
    much worse.

    Another advantage in this algorithm is that you can supply it with
    'guess' vectors in the `outer_v` argument that augment the Krylov
    subspace. If the solution lies close to the span of these vectors,
    the algorithm converges faster. This can be useful if several very
    similar matrices need to be inverted one after another, such as in
    Newton-Krylov iteration where the Jacobian matrix often changes
    little in the nonlinear steps.

    References
    ----------
    .. [1] A.H. Baker and E.R. Jessup and T. Manteuffel, "A Technique for
             Accelerating the Convergence of Restarted GMRES", SIAM J. Matrix
             Anal. Appl. 26, 962 (2005).
    .. [2] A.H. Baker, "On Improving the Performance of the Linear Solver
             restarted GMRES", PhD thesis, University of Colorado (2003).

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.sparse import csc_array
    >>> from scipy.sparse.linalg import lgmres
    >>> A = csc_array([[3, 2, 0], [1, -1, 0], [0, 5, 1]], dtype=float)
    >>> b = np.array([2, 4, -1], dtype=float)
    >>> x, exitCode = lgmres(A, b, atol=1e-5)
    >>> print(exitCode)            # 0 indicates successful convergence
    0
    >>> np.allclose(A.dot(x), b)
    True
    """
    A,M,x,b = make_system(A,M,x0,b)

    if not np.isfinite(b).all():
        raise ValueError("RHS must contain only finite numbers")

    matvec = A.matvec
    psolve = M.matvec

    if outer_v is None:
        outer_v = []

    axpy, dot, scal = None, None, None
    nrm2 = get_blas_funcs('nrm2', [b])

    b_norm = nrm2(b)

    # we call this to get the right atol/rtol and raise errors as necessary
    atol, rtol = _get_atol_rtol('lgmres', b_norm, atol, rtol)

    if b_norm == 0:
        x = b
        return (x, 0)

    ptol_max_factor = 1.0

    for k_outer in range(maxiter):
        r_outer = matvec(x) - b

        # -- callback
        if callback is not None:
            callback(x)

        # -- determine input type routines
        if axpy is None:
            if np.iscomplexobj(r_outer) and not np.iscomplexobj(x):
                x = x.astype(r_outer.dtype)
            axpy, dot, scal, nrm2 = get_blas_funcs(['axpy', 'dot', 'scal', 'nrm2'],
                                                   (x, r_outer))

        # -- check stopping condition
        r_norm = nrm2(r_outer)
        if r_norm <= max(atol, rtol * b_norm):
            break

        # -- inner LGMRES iteration
        v0 = -psolve(r_outer)
        inner_res_0 = nrm2(v0)

        if inner_res_0 == 0:
            rnorm = nrm2(r_outer)
            raise RuntimeError("Preconditioner returned a zero vector; "
                               f"|v| ~ {rnorm:.1g}, |M v| = 0")

        v0 = scal(1.0/inner_res_0, v0)

        ptol = min(ptol_max_factor, max(atol, rtol*b_norm)/r_norm)

        try:
            Q, R, B, vs, zs, y, pres = _fgmres(matvec,
                                               v0,
                                               inner_m,
                                               lpsolve=psolve,
                                               atol=ptol,
                                               outer_v=outer_v,
                                               prepend_outer_v=prepend_outer_v)
            y *= inner_res_0
            if not np.isfinite(y).all():
                # Overflow etc. in computation. There's no way to
                # recover from this, so we have to bail out.
                raise LinAlgError()
        except LinAlgError:
            # Floating point over/underflow, non-finite result from
            # matmul etc. -- report failure.
            return x, k_outer + 1

        # Inner loop tolerance control
        if pres > ptol:
            ptol_max_factor = min(1.0, 1.5 * ptol_max_factor)
        else:
            ptol_max_factor = max(1e-16, 0.25 * ptol_max_factor)

        # -- GMRES terminated: eval solution
        dx = zs[0]*y[0]
        for w, yc in zip(zs[1:], y[1:]):
            dx = axpy(w, dx, dx.shape[0], yc)  # dx += w*yc

        # -- Store LGMRES augmentation vectors
        nx = nrm2(dx)
        if nx > 0:
            if store_outer_Av:
                q = Q.dot(R.dot(y))
                ax = vs[0]*q[0]
                for v, qc in zip(vs[1:], q[1:]):
                    ax = axpy(v, ax, ax.shape[0], qc)
                outer_v.append((dx/nx, ax/nx))
            else:
                outer_v.append((dx/nx, None))

        # -- Retain only a finite number of augmentation vectors
        while len(outer_v) > outer_k:
            del outer_v[0]

        # -- Apply step
        x += dx
    else:
        # didn't converge ...
        return x, maxiter

    return x, 0

################################################################################
# Method Path: scipy.linalg._matfuncs_inv_ssq._remainder_matrix_power
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs_inv_ssq.py
# Target:      dot
# Call Count:  890
# Total Exec:  560
# Func Size:   73
################################################################################

def _remainder_matrix_power(A, t):
    """
    Compute the fractional power of a matrix, for fractions -1 < t < 1.

    This uses algorithm (3.1) of [1]_.
    The Pade approximation itself uses algorithm (4.1) of [2]_.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix whose fractional power to evaluate.
    t : float
        Fractional power between -1 and 1 exclusive.

    Returns
    -------
    X : (N, N) array_like
        The fractional power of the matrix.

    References
    ----------
    .. [1] Nicholas J. Higham and Lijing Lin (2013)
           "An Improved Schur-Pade Algorithm for Fractional Powers
           of a Matrix and their Frechet Derivatives."

    .. [2] Nicholas J. Higham and Lijing lin (2011)
           "A Schur-Pade Algorithm for Fractional Powers of a Matrix."
           SIAM Journal on Matrix Analysis and Applications,
           32 (3). pp. 1056-1078. ISSN 0895-4798

    """
    # This code block is copied from numpy.matrix_power().
    A = np.asarray(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('input must be a square array')

    # Get the number of rows and columns.
    n, n = A.shape

    # Triangularize the matrix if necessary,
    # attempting to preserve dtype if possible.
    if np.array_equal(A, np.triu(A)):
        Z = None
        T = A
    else:
        if np.isrealobj(A):
            T, Z = schur(A)
            if not np.array_equal(T, np.triu(T)):
                T, Z = rsf2csf(T, Z)
        else:
            T, Z = schur(A, output='complex')

    # Zeros on the diagonal of the triangular matrix are forbidden,
    # because the inverse scaling and squaring cannot deal with it.
    T_diag = np.diag(T)
    if np.count_nonzero(T_diag) != n:
        raise FractionalMatrixPowerError(
                'cannot use inverse scaling and squaring to find '
                'the fractional matrix power of a singular matrix')

    # If the triangular matrix is real and has a negative
    # entry on the diagonal, then force the matrix to be complex.
    if np.isrealobj(T) and np.min(T_diag) < 0:
        T = T.astype(complex)

    # Get the fractional power of the triangular matrix,
    # and de-triangularize it if necessary.
    U = _remainder_matrix_power_triu(T, t)
    if Z is not None:
        ZH = np.conjugate(Z).T
        return Z.dot(U).dot(ZH)
    else:
        return U

################################################################################
# Method Path: scipy.sparse.linalg._expm_multiply.traceest
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_expm_multiply.py
# Target:      trace
# Call Count:  862
# Total Exec:  431
# Func Size:   45
################################################################################

def traceest(A, m3, seed=None):
    """Estimate `np.trace(A)` using `3*m3` matrix-vector products.

    The result is not deterministic.

    Parameters
    ----------
    A : LinearOperator
        Linear operator whose trace will be estimated. Has to be square.
    m3 : int
        Number of matrix-vector products divided by 3 used to estimate the
        trace.
    seed : optional
        Seed for `numpy.random.default_rng`.
        Can be provided to obtain deterministic results.

    Returns
    -------
    trace : LinearOperator.dtype
        Estimate of the trace

    Notes
    -----
    This is the Hutch++ algorithm given in [1]_.

    References
    ----------
    .. [1] Meyer, Raphael A., Cameron Musco, Christopher Musco, and David P.
       Woodruff. "Hutch++: Optimal Stochastic Trace Estimation." In Symposium
       on Simplicity in Algorithms (SOSA), pp. 142-155. Society for Industrial
       and Applied Mathematics, 2021.
       :doi:`10.1137/1.9781611976496.16`.

    """
    rng = np.random.default_rng(seed)
    if len(A.shape) != 2 or A.shape[-1] != A.shape[-2]:
        raise ValueError("Expected A to be like a square matrix.")
    n = A.shape[-1]
    S = rng.choice([-1.0, +1.0], [n, m3])
    Q, _ = qr(A.matmat(S), overwrite_a=True, mode='economic')
    trQAQ = np.trace(Q.conj().T @ A.matmat(Q))
    G = rng.choice([-1, +1], [n, m3])
    right = G - Q@(Q.conj().T @ G)
    trGAG = np.trace(right.conj().T @ A.matmat(right))
    return trQAQ + trGAG/m3

################################################################################
# Method Path: scipy.optimize._lsq.least_squares.least_squares
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/least_squares.py
# Target:      dot
# Call Count:  790
# Total Exec:  993
# Func Size:   779
################################################################################

def least_squares(
        fun, x0, jac='2-point', bounds=(-np.inf, np.inf), method='trf',
        ftol=1e-8, xtol=1e-8, gtol=1e-8, x_scale=None, loss='linear',
        f_scale=1.0, diff_step=None, tr_solver=None, tr_options=None,
        jac_sparsity=None, max_nfev=None, verbose=0, args=(), kwargs=None,
        callback=None, workers=None
):
    """Solve a nonlinear least-squares problem with bounds on the variables.

    Given the residuals f(x) (an m-D real function of n real
    variables) and the loss function rho(s) (a scalar function), `least_squares`
    finds a local minimum of the cost function F(x)::

        minimize F(x) = 0.5 * sum(rho(f_i(x)**2), i = 0, ..., m - 1)
        subject to lb <= x <= ub

    The purpose of the loss function rho(s) is to reduce the influence of
    outliers on the solution.

    Parameters
    ----------
    fun : callable
        Function which computes the vector of residuals, with the signature
        ``fun(x, *args, **kwargs)``, i.e., the minimization proceeds with
        respect to its first argument. The argument ``x`` passed to this
        function is an ndarray of shape (n,) (never a scalar, even for n=1).
        It must allocate and return a 1-D array_like of shape (m,) or a scalar.
        If the argument ``x`` is complex or the function ``fun`` returns
        complex residuals, it must be wrapped in a real function of real
        arguments, as shown at the end of the Examples section.
    x0 : array_like with shape (n,) or float
        Initial guess on independent variables. If float, it will be treated
        as a 1-D array with one element. When `method` is 'trf', the initial
        guess might be slightly adjusted to lie sufficiently within the given
        `bounds`.
    jac : {'2-point', '3-point', 'cs', callable}, optional
        Method of computing the Jacobian matrix (an m-by-n matrix, where
        element (i, j) is the partial derivative of f[i] with respect to
        x[j]). The keywords select a finite difference scheme for numerical
        estimation. The scheme '3-point' is more accurate, but requires
        twice as many operations as '2-point' (default). The scheme 'cs'
        uses complex steps, and while potentially the most accurate, it is
        applicable only when `fun` correctly handles complex inputs and
        can be analytically continued to the complex plane. If callable, it is used as
        ``jac(x, *args, **kwargs)`` and should return a good approximation
        (or the exact value) for the Jacobian as an array_like (np.atleast_2d
        is applied), a sparse array (csr_array preferred for performance) or
        a `scipy.sparse.linalg.LinearOperator`.

        .. versionchanged:: 1.16.0
            An ability to use the '3-point', 'cs' keywords with the 'lm' method.
            Previously 'lm' was limited to '2-point' and callable.

    bounds : 2-tuple of array_like or `Bounds`, optional
        There are two ways to specify bounds:

        1. Instance of `Bounds` class
        2. Lower and upper bounds on independent variables. Defaults to no
           bounds. Each array must match the size of `x0` or be a scalar,
           in the latter case a bound will be the same for all variables.
           Use ``np.inf`` with an appropriate sign to disable bounds on all
           or some variables.

    method : {'trf', 'dogbox', 'lm'}, optional
        Algorithm to perform minimization.

        * 'trf' : Trust Region Reflective algorithm, particularly suitable
          for large sparse problems with bounds. Generally robust method.
        * 'dogbox' : dogleg algorithm with rectangular trust regions,
          typical use case is small problems with bounds. Not recommended
          for problems with rank-deficient Jacobian.
        * 'lm' : Levenberg-Marquardt algorithm as implemented in MINPACK.
          Doesn't handle bounds and sparse Jacobians. Usually the most
          efficient method for small unconstrained problems.

        Default is 'trf'. See Notes for more information.
    ftol : float or None, optional
        Tolerance for termination by the change of the cost function. Default
        is 1e-8. The optimization process is stopped when ``dF < ftol * F``,
        and there was an adequate agreement between a local quadratic model and
        the true model in the last step.

        If None and 'method' is not 'lm', the termination by this condition is
        disabled. If 'method' is 'lm', this tolerance must be higher than
        machine epsilon.
    xtol : float or None, optional
        Tolerance for termination by the change of the independent variables.
        Default is 1e-8. The exact condition depends on the `method` used:

        * For 'trf' and 'dogbox' : ``norm(dx) < xtol * (xtol + norm(x))``.
        * For 'lm' : ``Delta < xtol * norm(xs)``, where ``Delta`` is
          a trust-region radius and ``xs`` is the value of ``x``
          scaled according to `x_scale` parameter (see below).

        If None and 'method' is not 'lm', the termination by this condition is
        disabled. If 'method' is 'lm', this tolerance must be higher than
        machine epsilon.
    gtol : float or None, optional
        Tolerance for termination by the norm of the gradient. Default is 1e-8.
        The exact condition depends on a `method` used:

        * For 'trf' : ``norm(g_scaled, ord=np.inf) < gtol``, where
          ``g_scaled`` is the value of the gradient scaled to account for
          the presence of the bounds [STIR]_.
        * For 'dogbox' : ``norm(g_free, ord=np.inf) < gtol``, where
          ``g_free`` is the gradient with respect to the variables which
          are not in the optimal state on the boundary.
        * For 'lm' : the maximum absolute value of the cosine of angles
          between columns of the Jacobian and the residual vector is less
          than `gtol`, or the residual vector is zero.

        If None and 'method' is not 'lm', the termination by this condition is
        disabled. If 'method' is 'lm', this tolerance must be higher than
        machine epsilon.
    x_scale : {None, array_like, 'jac'}, optional
        Characteristic scale of each variable. Setting `x_scale` is equivalent
        to reformulating the problem in scaled variables ``xs = x / x_scale``.
        An alternative view is that the size of a trust region along jth
        dimension is proportional to ``x_scale[j]``. Improved convergence may
        be achieved by setting `x_scale` such that a step of a given size
        along any of the scaled variables has a similar effect on the cost
        function. If set to 'jac', the scale is iteratively updated using the
        inverse norms of the columns of the Jacobian matrix (as described in
        [JJMore]_). The default scaling for each method (i.e.
        if ``x_scale is None``) is as follows:

        * For 'trf'    : ``x_scale == 1``
        * For 'dogbox' : ``x_scale == 1``
        * For 'lm'     : ``x_scale == 'jac'``

        .. versionchanged:: 1.16.0
            The default keyword value is changed from 1 to None to indicate that
            a default approach to scaling is used.
            For the 'lm' method the default scaling is changed from 1 to 'jac'.
            This has been found to give better performance, and is the same
            scaling as performed by ``leastsq``.

    loss : str or callable, optional
        Determines the loss function. The following keyword values are allowed:

        * 'linear' (default) : ``rho(z) = z``. Gives a standard
          least-squares problem.
        * 'soft_l1' : ``rho(z) = 2 * ((1 + z)**0.5 - 1)``. The smooth
          approximation of l1 (absolute value) loss. Usually a good
          choice for robust least squares.
        * 'huber' : ``rho(z) = z if z <= 1 else 2*z**0.5 - 1``. Works
          similarly to 'soft_l1'.
        * 'cauchy' : ``rho(z) = ln(1 + z)``. Severely weakens outliers
          influence, but may cause difficulties in optimization process.
        * 'arctan' : ``rho(z) = arctan(z)``. Limits a maximum loss on
          a single residual, has properties similar to 'cauchy'.

        If callable, it must take a 1-D ndarray ``z=f**2`` and return an
        array_like with shape (3, m) where row 0 contains function values,
        row 1 contains first derivatives and row 2 contains second
        derivatives. Method 'lm' supports only 'linear' loss.
    f_scale : float, optional
        Value of soft margin between inlier and outlier residuals, default
        is 1.0. The loss function is evaluated as follows
        ``rho_(f**2) = C**2 * rho(f**2 / C**2)``, where ``C`` is `f_scale`,
        and ``rho`` is determined by `loss` parameter. This parameter has
        no effect with ``loss='linear'``, but for other `loss` values it is
        of crucial importance.
    diff_step : None or array_like, optional
        Determines the relative step size for the finite difference
        approximation of the Jacobian. The actual step is computed as
        ``x * diff_step``. If None (default), then `diff_step` is taken to be
        a conventional "optimal" power of machine epsilon for the finite
        difference scheme used [NR]_.
    tr_solver : {None, 'exact', 'lsmr'}, optional
        Method for solving trust-region subproblems, relevant only for 'trf'
        and 'dogbox' methods.

        * 'exact' is suitable for not very large problems with dense
          Jacobian matrices. The computational complexity per iteration is
          comparable to a singular value decomposition of the Jacobian
          matrix.
        * 'lsmr' is suitable for problems with sparse and large Jacobian
          matrices. It uses the iterative procedure
          `scipy.sparse.linalg.lsmr` for finding a solution of a linear
          least-squares problem and only requires matrix-vector product
          evaluations.

        If None (default), the solver is chosen based on the type of Jacobian
        returned on the first iteration.
    tr_options : dict, optional
        Keyword options passed to trust-region solver.

        * ``tr_solver='exact'``: `tr_options` are ignored.
        * ``tr_solver='lsmr'``: options for `scipy.sparse.linalg.lsmr`.
          Additionally,  ``method='trf'`` supports  'regularize' option
          (bool, default is True), which adds a regularization term to the
          normal equation, which improves convergence if the Jacobian is
          rank-deficient [Byrd]_ (eq. 3.4).

    jac_sparsity : {None, array_like, sparse array}, optional
        Defines the sparsity structure of the Jacobian matrix for finite
        difference estimation, its shape must be (m, n). If the Jacobian has
        only few non-zero elements in *each* row, providing the sparsity
        structure will greatly speed up the computations [Curtis]_. A zero
        entry means that a corresponding element in the Jacobian is identically
        zero. If provided, forces the use of 'lsmr' trust-region solver.
        If None (default), then dense differencing will be used. Has no effect
        for 'lm' method.
    max_nfev : None or int, optional
        For all methods this parameter controls the maximum number of function
        evaluations used by each method, separate to those used in numerical
        approximation of the jacobian.
        If None (default), the value is chosen automatically as 100 * n.

        .. versionchanged:: 1.16.0
            The default for the 'lm' method is changed to 100 * n, for both a callable
            and a numerically estimated jacobian. Previously the default when using an
            estimated jacobian was 100 * n * (n + 1), because the method included
            evaluations used in the estimation.

    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

        * 0 (default) : work silently.
        * 1 : display a termination report.
        * 2 : display progress during iterations (not supported by 'lm'
          method).

    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun` and `jac`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)`` and the same for
        `jac`.
    callback : None or callable, optional
        Callback function that is called by the algorithm on each iteration.
        This can be used to print or plot the optimization results at each
        step, and to stop the optimization algorithm based on some user-defined
        condition.  Only implemented for the `trf` and `dogbox` methods.

        The signature is ``callback(intermediate_result: OptimizeResult)``

        `intermediate_result is a `scipy.optimize.OptimizeResult`
        which contains the intermediate results of the optimization at the
        current iteration.

        The callback also supports a signature like: ``callback(x)``

        Introspection is used to determine which of the signatures is invoked.

        If the `callback` function raises `StopIteration` the optimization algorithm
        will stop and return with status code -2.

        .. versionadded:: 1.16.0
    workers : map-like callable, optional
        A map-like callable, such as `multiprocessing.Pool.map` for evaluating
        any numerical differentiation in parallel.
        This evaluation is carried out as ``workers(fun, iterable)``.

        .. versionadded:: 1.16.0

    Returns
    -------
    result : OptimizeResult
        `OptimizeResult` with the following fields defined:

        x : ndarray, shape (n,)
            Solution found.
        cost : float
            Value of the cost function at the solution.
        fun : ndarray, shape (m,)
            Vector of residuals at the solution.
        jac : ndarray, sparse array or LinearOperator, shape (m, n)
            Modified Jacobian matrix at the solution, in the sense that J^T J
            is a Gauss-Newton approximation of the Hessian of the cost function.
            The type is the same as the one used by the algorithm.
        grad : ndarray, shape (m,)
            Gradient of the cost function at the solution.
        optimality : float
            First-order optimality measure. In unconstrained problems, it is
            always the uniform norm of the gradient. In constrained problems,
            it is the quantity which was compared with `gtol` during iterations.
        active_mask : ndarray of int, shape (n,)
            Each component shows whether a corresponding constraint is active
            (that is, whether a variable is at the bound):

            *  0 : a constraint is not active.
            * -1 : a lower bound is active.
            *  1 : an upper bound is active.

            Might be somewhat arbitrary for 'trf' method as it generates a
            sequence of strictly feasible iterates and `active_mask` is
            determined within a tolerance threshold.
        nfev : int
            Number of function evaluations done. This number does not include
            the function calls used for numerical Jacobian approximation.

            .. versionchanged:: 1.16.0
                For the 'lm' method the number of function calls used in numerical
                Jacobian approximation is no longer included. This is to bring all
                methods into line.

        njev : int or None
            Number of Jacobian evaluations done. If numerical Jacobian
            approximation is used in 'lm' method, it is set to None.
        status : int
            The reason for algorithm termination:

            * -2 : terminated because callback raised StopIteration.
            * -1 : improper input parameters status returned from MINPACK.
            *  0 : the maximum number of function evaluations is exceeded.
            *  1 : `gtol` termination condition is satisfied.
            *  2 : `ftol` termination condition is satisfied.
            *  3 : `xtol` termination condition is satisfied.
            *  4 : Both `ftol` and `xtol` termination conditions are satisfied.

        message : str
            Verbal description of the termination reason.
        success : bool
            True if one of the convergence criteria is satisfied (`status` > 0).

    See Also
    --------
    leastsq : A legacy wrapper for the MINPACK implementation of the
              Levenberg-Marquadt algorithm.
    curve_fit : Least-squares minimization applied to a curve-fitting problem.

    Notes
    -----
    Method 'lm' (Levenberg-Marquardt) calls a wrapper over a least-squares
    algorithm implemented in MINPACK (lmder). It runs the
    Levenberg-Marquardt algorithm formulated as a trust-region type algorithm.
    The implementation is based on paper [JJMore]_, it is very robust and
    efficient with a lot of smart tricks. It should be your first choice
    for unconstrained problems. Note that it doesn't support bounds. Also,
    it doesn't work when m < n.

    Method 'trf' (Trust Region Reflective) is motivated by the process of
    solving a system of equations, which constitute the first-order optimality
    condition for a bound-constrained minimization problem as formulated in
    [STIR]_. The algorithm iteratively solves trust-region subproblems
    augmented by a special diagonal quadratic term and with trust-region shape
    determined by the distance from the bounds and the direction of the
    gradient. This enhancements help to avoid making steps directly into bounds
    and efficiently explore the whole space of variables. To further improve
    convergence, the algorithm considers search directions reflected from the
    bounds. To obey theoretical requirements, the algorithm keeps iterates
    strictly feasible. With dense Jacobians trust-region subproblems are
    solved by an exact method very similar to the one described in [JJMore]_
    (and implemented in MINPACK). The difference from the MINPACK
    implementation is that a singular value decomposition of a Jacobian
    matrix is done once per iteration, instead of a QR decomposition and series
    of Givens rotation eliminations. For large sparse Jacobians a 2-D subspace
    approach of solving trust-region subproblems is used [STIR]_, [Byrd]_.
    The subspace is spanned by a scaled gradient and an approximate
    Gauss-Newton solution delivered by `scipy.sparse.linalg.lsmr`. When no
    constraints are imposed the algorithm is very similar to MINPACK and has
    generally comparable performance. The algorithm works quite robust in
    unbounded and bounded problems, thus it is chosen as a default algorithm.

    Method 'dogbox' operates in a trust-region framework, but considers
    rectangular trust regions as opposed to conventional ellipsoids [Voglis]_.
    The intersection of a current trust region and initial bounds is again
    rectangular, so on each iteration a quadratic minimization problem subject
    to bound constraints is solved approximately by Powell's dogleg method
    [NumOpt]_. The required Gauss-Newton step can be computed exactly for
    dense Jacobians or approximately by `scipy.sparse.linalg.lsmr` for large
    sparse Jacobians. The algorithm is likely to exhibit slow convergence when
    the rank of Jacobian is less than the number of variables. The algorithm
    often outperforms 'trf' in bounded problems with a small number of
    variables.

    Robust loss functions are implemented as described in [BA]_. The idea
    is to modify a residual vector and a Jacobian matrix on each iteration
    such that computed gradient and Gauss-Newton Hessian approximation match
    the true gradient and Hessian approximation of the cost function. Then
    the algorithm proceeds in a normal way, i.e., robust loss functions are
    implemented as a simple wrapper over standard least-squares algorithms.

    .. versionadded:: 0.17.0

    References
    ----------
    .. [STIR] M. A. Branch, T. F. Coleman, and Y. Li, "A Subspace, Interior,
              and Conjugate Gradient Method for Large-Scale Bound-Constrained
              Minimization Problems," SIAM Journal on Scientific Computing,
              Vol. 21, Number 1, pp 1-23, 1999.
    .. [NR] William H. Press et. al., "Numerical Recipes. The Art of Scientific
            Computing. 3rd edition", Sec. 5.7.
    .. [Byrd] R. H. Byrd, R. B. Schnabel and G. A. Shultz, "Approximate
              solution of the trust region problem by minimization over
              two-dimensional subspaces", Math. Programming, 40, pp. 247-263,
              1988.
    .. [Curtis] A. Curtis, M. J. D. Powell, and J. Reid, "On the estimation of
                sparse Jacobian matrices", Journal of the Institute of
                Mathematics and its Applications, 13, pp. 117-120, 1974.
    .. [JJMore] J. J. More, "The Levenberg-Marquardt Algorithm: Implementation
                and Theory," Numerical Analysis, ed. G. A. Watson, Lecture
                Notes in Mathematics 630, Springer Verlag, pp. 105-116, 1977.
    .. [Voglis] C. Voglis and I. E. Lagaris, "A Rectangular Trust Region
                Dogleg Approach for Unconstrained and Bound Constrained
                Nonlinear Optimization", WSEAS International Conference on
                Applied Mathematics, Corfu, Greece, 2004.
    .. [NumOpt] J. Nocedal and S. J. Wright, "Numerical optimization,
                2nd edition", Chapter 4.
    .. [BA] B. Triggs et. al., "Bundle Adjustment - A Modern Synthesis",
            Proceedings of the International Workshop on Vision Algorithms:
            Theory and Practice, pp. 298-372, 1999.

    Examples
    --------
    In this example we find a minimum of the Rosenbrock function without bounds
    on independent variables.

    >>> import numpy as np
    >>> def fun_rosenbrock(x):
    ...     return np.array([10 * (x[1] - x[0]**2), (1 - x[0])])

    Notice that we only provide the vector of the residuals. The algorithm
    constructs the cost function as a sum of squares of the residuals, which
    gives the Rosenbrock function. The exact minimum is at ``x = [1.0, 1.0]``.

    >>> from scipy.optimize import least_squares
    >>> x0_rosenbrock = np.array([2, 2])
    >>> res_1 = least_squares(fun_rosenbrock, x0_rosenbrock)
    >>> res_1.x
    array([ 1.,  1.])
    >>> res_1.cost
    9.8669242910846867e-30
    >>> res_1.optimality
    8.8928864934219529e-14

    We now constrain the variables, in such a way that the previous solution
    becomes infeasible. Specifically, we require that ``x[1] >= 1.5``, and
    ``x[0]`` left unconstrained. To this end, we specify the `bounds` parameter
    to `least_squares` in the form ``bounds=([-np.inf, 1.5], np.inf)``.

    We also provide the analytic Jacobian:

    >>> def jac_rosenbrock(x):
    ...     return np.array([
    ...         [-20 * x[0], 10],
    ...         [-1, 0]])

    Putting this all together, we see that the new solution lies on the bound:

    >>> res_2 = least_squares(fun_rosenbrock, x0_rosenbrock, jac_rosenbrock,
    ...                       bounds=([-np.inf, 1.5], np.inf))
    >>> res_2.x
    array([ 1.22437075,  1.5       ])
    >>> res_2.cost
    0.025213093946805685
    >>> res_2.optimality
    1.5885401433157753e-07

    Now we solve a system of equations (i.e., the cost function should be zero
    at a minimum) for a Broyden tridiagonal vector-valued function of 100000
    variables:

    >>> def fun_broyden(x):
    ...     f = (3 - x) * x + 1
    ...     f[1:] -= x[:-1]
    ...     f[:-1] -= 2 * x[1:]
    ...     return f

    The corresponding Jacobian matrix is sparse. We tell the algorithm to
    estimate it by finite differences and provide the sparsity structure of
    Jacobian to significantly speed up this process.

    >>> from scipy.sparse import lil_array
    >>> def sparsity_broyden(n):
    ...     sparsity = lil_array((n, n), dtype=int)
    ...     i = np.arange(n)
    ...     sparsity[i, i] = 1
    ...     i = np.arange(1, n)
    ...     sparsity[i, i - 1] = 1
    ...     i = np.arange(n - 1)
    ...     sparsity[i, i + 1] = 1
    ...     return sparsity
    ...
    >>> n = 100000
    >>> x0_broyden = -np.ones(n)
    ...
    >>> res_3 = least_squares(fun_broyden, x0_broyden,
    ...                       jac_sparsity=sparsity_broyden(n))
    >>> res_3.cost
    4.5687069299604613e-23
    >>> res_3.optimality
    1.1650454296851518e-11

    Let's also solve a curve fitting problem using robust loss function to
    take care of outliers in the data. Define the model function as
    ``y = a + b * exp(c * t)``, where t is a predictor variable, y is an
    observation and a, b, c are parameters to estimate.

    First, define the function which generates the data with noise and
    outliers, define the model parameters, and generate data:

    >>> from numpy.random import default_rng
    >>> rng = default_rng()
    >>> def gen_data(t, a, b, c, noise=0., n_outliers=0, seed=None):
    ...     rng = default_rng(seed)
    ...
    ...     y = a + b * np.exp(t * c)
    ...
    ...     error = noise * rng.standard_normal(t.size)
    ...     outliers = rng.integers(0, t.size, n_outliers)
    ...     error[outliers] *= 10
    ...
    ...     return y + error
    ...
    >>> a = 0.5
    >>> b = 2.0
    >>> c = -1
    >>> t_min = 0
    >>> t_max = 10
    >>> n_points = 15
    ...
    >>> t_train = np.linspace(t_min, t_max, n_points)
    >>> y_train = gen_data(t_train, a, b, c, noise=0.1, n_outliers=3)

    Define function for computing residuals and initial estimate of
    parameters.

    >>> def fun(x, t, y):
    ...     return x[0] + x[1] * np.exp(x[2] * t) - y
    ...
    >>> x0 = np.array([1.0, 1.0, 0.0])

    Compute a standard least-squares solution:

    >>> res_lsq = least_squares(fun, x0, args=(t_train, y_train))

    Now compute two solutions with two different robust loss functions. The
    parameter `f_scale` is set to 0.1, meaning that inlier residuals should
    not significantly exceed 0.1 (the noise level used).

    >>> res_soft_l1 = least_squares(fun, x0, loss='soft_l1', f_scale=0.1,
    ...                             args=(t_train, y_train))
    >>> res_log = least_squares(fun, x0, loss='cauchy', f_scale=0.1,
    ...                         args=(t_train, y_train))

    And, finally, plot all the curves. We see that by selecting an appropriate
    `loss`  we can get estimates close to optimal even in the presence of
    strong outliers. But keep in mind that generally it is recommended to try
    'soft_l1' or 'huber' losses first (if at all necessary) as the other two
    options may cause difficulties in optimization process.

    >>> t_test = np.linspace(t_min, t_max, n_points * 10)
    >>> y_true = gen_data(t_test, a, b, c)
    >>> y_lsq = gen_data(t_test, *res_lsq.x)
    >>> y_soft_l1 = gen_data(t_test, *res_soft_l1.x)
    >>> y_log = gen_data(t_test, *res_log.x)
    ...
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(t_train, y_train, 'o')
    >>> plt.plot(t_test, y_true, 'k', linewidth=2, label='true')
    >>> plt.plot(t_test, y_lsq, label='linear loss')
    >>> plt.plot(t_test, y_soft_l1, label='soft_l1 loss')
    >>> plt.plot(t_test, y_log, label='cauchy loss')
    >>> plt.xlabel("t")
    >>> plt.ylabel("y")
    >>> plt.legend()
    >>> plt.show()

    In the next example, we show how complex-valued residual functions of
    complex variables can be optimized with ``least_squares()``. Consider the
    following function:

    >>> def f(z):
    ...     return z - (0.5 + 0.5j)

    We wrap it into a function of real variables that returns real residuals
    by simply handling the real and imaginary parts as independent variables:

    >>> def f_wrap(x):
    ...     fx = f(x[0] + 1j*x[1])
    ...     return np.array([fx.real, fx.imag])

    Thus, instead of the original m-D complex function of n complex
    variables we optimize a 2m-D real function of 2n real variables:

    >>> from scipy.optimize import least_squares
    >>> res_wrapped = least_squares(f_wrap, (0.1, 0.1), bounds=([0, 0], [1, 1]))
    >>> z = res_wrapped.x[0] + res_wrapped.x[1]*1j
    >>> z
    (0.49999999999925893+0.49999999999925893j)

    """
    if method not in ['trf', 'dogbox', 'lm']:
        raise ValueError("`method` must be 'trf', 'dogbox' or 'lm'.")

    if jac not in ['2-point', '3-point', 'cs'] and not callable(jac):
        raise ValueError("`jac` must be '2-point', '3-point', 'cs' or "
                         "callable.")

    if tr_solver not in [None, 'exact', 'lsmr']:
        raise ValueError("`tr_solver` must be None, 'exact' or 'lsmr'.")

    if loss not in IMPLEMENTED_LOSSES and not callable(loss):
        raise ValueError(f"`loss` must be one of {IMPLEMENTED_LOSSES.keys()}"
                         " or a callable.")

    if method == 'lm' and loss != 'linear':
        raise ValueError("method='lm' supports only 'linear' loss function.")

    if verbose not in [0, 1, 2]:
        raise ValueError("`verbose` must be in [0, 1, 2].")

    if max_nfev is not None and max_nfev <= 0:
        raise ValueError("`max_nfev` must be None or positive integer.")

    if np.iscomplexobj(x0):
        raise ValueError("`x0` must be real.")

    x0 = np.atleast_1d(x0).astype(float)

    if x0.ndim > 1:
        raise ValueError("`x0` must have at most 1 dimension.")

    if isinstance(bounds, Bounds):
        lb, ub = bounds.lb, bounds.ub
        bounds = (lb, ub)
    else:
        if len(bounds) == 2:
            lb, ub = prepare_bounds(bounds, x0.shape[0])
        else:
            raise ValueError("`bounds` must contain 2 elements.")

    if method == 'lm' and not np.all((lb == -np.inf) & (ub == np.inf)):
        raise ValueError("Method 'lm' doesn't support bounds.")

    if lb.shape != x0.shape or ub.shape != x0.shape:
        raise ValueError("Inconsistent shapes between bounds and `x0`.")

    if np.any(lb >= ub):
        raise ValueError("Each lower bound must be strictly less than each "
                         "upper bound.")

    if not in_bounds(x0, lb, ub):
        raise ValueError("Initial guess is outside of provided bounds")

    x_scale = check_x_scale(x_scale, x0, method)

    ftol, xtol, gtol = check_tolerance(ftol, xtol, gtol, method)

    if method == 'trf':
        x0 = make_strictly_feasible(x0, lb, ub)

    if tr_options is None:
        tr_options = {}

    ###########################################################################
    # assemble VectorFunction
    ###########################################################################
    # first wrap the args/kwargs
    fun_wrapped = _WrapArgsKwargs(fun, args=args, kwargs=kwargs)
    jac_wrapped = jac
    if callable(jac):
        jac_wrapped = _WrapArgsKwargs(jac, args=args, kwargs=kwargs)

    def _dummy_hess(x, *args):
        # we don't care about Hessian evaluations
        return x

    vector_fun = VectorFunction(
        fun_wrapped,
        x0,
        jac_wrapped,
        _dummy_hess,
        finite_diff_rel_step=diff_step,
        finite_diff_jac_sparsity=jac_sparsity,
        finite_diff_bounds=bounds,
        workers=workers
    )
    ###########################################################################

    f0 = vector_fun.fun(x0)
    J0 = vector_fun.jac(x0)

    if f0.ndim != 1:
        raise ValueError("`fun` must return at most 1-d array_like. "
                         f"f0.shape: {f0.shape}")

    if not np.all(np.isfinite(f0)):
        raise ValueError("Residuals are not finite in the initial point.")

    n = x0.size
    m = f0.size

    if method == 'lm' and m < n:
        raise ValueError("Method 'lm' doesn't work when the number of "
                         "residuals is less than the number of variables.")

    loss_function = construct_loss_function(m, loss, f_scale)
    if callable(loss):
        rho = loss_function(f0)
        if rho.shape != (3, m):
            raise ValueError("The return value of `loss` callable has wrong "
                             "shape.")
        initial_cost = 0.5 * np.sum(rho[0])
    elif loss_function is not None:
        initial_cost = loss_function(f0, cost_only=True)
    else:
        initial_cost = 0.5 * np.dot(f0, f0)

    if not callable(jac):
        # Estimate Jacobian by finite differences.
        if method == 'lm':
            if jac_sparsity is not None:
                raise ValueError("method='lm' does not support "
                                 "`jac_sparsity`.")
        else:
            # this will raise a ValueError if the jac_sparsity isn't correct
            _ = check_jac_sparsity(jac_sparsity, m, n)

            if jac_sparsity is not None and tr_solver == 'exact':
                raise ValueError("tr_solver='exact' is incompatible "
                                 "with `jac_sparsity`.")

    if J0.shape != (m, n):
        raise ValueError(
            f"The return value of `jac` has wrong shape: expected {(m, n)}, "
            f"actual {J0.shape}."
        )

    if not isinstance(J0, np.ndarray):
        if method == 'lm':
            raise ValueError("method='lm' works only with dense "
                             "Jacobian matrices.")

        if tr_solver == 'exact':
            raise ValueError(
                "tr_solver='exact' works only with dense "
                "Jacobian matrices.")

    jac_scale = isinstance(x_scale, str) and x_scale == 'jac'
    if isinstance(J0, LinearOperator) and jac_scale:
        raise ValueError("x_scale='jac' can't be used when `jac` "
                         "returns LinearOperator.")

    if tr_solver is None:
        if isinstance(J0, np.ndarray):
            tr_solver = 'exact'
        else:
            tr_solver = 'lsmr'

    # Wrap callback function.  If callback is None, callback_wrapped also is None
    callback_wrapped = _wrap_callback(callback)

    if method == 'lm':
        if callback is not None:
            warn("Callback function specified, but not supported with `lm` method.",
                 stacklevel=2)
        result = call_minpack(vector_fun.fun, x0, vector_fun.jac, ftol, xtol, gtol,
                              max_nfev, x_scale, jac_method=jac)

    elif method == 'trf':
        result = trf(vector_fun.fun, vector_fun.jac, x0, f0, J0, lb, ub, ftol, xtol,
                     gtol, max_nfev, x_scale, loss_function, tr_solver,
                     tr_options.copy(), verbose, callback=callback_wrapped)

    elif method == 'dogbox':
        if tr_solver == 'lsmr' and 'regularize' in tr_options:
            warn("The keyword 'regularize' in `tr_options` is not relevant "
                 "for 'dogbox' method.",
                 stacklevel=2)
            tr_options = tr_options.copy()
            del tr_options['regularize']

        result = dogbox(vector_fun.fun, vector_fun.jac, x0, f0, J0, lb, ub, ftol,
                        xtol, gtol, max_nfev, x_scale, loss_function,
                        tr_solver, tr_options, verbose, callback=callback_wrapped)

    result.message = TERMINATION_MESSAGES[result.status]
    result.success = result.status > 0

    if verbose >= 1:
        print(result.message)
        print(f"Function evaluations {result.nfev}, initial cost {initial_cost:.4e}, "
              f"final cost {result.cost:.4e}, "
              f"first-order optimality {result.optimality:.2e}.")

    return result

################################################################################
# Method Path: scipy.linalg._expm_frechet._diff_pade9
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_expm_frechet.py
# Target:      dot
# Call Count:  780
# Total Exec:  52
# Func Size:   17
################################################################################

def _diff_pade9(A, E, ident):
    b = (17643225600., 8821612800., 2075673600., 302702400., 30270240.,
            2162160., 110880., 3960., 90., 1.)
    A2 = A.dot(A)
    M2 = np.dot(A, E) + np.dot(E, A)
    A4 = np.dot(A2, A2)
    M4 = np.dot(A2, M2) + np.dot(M2, A2)
    A6 = np.dot(A2, A4)
    M6 = np.dot(A4, M2) + np.dot(M4, A2)
    A8 = np.dot(A4, A4)
    M8 = np.dot(A4, M4) + np.dot(M4, A4)
    U = A.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[8]*A8 + b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    Lu = (A.dot(b[9]*M8 + b[7]*M6 + b[5]*M4 + b[3]*M2) +
            E.dot(b[9]*A8 + b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident))
    Lv = b[8]*M8 + b[6]*M6 + b[4]*M4 + b[2]*M2
    return U, V, Lu, Lv

################################################################################
# Method Path: scipy.integrate._ivp.rk.RK45._dense_output_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  739
# Total Exec:  739
# Func Size:   3
################################################################################

def _dense_output_impl(self):
        Q = self.K.T.dot(self.P)
        return RkDenseOutput(self.t_old, self.t, self.y_old, Q)

################################################################################
# Method Path: scipy.optimize._remove_redundancy._remove_redundancy_pivot_sparse
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_remove_redundancy.py
# Target:      dot
# Call Count:  721
# Total Exec:  37
# Func Size:   125
################################################################################

def _remove_redundancy_pivot_sparse(A, rhs):
    """
    Eliminates redundant equations from system of equations defined by Ax = b
    and identifies infeasibilities.

    Parameters
    ----------
    A : 2-D sparse array
        An matrix representing the left-hand side of a system of equations
    rhs : 1-D array
        An array representing the right-hand side of a system of equations

    Returns
    -------
    A : 2-D sparse array
        A matrix representing the left-hand side of a system of equations
    rhs : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the system
        0: No infeasibility identified
        2: Trivially infeasible
    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [2] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.

    """

    tolapiv = 1e-8
    tolprimal = 1e-8
    status = 0
    message = ""
    inconsistent = ("There is a linear combination of rows of A_eq that "
                    "results in zero, suggesting a redundant constraint. "
                    "However the same linear combination of b_eq is "
                    "nonzero, suggesting that the constraints conflict "
                    "and the problem is infeasible.")
    A, rhs, status, message = _remove_zero_rows(A, rhs)

    if status != 0:
        return A, rhs, status, message

    m, n = A.shape

    v = list(range(m))      # Artificial column indices.
    b = list(v)             # Basis column indices.
    # This is better as a list than a set because column order of basis matrix
    # needs to be consistent.
    k = set(range(m, m+n))  # Structural column indices.
    d = []                  # Indices of dependent rows

    A_orig = A
    A = scipy.sparse.hstack((scipy.sparse.eye_array(m), A)).tocsc()
    e = np.zeros(m)

    # Implements basic algorithm from [2]
    # Uses only one of the suggested improvements (removing zero rows).
    # Removing column singletons would be easy, but it is not as important
    # because the procedure is performed only on the equality constraint
    # matrix from the original problem - not on the canonical form matrix,
    # which would have many more column singletons due to slack variables
    # from the inequality constraints.
    # The thoughts on "crashing" the initial basis sound useful, but the
    # description of the procedure seems to assume a lot of familiarity with
    # the subject; it is not very explicit. I already went through enough
    # trouble getting the basic algorithm working, so I was not interested in
    # trying to decipher this, too. (Overall, the paper is fraught with
    # mistakes and ambiguities - which is strange, because the rest of
    # Andersen's papers are quite good.)
    # I tried and tried and tried to improve performance using the
    # Bartels-Golub update. It works, but it's only practical if the LU
    # factorization can be specialized as described, and that is not possible
    # until the SciPy SuperLU interface permits control over column
    # permutation - see issue #7700.

    for i in v:
        B = A[:, b]

        e[i] = 1
        if i > 0:
            e[i-1] = 0

        pi = scipy.sparse.linalg.spsolve(B.transpose(), e).reshape(-1, 1)

        js = list(k-set(b))  # not efficient, but this is not the time sink...

        # Due to overhead, it tends to be faster (for problems tested) to
        # compute the full matrix-vector product rather than individual
        # vector-vector products (with the chance of terminating as soon
        # as any are nonzero). For very large matrices, it might be worth
        # it to compute, say, 100 or 1000 at a time and stop when a nonzero
        # is found.

        c = (np.abs(A[:, js].transpose().dot(pi)) > tolapiv).nonzero()[0]
        if len(c) > 0:  # independent
            j = js[c[0]]
            # in a previous commit, the previous line was changed to choose
            # index j corresponding with the maximum dot product.
            # While this avoided issues with almost
            # singular matrices, it slowed the routine in most NETLIB tests.
            # I think this is because these columns were denser than the
            # first column with nonzero dot product (c[0]).
            # It would be nice to have a heuristic that balances sparsity with
            # high dot product, but I don't think it's worth the time to
            # develop one right now. Bartels-Golub update is a much higher
            # priority.
            b[i] = j  # replace artificial column
        else:
            bibar = pi.T.dot(rhs.reshape(-1, 1))
            bnorm = np.linalg.norm(rhs)
            if abs(bibar)/(1 + bnorm) > tolprimal:
                status = 2
                message = inconsistent
                return A_orig, rhs, status, message
            else:  # dependent
                d.append(i)

    keep = set(range(m))
    keep = list(keep - set(d))
    return A_orig[keep, :], rhs[keep], status, message

################################################################################
# Method Path: scipy.optimize._differentiable_functions._VectorHessWrapper._fd_hess
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_differentiable_functions.py
# Target:      dot
# Call Count:  697
# Total Exec:  697
# Func Size:   11
################################################################################

def _fd_hess(self, x, v, J0=None):
        if J0 is None:
            J0 = self.jac(x)
            self.njev += 1

        # H will be a LinearOperator
        H = approx_derivative(self.jac_dot_v, x,
                              f0=J0.T.dot(v),
                              args=(v,),
                              **self.finite_diff_options)
        return H

################################################################################
# Method Path: scipy.optimize._lsq.trf_linear.trf_linear
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/trf_linear.py
# Target:      dot
# Call Count:  690
# Total Exec:  24
# Func Size:   108
################################################################################

def trf_linear(A, b, x_lsq, lb, ub, tol, lsq_solver, lsmr_tol,
               max_iter, verbose, *, lsmr_maxiter=None):
    m, n = A.shape
    x, _ = reflective_transformation(x_lsq, lb, ub)
    x = make_strictly_feasible(x, lb, ub, rstep=0.1)

    if lsq_solver == 'exact':
        QT, R, perm = qr(A, mode='economic', pivoting=True)
        QT = QT.T

        if m < n:
            R = np.vstack((R, np.zeros((n - m, n))))

        QTr = np.zeros(n)
        k = min(m, n)
    elif lsq_solver == 'lsmr':
        r_aug = np.zeros(m + n)
        auto_lsmr_tol = False
        if lsmr_tol is None:
            lsmr_tol = 1e-2 * tol
        elif lsmr_tol == 'auto':
            auto_lsmr_tol = True

    r = A.dot(x) - b
    g = compute_grad(A, r)
    cost = 0.5 * np.dot(r, r)
    initial_cost = cost

    termination_status = None
    step_norm = None
    cost_change = None

    if max_iter is None:
        max_iter = 100

    if verbose == 2:
        print_header_linear()

    for iteration in range(max_iter):
        v, dv = CL_scaling_vector(x, g, lb, ub)
        g_scaled = g * v
        g_norm = norm(g_scaled, ord=np.inf)
        if g_norm < tol:
            termination_status = 1

        if verbose == 2:
            print_iteration_linear(iteration, cost, cost_change,
                                   step_norm, g_norm)

        if termination_status is not None:
            break

        diag_h = g * dv
        diag_root_h = diag_h ** 0.5
        d = v ** 0.5
        g_h = d * g

        A_h = right_multiplied_operator(A, d)
        if lsq_solver == 'exact':
            QTr[:k] = QT.dot(r)
            p_h = -regularized_lsq_with_qr(m, n, R * d[perm], QTr, perm,
                                           diag_root_h, copy_R=False)
        elif lsq_solver == 'lsmr':
            lsmr_op = regularized_lsq_operator(A_h, diag_root_h)
            r_aug[:m] = r
            if auto_lsmr_tol:
                eta = 1e-2 * min(0.5, g_norm)
                lsmr_tol = max(EPS, min(0.1, eta * g_norm))
            p_h = -lsmr(lsmr_op, r_aug, maxiter=lsmr_maxiter,
                        atol=lsmr_tol, btol=lsmr_tol)[0]

        p = d * p_h

        p_dot_g = np.dot(p, g)
        if p_dot_g > 0:
            termination_status = -1

        theta = 1 - min(0.005, g_norm)
        step = select_step(x, A_h, g_h, diag_h, p, p_h, d, lb, ub, theta)
        cost_change = -evaluate_quadratic(A, g, step)

        # Perhaps almost never executed, the idea is that `p` is descent
        # direction thus we must find acceptable cost decrease using simple
        # "backtracking", otherwise the algorithm's logic would break.
        if cost_change < 0:
            x, step, cost_change = backtracking(
                A, g, x, p, theta, p_dot_g, lb, ub)
        else:
            x = make_strictly_feasible(x + step, lb, ub, rstep=0)

        step_norm = norm(step)
        r = A.dot(x) - b
        g = compute_grad(A, r)

        if cost_change < tol * cost:
            termination_status = 2

        cost = 0.5 * np.dot(r, r)

    if termination_status is None:
        termination_status = 0

    active_mask = find_active_constraints(x, lb, ub, rtol=tol)

    return OptimizeResult(
        x=x, fun=r, cost=cost, optimality=g_norm, active_mask=active_mask,
        nit=iteration + 1, status=termination_status,
        initial_cost=initial_cost)

################################################################################
# Method Path: scipy.integrate._ivp.bdf.change_D
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/bdf.py
# Target:      dot
# Call Count:  676
# Total Exec:  338
# Func Size:   6
################################################################################

def change_D(D, order, factor):
    """Change differences array in-place when step size is changed."""
    R = compute_R(order, factor)
    U = compute_R(order, 1)
    RU = R.dot(U)
    D[:order + 1] = np.dot(RU.T, D[:order + 1])

################################################################################
# Method Path: scipy.linalg._matfuncs.funm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs.py
# Target:      dot
# Call Count:  622
# Total Exec:  311
# Func Size:   93
################################################################################

def funm(A, func, disp=True):
    """
    Evaluate a matrix function specified by a callable.

    Returns the value of matrix-valued function ``f`` at `A`. The
    function ``f`` is an extension of the scalar-valued function `func`
    to matrices.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix at which to evaluate the function
    func : callable
        Callable object that evaluates a scalar function f.
        Must be vectorized (eg. using vectorize).
    disp : bool, optional
        Print warning if error in the result is estimated large
        instead of returning estimated error. (Default: True)

    Returns
    -------
    funm : (N, N) ndarray
        Value of the matrix function specified by func evaluated at `A`
    errest : float
        (if disp == False)

        1-norm of the estimated error, ||err||_1 / ||A||_1

    Notes
    -----
    This function implements the general algorithm based on Schur decomposition
    (Algorithm 9.1.1. in [1]_).

    If the input matrix is known to be diagonalizable, then relying on the
    eigendecomposition is likely to be faster. For example, if your matrix is
    Hermitian, you can do

    >>> from scipy.linalg import eigh
    >>> def funm_herm(a, func, check_finite=False):
    ...     w, v = eigh(a, check_finite=check_finite)
    ...     ## if you further know that your matrix is positive semidefinite,
    ...     ## you can optionally guard against precision errors by doing
    ...     # w = np.maximum(w, 0)
    ...     w = func(w)
    ...     return (v * w).dot(v.conj().T)

    References
    ----------
    .. [1] Gene H. Golub, Charles F. van Loan, Matrix Computations 4th ed.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import funm
    >>> a = np.array([[1.0, 3.0], [1.0, 4.0]])
    >>> funm(a, lambda x: x*x)
    array([[  4.,  15.],
           [  5.,  19.]])
    >>> a.dot(a)
    array([[  4.,  15.],
           [  5.,  19.]])

    """
    A = _asarray_square(A)
    # Perform Shur decomposition (lapack ?gees)
    T, Z = schur(A)
    T, Z = rsf2csf(T, Z)
    n, n = T.shape
    F = diag(func(diag(T)))  # apply function to diagonal elements
    F = F.astype(T.dtype.char)  # e.g., when F is real but T is complex

    minden = abs(T[0, 0])

    # implement Algorithm 11.1.1 from Golub and Van Loan
    #                 "matrix Computations."
    F, minden = _funm_loops(F, T, n, minden)

    F = dot(dot(Z, F), transpose(conjugate(Z)))
    F = _maybe_real(A, F)

    tol = {0: feps, 1: eps}[_array_precision[F.dtype.char]]
    if minden == 0.0:
        minden = tol
    err = min(1, max(tol, (tol/minden)*norm(triu(T, 1), 1)))
    if prod(ravel(logical_not(isfinite(F))), axis=0):
        err = np.inf
    if disp:
        if err > 1000*tol:
            print("funm result may be inaccurate, approximate err =", err)
        return F
    else:
        return F, err

################################################################################
# Method Path: scipy.optimize._linprog_util._presolve
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_util.py
# Target:      dot
# Call Count:  590
# Total Exec:  587
# Func Size:   440
################################################################################

def _presolve(lp, rr, rr_method, tol=1e-9):
    """
    Given inputs for a linear programming problem in preferred format,
    presolve the problem: identify trivial infeasibilities, redundancies,
    and unboundedness, tighten bounds where possible, and eliminate fixed
    variables.

    Parameters
    ----------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, as ``min`` and ``max`` pairs, one for each of the N
            elements of ``x``. The N x 2 array contains lower bounds in the first
            column and upper bounds in the 2nd. Unbounded variables have lower
            bound -np.inf and/or upper bound np.inf.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    rr : bool
        If ``True`` attempts to eliminate any redundant rows in ``A_eq``.
        Set False if ``A_eq`` is known to be of full row rank, or if you are
        looking for a potential speedup (at the expense of reliability).
    rr_method : str
        Method used to identify and remove redundant rows from the
        equality constraint matrix after presolve.
    tol : float
        The tolerance which determines when a solution is "close enough" to
        zero in Phase 1 to be considered a basic feasible solution or close
        enough to positive to serve as an optimal solution.

    Returns
    -------
    lp : A `scipy.optimize._linprog_util._LPProblem` consisting of the following fields:

        c : 1D array
            The coefficients of the linear objective function to be minimized.
        A_ub : 2D array, optional
            The inequality constraint matrix. Each row of ``A_ub`` specifies the
            coefficients of a linear inequality constraint on ``x``.
        b_ub : 1D array, optional
            The inequality constraint vector. Each element represents an
            upper bound on the corresponding value of ``A_ub @ x``.
        A_eq : 2D array, optional
            The equality constraint matrix. Each row of ``A_eq`` specifies the
            coefficients of a linear equality constraint on ``x``.
        b_eq : 1D array, optional
            The equality constraint vector. Each element of ``A_eq @ x`` must equal
            the corresponding element of ``b_eq``.
        bounds : 2D array
            The bounds of ``x``, as ``min`` and ``max`` pairs, possibly tightened.
        x0 : 1D array, optional
            Guess values of the decision variables, which will be refined by
            the optimization algorithm. This argument is currently used only by the
            'revised simplex' method, and can only be used if `x0` represents a
            basic feasible solution.

    c0 : 1D array
        Constant term in objective function due to fixed (and eliminated)
        variables.
    x : 1D array
        Solution vector (when the solution is trivial and can be determined
        in presolve)
    revstack: list of functions
        the functions in the list reverse the operations of _presolve()
        the function signature is x_org = f(x_mod), where x_mod is the result
        of a presolve step and x_org the value at the start of the step
        (currently, the revstack contains only one function)
    complete: bool
        Whether the solution is complete (solved or determined to be infeasible
        or unbounded in presolve)
    status : int
        An integer representing the exit status of the optimization::

         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded
         4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [5] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.
    .. [8] Andersen, Erling D., and Knud D. Andersen. "Presolving in linear
           programming." Mathematical Programming 71.2 (1995): 221-245.

    """
    # ideas from Reference [5] by Andersen and Andersen
    # however, unlike the reference, this is performed before converting
    # problem to standard form
    # There are a few advantages:
    #  * artificial variables have not been added, so matrices are smaller
    #  * bounds have not been converted to constraints yet. (It is better to
    #    do that after presolve because presolve may adjust the simple bounds.)
    # There are many improvements that can be made, namely:
    #  * implement remaining checks from [5]
    #  * loop presolve until no additional changes are made
    #  * implement additional efficiency improvements in redundancy removal [2]

    c, A_ub, b_ub, A_eq, b_eq, bounds, x0, _ = lp

    revstack = []               # record of variables eliminated from problem
    # constant term in cost function may be added if variables are eliminated
    c0 = 0
    complete = False        # complete is True if detected infeasible/unbounded
    x = np.zeros(c.shape)   # this is solution vector if completed in presolve

    status = 0              # all OK unless determined otherwise
    message = ""

    # Lower and upper bounds. Copy to prevent feedback.
    lb = bounds[:, 0].copy()
    ub = bounds[:, 1].copy()

    m_eq, n = A_eq.shape
    m_ub, n = A_ub.shape

    if (rr_method is not None
            and rr_method.lower() not in {"svd", "pivot", "id"}):
        message = ("'" + str(rr_method) + "' is not a valid option "
                   "for redundancy removal. Valid options are 'SVD', "
                   "'pivot', and 'ID'.")
        raise ValueError(message)

    if sps.issparse(A_eq):
        A_eq = A_eq.tocsr()
        A_ub = A_ub.tocsr()

        def where(A):
            return A.nonzero()

        vstack = sps.vstack
    else:
        where = np.where
        vstack = np.vstack

    # upper bounds > lower bounds
    if np.any(ub < lb) or np.any(lb == np.inf) or np.any(ub == -np.inf):
        status = 2
        message = ("The problem is (trivially) infeasible since one "
                   "or more upper bounds are smaller than the corresponding "
                   "lower bounds, a lower bound is np.inf or an upper bound "
                   "is -np.inf.")
        complete = True
        return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                c0, x, revstack, complete, status, message)

    # zero row in equality constraints
    zero_row = np.array(np.sum(A_eq != 0, axis=1) == 0).flatten()
    if np.any(zero_row):
        if np.any(
            np.logical_and(
                zero_row,
                np.abs(b_eq) > tol)):  # test_zero_row_1
            # infeasible if RHS is not zero
            status = 2
            message = ("The problem is (trivially) infeasible due to a row "
                       "of zeros in the equality constraint matrix with a "
                       "nonzero corresponding constraint value.")
            complete = True
            return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                    c0, x, revstack, complete, status, message)
        else:  # test_zero_row_2
            # if RHS is zero, we can eliminate this equation entirely
            A_eq = A_eq[np.logical_not(zero_row), :]
            b_eq = b_eq[np.logical_not(zero_row)]

    # zero row in inequality constraints
    zero_row = np.array(np.sum(A_ub != 0, axis=1) == 0).flatten()
    if np.any(zero_row):
        if np.any(np.logical_and(zero_row, b_ub < -tol)):  # test_zero_row_1
            # infeasible if RHS is less than zero (because LHS is zero)
            status = 2
            message = ("The problem is (trivially) infeasible due to a row "
                       "of zeros in the equality constraint matrix with a "
                       "nonzero corresponding  constraint value.")
            complete = True
            return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                    c0, x, revstack, complete, status, message)
        else:  # test_zero_row_2
            # if LHS is >= 0, we can eliminate this constraint entirely
            A_ub = A_ub[np.logical_not(zero_row), :]
            b_ub = b_ub[np.logical_not(zero_row)]

    # zero column in (both) constraints
    # this indicates that a variable isn't constrained and can be removed
    A = vstack((A_eq, A_ub))
    if A.shape[0] > 0:
        zero_col = np.array(np.sum(A != 0, axis=0) == 0).flatten()
        # variable will be at upper or lower bound, depending on objective
        x[np.logical_and(zero_col, c < 0)] = ub[
            np.logical_and(zero_col, c < 0)]
        x[np.logical_and(zero_col, c > 0)] = lb[
            np.logical_and(zero_col, c > 0)]
        if np.any(np.isinf(x)):  # if an unconstrained variable has no bound
            status = 3
            message = ("If feasible, the problem is (trivially) unbounded "
                       "due  to a zero column in the constraint matrices. If "
                       "you wish to check whether the problem is infeasible, "
                       "turn presolve off.")
            complete = True
            return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                    c0, x, revstack, complete, status, message)
        # variables will equal upper/lower bounds will be removed later
        lb[np.logical_and(zero_col, c < 0)] = ub[
            np.logical_and(zero_col, c < 0)]
        ub[np.logical_and(zero_col, c > 0)] = lb[
            np.logical_and(zero_col, c > 0)]

    # row singleton in equality constraints
    # this fixes a variable and removes the constraint
    singleton_row = np.array(np.sum(A_eq != 0, axis=1) == 1).flatten()
    rows = where(singleton_row)[0]
    cols = where(A_eq[rows, :])[1]
    if len(rows) > 0:
        for row, col in zip(rows, cols):
            val = b_eq[row] / A_eq[row, col]
            if not lb[col] - tol <= val <= ub[col] + tol:
                # infeasible if fixed value is not within bounds
                status = 2
                message = ("The problem is (trivially) infeasible because a "
                           "singleton row in the equality constraints is "
                           "inconsistent with the bounds.")
                complete = True
                return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                        c0, x, revstack, complete, status, message)
            else:
                # sets upper and lower bounds at that fixed value - variable
                # will be removed later
                lb[col] = val
                ub[col] = val
        A_eq = A_eq[np.logical_not(singleton_row), :]
        b_eq = b_eq[np.logical_not(singleton_row)]

    # row singleton in inequality constraints
    # this indicates a simple bound and the constraint can be removed
    # simple bounds may be adjusted here
    # After all of the simple bound information is combined here, get_Abc will
    # turn the simple bounds into constraints
    singleton_row = np.array(np.sum(A_ub != 0, axis=1) == 1).flatten()
    cols = where(A_ub[singleton_row, :])[1]
    rows = where(singleton_row)[0]
    if len(rows) > 0:
        for row, col in zip(rows, cols):
            val = b_ub[row] / A_ub[row, col]
            if A_ub[row, col] > 0:  # upper bound
                if val < lb[col] - tol:  # infeasible
                    complete = True
                elif val < ub[col]:  # new upper bound
                    ub[col] = val
            else:  # lower bound
                if val > ub[col] + tol:  # infeasible
                    complete = True
                elif val > lb[col]:  # new lower bound
                    lb[col] = val
            if complete:
                status = 2
                message = ("The problem is (trivially) infeasible because a "
                           "singleton row in the upper bound constraints is "
                           "inconsistent with the bounds.")
                return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                        c0, x, revstack, complete, status, message)
        A_ub = A_ub[np.logical_not(singleton_row), :]
        b_ub = b_ub[np.logical_not(singleton_row)]

    # identical bounds indicate that variable can be removed
    i_f = np.abs(lb - ub) < tol   # indices of "fixed" variables
    i_nf = np.logical_not(i_f)  # indices of "not fixed" variables

    # test_bounds_equal_but_infeasible
    if np.all(i_f):  # if bounds define solution, check for consistency
        residual = b_eq - A_eq.dot(lb)
        slack = b_ub - A_ub.dot(lb)
        if ((A_ub.size > 0 and np.any(slack < 0)) or
                (A_eq.size > 0 and not np.allclose(residual, 0))):
            status = 2
            message = ("The problem is (trivially) infeasible because the "
                       "bounds fix all variables to values inconsistent with "
                       "the constraints")
            complete = True
            return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                    c0, x, revstack, complete, status, message)

    ub_mod = ub
    lb_mod = lb
    if np.any(i_f):
        c0 += c[i_f].dot(lb[i_f])
        b_eq = b_eq - A_eq[:, i_f].dot(lb[i_f])
        b_ub = b_ub - A_ub[:, i_f].dot(lb[i_f])
        c = c[i_nf]
        x_undo = lb[i_f]  # not x[i_f], x is just zeroes
        x = x[i_nf]
        # user guess x0 stays separate from presolve solution x
        if x0 is not None:
            x0 = x0[i_nf]
        A_eq = A_eq[:, i_nf]
        A_ub = A_ub[:, i_nf]
        # modify bounds
        lb_mod = lb[i_nf]
        ub_mod = ub[i_nf]

        def rev(x_mod):
            # Function to restore x: insert x_undo into x_mod.
            # When elements have been removed at positions k1, k2, k3, ...
            # then these must be replaced at (after) positions k1-1, k2-2,
            # k3-3, ... in the modified array to recreate the original
            i = np.flatnonzero(i_f)
            # Number of variables to restore
            N = len(i)
            index_offset = np.arange(N)
            # Create insert indices
            insert_indices = i - index_offset
            x_rev = np.insert(x_mod.astype(float), insert_indices, x_undo)
            return x_rev

        # Use revstack as a list of functions, currently just this one.
        revstack.append(rev)

    # no constraints indicates that problem is trivial
    if A_eq.size == 0 and A_ub.size == 0:
        b_eq = np.array([])
        b_ub = np.array([])
        # test_empty_constraint_1
        if c.size == 0:
            status = 0
            message = ("The solution was determined in presolve as there are "
                       "no non-trivial constraints.")
        elif (np.any(np.logical_and(c < 0, ub_mod == np.inf)) or
              np.any(np.logical_and(c > 0, lb_mod == -np.inf))):
            # test_no_constraints()
            # test_unbounded_no_nontrivial_constraints_1
            # test_unbounded_no_nontrivial_constraints_2
            status = 3
            message = ("The problem is (trivially) unbounded "
                       "because there are no non-trivial constraints and "
                       "a) at least one decision variable is unbounded "
                       "above and its corresponding cost is negative, or "
                       "b) at least one decision variable is unbounded below "
                       "and its corresponding cost is positive. ")
        else:  # test_empty_constraint_2
            status = 0
            message = ("The solution was determined in presolve as there are "
                       "no non-trivial constraints.")
        complete = True
        x[c < 0] = ub_mod[c < 0]
        x[c > 0] = lb_mod[c > 0]
        # where c is zero, set x to a finite bound or zero
        x_zero_c = ub_mod[c == 0]
        x_zero_c[np.isinf(x_zero_c)] = ub_mod[c == 0][np.isinf(x_zero_c)]
        x_zero_c[np.isinf(x_zero_c)] = 0
        x[c == 0] = x_zero_c
        # if this is not the last step of presolve, should convert bounds back
        # to array and return here

    # Convert modified lb and ub back into N x 2 bounds
    bounds = np.hstack((lb_mod[:, np.newaxis], ub_mod[:, np.newaxis]))

    # remove redundant (linearly dependent) rows from equality constraints
    n_rows_A = A_eq.shape[0]
    redundancy_warning = ("A_eq does not appear to be of full row rank. To "
                          "improve performance, check the problem formulation "
                          "for redundant equality constraints.")
    if (sps.issparse(A_eq)):
        if rr and A_eq.size > 0:  # TODO: Fast sparse rank check?
            rr_res = _remove_redundancy_pivot_sparse(A_eq, b_eq)
            A_eq, b_eq, status, message = rr_res
            if A_eq.shape[0] < n_rows_A:
                warn(redundancy_warning, OptimizeWarning, stacklevel=1)
            if status != 0:
                complete = True
        return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
                c0, x, revstack, complete, status, message)

    # This is a wild guess for which redundancy removal algorithm will be
    # faster. More testing would be good.
    small_nullspace = 5
    if rr and A_eq.size > 0:
        try:  # TODO: use results of first SVD in _remove_redundancy_svd
            rank = np.linalg.matrix_rank(A_eq)
        # oh well, we'll have to go with _remove_redundancy_pivot_dense
        except Exception:
            rank = 0
    if rr and A_eq.size > 0 and rank < A_eq.shape[0]:
        warn(redundancy_warning, OptimizeWarning, stacklevel=3)
        dim_row_nullspace = A_eq.shape[0]-rank
        if rr_method is None:
            if dim_row_nullspace <= small_nullspace:
                rr_res = _remove_redundancy_svd(A_eq, b_eq)
                A_eq, b_eq, status, message = rr_res
            if dim_row_nullspace > small_nullspace or status == 4:
                rr_res = _remove_redundancy_pivot_dense(A_eq, b_eq)
                A_eq, b_eq, status, message = rr_res

        else:
            rr_method = rr_method.lower()
            if rr_method == "svd":
                rr_res = _remove_redundancy_svd(A_eq, b_eq)
                A_eq, b_eq, status, message = rr_res
            elif rr_method == "pivot":
                rr_res = _remove_redundancy_pivot_dense(A_eq, b_eq)
                A_eq, b_eq, status, message = rr_res
            elif rr_method == "id":
                rr_res = _remove_redundancy_id(A_eq, b_eq, rank)
                A_eq, b_eq, status, message = rr_res
            else:  # shouldn't get here; option validity checked above
                pass
        if A_eq.shape[0] < rank:
            message = ("Due to numerical issues, redundant equality "
                       "constraints could not be removed automatically. "
                       "Try providing your constraint matrices as sparse "
                       "matrices to activate sparse presolve, try turning "
                       "off redundancy removal, or try turning off presolve "
                       "altogether.")
            status = 4
        if status != 0:
            complete = True
    return (_LPProblem(c, A_ub, b_ub, A_eq, b_eq, bounds, x0),
            c0, x, revstack, complete, status, message)

################################################################################
# Method Path: scipy.sparse.linalg._expm_multiply._trace
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_expm_multiply.py
# Target:      trace
# Call Count:  564
# Total Exec:  564
# Func Size:   6
################################################################################

def _trace(A):
    # A compatibility function which should eventually disappear.
    if is_pydata_spmatrix(A):
        return A.to_scipy_sparse().trace()
    else:
        return A.trace()

################################################################################
# Method Path: scipy.linalg._decomp_polar.polar
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_decomp_polar.py
# Target:      dot
# Call Count:  560
# Total Exec:  280
# Func Size:   105
################################################################################

def polar(a, side="right"):
    """
    Compute the polar decomposition.

    Returns the factors of the polar decomposition [1]_ `u` and `p` such
    that ``a = up`` (if `side` is "right") or ``a = pu`` (if `side` is
    "left"), where `p` is positive semidefinite. Depending on the shape
    of `a`, either the rows or columns of `u` are orthonormal. When `a`
    is a square array, `u` is a square unitary array. When `a` is not
    square, the "canonical polar decomposition" [2]_ is computed.

    Parameters
    ----------
    a : (m, n) array_like
        The array to be factored.
    side : {'left', 'right'}, optional
        Determines whether a right or left polar decomposition is computed.
        If `side` is "right", then ``a = up``.  If `side` is "left",  then
        ``a = pu``.  The default is "right".

    Returns
    -------
    u : (m, n) ndarray
        If `a` is square, then `u` is unitary. If m > n, then the columns
        of `u` are orthonormal, and if m < n, then the rows of `u` are
        orthonormal.
    p : ndarray
        `p` is Hermitian positive semidefinite. If `a` is nonsingular, `p`
        is positive definite. The shape of `p` is (n, n) or (m, m), depending
        on whether `side` is "right" or "left", respectively.

    References
    ----------
    .. [1] R. A. Horn and C. R. Johnson, "Matrix Analysis", Cambridge
           University Press, 1985.
    .. [2] N. J. Higham, "Functions of Matrices: Theory and Computation",
           SIAM, 2008.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.linalg import polar
    >>> a = np.array([[1, -1], [2, 4]])
    >>> u, p = polar(a)
    >>> u
    array([[ 0.85749293, -0.51449576],
           [ 0.51449576,  0.85749293]])
    >>> p
    array([[ 1.88648444,  1.2004901 ],
           [ 1.2004901 ,  3.94446746]])

    A non-square example, with m < n:

    >>> b = np.array([[0.5, 1, 2], [1.5, 3, 4]])
    >>> u, p = polar(b)
    >>> u
    array([[-0.21196618, -0.42393237,  0.88054056],
           [ 0.39378971,  0.78757942,  0.4739708 ]])
    >>> p
    array([[ 0.48470147,  0.96940295,  1.15122648],
           [ 0.96940295,  1.9388059 ,  2.30245295],
           [ 1.15122648,  2.30245295,  3.65696431]])
    >>> u.dot(p)   # Verify the decomposition.
    array([[ 0.5,  1. ,  2. ],
           [ 1.5,  3. ,  4. ]])
    >>> u.dot(u.T)   # The rows of u are orthonormal.
    array([[  1.00000000e+00,  -2.07353665e-17],
           [ -2.07353665e-17,   1.00000000e+00]])

    Another non-square example, with m > n:

    >>> c = b.T
    >>> u, p = polar(c)
    >>> u
    array([[-0.21196618,  0.39378971],
           [-0.42393237,  0.78757942],
           [ 0.88054056,  0.4739708 ]])
    >>> p
    array([[ 1.23116567,  1.93241587],
           [ 1.93241587,  4.84930602]])
    >>> u.dot(p)   # Verify the decomposition.
    array([[ 0.5,  1.5],
           [ 1. ,  3. ],
           [ 2. ,  4. ]])
    >>> u.T.dot(u)  # The columns of u are orthonormal.
    array([[  1.00000000e+00,  -1.26363763e-16],
           [ -1.26363763e-16,   1.00000000e+00]])

    """
    if side not in ['right', 'left']:
        raise ValueError("`side` must be either 'right' or 'left'")
    a = np.asarray(a)
    if a.ndim != 2:
        raise ValueError("`a` must be a 2-D array.")

    w, s, vh = svd(a, full_matrices=False)
    u = w.dot(vh)
    if side == 'right':
        # a = up
        p = (vh.T.conj() * s).dot(vh)
    else:
        # a = pu
        p = (w * s).dot(w.T.conj())
    return u, p

################################################################################
# Method Path: scipy.linalg._matfuncs_inv_ssq._fractional_matrix_power
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs_inv_ssq.py
# Target:      dot
# Call Count:  555
# Total Exec:  670
# Func Size:   48
################################################################################

def _fractional_matrix_power(A, p):
    """
    Compute the fractional power of a matrix.

    See the fractional_matrix_power docstring in matfuncs.py for more info.

    """
    A = np.asarray(A)
    if len(A.shape) != 2 or A.shape[0] != A.shape[1]:
        raise ValueError('expected a square matrix')
    if p == int(p):
        return np.linalg.matrix_power(A, int(p))
    # Compute singular values.
    s = svdvals(A)
    # Inverse scaling and squaring cannot deal with a singular matrix,
    # because the process of repeatedly taking square roots
    # would not converge to the identity matrix.
    if s[-1]:
        # Compute the condition number relative to matrix inversion,
        # and use this to decide between floor(p) and ceil(p).
        k2 = s[0] / s[-1]
        p1 = p - np.floor(p)
        p2 = p - np.ceil(p)
        if p1 * k2 ** (1 - p1) <= -p2 * k2:
            a = int(np.floor(p))
            b = p1
        else:
            a = int(np.ceil(p))
            b = p2
        try:
            R = _remainder_matrix_power(A, b)
            Q = np.linalg.matrix_power(A, a)
            return Q.dot(R)
        except np.linalg.LinAlgError:
            pass
    # If p is negative then we are going to give up.
    # If p is non-negative then we can fall back to generic funm.
    if p < 0:
        X = np.empty_like(A)
        X.fill(np.nan)
        return X
    else:
        p1 = p - np.floor(p)
        a = int(np.floor(p))
        b = p1
        R, info = funm(A, lambda x: pow(x, b), disp=False)
        Q = np.linalg.matrix_power(A, a)
        return Q.dot(R)

################################################################################
# Method Path: scipy.integrate._ivp.bdf.BdfDenseOutput._call_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/bdf.py
# Target:      dot
# Call Count:  555
# Total Exec:  555
# Func Size:   15
################################################################################

def _call_impl(self, t):
        if t.ndim == 0:
            x = (t - self.t_shift) / self.denom
            p = np.cumprod(x)
        else:
            x = (t - self.t_shift[:, None]) / self.denom[:, None]
            p = np.cumprod(x, axis=0)

        y = np.dot(self.D[1:].T, p)
        if y.ndim == 1:
            y += self.D[0]
        else:
            y += self.D[0, :, None]

        return y

################################################################################
# Method Path: scipy.optimize._optimize._hessp
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  508
# Total Exec:  508
# Func Size:   2
################################################################################

def _hessp(x, p, *args):
            return sf.hess(x).dot(p)

################################################################################
# Method Path: scipy.linalg._solvers.solve_continuous_lyapunov
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_solvers.py
# Target:      dot
# Call Count:  488
# Total Exec:  152
# Func Size:   98
################################################################################

def solve_continuous_lyapunov(a, q):
    """
    Solves the continuous Lyapunov equation :math:`AX + XA^H = Q`.

    Uses the Bartels-Stewart algorithm to find :math:`X`.

    Parameters
    ----------
    a : array_like
        A square matrix

    q : array_like
        Right-hand side square matrix

    Returns
    -------
    x : ndarray
        Solution to the continuous Lyapunov equation

    See Also
    --------
    solve_discrete_lyapunov : computes the solution to the discrete-time
        Lyapunov equation
    solve_sylvester : computes the solution to the Sylvester equation

    Notes
    -----
    The continuous Lyapunov equation is a special form of the Sylvester
    equation, hence this solver relies on LAPACK routine ?TRSYL.

    .. versionadded:: 0.11.0

    Examples
    --------
    Given `a` and `q` solve for `x`:

    >>> import numpy as np
    >>> from scipy import linalg
    >>> a = np.array([[-3, -2, 0], [-1, -1, 0], [0, -5, -1]])
    >>> b = np.array([2, 4, -1])
    >>> q = np.eye(3)
    >>> x = linalg.solve_continuous_lyapunov(a, q)
    >>> x
    array([[ -0.75  ,   0.875 ,  -3.75  ],
           [  0.875 ,  -1.375 ,   5.3125],
           [ -3.75  ,   5.3125, -27.0625]])
    >>> np.allclose(a.dot(x) + x.dot(a.T), q)
    True
    """

    a = np.atleast_2d(_asarray_validated(a, check_finite=True))
    q = np.atleast_2d(_asarray_validated(q, check_finite=True))

    r_or_c = float

    for ind, _ in enumerate((a, q)):
        if np.iscomplexobj(_):
            r_or_c = complex

        if not np.equal(*_.shape):
            raise ValueError(f"Matrix {'aq'[ind]} should be square.")

    # Shape consistency check
    if a.shape != q.shape:
        raise ValueError("Matrix a and q should have the same shape.")

    # Accommodate empty array
    if a.size == 0:
        tdict = {'s': np.float32, 'd': np.float64,
                 'c': np.complex64, 'z': np.complex128}
        func, = get_lapack_funcs(('trsyl',), arrays=(a, q))
        return np.empty(a.shape, dtype=tdict[func.typecode])

    # Compute the Schur decomposition form of a
    r, u = schur(a, output='real')

    # Construct f = u'*q*u
    f = u.conj().T.dot(q.dot(u))

    # Call the Sylvester equation solver
    trsyl = get_lapack_funcs('trsyl', (r, f))

    dtype_string = 'T' if r_or_c is float else 'C'
    y, scale, info = trsyl(r, r, f, tranb=dtype_string)

    if info < 0:
        raise ValueError('?TRSYL exited with the internal error '
                         f'"illegal value in argument number {-info}.". See '
                         'LAPACK documentation for the ?TRSYL error codes.')
    elif info == 1:
        warnings.warn('Input "a" has an eigenvalue pair whose sum is '
                      'very close to or exactly zero. The solution is '
                      'obtained via perturbing the coefficients.',
                      RuntimeWarning, stacklevel=2)
    y *= scale

    return u.dot(y).dot(u.conj().T)

################################################################################
# Method Path: scipy.optimize._nonlin.LowRankMatrix.rmatvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_nonlin.py
# Target:      dot
# Call Count:  470
# Total Exec:  636
# Func Size:   5
################################################################################

def rmatvec(self, v):
        """Evaluate w = M^H v"""
        if self.collapsed is not None:
            return np.dot(self.collapsed.T.conj(), v)
        return LowRankMatrix._matvec(v, np.conj(self.alpha), self.ds, self.cs)

################################################################################
# Method Path: scipy.optimize._linesearch.line_search_wolfe2
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linesearch.py
# Target:      dot
# Call Count:  458
# Total Exec:  458
# Func Size:   149
################################################################################

def line_search_wolfe2(f, myfprime, xk, pk, gfk=None, old_fval=None,
                       old_old_fval=None, args=(), c1=1e-4, c2=0.9, amax=None,
                       extra_condition=None, maxiter=10):
    """Find alpha that satisfies strong Wolfe conditions.

    Parameters
    ----------
    f : callable f(x,*args)
        Objective function.
    myfprime : callable f'(x,*args)
        Objective function gradient.
    xk : ndarray
        Starting point.
    pk : ndarray
        Search direction. The search direction must be a descent direction
        for the algorithm to converge.
    gfk : ndarray, optional
        Gradient value for x=xk (xk being the current parameter
        estimate). Will be recomputed if omitted.
    old_fval : float, optional
        Function value for x=xk. Will be recomputed if omitted.
    old_old_fval : float, optional
        Function value for the point preceding x=xk.
    args : tuple, optional
        Additional arguments passed to objective function.
    c1 : float, optional
        Parameter for Armijo condition rule.
    c2 : float, optional
        Parameter for curvature condition rule.
    amax : float, optional
        Maximum step size
    extra_condition : callable, optional
        A callable of the form ``extra_condition(alpha, x, f, g)``
        returning a boolean. Arguments are the proposed step ``alpha``
        and the corresponding ``x``, ``f`` and ``g`` values. The line search
        accepts the value of ``alpha`` only if this
        callable returns ``True``. If the callable returns ``False``
        for the step length, the algorithm will continue with
        new iterates. The callable is only called for iterates
        satisfying the strong Wolfe conditions.
    maxiter : int, optional
        Maximum number of iterations to perform.

    Returns
    -------
    alpha : float or None
        Alpha for which ``x_new = x0 + alpha * pk``,
        or None if the line search algorithm did not converge.
    fc : int
        Number of function evaluations made.
    gc : int
        Number of gradient evaluations made.
    new_fval : float or None
        New function value ``f(x_new)=f(x0+alpha*pk)``,
        or None if the line search algorithm did not converge.
    old_fval : float
        Old function value ``f(x0)``.
    new_slope : float or None
        The local slope along the search direction at the
        new value ``<myfprime(x_new), pk>``,
        or None if the line search algorithm did not converge.


    Notes
    -----
    Uses the line search algorithm to enforce strong Wolfe conditions. See algorithms
    3.5 and 3.6 on pp. 60-61 in [1]_ (first edition, algorithms 3.2 and 3.3 on
    pp. 59-61), see also [2]_.

    The search direction `pk` must be a descent direction (e.g.
    ``-myfprime(xk)``) to find a step length that satisfies the strong Wolfe
    conditions. If the search direction is not a descent direction (e.g.
    ``myfprime(xk)``), then `alpha`, `new_fval`, and `new_slope` will be None.

    References
    ----------
    .. [1] J. Nocedal and S. J. Wright. "Numerical Optimization". Springer Ser.
       Oper. Res. Financ. Eng. Springer, New York, NY, USA, second edition,
       2006. doi:`10.1007/978-0-387-40065-5`
    .. [2] R. Fletcher, "Practical Methods of Optimization." Wiley, May 23, 2000.
       doi:`10.1002/9781118723203`

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.optimize import line_search

    An objective function and its gradient are defined.

    >>> def obj_func(x):
    ...     return (x[0])**2+(x[1])**2
    >>> def obj_grad(x):
    ...     return [2*x[0], 2*x[1]]

    We can find alpha that satisfies strong Wolfe conditions.

    >>> start_point = np.array([1.8, 1.7])
    >>> search_gradient = np.array([-1.0, -1.0])
    >>> line_search(obj_func, obj_grad, start_point, search_gradient)
    (1.0, 2, 1, 1.1300000000000001, 6.13, [1.6, 1.4])

    """
    fc = [0]
    gc = [0]
    gval = [None]
    gval_alpha = [None]

    def phi(alpha):
        fc[0] += 1
        return f(xk + alpha * pk, *args)

    fprime = myfprime

    def derphi(alpha):
        gc[0] += 1
        gval[0] = fprime(xk + alpha * pk, *args)  # store for later use
        gval_alpha[0] = alpha
        return np.dot(gval[0], pk)

    if gfk is None:
        gfk = fprime(xk, *args)
    derphi0 = np.dot(gfk, pk)

    if extra_condition is not None:
        # Add the current gradient as argument, to avoid needless
        # re-evaluation
        def extra_condition2(alpha, phi):
            if gval_alpha[0] != alpha:
                derphi(alpha)
            x = xk + alpha * pk
            return extra_condition(alpha, x, phi, gval[0])
    else:
        extra_condition2 = None

    alpha_star, phi_star, old_fval, derphi_star = scalar_search_wolfe2(
            phi, derphi, old_fval, old_old_fval, derphi0, c1, c2, amax,
            extra_condition2, maxiter=maxiter)

    if derphi_star is None:
        warn('The line search algorithm did not converge',
             LineSearchWarning, stacklevel=2)
    else:
        # derphi_star is a number (derphi) -- so use the most recently
        # calculated gradient used in computing it derphi = gfk*pk
        # this is the gradient at the next step no need to compute it
        # again in the outer loop.
        derphi_star = gval[0]

    return alpha_star, fc[0], gc[0], phi_star, old_fval, derphi_star

################################################################################
# Method Path: scipy.linalg._expm_frechet._diff_pade7
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_expm_frechet.py
# Target:      dot
# Call Count:  444
# Total Exec:  37
# Func Size:   14
################################################################################

def _diff_pade7(A, E, ident):
    b = (17297280., 8648640., 1995840., 277200., 25200., 1512., 56., 1.)
    A2 = A.dot(A)
    M2 = np.dot(A, E) + np.dot(E, A)
    A4 = np.dot(A2, A2)
    M4 = np.dot(A2, M2) + np.dot(M2, A2)
    A6 = np.dot(A2, A4)
    M6 = np.dot(A4, M2) + np.dot(M4, A2)
    U = A.dot(b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[6]*A6 + b[4]*A4 + b[2]*A2 + b[0]*ident
    Lu = (A.dot(b[7]*M6 + b[5]*M4 + b[3]*M2) +
            E.dot(b[7]*A6 + b[5]*A4 + b[3]*A2 + b[1]*ident))
    Lv = b[6]*M6 + b[4]*M4 + b[2]*M2
    return U, V, Lu, Lv

################################################################################
# Method Path: scipy.optimize._remove_redundancy._remove_redundancy_pivot_dense
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_remove_redundancy.py
# Target:      dot
# Call Count:  441
# Total Exec:  22
# Func Size:   125
################################################################################

def _remove_redundancy_pivot_dense(A, rhs, true_rank=None):
    """
    Eliminates redundant equations from system of equations defined by Ax = b
    and identifies infeasibilities.

    Parameters
    ----------
    A : 2-D array
        An matrix representing the left-hand side of a system of equations
    rhs : 1-D array
        An array representing the right-hand side of a system of equations

    Returns
    -------
    A : 2-D array
        A matrix representing the left-hand side of a system of equations
    rhs : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the system
        0: No infeasibility identified
        2: Trivially infeasible
    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [2] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.

    """
    tolapiv = 1e-8
    tolprimal = 1e-8
    status = 0
    message = ""
    inconsistent = ("There is a linear combination of rows of A_eq that "
                    "results in zero, suggesting a redundant constraint. "
                    "However the same linear combination of b_eq is "
                    "nonzero, suggesting that the constraints conflict "
                    "and the problem is infeasible.")
    A, rhs, status, message = _remove_zero_rows(A, rhs)

    if status != 0:
        return A, rhs, status, message

    m, n = A.shape

    v = list(range(m))      # Artificial column indices.
    b = list(v)             # Basis column indices.
    # This is better as a list than a set because column order of basis matrix
    # needs to be consistent.
    d = []                  # Indices of dependent rows
    perm_r = None

    A_orig = A
    A = np.zeros((m, m + n), order='F')
    np.fill_diagonal(A, 1)
    A[:, m:] = A_orig
    e = np.zeros(m)

    js_candidates = np.arange(m, m+n, dtype=int)  # candidate columns for basis
    # manual masking was faster than masked array
    js_mask = np.ones(js_candidates.shape, dtype=bool)

    # Implements basic algorithm from [2]
    # Uses some of the suggested improvements (removing zero rows and
    # Bartels-Golub update idea).
    # Removing column singletons would be easy, but it is not as important
    # because the procedure is performed only on the equality constraint
    # matrix from the original problem - not on the canonical form matrix,
    # which would have many more column singletons due to slack variables
    # from the inequality constraints.
    # The thoughts on "crashing" the initial basis are only really useful if
    # the matrix is sparse.

    lu = np.eye(m, order='F'), np.arange(m)  # initial LU is trivial
    perm_r = lu[1]
    for i in v:

        e[i] = 1
        if i > 0:
            e[i-1] = 0

        try:  # fails for i==0 and any time it gets ill-conditioned
            j = b[i-1]
            lu = bg_update_dense(lu, perm_r, A[:, j], i-1)
        except Exception:
            lu = scipy.linalg.lu_factor(A[:, b])
            LU, p = lu
            perm_r = list(range(m))
            for i1, i2 in enumerate(p):
                perm_r[i1], perm_r[i2] = perm_r[i2], perm_r[i1]

        pi = scipy.linalg.lu_solve(lu, e, trans=1)

        js = js_candidates[js_mask]
        batch = 50

        # This is a tiny bit faster than looping over columns individually,
        # like for j in js: if abs(A[:,j].transpose().dot(pi)) > tolapiv:
        for j_index in range(0, len(js), batch):
            j_indices = js[j_index: min(j_index+batch, len(js))]

            c = abs(A[:, j_indices].transpose().dot(pi))
            if (c > tolapiv).any():
                j = js[j_index + np.argmax(c)]  # very independent column
                b[i] = j
                js_mask[j-m] = False
                break
        else:
            bibar = pi.T.dot(rhs.reshape(-1, 1))
            bnorm = np.linalg.norm(rhs)
            if abs(bibar)/(1+bnorm) > tolprimal:  # inconsistent
                status = 2
                message = inconsistent
                return A_orig, rhs, status, message
            else:  # dependent
                d.append(i)
                if true_rank is not None and len(d) == m - true_rank:
                    break   # found all redundancies

    keep = set(range(m))
    keep = list(keep - set(d))
    return A_orig[keep, :], rhs[keep], status, message

################################################################################
# Method Path: scipy.integrate._ivp.lsoda.LsodaDenseOutput._call_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/lsoda.py
# Target:      dot
# Call Count:  426
# Total Exec:  426
# Func Size:   7
################################################################################

def _call_impl(self, t):
        if t.ndim == 0:
            x = ((t - self.t) / self.h) ** self.p
        else:
            x = ((t - self.t) / self.h) ** self.p[:, None]

        return np.dot(self.yh, x)

################################################################################
# Method Path: scipy.interpolate._polyint.BarycentricInterpolator._evaluate
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/interpolate/_polyint.py
# Target:      dot
# Call Count:  421
# Total Exec:  439
# Func Size:   18
################################################################################

def _evaluate(self, x):
        if x.size == 0:
            p = np.zeros((0, self.r), dtype=self.dtype)
        else:
            c = x[..., np.newaxis] - self.xi
            z = c == 0
            c[z] = 1
            c = self.wi / c
            with np.errstate(divide='ignore'):
                p = np.dot(c, self.yi) / np.sum(c, axis=-1)[..., np.newaxis]
            # Now fix where x==some xi
            r = np.nonzero(z)
            if len(r) == 1:  # evaluation at a scalar
                if len(r[0]) > 0:  # equals one of the points
                    p = self.yi[r[0][0]]
            else:
                p[r[:-1]] = self.yi[r[-1]]
        return p

################################################################################
# Method Path: scipy.integrate._bvp.stacked_matmul
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_bvp.py
# Target:      dot
# Call Count:  392
# Total Exec:  162
# Func Size:   13
################################################################################

def stacked_matmul(a, b):
    """Stacked matrix multiply: out[i,:,:] = np.dot(a[i,:,:], b[i,:,:]).

    Empirical optimization. Use outer Python loop and BLAS for large
    matrices, otherwise use a single einsum call.
    """
    if a.shape[1] > 50:
        out = np.empty((a.shape[0], a.shape[1], b.shape[2]))
        for i in range(a.shape[0]):
            out[i] = np.dot(a[i], b[i])
        return out
    else:
        return np.einsum('...ij,...jk->...ik', a, b)

################################################################################
# Method Path: scipy.integrate._ivp.rk.DOP853._dense_output_impl
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  368
# Total Exec:  92
# Func Size:   20
################################################################################

def _dense_output_impl(self):
        K = self.K_extended
        h = self.h_previous
        for s, (a, c) in enumerate(zip(self.A_EXTRA, self.C_EXTRA),
                                   start=self.n_stages + 1):
            dy = np.dot(K[:s].T, a[:s]) * h
            K[s] = self.fun(self.t_old + c * h, self.y_old + dy)

        F = np.empty((dop853_coefficients.INTERPOLATOR_POWER, self.n),
                     dtype=self.y_old.dtype)

        f_old = K[0]
        delta_y = self.y - self.y_old

        F[0] = delta_y
        F[1] = h * f_old - delta_y
        F[2] = 2 * delta_y - h * (self.f + f_old)
        F[3:] = h * np.dot(self.D, K)

        return Dop853DenseOutput(self.t_old, self.t, self.y_old, F)

################################################################################
# Method Path: scipy._external.array_api_compat._internal.wrapped_f
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/subprojects/array_api_compat/array_api_compat/array_api_compat/_internal.py
# Target:      @
# Call Count:  362
# Total Exec:  2130950
# Func Size:   3
################################################################################

def wrapped_f(*args: object, **kwargs: object) -> object:
            return f(*args, xp=xp, **kwargs)

################################################################################
# Method Path: scipy.optimize._differentiable_functions.IdentityVectorFunction.__init__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_differentiable_functions.py
# Target:      dot
# Call Count:  353
# Total Exec:  353
# Func Size:   29
################################################################################

def __init__(self, x0, sparse_jacobian):
        n = len(x0)
        if sparse_jacobian or sparse_jacobian is None:
            A = sps.eye_array(n, format='csr')
            sparse_jacobian = True
        else:
            A = np.eye(n)
            sparse_jacobian = False
        super().__init__(A, x0, sparse_jacobian)

################################################################################
# Method Path: scipy.sparse.linalg._eigen.arpack.arpack._UnsymmetricArpackParams.extract
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_eigen/arpack/arpack.py
# Target:      dot
# Call Count:  312
# Total Exec:  803
# Func Size:   143
################################################################################

def extract(self, return_eigenvectors):
        k, n = self.k, self.n

        ierr = 0
        self.arpack_dict['info'] = 0  # Clear, if any, previous error from naupd
        howmny = HOWMNY_DICT['A']  # return all eigenvectors
        sselect = np.zeros(self.ncv, dtype=_int_dtype)
        sigmar = float(np.real(self.sigma))
        sigmai = float(np.imag(self.sigma))
        workev = np.zeros(3 * self.ncv, self.tp)

        if self.tp in 'fd':
            dr = np.zeros([k + 1], dtype=self.tp)
            di = np.zeros([k + 1], dtype=self.tp)
            # Using a Fortran ordered array for NumPy parse the result correctly
            zr = np.zeros([n, k + 1], dtype=self.tp, order='F')

            # ARPACK _neupd call
            self._arpack_extract(
                self.arpack_dict, return_eigenvectors, howmny, sselect, dr, di,
                zr, sigmar, sigmai, workev, self.resid, self.v, self.ipntr,
                self.workd, self.workl
            )

            ierr = self.arpack_dict['info']

            if ierr != 0:
                raise ArpackError(ierr, infodict=self.extract_infodict)
            nreturned = self.arpack_dict['nconv']  # number of good eigs returned

            # Build complex eigenvalues from real and imaginary parts
            d = dr + 1.0j * di

            # Arrange the eigenvectors: complex eigenvectors are stored as
            # real,imaginary in consecutive columns
            z = zr.astype(self.tp.upper())

            # The ARPACK nonsymmetric real and double interface (s,d)naupd
            # return eigenvalues and eigenvectors in real (float,double)
            # arrays.

            # Efficiency: this should check that return_eigenvectors == True
            #  before going through this construction.
            if sigmai == 0:
                i = 0
                while i <= k:
                    # check if complex
                    if abs(d[i].imag) != 0:
                        # this is a complex conjugate pair with eigenvalues
                        # in consecutive columns
                        if i < k:
                            z[:, i] = zr[:, i] + 1.0j * zr[:, i + 1]
                            z[:, i + 1] = z[:, i].conjugate()
                            i += 1
                        else:
                            #last eigenvalue is complex: the imaginary part of
                            # the eigenvector has not been returned
                            #this can only happen if nreturned > k, so we'll
                            # throw out this case.
                            nreturned -= 1
                    i += 1

            else:
                # real matrix, mode 3 or 4, imag(sigma) is nonzero:
                # see remark 3 in <s,d>neupd.f
                # Build complex eigenvalues from real and imaginary parts
                i = 0
                while i <= k:
                    if abs(d[i].imag) == 0:
                        d[i] = np.dot(zr[:, i], self.matvec(zr[:, i]))
                    else:
                        if i < k:
                            z[:, i] = zr[:, i] + 1.0j * zr[:, i + 1]
                            z[:, i + 1] = z[:, i].conjugate()
                            d[i] = ((np.dot(zr[:, i],
                                            self.matvec(zr[:, i]))
                                     + np.dot(zr[:, i + 1],
                                              self.matvec(zr[:, i + 1])))
                                    + 1j * (np.dot(zr[:, i],
                                                   self.matvec(zr[:, i + 1]))
                                            - np.dot(zr[:, i + 1],
                                                     self.matvec(zr[:, i]))))
                            d[i + 1] = d[i].conj()
                            i += 1
                        else:
                            #last eigenvalue is complex: the imaginary part of
                            # the eigenvector has not been returned
                            #this can only happen if nreturned > k, so we'll
                            # throw out this case.
                            nreturned -= 1
                    i += 1

            # Now we have k+1 possible eigenvalues and eigenvectors
            # Return the ones specified by the keyword "which"

            if nreturned <= k:
                # we got less or equal as many eigenvalues we wanted
                d = d[:nreturned]
                z = z[:, :nreturned]
            else:
                # we got one extra eigenvalue (likely a cc pair, but which?)
                if self.mode in (1, 2):
                    rd = d
                elif self.mode in (3, 4):
                    rd = 1 / (d - self.sigma)

                if self.which in ['LR', 'SR']:
                    ind = np.argsort(rd.real)
                elif self.which in ['LI', 'SI']:
                    # for LI,SI ARPACK returns largest,smallest
                    # abs(imaginary) (complex pairs come together)
                    ind = np.argsort(abs(rd.imag))
                else:
                    ind = np.argsort(abs(rd))

                if self.which in ['LR', 'LM', 'LI']:
                    ind = ind[-k:][::-1]
                elif self.which in ['SR', 'SM', 'SI']:
                    ind = ind[:k]

                d = d[ind]
                z = z[:, ind]
        else:
            d = np.zeros([k], dtype=self.tp)
            z = np.zeros([n, k], dtype=self.tp, order='F')
            self._arpack_extract(
                self.arpack_dict, return_eigenvectors, howmny, sselect, d, z,
                self.sigma, workev, self.resid, self.v, self.ipntr, self.workd,
                self.workl, self.rwork)

            ierr = self.arpack_dict['info']

            if ierr != 0:
                raise ArpackError(ierr, infodict=self.extract_infodict)

            k_ok = self.arpack_dict['nconv']
            d = d[:k_ok]
            z = z[:, :k_ok]

        if return_eigenvectors:
            return d, z
        else:
            return d

################################################################################
# Method Path: scipy.optimize._remove_redundancy._remove_redundancy_svd
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_remove_redundancy.py
# Target:      dot
# Call Count:  268
# Total Exec:  66
# Func Size:   90
################################################################################

def _remove_redundancy_svd(A, b):
    """
    Eliminates redundant equations from system of equations defined by Ax = b
    and identifies infeasibilities.

    Parameters
    ----------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations

    Returns
    -------
    A : 2-D array
        An array representing the left-hand side of a system of equations
    b : 1-D array
        An array representing the right-hand side of a system of equations
    status: int
        An integer indicating the status of the system
        0: No infeasibility identified
        2: Trivially infeasible
    message : str
        A string descriptor of the exit status of the optimization.

    References
    ----------
    .. [2] Andersen, Erling D. "Finding all linearly dependent rows in
           large-scale linear programming." Optimization Methods and Software
           6.3 (1995): 219-227.

    """

    A, b, status, message = _remove_zero_rows(A, b)

    if status != 0:
        return A, b, status, message

    U, s, Vh = svd(A)
    eps = np.finfo(float).eps
    tol = s.max() * max(A.shape) * eps

    m, n = A.shape
    s_min = s[-1] if m <= n else 0

    # this algorithm is faster than that of [2] when the nullspace is small
    # but it could probably be improvement by randomized algorithms and with
    # a sparse implementation.
    # it relies on repeated singular value decomposition to find linearly
    # dependent rows (as identified by columns of U that correspond with zero
    # singular values). Unfortunately, only one row can be removed per
    # decomposition (I tried otherwise; doing so can cause problems.)
    # It would be nice if we could do truncated SVD like sp.sparse.linalg.svds
    # but that function is unreliable at finding singular values near zero.
    # Finding max eigenvalue L of A A^T, then largest eigenvalue (and
    # associated eigenvector) of -A A^T + L I (I is identity) via power
    # iteration would also work in theory, but is only efficient if the
    # smallest nonzero eigenvalue of A A^T is close to the largest nonzero
    # eigenvalue.

    while abs(s_min) < tol:
        v = U[:, -1]  # TODO: return these so user can eliminate from problem?
        # rows need to be represented in significant amount
        eligibleRows = np.abs(v) > tol * 10e6
        if not np.any(eligibleRows) or np.any(np.abs(v.dot(A)) > tol):
            status = 4
            message = ("Due to numerical issues, redundant equality "
                       "constraints could not be removed automatically. "
                       "Try providing your constraint matrices as sparse "
                       "matrices to activate sparse presolve, try turning "
                       "off redundancy removal, or try turning off presolve "
                       "altogether.")
            break
        if np.any(np.abs(v.dot(b)) > tol * 100):  # factor of 100 to fix 10038 and 10349
            status = 2
            message = ("There is a linear combination of rows of A_eq that "
                       "results in zero, suggesting a redundant constraint. "
                       "However the same linear combination of b_eq is "
                       "nonzero, suggesting that the constraints conflict "
                       "and the problem is infeasible.")
            break

        i_remove = _get_densest(A, eligibleRows)
        A = np.delete(A, i_remove, axis=0)
        b = np.delete(b, i_remove)
        U, s, Vh = svd(A)
        m, n = A.shape
        s_min = s[-1] if m <= n else 0

    return A, b, status, message

################################################################################
# Method Path: scipy.linalg._decomp_svd.subspace_angles
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_decomp_svd.py
# Target:      dot
# Call Count:  236
# Total Exec:  121
# Func Size:   106
################################################################################

def subspace_angles(A, B):
    r"""
    Compute the subspace angles between two matrices.

    Parameters
    ----------
    A : (M, N) array_like
        The first input array.
    B : (M, K) array_like
        The second input array.

    Returns
    -------
    angles : ndarray, shape (min(N, K),)
        The subspace angles between the column spaces of `A` and `B` in
        descending order.

    See Also
    --------
    orth
    svd

    Notes
    -----
    This computes the subspace angles according to the formula
    provided in [1]_. For equivalence with MATLAB and Octave behavior,
    use ``angles[0]``.

    .. versionadded:: 1.0

    References
    ----------
    .. [1] Knyazev A, Argentati M (2002) Principal Angles between Subspaces
           in an A-Based Scalar Product: Algorithms and Perturbation
           Estimates. SIAM J. Sci. Comput. 23:2008-2040.

    Examples
    --------
    A Hadamard matrix, which has orthogonal columns, so we expect that
    the subspace angle to be :math:`\frac{\pi}{2}`:

    >>> import numpy as np
    >>> from scipy.linalg import hadamard, subspace_angles
    >>> rng = np.random.default_rng()
    >>> H = hadamard(4)
    >>> print(H)
    [[ 1  1  1  1]
     [ 1 -1  1 -1]
     [ 1  1 -1 -1]
     [ 1 -1 -1  1]]
    >>> np.rad2deg(subspace_angles(H[:, :2], H[:, 2:]))
    array([ 90.,  90.])

    And the subspace angle of a matrix to itself should be zero:

    >>> subspace_angles(H[:, :2], H[:, :2]) <= 2 * np.finfo(float).eps
    array([ True,  True], dtype=bool)

    The angles between non-orthogonal subspaces are in between these extremes:

    >>> x = rng.standard_normal((4, 3))
    >>> np.rad2deg(subspace_angles(x[:, :2], x[:, [2]]))
    array([ 55.832])  # random
    """
    # Steps here omit the U and V calculation steps from the paper

    # 1. Compute orthonormal bases of column-spaces
    A = _asarray_validated(A, check_finite=True)
    if len(A.shape) != 2:
        raise ValueError(f'expected 2D array, got shape {A.shape}')
    QA = orth(A)
    del A

    B = _asarray_validated(B, check_finite=True)
    if len(B.shape) != 2:
        raise ValueError(f'expected 2D array, got shape {B.shape}')
    if len(B) != len(QA):
        raise ValueError('A and B must have the same number of rows, got '
                         f'{QA.shape[0]} and {B.shape[0]}')
    QB = orth(B)
    del B

    # 2. Compute SVD for cosine
    QA_H_QB = np.dot(QA.T.conj(), QB)
    sigma = svdvals(QA_H_QB)

    # 3. Compute matrix B
    if QA.shape[1] >= QB.shape[1]:
        B = QB - np.dot(QA, QA_H_QB)
    else:
        B = QA - np.dot(QB, QA_H_QB.T.conj())
    del QA, QB, QA_H_QB

    # 4. Compute SVD for sine
    mask = sigma ** 2 >= 0.5
    if mask.any():
        mu_arcsin = np.arcsin(np.clip(svdvals(B, overwrite_a=True), -1., 1.))
    else:
        mu_arcsin = 0.

    # 5. Compute the principal angles
    # with reverse ordering of sigma because smallest sigma belongs to largest
    # angle theta
    theta = np.where(mask, mu_arcsin, np.arccos(np.clip(sigma[::-1], -1., 1.)))
    return theta

################################################################################
# Method Path: scipy.integrate._bvp.solve_newton
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_bvp.py
# Target:      dot
# Call Count:  234
# Total Exec:  51
# Func Size:   155
################################################################################

def solve_newton(n, m, h, col_fun, bc, jac, y, p, B, bvp_tol, bc_tol):
    """Solve the nonlinear collocation system by a Newton method.

    This is a simple Newton method with a backtracking line search. As
    advised in [1]_, an affine-invariant criterion function F = ||J^-1 r||^2
    is used, where J is the Jacobian matrix at the current iteration and r is
    the vector or collocation residuals (values of the system lhs).

    The method alters between full Newton iterations and the fixed-Jacobian
    iterations based

    There are other tricks proposed in [1]_, but they are not used as they
    don't seem to improve anything significantly, and even break the
    convergence on some test problems I tried.

    All important parameters of the algorithm are defined inside the function.

    Parameters
    ----------
    n : int
        Number of equations in the ODE system.
    m : int
        Number of nodes in the mesh.
    h : ndarray, shape (m-1,)
        Mesh intervals.
    col_fun : callable
        Function computing collocation residuals.
    bc : callable
        Function computing boundary condition residuals.
    jac : callable
        Function computing the Jacobian of the whole system (including
        collocation and boundary condition residuals). It must return
        scipy.sparse in CSC format ready for use with `scipy.sparse.linalg.splu`.
    y : ndarray, shape (n, m)
        Initial guess for the function values at the mesh nodes.
    p : ndarray, shape (k,)
        Initial guess for the unknown parameters.
    B : ndarray with shape (n, n) or None
        Matrix to force the S y(a) = 0 condition for a problems with the
        singular term. If None, the singular term is assumed to be absent.
    bvp_tol : float
        Tolerance to which we want to solve a BVP.
    bc_tol : float
        Tolerance to which we want to satisfy the boundary conditions.

    Returns
    -------
    y : ndarray, shape (n, m)
        Final iterate for the function values at the mesh nodes.
    p : ndarray, shape (k,)
        Final iterate for the unknown parameters.
    singular : bool
        True, if the LU decomposition failed because Jacobian turned out
        to be singular.

    References
    ----------
    .. [1]  U. Ascher, R. Mattheij and R. Russell "Numerical Solution of
       Boundary Value Problems for Ordinary Differential Equations",
       Philidelphia, PA: Society for Industrial and Applied Mathematics,
       1995.
    """
    # We know that the solution residuals at the middle points of the mesh
    # are connected with collocation residuals  r_middle = 1.5 * col_res / h.
    # As our BVP solver tries to decrease relative residuals below a certain
    # tolerance, it seems reasonable to terminated Newton iterations by
    # comparison of r_middle / (1 + np.abs(f_middle)) with a certain threshold,
    # which we choose to be 1.5 orders lower than the BVP tolerance. We rewrite
    # the condition as col_res < tol_r * (1 + np.abs(f_middle)), then tol_r
    # should be computed as follows:
    tol_r = 2/3 * h * 5e-2 * bvp_tol

    # Maximum allowed number of Jacobian evaluation and factorization, in
    # other words, the maximum number of full Newton iterations. A small value
    # is recommended in the literature.
    max_njev = 4

    # Maximum number of iterations, considering that some of them can be
    # performed with the fixed Jacobian. In theory, such iterations are cheap,
    # but it's not that simple in Python.
    max_iter = 8

    # Minimum relative improvement of the criterion function to accept the
    # step (Armijo constant).
    sigma = 0.2

    # Step size decrease factor for backtracking.
    tau = 0.5

    # Maximum number of backtracking steps, the minimum step is then
    # tau ** n_trial.
    n_trial = 4

    col_res, y_middle, f, f_middle = col_fun(y, p)
    bc_res = bc(y[:, 0], y[:, -1], p)
    res = np.hstack((col_res.ravel(order='F'), bc_res))

    njev = 0
    singular = False
    recompute_jac = True
    for iteration in range(max_iter):
        if recompute_jac:
            J = jac(y, p, y_middle, f, f_middle, bc_res)
            njev += 1
            try:
                LU = splu(J)
            except RuntimeError:
                singular = True
                break

            step = LU.solve(res)
            cost = np.dot(step, step)

        y_step = step[:m * n].reshape((n, m), order='F')
        p_step = step[m * n:]

        alpha = 1
        for trial in range(n_trial + 1):
            y_new = y - alpha * y_step
            if B is not None:
                y_new[:, 0] = np.dot(B, y_new[:, 0])
            p_new = p - alpha * p_step

            col_res, y_middle, f, f_middle = col_fun(y_new, p_new)
            bc_res = bc(y_new[:, 0], y_new[:, -1], p_new)
            res = np.hstack((col_res.ravel(order='F'), bc_res))

            step_new = LU.solve(res)
            cost_new = np.dot(step_new, step_new)
            if cost_new < (1 - 2 * alpha * sigma) * cost:
                break

            if trial < n_trial:
                alpha *= tau

        y = y_new
        p = p_new

        if njev == max_njev:
            break

        if (np.all(np.abs(col_res) < tol_r * (1 + np.abs(f_middle))) and
                np.all(np.abs(bc_res) < bc_tol)):
            break

        # If the full step was taken, then we are going to continue with
        # the same Jacobian. This is the approach of BVP_SOLVER.
        if alpha == 1:
            step = step_new
            cost = cost_new
            recompute_jac = False
        else:
            recompute_jac = True

    return y, p, singular

################################################################################
# Method Path: scipy.integrate._ivp.rk.DOP853._estimate_error_norm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  226
# Total Exec:  113
# Func Size:   9
################################################################################

def _estimate_error_norm(self, K, h, scale):
        err5 = np.dot(K.T, self.E5) / scale
        err3 = np.dot(K.T, self.E3) / scale
        err5_norm_2 = np.linalg.norm(err5)**2
        err3_norm_2 = np.linalg.norm(err3)**2
        if err5_norm_2 == 0 and err3_norm_2 == 0:
            return 0.0
        denom = err5_norm_2 + 0.01 * err3_norm_2
        return np.abs(h) * err5_norm_2 / np.sqrt(denom * len(scale))

################################################################################
# Method Path: scipy.signal._spectral_py.lombscargle
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_spectral_py.py
# Target:      dot
# Call Count:  199
# Total Exec:  38
# Func Size:   344
################################################################################

def lombscargle(
    x: npt.ArrayLike,
    y: npt.ArrayLike,
    freqs: npt.ArrayLike,
    *,
    precenter: bool = _NoValue,  # type:ignore[assignment]
    normalize: bool | Literal["power", "normalize", "amplitude"] = False,
    weights: npt.NDArray | None = None,
    floating_mean: bool = False,
) -> npt.NDArray:
    """
    Compute the generalized Lomb-Scargle periodogram.

    The Lomb-Scargle periodogram was developed by Lomb [1]_ and further
    extended by Scargle [2]_ to find, and test the significance of weak
    periodic signals with uneven temporal sampling. The algorithm used
    here is based on a weighted least-squares fit of the form
    ``y(ω) = a*cos(ω*x) + b*sin(ω*x) + c``, where the fit is calculated for
    each frequency independently. This algorithm was developed by Zechmeister
    and Kürster which improves the Lomb-Scargle periodogram by enabling
    the weighting of individual samples and calculating an unknown y offset
    (also called a "floating-mean" model) [3]_. For more details, and practical
    considerations, see the excellent reference on the Lomb-Scargle periodogram [4]_.

    When *normalize* is False (or "power") (default) the computed periodogram
    is unnormalized, it takes the value ``(A**2) * N/4`` for a harmonic
    signal with amplitude A for sufficiently large N. Where N is the length of x or y.

    When *normalize* is True (or "normalize") the computed periodogram is normalized
    by the residuals of the data around a constant reference model (at zero).

    When *normalize* is "amplitude" the computed periodogram is the complex
    representation of the amplitude and phase.

    Input arrays should be 1-D of a real floating data type, which are converted into
    float64 arrays before processing.

    Parameters
    ----------
    x : array_like
        Sample times.
    y : array_like
        Measurement values. Values are assumed to have a baseline of ``y = 0``. If
        there is a possibility of a y offset, it is recommended to set `floating_mean`
        to True.
    freqs : array_like
        Angular frequencies (e.g., having unit rad/s=2π/s for `x` having unit s) for
        output periodogram. Frequencies are normally >= 0, as any peak at ``-freq`` will
        also exist at ``+freq``.
    precenter : bool, optional
        Pre-center measurement values by subtracting the mean, if True. This is
        a legacy parameter and unnecessary if `floating_mean` is True.

        .. deprecated:: 1.17.0
            The `precenter` argument is deprecated and will be removed in SciPy 1.19.0.
            The functionality can be substituted by passing ``y - y.mean()`` to `y`.

    normalize : bool | str, optional
        Compute normalized or complex (amplitude + phase) periodogram.
        Valid options are: ``False``/``"power"``, ``True``/``"normalize"``, or
        ``"amplitude"``.
    weights : array_like, optional
        Weights for each sample. Weights must be nonnegative.
    floating_mean : bool, optional
        Determines a y offset for each frequency independently, if True.
        Else the y offset is assumed to be `0`.

    Returns
    -------
    pgram : array_like
        Lomb-Scargle periodogram.

    Raises
    ------
    ValueError
        If any of the input arrays x, y, freqs, or weights are not 1D, or if any are
        zero length. Or, if the input arrays x, y, and weights do not have the same
        shape as each other.
    ValueError
        If any weight is < 0, or the sum of the weights is <= 0.
    ValueError
        If the normalize parameter is not one of the allowed options.

    See Also
    --------
    periodogram: Power spectral density using a periodogram
    welch: Power spectral density by Welch's method
    csd: Cross spectral density by Welch's method

    Notes
    -----
    The algorithm used will not automatically account for any unknown y offset, unless
    `floating_mean` is ``True``. Therefore, for most use cases, if there is a
    possibility of a y offset, it is recommended to set `floating_mean` to ``True``.
    Furthermore, `floating_mean` accounts for sample weights, and will also correct for
    any bias due to consistently missing observations at peaks and/or troughs.

    The legacy concept of "pre-centering" entails removing the mean from parameter `y`
    before processing, i.e., passing ``y - y.mean()`` instead of setting the parameter
    `floating_mean` to ``True``.

    When the normalize parameter is "amplitude", for any frequency in freqs that is
    below ``(2*pi)/(x.max() - x.min())``, the predicted amplitude will tend towards
    infinity. The concept of a "Nyquist frequency" limit (see Nyquist-Shannon sampling
    theorem) is not generally applicable to unevenly sampled data. Therefore, with
    unevenly sampled data, valid frequencies in freqs can often be much higher than
    expected for those familiar with methods like FFT.

    References
    ----------
    .. [1] N.R. Lomb "Least-squares frequency analysis of unequally spaced
           data", Astrophysics and Space Science, vol 39, pp. 447-462, 1976
           :doi:`10.1007/bf00648343`

    .. [2] J.D. Scargle "Studies in astronomical time series analysis. II -
           Statistical aspects of spectral analysis of unevenly spaced data",
           The Astrophysical Journal, vol 263, pp. 835-853, 1982
           :doi:`10.1086/160554`

    .. [3] M. Zechmeister and M. Kürster, "The generalised Lomb-Scargle periodogram.
           A new formalism for the floating-mean and Keplerian periodograms,"
           Astronomy and Astrophysics, vol. 496, pp. 577-584, 2009
           :doi:`10.1051/0004-6361:200811296`

    .. [4] J.T. VanderPlas, "Understanding the Lomb-Scargle Periodogram,"
           The Astrophysical Journal Supplement Series, vol. 236, no. 1, p. 16,
           May 2018
           :doi:`10.3847/1538-4365/aab766`


    Examples
    --------
    >>> import numpy as np
    >>> rng = np.random.default_rng()

    First define some input parameters for the signal:

    >>> A = 2.  # amplitude
    >>> c = 2.  # offset
    >>> w0 = 1.  # rad/sec
    >>> nin = 150
    >>> nout = 1002

    Randomly generate sample times:

    >>> x = rng.uniform(0, 10*np.pi, nin)

    Plot a sine wave for the selected times:

    >>> y = A * np.cos(w0*x) + c

    Define the array of frequencies for which to compute the periodogram:

    >>> w = np.linspace(0.25, 10, nout)

    Calculate Lomb-Scargle periodogram for each of the normalize options:

    >>> from scipy.signal import lombscargle
    >>> pgram_power = lombscargle(x, y, w, normalize=False)
    >>> pgram_norm = lombscargle(x, y, w, normalize=True)
    >>> pgram_amp = lombscargle(x, y, w, normalize='amplitude')
    ...
    >>> pgram_power_f = lombscargle(x, y, w, normalize=False, floating_mean=True)
    >>> pgram_norm_f = lombscargle(x, y, w, normalize=True, floating_mean=True)
    >>> pgram_amp_f = lombscargle(x, y, w, normalize='amplitude', floating_mean=True)

    Now make a plot of the input data:

    >>> import matplotlib.pyplot as plt
    >>> fig, (ax_t, ax_p, ax_n, ax_a) = plt.subplots(4, 1, figsize=(5, 6))
    >>> ax_t.plot(x, y, 'b+')
    >>> ax_t.set_xlabel('Time [s]')
    >>> ax_t.set_ylabel('Amplitude')

    Then plot the periodogram for each of the normalize options, as well as with and
    without floating_mean=True:

    >>> ax_p.plot(w, pgram_power, label='default')
    >>> ax_p.plot(w, pgram_power_f, label='floating_mean=True')
    >>> ax_p.set_xlabel('Angular frequency [rad/s]')
    >>> ax_p.set_ylabel('Power')
    >>> ax_p.legend(prop={'size': 7})
    ...
    >>> ax_n.plot(w, pgram_norm, label='default')
    >>> ax_n.plot(w, pgram_norm_f, label='floating_mean=True')
    >>> ax_n.set_xlabel('Angular frequency [rad/s]')
    >>> ax_n.set_ylabel('Normalized')
    >>> ax_n.legend(prop={'size': 7})
    ...
    >>> ax_a.plot(w, np.abs(pgram_amp), label='default')
    >>> ax_a.plot(w, np.abs(pgram_amp_f), label='floating_mean=True')
    >>> ax_a.set_xlabel('Angular frequency [rad/s]')
    >>> ax_a.set_ylabel('Amplitude')
    >>> ax_a.legend(prop={'size': 7})
    ...
    >>> plt.tight_layout()
    >>> plt.show()

    """

    # if no weights are provided, assume all data points are equally important
    if weights is None:
        weights = np.ones_like(y, dtype=np.float64)
    else:
        # if provided, make sure weights is an array and cast to float64
        weights = np.asarray(weights, dtype=np.float64)

    # make sure other inputs are arrays and cast to float64
    # done before validation, in case they were not arrays
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    freqs = np.asarray(freqs, dtype=np.float64)

    # validate input shapes
    if not (x.ndim == 1 and x.size > 0 and x.shape == y.shape == weights.shape):
        raise ValueError("Parameters x, y, weights must be 1-D arrays of "
                         "equal non-zero length!")
    if not (freqs.ndim == 1 and freqs.size > 0):
        raise ValueError("Parameter freqs must be a 1-D array of non-zero length!")

    # validate weights
    if not (np.all(weights >= 0) and np.sum(weights) > 0):
        raise ValueError("Parameter weights must have only non-negative entries "
                         "which sum to a positive value!")

    # validate normalize parameter
    if isinstance(normalize, bool):
        # if bool, convert to str literal
        normalize = "normalize" if normalize else "power"

    if normalize not in ["power", "normalize", "amplitude"]:
        raise ValueError(
            "Normalize must be: False (or 'power'), True (or 'normalize'), "
            "or 'amplitude'."
        )

    # weight vector must sum to 1
    weights = weights * (1.0 / weights.sum())

    # if requested, perform precenter
    if precenter is not _NoValue:
        msg = ("Use of parameter 'precenter' is deprecated as of SciPy 1.17.0 and "
               "will be removed in 1.19.0. Please leave 'precenter' unspecified. "
               "Passing True to 'precenter' "
               "can be exactly substituted by passing 'y = (y - y.mean())' into "
               "the input. Consider setting `floating_mean` to True instead.")
        warnings.warn(msg, DeprecationWarning, stacklevel=2)
        if precenter:
            y = y - y.mean()

    # transform arrays
    # row vector
    freqs = freqs.reshape(1, -1)
    # column vectors
    x = x.reshape(-1, 1)
    y = y.reshape(-1, 1)  # type:ignore[union-attr]
    weights = weights.reshape(-1, 1)

    # store frequent intermediates
    weights_y = weights * y
    freqst = freqs * x
    coswt = np.cos(freqst)
    sinwt = np.sin(freqst)

    Y = np.dot(weights.T, y)  # Eq. 7
    CC = np.dot(weights.T, coswt * coswt)  # Eq. 13
    SS = 1.0 - CC  # trig identity: S^2 = 1 - C^2  Eq.14
    CS = np.dot(weights.T, coswt * sinwt)  # Eq. 15

    if floating_mean:
        C = np.dot(weights.T, coswt)  # Eq. 8
        S = np.dot(weights.T, sinwt)  # Eq. 9
        CC -= C * C  # Eq. 13
        SS -= S * S  # Eq. 14
        CS -= C * S  # Eq. 15

    # calculate tau (phase offset to eliminate CS variable)
    tau = 0.5 * np.arctan2(2.0 * CS, CC - SS)  # Eq. 19
    freqst_tau = freqst - tau

    # coswt and sinwt are now offset by tau, which eliminates CS
    coswt_tau = np.cos(freqst_tau)
    sinwt_tau = np.sin(freqst_tau)

    YC = np.dot(weights_y.T, coswt_tau)  # Eq. 11
    YS = np.dot(weights_y.T, sinwt_tau)  # Eq. 12
    CC = np.dot(weights.T, coswt_tau * coswt_tau)  # Eq. 13, CC range is [0.5, 1.0]
    SS = 1.0 - CC  # trig identity: S^2 = 1 - C^2    Eq. 14, SS range is [0.0, 0.5]

    if floating_mean:
        C = np.dot(weights.T, coswt_tau)  # Eq. 8
        S = np.dot(weights.T, sinwt_tau)  # Eq. 9
        YC -= Y * C  # Eq. 11
        YS -= Y * S  # Eq. 12
        CC -= C * C  # Eq. 13, CC range is now [0.0, 1.0]
        SS -= S * S  # Eq. 14, SS range is now [0.0, 0.5]

    # to prevent division by zero errors with a and b, as well as correcting for
    # numerical precision errors that lead to CC or SS being approximately -0.0,
    # make sure CC and SS are both > 0
    epsneg = np.finfo(dtype=y.dtype).epsneg  # type:ignore[union-attr]
    CC[CC < epsneg] = epsneg
    SS[SS < epsneg] = epsneg

    # calculate a and b
    # where: y(w) = a*cos(w) + b*sin(w) + c
    a = YC / CC  # Eq. A.4 and 6, eliminating CS
    b = YS / SS  # Eq. A.4 and 6, eliminating CS
    # c = Y - a * C - b * S

    # store final value as power in A^2 (i.e., (y units)^2)
    pgram = 2.0 * (a * YC + b * YS)

    # squeeze back to a vector
    pgram = np.squeeze(pgram)

    if normalize == "power":  # (default)
        # return the legacy power units ((A**2) * N/4)

        pgram *= float(x.shape[0]) / 4.0

    elif normalize == "normalize":
        # return the normalized power (power at current frequency wrt the entire signal)
        # range will be [0, 1]

        YY = np.dot(weights_y.T, y)  # Eq. 10
        if floating_mean:
            YY -= Y * Y  # Eq. 10

        pgram *= 0.5 / np.squeeze(YY)  # Eq. 20

    else:  # normalize == "amplitude":
        # return the complex representation of the best-fit amplitude and phase

        # squeeze back to vectors
        a = np.squeeze(a)
        b = np.squeeze(b)
        tau = np.squeeze(tau)

        # calculate the complex representation, and correct for tau rotation
        pgram = (a + 1j * b) * np.exp(1j * tau)

    return pgram

################################################################################
# Method Path: scipy.linalg._expm_frechet._diff_pade5
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_expm_frechet.py
# Target:      dot
# Call Count:  189
# Total Exec:  21
# Func Size:   12
################################################################################

def _diff_pade5(A, E, ident):
    b = (30240., 15120., 3360., 420., 30., 1.)
    A2 = A.dot(A)
    M2 = np.dot(A, E) + np.dot(E, A)
    A4 = np.dot(A2, A2)
    M4 = np.dot(A2, M2) + np.dot(M2, A2)
    U = A.dot(b[5]*A4 + b[3]*A2 + b[1]*ident)
    V = b[4]*A4 + b[2]*A2 + b[0]*ident
    Lu = (A.dot(b[5]*M4 + b[3]*M2) +
            E.dot(b[5]*A4 + b[3]*A2 + b[1]*ident))
    Lv = b[4]*M4 + b[2]*M2
    return U, V, Lu, Lv

################################################################################
# Method Path: scipy.optimize._lsq.lsq_linear.lsq_linear
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/lsq_linear.py
# Target:      dot
# Call Count:  162
# Total Exec:  214
# Func Size:   327
################################################################################

def lsq_linear(A, b, bounds=(-np.inf, np.inf), method='trf', tol=1e-10,
               lsq_solver=None, lsmr_tol=None, max_iter=None,
               verbose=0, *, lsmr_maxiter=None,):
    r"""Solve a linear least-squares problem with bounds on the variables.

    Given an m-by-n design matrix A and a target vector b with m elements,
    `lsq_linear` solves the following optimization problem::

        minimize 0.5 * ||A x - b||**2
        subject to lb <= x <= ub

    This optimization problem is convex, hence a found minimum (if iterations
    have converged) is guaranteed to be global.

    Parameters
    ----------
    A : array_like, sparse array or LinearOperator, shape (m, n)
        Design matrix. Can be `scipy.sparse.linalg.LinearOperator`.
    b : array_like, shape (m,)
        Target vector.
    bounds : 2-tuple of array_like or `Bounds`, optional
        Lower and upper bounds on parameters. Defaults to no bounds.
        There are two ways to specify the bounds:

        - Instance of `Bounds` class.
        - 2-tuple of array_like: Each element of the tuple must be either
          an array with the length equal to the number of parameters, or a
          scalar (in which case the bound is taken to be the same for all
          parameters). Use ``np.inf`` with an appropriate sign to disable
          bounds on all or some parameters.

    method : 'trf' or 'bvls', optional
        Method to perform minimization.

        * 'trf' : Trust Region Reflective algorithm adapted for a linear
          least-squares problem. This is an interior-point-like method
          and the required number of iterations is weakly correlated with
          the number of variables.
        * 'bvls' : Bounded-variable least-squares algorithm. This is
          an active set method, which requires the number of iterations
          comparable to the number of variables. Can't be used when `A` is
          sparse or LinearOperator.

        Default is 'trf'.
    tol : float, optional
        Tolerance parameter. The algorithm terminates if a relative change
        of the cost function is less than `tol` on the last iteration.
        Additionally, the first-order optimality measure is considered:

        * ``method='trf'`` terminates if the uniform norm of the gradient,
          scaled to account for the presence of the bounds, is less than
          `tol`.
        * ``method='bvls'`` terminates if Karush-Kuhn-Tucker conditions
          are satisfied within `tol` tolerance.

    lsq_solver : {None, 'exact', 'lsmr'}, optional
        Method of solving unbounded least-squares problems throughout
        iterations:

        * 'exact' : Use dense QR or SVD decomposition approach. Can't be
          used when `A` is sparse or LinearOperator.
        * 'lsmr' : Use `scipy.sparse.linalg.lsmr` iterative procedure
          which requires only matrix-vector product evaluations. Can't
          be used with ``method='bvls'``.

        If None (default), the solver is chosen based on type of `A`.
    lsmr_tol : None, float or 'auto', optional
        Tolerance parameters 'atol' and 'btol' for `scipy.sparse.linalg.lsmr`
        If None (default), it is set to ``1e-2 * tol``. If 'auto', the
        tolerance will be adjusted based on the optimality of the current
        iterate, which can speed up the optimization process, but is not always
        reliable.
    max_iter : None or int, optional
        Maximum number of iterations before termination. If None (default), it
        is set to 100 for ``method='trf'`` or to the number of variables for
        ``method='bvls'`` (not counting iterations for 'bvls' initialization).
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

        * 0 : work silently (default).
        * 1 : display a termination report.
        * 2 : display progress during iterations.

    lsmr_maxiter : None or int, optional
        Maximum number of iterations for the lsmr least squares solver,
        if it is used (by setting ``lsq_solver='lsmr'``). If None (default), it
        uses lsmr's default of ``min(m, n)`` where ``m`` and ``n`` are the
        number of rows and columns of `A`, respectively. Has no effect if
        ``lsq_solver='exact'``.

    Returns
    -------
    result : OptimizeResult
        Result object with the following fields defined:

        x : ndarray, shape (n,)
            Solution found.
        cost : float
            Value of the cost function at the solution.
        fun : ndarray, shape (m,)
            Vector of residuals at the solution.
        optimality : float
            First-order optimality measure. The exact meaning depends on `method`,
            refer to the description of `tol` parameter.
        active_mask : ndarray of int, shape (n,)
            Each component shows whether a corresponding constraint is active
            (that is, whether a variable is at the bound):

            *  0 : a constraint is not active.
            * -1 : a lower bound is active.
            *  1 : an upper bound is active.

            Might be somewhat arbitrary for the `trf` method as it generates a
            sequence of strictly feasible iterates and active_mask is determined
            within a tolerance threshold.
        unbounded_sol : tuple
            Unbounded least squares solution tuple returned by the least squares
            solver (set with `lsq_solver` option). If `lsq_solver` is not set or is
            set to ``'exact'``, the tuple contains an ndarray of shape (n,) with
            the unbounded solution, an ndarray with the sum of squared residuals,
            an int with the rank of `A`, and an ndarray with the singular values
            of `A` (see NumPy's ``linalg.lstsq`` for more information). If
            `lsq_solver` is set to ``'lsmr'``, the tuple contains an ndarray of
            shape (n,) with the unbounded solution, an int with the exit code,
            an int with the number of iterations, and five floats with
            various norms and the condition number of `A` (see SciPy's
            ``sparse.linalg.lsmr`` for more information). This output can be
            useful for determining the convergence of the least squares solver,
            particularly the iterative ``'lsmr'`` solver. The unbounded least
            squares problem is to minimize ``0.5 * ||A x - b||**2``.
        nit : int
            Number of iterations. Zero if the unconstrained solution is optimal.
        status : int
            Reason for algorithm termination:

            * -1 : the algorithm was not able to make progress on the last
              iteration.
            *  0 : the maximum number of iterations is exceeded.
            *  1 : the first-order optimality measure is less than `tol`.
            *  2 : the relative change of the cost function is less than `tol`.
            *  3 : the unconstrained solution is optimal.

        message : str
            Verbal description of the termination reason.
        success : bool
            True if one of the convergence criteria is satisfied (`status` > 0).

    See Also
    --------
    nnls : Linear least squares with non-negativity constraint.
    least_squares : Nonlinear least squares with bounds on the variables.

    Notes
    -----
    The algorithm first computes the unconstrained least-squares solution by
    `numpy.linalg.lstsq` or `scipy.sparse.linalg.lsmr` depending on
    `lsq_solver`. This solution is returned as optimal if it lies within the
    bounds.

    Method 'trf' runs the adaptation of the algorithm described in [STIR]_ for
    a linear least-squares problem. The iterations are essentially the same as
    in the nonlinear least-squares algorithm, but as the quadratic function
    model is always accurate, we don't need to track or modify the radius of
    a trust region. The line search (backtracking) is used as a safety net
    when a selected step does not decrease the cost function. Read more
    detailed description of the algorithm in `scipy.optimize.least_squares`.

    Method 'bvls' runs a Python implementation of the algorithm described in
    [BVLS]_. The algorithm maintains active and free sets of variables, on
    each iteration chooses a new variable to move from the active set to the
    free set and then solves the unconstrained least-squares problem on free
    variables. This algorithm is guaranteed to give an accurate solution
    eventually, but may require up to n iterations for a problem with n
    variables. Additionally, an ad-hoc initialization procedure is
    implemented, that determines which variables to set free or active
    initially. It takes some number of iterations before actual BVLS starts,
    but can significantly reduce the number of further iterations.

    References
    ----------
    .. [STIR] M. A. Branch, T. F. Coleman, and Y. Li, "A Subspace, Interior,
              and Conjugate Gradient Method for Large-Scale Bound-Constrained
              Minimization Problems," SIAM Journal on Scientific Computing,
              Vol. 21, Number 1, pp 1-23, 1999.
    .. [BVLS] P. B. Start and R. L. Parker, "Bounded-Variable Least-Squares:
              an Algorithm and Applications", Computational Statistics, 10,
              129-141, 1995.

    Examples
    --------
    In this example, a problem with a large sparse arrays and bounds on the
    variables is solved.

    >>> import numpy as np
    >>> from scipy.sparse import random_array
    >>> from scipy.optimize import lsq_linear
    >>> rng = np.random.default_rng()
    ...
    >>> m = 2000
    >>> n = 1000
    ...
    >>> A = random_array((m, n), density=1e-4, random_state=rng)
    >>> b = rng.standard_normal(m)
    ...
    >>> lb = rng.standard_normal(n)
    >>> ub = lb + 1
    ...
    >>> res = lsq_linear(A, b, bounds=(lb, ub), lsmr_tol='auto', verbose=1)
    The relative change of the cost function is less than `tol`.
    Number of iterations 10, initial cost 1.0070e+03, final cost 9.6602e+02,
    first-order optimality 2.21e-09.        # may vary
    """
    if method not in ['trf', 'bvls']:
        raise ValueError("`method` must be 'trf' or 'bvls'")

    if lsq_solver not in [None, 'exact', 'lsmr']:
        raise ValueError("`solver` must be None, 'exact' or 'lsmr'.")

    if verbose not in [0, 1, 2]:
        raise ValueError("`verbose` must be in [0, 1, 2].")

    if issparse(A):
        A = csr_array(A)
    elif not isinstance(A, LinearOperator):
        A = np.atleast_2d(np.asarray(A))

    if method == 'bvls':
        if lsq_solver == 'lsmr':
            raise ValueError("method='bvls' can't be used with "
                             "lsq_solver='lsmr'")

        if not isinstance(A, np.ndarray):
            raise ValueError("method='bvls' can't be used with `A` being "
                             "sparse or LinearOperator.")

    if lsq_solver is None:
        if isinstance(A, np.ndarray):
            lsq_solver = 'exact'
        else:
            lsq_solver = 'lsmr'
    elif lsq_solver == 'exact' and not isinstance(A, np.ndarray):
        raise ValueError("`exact` solver can't be used when `A` is "
                         "sparse or LinearOperator.")

    if len(A.shape) != 2:  # No ndim for LinearOperator.
        raise ValueError("`A` must have at most 2 dimensions.")

    if max_iter is not None and max_iter <= 0:
        raise ValueError("`max_iter` must be None or positive integer.")

    m, n = A.shape

    b = np.atleast_1d(b)
    if b.ndim != 1:
        raise ValueError("`b` must have at most 1 dimension.")

    if b.size != m:
        raise ValueError("Inconsistent shapes between `A` and `b`.")

    if isinstance(bounds, Bounds):
        lb = bounds.lb
        ub = bounds.ub
    else:
        lb, ub = prepare_bounds(bounds, n)

    if lb.shape != (n,) and ub.shape != (n,):
        raise ValueError("Bounds have wrong shape.")

    if np.any(lb >= ub):
        raise ValueError("Each lower bound must be strictly less than each "
                         "upper bound.")

    if lsmr_maxiter is not None and lsmr_maxiter < 1:
        raise ValueError("`lsmr_maxiter` must be None or positive integer.")

    if not ((isinstance(lsmr_tol, float) and lsmr_tol > 0) or
            lsmr_tol in ('auto', None)):
        raise ValueError("`lsmr_tol` must be None, 'auto', or positive float.")

    if lsq_solver == 'exact':
        unbd_lsq = np.linalg.lstsq(A, b, rcond=-1)
    elif lsq_solver == 'lsmr':
        first_lsmr_tol = lsmr_tol  # tol of first call to lsmr
        if lsmr_tol is None or lsmr_tol == 'auto':
            first_lsmr_tol = 1e-2 * tol  # default if lsmr_tol not defined
        unbd_lsq = lsmr(A, b, maxiter=lsmr_maxiter,
                        atol=first_lsmr_tol, btol=first_lsmr_tol)
    x_lsq = unbd_lsq[0]  # extract the solution from the least squares solver

    if in_bounds(x_lsq, lb, ub):
        r = A @ x_lsq - b
        cost = 0.5 * np.dot(r, r)
        termination_status = 3
        termination_message = TERMINATION_MESSAGES[termination_status]
        g = compute_grad(A, r)
        g_norm = norm(g, ord=np.inf)

        if verbose > 0:
            print(termination_message)
            print(f"Final cost {cost:.4e}, first-order optimality {g_norm:.2e}")

        return OptimizeResult(
            x=x_lsq, fun=r, cost=cost, optimality=g_norm,
            active_mask=np.zeros(n), unbounded_sol=unbd_lsq,
            nit=0, status=termination_status,
            message=termination_message, success=True)

    if method == 'trf':
        res = trf_linear(A, b, x_lsq, lb, ub, tol, lsq_solver, lsmr_tol,
                         max_iter, verbose, lsmr_maxiter=lsmr_maxiter)
    elif method == 'bvls':
        res = bvls(A, b, x_lsq, lb, ub, tol, max_iter, verbose)

    res.unbounded_sol = unbd_lsq
    res.message = TERMINATION_MESSAGES[res.status]
    res.success = res.status > 0

    if verbose > 0:
        print(res.message)
        print(
            f"Number of iterations {res.nit}, initial cost {res.initial_cost:.4e}, "
            f"final cost {res.cost:.4e}, first-order optimality {res.optimality:.2e}."
        )

    del res.initial_cost

    return res

################################################################################
# Method Path: scipy.signal._ltisys._YT_loop
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      det
# Call Count:  138
# Total Exec:  11
# Func Size:   121
################################################################################

def _YT_loop(ker_pole, transfer_matrix, poles, B, maxiter, rtol):
    """
    Algorithm "YT" Tits, Yang. Globally Convergent
    Algorithms for Robust Pole Assignment by State Feedback
    https://hdl.handle.net/1903/5598
    The poles P have to be sorted accordingly to section 6.2 page 20

    """
    # The IEEE edition of the YT paper gives useful information on the
    # optimal update order for the real poles in order to minimize the number
    # of times we have to loop over all poles, see page 1442
    nb_real = poles[np.isreal(poles)].shape[0]
    # hnb => Half Nb Real
    hnb = nb_real // 2

    # Stick to the indices in the paper and then remove one to get numpy array
    # index it is a bit easier to link the code to the paper this way even if it
    # is not very clean. The paper is unclear about what should be done when
    # there is only one real pole => use KNV0 on this real pole seem to work
    if nb_real > 0:
        #update the biggest real pole with the smallest one
        update_order = [[nb_real], [1]]
    else:
        update_order = [[],[]]

    r_comp = np.arange(nb_real+1, len(poles)+1, 2)
    # step 1.a
    r_p = np.arange(1, hnb+nb_real % 2)
    update_order[0].extend(2*r_p)
    update_order[1].extend(2*r_p+1)
    # step 1.b
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 1.c
    r_p = np.arange(1, hnb+1)
    update_order[0].extend(2*r_p-1)
    update_order[1].extend(2*r_p)
    # step 1.d
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 2.a
    r_j = np.arange(2, hnb+nb_real % 2)
    for j in r_j:
        for i in range(1, hnb+1):
            update_order[0].append(i)
            update_order[1].append(i+j)
    # step 2.b
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 2.c
    r_j = np.arange(2, hnb+nb_real % 2)
    for j in r_j:
        for i in range(hnb+1, nb_real+1):
            idx_1 = i+j
            if idx_1 > nb_real:
                idx_1 = i+j-nb_real
            update_order[0].append(i)
            update_order[1].append(idx_1)
    # step 2.d
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)
    # step 3.a
    for i in range(1, hnb+1):
        update_order[0].append(i)
        update_order[1].append(i+hnb)
    # step 3.b
    if hnb == 0 and np.isreal(poles[0]):
        update_order[0].append(1)
        update_order[1].append(1)
    update_order[0].extend(r_comp)
    update_order[1].extend(r_comp+1)

    update_order = np.array(update_order).T-1
    stop = False
    nb_try = 0
    while nb_try < maxiter and not stop:
        det_transfer_matrixb = np.abs(np.linalg.det(transfer_matrix))
        for i, j in update_order:
            if i == j:
                assert i == 0, "i!=0 for KNV call in YT"
                assert np.isreal(poles[i]), "calling KNV on a complex pole"
                _KNV0(B, ker_pole, transfer_matrix, i, poles)
            else:
                transfer_matrix_not_i_j = np.delete(transfer_matrix, (i, j),
                                                    axis=1)
                # after merge of gh-4249 great speed improvements could be
                # achieved using QR updates instead of full QR in the line below

                #to debug with numpy qr uncomment the line below
                #Q, _ = np.linalg.qr(transfer_matrix_not_i_j, mode="complete")
                Q, _ = s_qr(transfer_matrix_not_i_j, mode="full")

                if np.isreal(poles[i]):
                    assert np.isreal(poles[j]), "mixing real and complex " + \
                        "in YT_real" + str(poles)
                    _YT_real(ker_pole, Q, transfer_matrix, i, j)
                else:
                    assert ~np.isreal(poles[i]), "mixing real and complex " + \
                        "in YT_real" + str(poles)
                    _YT_complex(ker_pole, Q, transfer_matrix, i, j)

        det_transfer_matrix = np.max((np.sqrt(np.spacing(1)),
                                  np.abs(np.linalg.det(transfer_matrix))))
        cur_rtol = np.abs(
            (det_transfer_matrix -
             det_transfer_matrixb) /
            det_transfer_matrix)
        if cur_rtol < rtol and det_transfer_matrix > np.sqrt(np.spacing(1)):
            # Convergence test from YT page 21
            stop = True
        nb_try += 1
    return stop, cur_rtol, nb_try

################################################################################
# Method Path: scipy.integrate._bvp.fun_wrapped
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_bvp.py
# Target:      dot
# Call Count:  136
# Total Exec:  92
# Func Size:   8
################################################################################

def fun_wrapped(x, y, p):
            f = fun_p(x, y, p)
            if x[0] == a:
                f[:, 0] = np.dot(D, f[:, 0])
                f[:, 1:] += np.dot(S, y[:, 1:]) / (x[1:] - a)
            else:
                f += np.dot(S, y) / (x - a)
            return f

################################################################################
# Method Path: scipy.optimize._lsq.bvls.bvls
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/bvls.py
# Target:      dot
# Call Count:  135
# Total Exec:  23
# Func Size:   167
################################################################################

def bvls(A, b, x_lsq, lb, ub, tol, max_iter, verbose, rcond=None):
    m, n = A.shape

    x = x_lsq.copy()
    on_bound = np.zeros(n)

    mask = x <= lb
    x[mask] = lb[mask]
    on_bound[mask] = -1

    mask = x >= ub
    x[mask] = ub[mask]
    on_bound[mask] = 1

    free_set = on_bound == 0
    active_set = ~free_set
    free_set, = np.nonzero(free_set)

    r = A.dot(x) - b
    cost = 0.5 * np.dot(r, r)
    initial_cost = cost
    g = A.T.dot(r)

    cost_change = None
    step_norm = None
    iteration = 0

    if verbose == 2:
        print_header_linear()

    # This is the initialization loop. The requirement is that the
    # least-squares solution on free variables is feasible before BVLS starts.
    # One possible initialization is to set all variables to lower or upper
    # bounds, but many iterations may be required from this state later on.
    # The implemented ad-hoc procedure which intuitively should give a better
    # initial state: find the least-squares solution on current free variables,
    # if its feasible then stop, otherwise, set violating variables to
    # corresponding bounds and continue on the reduced set of free variables.

    while free_set.size > 0:
        if verbose == 2:
            optimality = compute_kkt_optimality(g, on_bound)
            print_iteration_linear(iteration, cost, cost_change, step_norm,
                                   optimality)

        iteration += 1
        x_free_old = x[free_set].copy()

        A_free = A[:, free_set]
        b_free = b - A.dot(x * active_set)
        z = lstsq(A_free, b_free, rcond=rcond)[0]

        lbv = z < lb[free_set]
        ubv = z > ub[free_set]
        v = lbv | ubv

        if np.any(lbv):
            ind = free_set[lbv]
            x[ind] = lb[ind]
            active_set[ind] = True
            on_bound[ind] = -1

        if np.any(ubv):
            ind = free_set[ubv]
            x[ind] = ub[ind]
            active_set[ind] = True
            on_bound[ind] = 1

        ind = free_set[~v]
        x[ind] = z[~v]

        r = A.dot(x) - b
        cost_new = 0.5 * np.dot(r, r)
        cost_change = cost - cost_new
        cost = cost_new
        g = A.T.dot(r)
        step_norm = norm(x[free_set] - x_free_old)

        if np.any(v):
            free_set = free_set[~v]
        else:
            break

    if max_iter is None:
        max_iter = n
    max_iter += iteration

    termination_status = None

    # Main BVLS loop.

    optimality = compute_kkt_optimality(g, on_bound)
    for iteration in range(iteration, max_iter):  # BVLS Loop A
        if verbose == 2:
            print_iteration_linear(iteration, cost, cost_change,
                                   step_norm, optimality)

        if optimality < tol:
            termination_status = 1

        if termination_status is not None:
            break

        move_to_free = np.argmax(g * on_bound)
        on_bound[move_to_free] = 0
        
        while True:   # BVLS Loop B

            free_set = on_bound == 0
            active_set = ~free_set
            free_set, = np.nonzero(free_set)
    
            x_free = x[free_set]
            x_free_old = x_free.copy()
            lb_free = lb[free_set]
            ub_free = ub[free_set]

            A_free = A[:, free_set]
            b_free = b - A.dot(x * active_set)
            z = lstsq(A_free, b_free, rcond=rcond)[0]

            lbv, = np.nonzero(z < lb_free)
            ubv, = np.nonzero(z > ub_free)
            v = np.hstack((lbv, ubv))

            if v.size > 0:
                alphas = np.hstack((
                    lb_free[lbv] - x_free[lbv],
                    ub_free[ubv] - x_free[ubv])) / (z[v] - x_free[v])

                i = np.argmin(alphas)
                i_free = v[i]
                alpha = alphas[i]

                x_free *= 1 - alpha
                x_free += alpha * z
                x[free_set] = x_free

                if i < lbv.size:
                    on_bound[free_set[i_free]] = -1
                else:
                    on_bound[free_set[i_free]] = 1
            else:
                x_free = z
                x[free_set] = x_free
                break

        step_norm = norm(x_free - x_free_old)

        r = A.dot(x) - b
        cost_new = 0.5 * np.dot(r, r)
        cost_change = cost - cost_new

        if cost_change < tol * cost:
            termination_status = 2
        cost = cost_new

        g = A.T.dot(r)
        optimality = compute_kkt_optimality(g, on_bound)

    if termination_status is None:
        termination_status = 0

    return OptimizeResult(
        x=x, fun=r, cost=cost, optimality=optimality, active_mask=on_bound,
        nit=iteration + 1, status=termination_status,
        initial_cost=initial_cost)

################################################################################
# Method Path: scipy.stats._multivariate._PSD.pinv
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  124
# Total Exec:  124
# Func Size:   5
################################################################################

def pinv(self):
        if self._pinv is None:
            self._pinv = np.dot(self.U, self.U.T)
        return self._pinv

################################################################################
# Method Path: scipy.stats._multivariate.matrix_normal_gen._logpdf
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  117
# Total Exec:  117
# Func Size:   37
################################################################################

def _logpdf(self, dims, X, mean, row_prec_rt, log_det_rowcov,
                col_prec_rt, log_det_colcov):
        """Log of the matrix normal probability density function.

        Parameters
        ----------
        dims : tuple
            Dimensions of the matrix variates
        X : ndarray
            Points at which to evaluate the log of the probability
            density function
        mean : ndarray
            Mean of the distribution
        row_prec_rt : ndarray
            A decomposition such that np.dot(row_prec_rt, row_prec_rt.T)
            is the inverse of the among-row covariance matrix
        log_det_rowcov : float
            Logarithm of the determinant of the among-row covariance matrix
        col_prec_rt : ndarray
            A decomposition such that np.dot(col_prec_rt, col_prec_rt.T)
            is the inverse of the among-column covariance matrix
        log_det_colcov : float
            Logarithm of the determinant of the among-column covariance matrix

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'logpdf' instead.

        """
        numrows, numcols = dims
        roll_dev = np.moveaxis(X-mean, -1, 0)
        scale_dev = np.tensordot(col_prec_rt.T,
                                 np.dot(roll_dev, row_prec_rt), 1)
        maha = np.sum(np.sum(np.square(scale_dev), axis=-1), axis=0)
        return -0.5 * (numrows*numcols*_LOG_2PI + numcols*log_det_rowcov
                       + numrows*log_det_colcov + maha)

################################################################################
# Method Path: scipy.optimize._constraints.fun
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_constraints.py
# Target:      dot
# Call Count:  110
# Total Exec:  110
# Func Size:   2
################################################################################

def fun(x):
            return np.dot(A, x)

################################################################################
# Method Path: scipy.stats._multivariate.wishart_gen._logpdf
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      trace
# Call Count:  106
# Total Exec:  28
# Func Size:   45
################################################################################

def _logpdf(self, x, dim, df, scale, log_det_scale, C):
        """Log of the Wishart probability density function.

        Parameters
        ----------
        x : ndarray
            Points at which to evaluate the log of the probability
            density function
        dim : int
            Dimension of the scale matrix
        df : int
            Degrees of freedom
        scale : ndarray
            Scale matrix
        log_det_scale : float
            Logarithm of the determinant of the scale matrix
        C : ndarray
            Cholesky factorization of the scale matrix, lower triangular.

        Notes
        -----
        As this function does no argument checking, it should not be
        called directly; use 'logpdf' instead.

        """
        # log determinant of x
        # Note: x has components along the last axis, so that x.T has
        # components alone the 0-th axis. Then since det(A) = det(A'), this
        # gives us a 1-dim vector of determinants

        # Retrieve tr(scale^{-1} x)
        log_det_x = np.empty(x.shape[-1])
        scale_inv_x = np.empty(x.shape)
        tr_scale_inv_x = np.empty(x.shape[-1])
        for i in range(x.shape[-1]):
            _, log_det_x[i] = self._cholesky_logdet(x[:, :, i])
            scale_inv_x[:, :, i] = scipy.linalg.cho_solve((C, True), x[:, :, i])
            tr_scale_inv_x[i] = scale_inv_x[:, :, i].trace()

        # Log PDF
        out = ((0.5 * (df - dim - 1) * log_det_x - 0.5 * tr_scale_inv_x) -
               (0.5 * df * dim * _LOG_2 + 0.5 * df * log_det_scale +
                multigammaln(0.5*df, dim)))

        return out

################################################################################
# Method Path: scipy.signal._ltisys.place_poles
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  105
# Total Exec:  29
# Func Size:   349
################################################################################

def place_poles(A, B, poles, method="YT", rtol=1e-3, maxiter=30):
    """
    Compute K such that eigenvalues (A - dot(B, K))=poles.

    K is the gain matrix such as the plant described by the linear system
    ``AX+BU`` will have its closed-loop poles, i.e the eigenvalues ``A - B*K``,
    as close as possible to those asked for in poles.

    SISO, MISO and MIMO systems are supported.

    Parameters
    ----------
    A, B : ndarray
        State-space representation of linear system ``AX + BU``.
    poles : array_like
        Desired real poles and/or complex conjugates poles.
        Complex poles are only supported with ``method="YT"`` (default).
    method : {'YT', 'KNV0'}, optional
        Which method to choose to find the gain matrix K. One of:

        - 'YT': Yang Tits
        - 'KNV0': Kautsky, Nichols, Van Dooren update method 0

        See References and Notes for details on the algorithms.
    rtol : float, optional
        After each iteration the determinant of the eigenvectors of
        ``A - B*K`` is compared to its previous value, when the relative
        error between these two values becomes lower than `rtol` the algorithm
        stops.  Default is 1e-3.
    maxiter : int, optional
        Maximum number of iterations to compute the gain matrix.
        Default is 30.

    Returns
    -------
    full_state_feedback : Bunch object
        full_state_feedback is composed of:
            gain_matrix : 2-D ndarray
                The closed loop matrix K such as the eigenvalues of ``A-BK``
                are as close as possible to the requested poles.
            computed_poles : 1-D ndarray
                The poles corresponding to ``A-BK`` sorted as first the real
                poles in increasing order, then the complex conjugates in
                lexicographic order.
            requested_poles : 1-D ndarray
                The poles the algorithm was asked to place sorted as above,
                they may differ from what was achieved.
            X : 2-D ndarray
                The transfer matrix such as ``X * diag(poles) = (A - B*K)*X``
                (see Notes)
            rtol : float
                The relative tolerance achieved on ``det(X)`` (see Notes).
                `rtol` will be NaN if it is possible to solve the system
                ``diag(poles) = (A - B*K)``, or 0 when the optimization
                algorithms can't do anything i.e when ``B.shape[1] == 1``.
            nb_iter : int
                The number of iterations performed before converging.
                `nb_iter` will be NaN if it is possible to solve the system
                ``diag(poles) = (A - B*K)``, or 0 when the optimization
                algorithms can't do anything i.e when ``B.shape[1] == 1``.

    Notes
    -----
    The Tits and Yang (YT), [2]_ paper is an update of the original Kautsky et
    al. (KNV) paper [1]_.  KNV relies on rank-1 updates to find the transfer
    matrix X such that ``X * diag(poles) = (A - B*K)*X``, whereas YT uses
    rank-2 updates. This yields on average more robust solutions (see [2]_
    pp 21-22), furthermore the YT algorithm supports complex poles whereas KNV
    does not in its original version.  Only update method 0 proposed by KNV has
    been implemented here, hence the name ``'KNV0'``.

    KNV extended to complex poles is used in Matlab's ``place`` function, YT is
    distributed under a non-free licence by Slicot under the name ``robpole``.
    It is unclear and undocumented how KNV0 has been extended to complex poles
    (Tits and Yang claim on page 14 of their paper that their method can not be
    used to extend KNV to complex poles), therefore only YT supports them in
    this implementation.

    As the solution to the problem of pole placement is not unique for MIMO
    systems, both methods start with a tentative transfer matrix which is
    altered in various way to increase its determinant.  Both methods have been
    proven to converge to a stable solution, however depending on the way the
    initial transfer matrix is chosen they will converge to different
    solutions and therefore there is absolutely no guarantee that using
    ``'KNV0'`` will yield results similar to Matlab's or any other
    implementation of these algorithms.

    Using the default method ``'YT'`` should be fine in most cases; ``'KNV0'``
    is only provided because it is needed by ``'YT'`` in some specific cases.
    Furthermore ``'YT'`` gives on average more robust results than ``'KNV0'``
    when ``abs(det(X))`` is used as a robustness indicator.

    [2]_ is available as a technical report on the following URL:
    https://hdl.handle.net/1903/5598

    References
    ----------
    .. [1] J. Kautsky, N.K. Nichols and P. van Dooren, "Robust pole assignment
           in linear state feedback", International Journal of Control, Vol. 41
           pp. 1129-1155, 1985.
    .. [2] A.L. Tits and Y. Yang, "Globally convergent algorithms for robust
           pole assignment by state feedback", IEEE Transactions on Automatic
           Control, Vol. 41, pp. 1432-1452, 1996.

    Examples
    --------
    A simple example demonstrating real pole placement using both KNV and YT
    algorithms.  This is example number 1 from section 4 of the reference KNV
    publication ([1]_):

    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt

    >>> A = np.array([[ 1.380,  -0.2077,  6.715, -5.676  ],
    ...               [-0.5814, -4.290,   0,      0.6750 ],
    ...               [ 1.067,   4.273,  -6.654,  5.893  ],
    ...               [ 0.0480,  4.273,   1.343, -2.104  ]])
    >>> B = np.array([[ 0,      5.679 ],
    ...               [ 1.136,  1.136 ],
    ...               [ 0,      0,    ],
    ...               [-3.146,  0     ]])
    >>> P = np.array([-0.2, -0.5, -5.0566, -8.6659])

    Now compute K with KNV method 0, with the default YT method and with the YT
    method while forcing 100 iterations of the algorithm and print some results
    after each call.

    >>> fsf1 = signal.place_poles(A, B, P, method='KNV0')
    >>> fsf1.gain_matrix
    array([[ 0.20071427, -0.96665799,  0.24066128, -0.10279785],
           [ 0.50587268,  0.57779091,  0.51795763, -0.41991442]])

    >>> fsf2 = signal.place_poles(A, B, P)  # uses YT method
    >>> fsf2.computed_poles
    array([-8.6659, -5.0566, -0.5   , -0.2   ])

    >>> fsf3 = signal.place_poles(A, B, P, rtol=-1, maxiter=100)
    >>> fsf3.X
    array([[ 0.52072442+0.j, -0.08409372+0.j, -0.56847937+0.j,  0.74823657+0.j],
           [-0.04977751+0.j, -0.80872954+0.j,  0.13566234+0.j, -0.29322906+0.j],
           [-0.82266932+0.j, -0.19168026+0.j, -0.56348322+0.j, -0.43815060+0.j],
           [ 0.22267347+0.j,  0.54967577+0.j, -0.58387806+0.j, -0.40271926+0.j]])

    The absolute value of the determinant of X is a good indicator to check the
    robustness of the results, both ``'KNV0'`` and ``'YT'`` aim at maximizing
    it.  Below a comparison of the robustness of the results above:

    >>> abs(np.linalg.det(fsf1.X)) < abs(np.linalg.det(fsf2.X))
    True
    >>> abs(np.linalg.det(fsf2.X)) < abs(np.linalg.det(fsf3.X))
    True

    Now a simple example for complex poles:

    >>> A = np.array([[ 0,  7/3.,  0,   0   ],
    ...               [ 0,   0,    0,  7/9. ],
    ...               [ 0,   0,    0,   0   ],
    ...               [ 0,   0,    0,   0   ]])
    >>> B = np.array([[ 0,  0 ],
    ...               [ 0,  0 ],
    ...               [ 1,  0 ],
    ...               [ 0,  1 ]])
    >>> P = np.array([-3, -1, -2-1j, -2+1j]) / 3.
    >>> fsf = signal.place_poles(A, B, P, method='YT')

    We can plot the desired and computed poles in the complex plane:

    >>> t = np.linspace(0, 2*np.pi, 401)
    >>> plt.plot(np.cos(t), np.sin(t), 'k--')  # unit circle
    >>> plt.plot(fsf.requested_poles.real, fsf.requested_poles.imag,
    ...          'wo', label='Desired')
    >>> plt.plot(fsf.computed_poles.real, fsf.computed_poles.imag, 'bx',
    ...          label='Placed')
    >>> plt.grid()
    >>> plt.axis('image')
    >>> plt.axis([-1.1, 1.1, -1.1, 1.1])
    >>> plt.legend(bbox_to_anchor=(1.05, 1), loc=2, numpoints=1)

    """
    # Move away all the inputs checking, it only adds noise to the code
    update_loop, poles = _valid_inputs(A, B, poles, method, rtol, maxiter)

    # The current value of the relative tolerance we achieved
    cur_rtol = 0
    # The number of iterations needed before converging
    nb_iter = 0

    # Step A: QR decomposition of B page 1132 KN
    # to debug with numpy qr uncomment the line below
    # u, z = np.linalg.qr(B, mode="complete")
    u, z = s_qr(B, mode="full")
    rankB = np.linalg.matrix_rank(B)
    u0 = u[:, :rankB]
    u1 = u[:, rankB:]
    z = z[:rankB, :]

    # If we can use the identity matrix as X the solution is obvious
    if B.shape[0] == rankB:
        # if B is square and full rank there is only one solution
        # such as (A+BK)=inv(X)*diag(P)*X with X=eye(A.shape[0])
        # i.e K=inv(B)*(diag(P)-A)
        # if B has as many lines as its rank (but not square) there are many
        # solutions and we can choose one using least squares
        # => use lstsq in both cases.
        # In both cases the transfer matrix X will be eye(A.shape[0]) and I
        # can hardly think of a better one so there is nothing to optimize
        #
        # for complex poles we use the following trick
        #
        # |a -b| has for eigenvalues a+b and a-b
        # |b a|
        #
        # |a+bi 0| has the obvious eigenvalues a+bi and a-bi
        # |0 a-bi|
        #
        # e.g solving the first one in R gives the solution
        # for the second one in C
        diag_poles = np.zeros(A.shape)
        idx = 0
        while idx < poles.shape[0]:
            p = poles[idx]
            diag_poles[idx, idx] = np.real(p)
            if ~np.isreal(p):
                diag_poles[idx, idx+1] = -np.imag(p)
                diag_poles[idx+1, idx+1] = np.real(p)
                diag_poles[idx+1, idx] = np.imag(p)
                idx += 1  # skip next one
            idx += 1
        gain_matrix = np.linalg.lstsq(B, diag_poles-A, rcond=-1)[0]
        transfer_matrix = np.eye(A.shape[0])
        cur_rtol = np.nan
        nb_iter = np.nan
    else:
        # step A (p1144 KNV) and beginning of step F: decompose
        # dot(U1.T, A-P[i]*I).T and build our set of transfer_matrix vectors
        # in the same loop
        ker_pole = []

        # flag to skip the conjugate of a complex pole
        skip_conjugate = False
        # select orthonormal base ker_pole for each Pole and vectors for
        # transfer_matrix
        for j in range(B.shape[0]):
            if skip_conjugate:
                skip_conjugate = False
                continue
            pole_space_j = np.dot(u1.T, A-poles[j]*np.eye(B.shape[0])).T

            # after QR Q=Q0|Q1
            # only Q0 is used to reconstruct  the qr'ed (dot Q, R) matrix.
            # Q1 is orthogonal to Q0 and will be multiplied by the zeros in
            # R when using mode "complete". In default mode Q1 and the zeros
            # in R are not computed

            # To debug with numpy qr uncomment the line below
            # Q, _ = np.linalg.qr(pole_space_j, mode="complete")
            Q, _ = s_qr(pole_space_j, mode="full")

            ker_pole_j = Q[:, pole_space_j.shape[1]:]

            # We want to select one vector in ker_pole_j to build the transfer
            # matrix, however qr returns sometimes vectors with zeros on the
            # same line for each pole and this yields very long convergence
            # times.
            # Or some other times a set of vectors, one with zero imaginary
            # part and one (or several) with imaginary parts. After trying
            # many ways to select the best possible one (eg ditch vectors
            # with zero imaginary part for complex poles) I ended up summing
            # all vectors in ker_pole_j, this solves 100% of the problems and
            # is a valid choice for transfer_matrix.
            # This way for complex poles we are sure to have a non zero
            # imaginary part that way, and the problem of lines full of zeros
            # in transfer_matrix is solved too as when a vector from
            # ker_pole_j has a zero the other one(s) when
            # ker_pole_j.shape[1]>1) for sure won't have a zero there.

            transfer_matrix_j = np.sum(ker_pole_j, axis=1)[:, np.newaxis]
            transfer_matrix_j = (transfer_matrix_j /
                                 np.linalg.norm(transfer_matrix_j))
            if ~np.isreal(poles[j]):  # complex pole
                transfer_matrix_j = np.hstack([np.real(transfer_matrix_j),
                                               np.imag(transfer_matrix_j)])
                ker_pole.extend([ker_pole_j, ker_pole_j])

                # Skip next pole as it is the conjugate
                skip_conjugate = True
            else:  # real pole, nothing to do
                ker_pole.append(ker_pole_j)

            if j == 0:
                transfer_matrix = transfer_matrix_j
            else:
                transfer_matrix = np.hstack((transfer_matrix, transfer_matrix_j))

        if rankB > 1:  # otherwise there is nothing we can optimize
            stop, cur_rtol, nb_iter = update_loop(ker_pole, transfer_matrix,
                                                  poles, B, maxiter, rtol)
            if not stop and rtol > 0:
                # if rtol<=0 the user has probably done that on purpose,
                # don't annoy them
                err_msg = (
                    "Convergence was not reached after maxiter iterations.\n"
                    f"You asked for a tolerance of {rtol}, we got {cur_rtol}."
                    )
                warnings.warn(err_msg, stacklevel=2)

        # reconstruct transfer_matrix to match complex conjugate pairs,
        # ie transfer_matrix_j/transfer_matrix_j+1 are
        # Re(Complex_pole), Im(Complex_pole) now and will be Re-Im/Re+Im after
        transfer_matrix = transfer_matrix.astype(complex)
        idx = 0
        while idx < poles.shape[0]-1:
            if ~np.isreal(poles[idx]):
                rel = transfer_matrix[:, idx].copy()
                img = transfer_matrix[:, idx+1]
                # rel will be an array referencing a column of transfer_matrix
                # if we don't copy() it will changer after the next line and
                # and the line after will not yield the correct value
                transfer_matrix[:, idx] = rel-1j*img
                transfer_matrix[:, idx+1] = rel+1j*img
                idx += 1  # skip next one
            idx += 1

        try:
            m = np.linalg.solve(transfer_matrix.T, np.dot(np.diag(poles),
                                                          transfer_matrix.T)).T
            gain_matrix = np.linalg.solve(z, np.dot(u0.T, m-A))
        except np.linalg.LinAlgError as e:
            raise ValueError("The poles you've chosen can't be placed. "
                             "Check the controllability matrix and try "
                             "another set of poles") from e

    # Beware: Kautsky solves A+BK but the usual form is A-BK
    gain_matrix = -gain_matrix
    # K still contains complex with ~=0j imaginary parts, get rid of them
    gain_matrix = np.real(gain_matrix)

    full_state_feedback = Bunch()
    full_state_feedback.gain_matrix = gain_matrix
    full_state_feedback.computed_poles = _order_complex_poles(
        np.linalg.eig(A - np.dot(B, gain_matrix))[0]
        )
    full_state_feedback.requested_poles = poles
    full_state_feedback.X = transfer_matrix
    full_state_feedback.rtol = cur_rtol
    full_state_feedback.nb_iter = nb_iter

    return full_state_feedback

################################################################################
# Method Path: scipy.optimize._lsq.least_squares.call_minpack
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_lsq/least_squares.py
# Target:      dot
# Call Count:  102
# Total Exec:  51
# Func Size:   54
################################################################################

def call_minpack(fun, x0, jac, ftol, xtol, gtol, max_nfev, x_scale, jac_method=None):
    n = x0.size

    # Compute MINPACK's `diag`, which is inverse of our `x_scale` and
    # ``x_scale='jac'`` corresponds to ``diag=None``.

    # 1.16.0 - default x_scale changed to 'jac', with diag=None
    if x_scale is None or (isinstance(x_scale, str) and x_scale == 'jac'):
        diag = None
    else:
        # x_scale specified, so use that
        diag = 1 / x_scale

    full_output = True
    col_deriv = False
    factor = 100.0

    if max_nfev is None:
        max_nfev = 100 * n

    # lmder is typically used for systems with analytic jacobians, with lmdif being
    # used if there is only an objective fun (lmdif uses finite differences to estimate
    # jacobian). Otherwise they're very similar internally.
    # We now do all the finite differencing in VectorFunction, which means we can drop
    # lmdif and just use lmder.

    # for sending a copy of x0 into _lmder
    xp = array_namespace(x0)

    x, info, status = _minpack._lmder(
        fun, jac, xp.astype(x0, x0.dtype), (), full_output, col_deriv,
        ftol, xtol, gtol, max_nfev, factor, diag)

    f = info['fvec']
    J = jac(x)

    cost = 0.5 * np.dot(f, f)
    g = J.T.dot(f)
    g_norm = norm(g, ord=np.inf)

    nfev = info['nfev']
    if callable(jac_method):
        # user supplied a callable ("analytic") jac
        njev = info.get('njev', None)
    else:
        # If there are no analytic jacobian evaluations we need to set `njev=None`.
        njev = None

    status = FROM_MINPACK_TO_COMMON[status]
    active_mask = np.zeros_like(x0, dtype=int)

    return OptimizeResult(
        x=x, cost=cost, fun=f, jac=J, grad=g, optimality=g_norm,
        active_mask=active_mask, nfev=nfev, njev=njev, status=status)

################################################################################
# Method Path: scipy.stats._multivariate.multivariate_normal_gen.rvs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  97
# Total Exec:  150
# Func Size:   36
################################################################################

def rvs(self, mean=None, cov=1, size=1, random_state=None):
        """Draw random samples from a multivariate normal distribution.

        Parameters
        ----------
        %(_mvn_doc_default_callparams)s
        size : int, optional
            Number of samples to draw (default 1).
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of size (`size`, `N`), where `N` is the
            dimension of the random variable.

        Notes
        -----
        %(_mvn_doc_callparams_note)s

        """
        dim, mean, cov_object = self._process_parameters(mean, cov)
        random_state = self._get_random_state(random_state)

        if isinstance(cov_object, _covariance.CovViaPSD):
            cov = cov_object.covariance
            out = random_state.multivariate_normal(mean, cov, size)
            out = _squeeze_output(out)
        else:
            size = size or tuple()
            if not np.iterable(size):
                size = (size,)
            shape = tuple(size) + (cov_object.shape[-1],)
            x = random_state.normal(size=shape)
            out = mean + cov_object.colorize(x)
        return out

################################################################################
# Method Path: scipy.optimize._linprog_rs._phase_one
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_rs.py
# Target:      dot
# Call Count:  88
# Total Exec:  88
# Func Size:   68
################################################################################

def _phase_one(A, b, x0, callback, postsolve_args, maxiter, tol, disp,
               maxupdate, mast, pivot):
    """
    The purpose of phase one is to find an initial basic feasible solution
    (BFS) to the original problem.

    Generates an auxiliary problem with a trivial BFS and an objective that
    minimizes infeasibility of the original problem. Solves the auxiliary
    problem using the main simplex routine (phase two). This either yields
    a BFS to the original problem or determines that the original problem is
    infeasible. If feasible, phase one detects redundant rows in the original
    constraint matrix and removes them, then chooses additional indices as
    necessary to complete a basis/BFS for the original problem.
    """

    m, n = A.shape
    status = 0

    # generate auxiliary problem to get initial BFS
    A, b, c, basis, x, status = _generate_auxiliary_problem(A, b, x0, tol)

    if status == 6:
        residual = c.dot(x)
        iter_k = 0
        return x, basis, A, b, residual, status, iter_k

    # solve auxiliary problem
    phase_one_n = n
    iter_k = 0
    x, basis, status, iter_k = _phase_two(c, A, x, basis, callback,
                                          postsolve_args,
                                          maxiter, tol, disp,
                                          maxupdate, mast, pivot,
                                          iter_k, phase_one_n)

    # check for infeasibility
    residual = c.dot(x)
    if status == 0 and residual > tol:
        status = 2

    # drive artificial variables out of basis
    # TODO: test redundant row removal better
    # TODO: make solve more efficient with BGLU? This could take a while.
    keep_rows = np.ones(m, dtype=bool)
    for basis_column in basis[basis >= n]:
        B = A[:, basis]
        try:
            basis_finder = np.abs(solve(B, A))  # inefficient
            pertinent_row = np.argmax(basis_finder[:, basis_column])
            eligible_columns = np.ones(n, dtype=bool)
            eligible_columns[basis[basis < n]] = 0
            eligible_column_indices = np.where(eligible_columns)[0]
            index = np.argmax(basis_finder[:, :n]
                              [pertinent_row, eligible_columns])
            new_basis_column = eligible_column_indices[index]
            if basis_finder[pertinent_row, new_basis_column] < tol:
                keep_rows[pertinent_row] = False
            else:
                basis[basis == basis_column] = new_basis_column
        except LinAlgError:
            status = 4

    # form solution to original problem
    A = A[keep_rows, :n]
    basis = basis[keep_rows]
    x = x[:n]
    m = A.shape[0]
    return x, basis, A, b, residual, status, iter_k

################################################################################
# Method Path: scipy.sparse.linalg._matfuncs.MatrixPowerOperator._rmatvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_matfuncs.py
# Target:      dot
# Call Count:  82
# Total Exec:  34
# Func Size:   6
################################################################################

def _rmatvec(self, x):
        A_T = self._A.T
        x = x.ravel()
        for i in range(self._p):
            x = A_T.dot(x)
        return x

################################################################################
# Method Path: scipy.interpolate._rbf.Rbf.__call__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/interpolate/_rbf.py
# Target:      dot
# Call Count:  80
# Total Exec:  80
# Func Size:   11
################################################################################

def __call__(self, *args):
        args = [np.asarray(x) for x in args]
        if not all([x.shape == y.shape for x in args for y in args]):
            raise ValueError("Array lengths must be equal")
        if self._target_dim > 1:
            shp = args[0].shape + (self._target_dim,)
        else:
            shp = args[0].shape
        xa = np.asarray([a.flatten() for a in args], dtype=np.float64)
        r = self._call_norm(xa, self.xi)
        return np.dot(self._function(r), self.nodes).reshape(shp)

################################################################################
# Method Path: scipy.optimize._trustregion_dogleg.DoglegSubproblem.cauchy_point
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_trustregion_dogleg.py
# Target:      dot
# Call Count:  76
# Total Exec:  86
# Func Size:   9
################################################################################

def cauchy_point(self):
        """
        The Cauchy point is minimal along the direction of steepest descent.
        """
        if self._cauchy_point is None:
            g = self.jac
            Bg = self.hessp(g)
            self._cauchy_point = -(np.dot(g, g) / np.dot(g, Bg)) * g
        return self._cauchy_point

################################################################################
# Method Path: scipy.signal._ltisys._KNV0
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  76
# Total Exec:  38
# Func Size:   45
################################################################################

def _KNV0(B, ker_pole, transfer_matrix, j, poles):
    """
    Algorithm "KNV0" Kautsky et Al. Robust pole
    assignment in linear state feedback, Int journal of Control
    1985, vol 41 p 1129->1155
    https://la.epfl.ch/files/content/sites/la/files/
        users/105941/public/KautskyNicholsDooren

    """
    # Remove xj form the base
    transfer_matrix_not_j = np.delete(transfer_matrix, j, axis=1)
    # If we QR this matrix in full mode Q=Q0|Q1
    # then Q1 will be a single column orthogonal to
    # Q0, that's what we are looking for !

    # After merge of gh-4249 great speed improvements could be achieved
    # using QR updates instead of full QR in the line below

    # To debug with numpy qr uncomment the line below
    # Q, R = np.linalg.qr(transfer_matrix_not_j, mode="complete")
    Q, R = s_qr(transfer_matrix_not_j, mode="full")

    mat_ker_pj = np.dot(ker_pole[j], ker_pole[j].T)
    yj = np.dot(mat_ker_pj, Q[:, -1])

    # If Q[:, -1] is "almost" orthogonal to ker_pole[j] its
    # projection into ker_pole[j] will yield a vector
    # close to 0.  As we are looking for a vector in ker_pole[j]
    # simply stick with transfer_matrix[:, j] (unless someone provides me with
    # a better choice ?)

    if not np.allclose(yj, 0):
        xj = yj/np.linalg.norm(yj)
        transfer_matrix[:, j] = xj

################################################################################
# Method Path: scipy.stats._qmc.Sobol._scramble
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_qmc.py
# Target:      dot
# Call Count:  69
# Total Exec:  69
# Func Size:   17
################################################################################

def _scramble(self) -> None:
        """Scramble the sequence using LMS+shift."""
        # Generate shift vector
        self._shift = np.dot(
            rng_integers(self.rng, 2, size=(self.d, self.bits),
                         dtype=self.dtype_i),
            # pyrefly:ignore[no-matching-overload]
            2 ** np.arange(self.bits, dtype=self.dtype_i),
        )
        # Generate lower triangular matrices (stacked across dimensions)
        ltm = np.tril(rng_integers(self.rng, 2,
                                   size=(self.d, self.bits, self.bits),
                                   dtype=self.dtype_i))
        _cscramble(
            dim=self.d, bits=self.bits,  # type: ignore[arg-type]
            ltm=ltm, sv=self._sv
        )

################################################################################
# Method Path: scipy.optimize._minpack_py.curve_fit
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_minpack_py.py
# Target:      dot
# Call Count:  59
# Total Exec:  172
# Func Size:   488
################################################################################

def curve_fit(f, xdata, ydata, p0=None, sigma=None, absolute_sigma=False,
              check_finite=None, bounds=(-np.inf, np.inf), method=None,
              jac=None, *, full_output=False, nan_policy=None,
              **kwargs):
    """
    Use non-linear least squares to fit a function, f, to data.

    Assumes ``ydata = f(xdata, *params) + eps``.

    Parameters
    ----------
    f : callable
        The model function, f(x, ...). It must take the independent
        variable as the first argument and the parameters to fit as
        separate remaining arguments.
    xdata : array_like
        The independent variable where the data is measured.
        Should usually be an M-length sequence or an (k,M)-shaped array for
        functions with k predictors, and each element should be float
        convertible if it is an array like object.
    ydata : array_like
        The dependent data, a length M array - nominally ``f(xdata, ...)``.
    p0 : array_like, optional
        Initial guess for the parameters (length N). If None, then the
        initial values will all be 1 (if the number of parameters for the
        function can be determined using introspection, otherwise a
        ValueError is raised).
    sigma : None or scalar or M-length sequence or MxM array, optional
        Determines the uncertainty in `ydata`. If we define residuals as
        ``r = ydata - f(xdata, *popt)``, then the interpretation of `sigma`
        depends on its number of dimensions:

        - A scalar or 1-D `sigma` should contain values of standard deviations of
          errors in `ydata`. In this case, the optimized function is
          ``chisq = sum((r / sigma) ** 2)``.

        - A 2-D `sigma` should contain the covariance matrix of
          errors in `ydata`. In this case, the optimized function is
          ``chisq = r.T @ inv(sigma) @ r``.

          .. versionadded:: 0.19

        None (default) is equivalent of 1-D `sigma` filled with ones.
    absolute_sigma : bool, optional
        If True, `sigma` is used in an absolute sense and the estimated parameter
        covariance `pcov` reflects these absolute values.

        If False (default), only the relative magnitudes of the `sigma` values matter.
        The returned parameter covariance matrix `pcov` is based on scaling
        `sigma` by a constant factor. This constant is set by demanding that the
        reduced `chisq` for the optimal parameters `popt` when using the
        *scaled* `sigma` equals unity. In other words, `sigma` is scaled to
        match the sample variance of the residuals after the fit. Default is False.
        Mathematically,
        ``pcov(absolute_sigma=False) = pcov(absolute_sigma=True) * chisq(popt)/(M-N)``
    check_finite : bool, optional
        If True, check that the input arrays do not contain nans of infs,
        and raise a ValueError if they do. Setting this parameter to
        False may silently produce nonsensical results if the input arrays
        do contain nans. Default is True if `nan_policy` is not specified
        explicitly and False otherwise.
    bounds : 2-tuple of array_like or `Bounds`, optional
        Lower and upper bounds on parameters. Defaults to no bounds.
        There are two ways to specify the bounds:

        - Instance of `Bounds` class.

        - 2-tuple of array_like: Each element of the tuple must be either
          an array with the length equal to the number of parameters, or a
          scalar (in which case the bound is taken to be the same for all
          parameters). Use ``np.inf`` with an appropriate sign to disable
          bounds on all or some parameters.

    method : {'lm', 'trf', 'dogbox'}, optional
        Method to use for optimization. See `least_squares` for more details.
        Default is 'lm' for unconstrained problems and 'trf' if `bounds` are
        provided. The method 'lm' won't work when the number of observations
        is less than the number of variables, use 'trf' or 'dogbox' in this
        case.

        .. versionadded:: 0.17
    jac : callable, str or None, optional
        Function with signature ``jac(x, ...)`` which computes the Jacobian
        matrix of the model function with respect to parameters as a dense
        array_like structure. It will be scaled according to provided `sigma`.
        If None (default), the Jacobian will be estimated numerically.
        String keywords for 'trf' and 'dogbox' methods can be used to select
        a finite difference scheme, see `least_squares`.

        .. versionadded:: 0.18
    full_output : bool, optional
        If True, this function returns additional information: `infodict`,
        `mesg`, and `ier`.

        .. versionadded:: 1.9
    nan_policy : {'raise', 'omit', None}, optional
        Defines how to handle when input contains nan.
        The following options are available (default is None):

        * 'raise': throws an error
        * 'omit': performs the calculations ignoring nan values
        * None: no special handling of NaNs is performed
          (except what is done by check_finite); the behavior when NaNs
          are present is implementation-dependent and may change.

        Note that if this value is specified explicitly (not None),
        `check_finite` will be set as False.

        .. versionadded:: 1.11
    **kwargs
        Keyword arguments passed to `leastsq` for ``method='lm'`` or
        `least_squares` otherwise.

    Returns
    -------
    popt : array
        Optimal values for the parameters so that the sum of the squared
        residuals of ``f(xdata, *popt) - ydata`` is minimized.
    pcov : 2-D array
        The estimated approximate covariance of popt. The diagonals provide
        the variance of the parameter estimate. To compute one standard
        deviation errors on the parameters, use
        ``perr = np.sqrt(np.diag(pcov))``. Note that the relationship between
        `cov` and parameter error estimates is derived based on a linear
        approximation to the model function around the optimum [1]_.
        When this approximation becomes inaccurate, `cov` may not provide an
        accurate measure of uncertainty.

        How the `sigma` parameter affects the estimated covariance
        depends on `absolute_sigma` argument, as described above.

        If the Jacobian matrix at the solution doesn't have a full rank, then
        'lm' method returns a matrix filled with ``np.inf``, on the other hand
        'trf'  and 'dogbox' methods use Moore-Penrose pseudoinverse to compute
        the covariance matrix. Covariance matrices with large condition numbers
        (e.g. computed with `numpy.linalg.cond`) may indicate that results are
        unreliable.
    infodict : dict (returned only if `full_output` is True)
        a dictionary of optional outputs with the keys:

        ``nfev``
            The number of function calls. Methods 'trf' and 'dogbox' do not
            count function calls for numerical Jacobian approximation,
            as opposed to 'lm' method.
        ``fvec``
            The residual values evaluated at the solution, for a 1-D `sigma`
            this is ``(f(x, *popt) - ydata)/sigma``.
        ``fjac``
            A permutation of the R matrix of a QR
            factorization of the final approximate
            Jacobian matrix, stored column wise.
            Together with ipvt, the covariance of the
            estimate can be approximated.
            Method 'lm' only provides this information.
        ``ipvt``
            An integer array of length N which defines
            a permutation matrix, p, such that
            fjac*p = q*r, where r is upper triangular
            with diagonal elements of nonincreasing
            magnitude. Column j of p is column ipvt(j)
            of the identity matrix.
            Method 'lm' only provides this information.
        ``qtf``
            The vector (transpose(q) * fvec).
            Method 'lm' only provides this information.

        .. versionadded:: 1.9
    mesg : str (returned only if `full_output` is True)
        A string message giving information about the solution.

        .. versionadded:: 1.9
    ier : int (returned only if `full_output` is True)
        An integer flag. If it is equal to 1, 2, 3 or 4, the solution was
        found. Otherwise, the solution was not found. In either case, the
        optional output variable `mesg` gives more information.

        .. versionadded:: 1.9

    Raises
    ------
    ValueError
        if either `ydata` or `xdata` contain NaNs, or if incompatible options
        are used.

    RuntimeError
        if the least-squares minimization fails.

    TypeError
        if the number of data points is fewer than the number of parameters

    OptimizeWarning
        if covariance of the parameters can not be estimated.

    See Also
    --------
    least_squares : Minimize the sum of squares of nonlinear functions.
    scipy.stats.linregress : Calculate a linear least squares regression for
                             two sets of measurements.

    Notes
    -----
    Users should ensure that inputs `xdata`, `ydata`, and the output of `f`
    are ``float64``, or else the optimization may return incorrect results.

    With ``method='lm'``, the algorithm uses the Levenberg-Marquardt algorithm
    through `leastsq`. Note that this algorithm can only deal with
    unconstrained problems.

    Box constraints can be handled by methods 'trf' and 'dogbox'. Refer to
    the docstring of `least_squares` for more information.

    Parameters to be fitted must have similar scale. Differences of multiple
    orders of magnitude can lead to incorrect results. For the 'trf' and
    'dogbox' methods, the `x_scale` keyword argument can be used to scale
    the parameters.

    `curve_fit` is for local optimization of parameters to minimize the sum of squares
    of residuals. For global optimization, other choices of objective function, and
    other advanced features, consider using SciPy's :ref:`tutorial_optimize_global`
    tools or the `LMFIT <https://lmfit.github.io/lmfit-py/index.html>`_ package.

    References
    ----------
    .. [1] K. Vugrin et al. Confidence region estimation techniques for nonlinear
           regression in groundwater flow: Three case studies. Water Resources
           Research, Vol. 43, W03423, :doi:`10.1029/2005WR004804`

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.optimize import curve_fit

    >>> def func(x, a, b, c):
    ...     return a * np.exp(-b * x) + c

    Define the data to be fit with some noise:

    >>> xdata = np.linspace(0, 4, 50)
    >>> y = func(xdata, 2.5, 1.3, 0.5)
    >>> rng = np.random.default_rng()
    >>> y_noise = 0.2 * rng.normal(size=xdata.size)
    >>> ydata = y + y_noise
    >>> plt.plot(xdata, ydata, 'b-', label='data')

    Fit for the parameters a, b, c of the function `func`:

    >>> popt, pcov = curve_fit(func, xdata, ydata)
    >>> popt
    array([2.56274217, 1.37268521, 0.47427475])
    >>> plt.plot(xdata, func(xdata, *popt), 'r-',
    ...          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

    Constrain the optimization to the region of ``0 <= a <= 3``,
    ``0 <= b <= 1`` and ``0 <= c <= 0.5``:

    >>> popt, pcov = curve_fit(func, xdata, ydata, bounds=(0, [3., 1., 0.5]))
    >>> popt
    array([2.43736712, 1.        , 0.34463856])
    >>> plt.plot(xdata, func(xdata, *popt), 'g--',
    ...          label='fit: a=%5.3f, b=%5.3f, c=%5.3f' % tuple(popt))

    >>> plt.xlabel('x')
    >>> plt.ylabel('y')
    >>> plt.legend()
    >>> plt.show()

    For reliable results, the model `func` should not be overparametrized;
    redundant parameters can cause unreliable covariance matrices and, in some
    cases, poorer quality fits. As a quick check of whether the model may be
    overparameterized, calculate the condition number of the covariance matrix:

    >>> np.linalg.cond(pcov)
    34.571092161547405  # may vary

    The value is small, so it does not raise much concern. If, however, we were
    to add a fourth parameter ``d`` to `func` with the same effect as ``a``:

    >>> def func2(x, a, b, c, d):
    ...     return a * d * np.exp(-b * x) + c  # a and d are redundant
    >>> popt, pcov = curve_fit(func2, xdata, ydata)
    >>> np.linalg.cond(pcov)
    1.13250718925596e+32  # may vary

    Such a large value is cause for concern. The diagonal elements of the
    covariance matrix, which is related to uncertainty of the fit, gives more
    information:

    >>> np.diag(pcov)
    array([1.48814742e+29, 3.78596560e-02, 5.39253738e-03, 2.76417220e+28])  # may vary

    Note that the first and last terms are much larger than the other elements,
    suggesting that the optimal values of these parameters are ambiguous and
    that only one of these parameters is needed in the model.

    If the optimal parameters of `f` differ by multiple orders of magnitude, the
    resulting fit can be inaccurate. Sometimes, `curve_fit` can fail to find any
    results:

    >>> ydata = func(xdata, 500000, 0.01, 15)
    >>> try:
    ...     popt, pcov = curve_fit(func, xdata, ydata, method = 'trf')
    ... except RuntimeError as e:
    ...     print(e)
    Optimal parameters not found: The maximum number of function evaluations is
    exceeded.

    If parameter scale is roughly known beforehand, it can be defined in
    `x_scale` argument:

    >>> popt, pcov = curve_fit(func, xdata, ydata, method = 'trf',
    ...                        x_scale = [1000, 1, 1])
    >>> popt
    array([5.00000000e+05, 1.00000000e-02, 1.49999999e+01])
    """
    if p0 is None:
        # determine number of parameters by inspecting the function
        sig = _getfullargspec(f)
        args = sig.args
        if len(args) < 2:
            raise ValueError("Unable to determine number of fit parameters.")
        n = len(args) - 1
    else:
        p0 = np.atleast_1d(p0)
        n = p0.size

    if isinstance(bounds, Bounds):
        lb, ub = bounds.lb, bounds.ub
    else:
        lb, ub = prepare_bounds(bounds, n)
    if p0 is None:
        p0 = _initialize_feasible(lb, ub)

    bounded_problem = np.any((lb > -np.inf) | (ub < np.inf))
    if method is None:
        if bounded_problem:
            method = 'trf'
        else:
            method = 'lm'

    if method == 'lm' and bounded_problem:
        raise ValueError("Method 'lm' only works for unconstrained problems. "
                         "Use 'trf' or 'dogbox' instead.")

    if check_finite is None:
        check_finite = True if nan_policy is None else False

    # optimization may produce garbage for float32 inputs, cast them to float64
    if check_finite:
        ydata = np.asarray_chkfinite(ydata, float)
    else:
        ydata = np.asarray(ydata, float)

    if isinstance(xdata, list | tuple | np.ndarray):
        # `xdata` is passed straight to the user-defined `f`, so allow
        # non-array_like `xdata`.
        if check_finite:
            xdata = np.asarray_chkfinite(xdata, float)
        else:
            xdata = np.asarray(xdata, float)

    if ydata.size == 0:
        raise ValueError("`ydata` must not be empty!")

    # nan handling is needed only if check_finite is False because if True,
    # the x-y data are already checked, and they don't contain nans.
    if not check_finite and nan_policy is not None:
        if nan_policy == "propagate":
            msg = "`nan_policy='propagate'` is not supported by this function."
            raise ValueError(msg)
        if nan_policy not in ("raise", "omit"):
            # Override error message raised by _contains_nan
            msg = "nan_policy must be one of {None, 'raise', 'omit'}"
            raise ValueError(msg)

        x_contains_nan = _contains_nan(xdata, nan_policy)
        y_contains_nan = _contains_nan(ydata, nan_policy)

        if (x_contains_nan or y_contains_nan) and nan_policy == 'omit':
            # ignore NaNs for N dimensional arrays
            has_nan = np.isnan(xdata)
            has_nan = has_nan.any(axis=tuple(range(has_nan.ndim-1)))
            has_nan |= np.isnan(ydata)

            xdata = xdata[..., ~has_nan]
            ydata = ydata[~has_nan]

            # Also omit the corresponding entries from sigma
            if sigma is not None:
                sigma = np.asarray(sigma)
                if sigma.ndim == 1:
                    sigma = sigma[~has_nan]
                elif sigma.ndim == 2:
                    sigma = sigma[~has_nan, :]
                    sigma = sigma[:, ~has_nan]

    # Determine type of sigma
    if sigma is not None:
        sigma = np.asarray(sigma)

        # if 1-D or a scalar, sigma are errors, define transform = 1/sigma
        if sigma.size == 1 or sigma.shape == (ydata.size,):
            transform = 1.0 / sigma
        # if 2-D, sigma is the covariance matrix,
        # define transform = L such that L L^T = C
        elif sigma.shape == (ydata.size, ydata.size):
            try:
                # scipy.linalg.cholesky requires lower=True to return L L^T = A
                transform = cholesky(sigma, lower=True)
            except LinAlgError as e:
                raise ValueError("`sigma` must be positive definite.") from e
        else:
            raise ValueError("`sigma` has incorrect shape.")
    else:
        transform = None

    func = _lightweight_memoizer(_wrap_func(f, xdata, ydata, transform))

    if callable(jac):
        jac = _lightweight_memoizer(_wrap_jac(jac, xdata, transform))
    elif jac is None and method != 'lm':
        jac = '2-point'

    if 'args' in kwargs:
        # The specification for the model function `f` does not support
        # additional arguments. Refer to the `curve_fit` docstring for
        # acceptable call signatures of `f`.
        raise ValueError("'args' is not a supported keyword argument.")

    if method == 'lm':
        # if ydata.size == 1, this might be used for broadcast.
        if ydata.size != 1 and n > ydata.size:
            raise TypeError(f"The number of func parameters={n} must not"
                            f" exceed the number of data points={ydata.size}")
        res = leastsq(func, p0, Dfun=jac, full_output=True, **kwargs)
        popt, pcov, infodict, errmsg, ier = res
        ysize = len(infodict['fvec'])
        cost = np.sum(infodict['fvec'] ** 2)
        if ier not in [1, 2, 3, 4]:
            raise RuntimeError("Optimal parameters not found: " + errmsg)
    else:
        # Rename maxfev (leastsq) to max_nfev (least_squares), if specified.
        if 'max_nfev' not in kwargs:
            kwargs['max_nfev'] = kwargs.pop('maxfev', None)

        res = least_squares(func, p0, jac=jac, bounds=bounds, method=method,
                            **kwargs)

        if not res.success:
            raise RuntimeError("Optimal parameters not found: " + res.message)

        infodict = dict(nfev=res.nfev, fvec=res.fun)
        ier = res.status
        errmsg = res.message

        ysize = len(res.fun)
        cost = 2 * res.cost  # res.cost is half sum of squares!
        popt = res.x

        # Do Moore-Penrose inverse discarding zero singular values.
        _, s, VT = svd(res.jac, full_matrices=False)
        threshold = np.finfo(float).eps * max(res.jac.shape) * s[0]
        s = s[s > threshold]
        VT = VT[:s.size]
        pcov = np.dot(VT.T / s**2, VT)

    warn_cov = False
    if pcov is None or np.isnan(pcov).any():
        # indeterminate covariance
        pcov = zeros((len(popt), len(popt)), dtype=float)
        pcov.fill(inf)
        warn_cov = True
    elif not absolute_sigma:
        if ysize > p0.size:
            s_sq = cost / (ysize - p0.size)
            pcov = pcov * s_sq
        else:
            pcov.fill(inf)
            warn_cov = True

    if warn_cov:
        warnings.warn('Covariance of the parameters could not be estimated',
                      category=OptimizeWarning, stacklevel=2)

    if full_output:
        return popt, pcov, infodict, errmsg, ier
    else:
        return popt, pcov

################################################################################
# Method Path: scipy._external.pyprima.cobyla.initialize.initxfc
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/subprojects/pyprima/pyprima/pyprima/src/pyprima/cobyla/initialize.py
# Target:      inv
# Call Count:  56
# Total Exec:  56
# Func Size:   126
################################################################################

def initxfc(calcfc, iprint, maxfun, constr0, amat, bvec, ctol, f0, ftarget, rhobeg, x0,
            xhist, fhist, chist, conhist, maxhist):
    '''
    This subroutine does the initialization concerning X, function values, and
    constraints.
    '''

    # Local variables
    solver = 'COBYLA'
    srname = "INITIALIZE"

    # Sizes
    num_constraints = np.size(constr0)
    m_lcon = np.size(bvec) if bvec is not None else 0
    m_nlcon = num_constraints - m_lcon
    num_vars = np.size(x0)

    # Preconditions
    if DEBUGGING:
        assert num_constraints >= 0, f'M >= 0 {srname}'
        assert num_vars >= 1, f'N >= 1 {srname}'
        assert abs(iprint) <= 3, f'IPRINT is 0, 1, -1, 2, -2, 3, or -3 {srname}'
        # assert conmat.shape == (num_constraints , num_vars + 1), f'CONMAT.shape = [M, N+1] {srname}'
        # assert cval.size == num_vars + 1, f'CVAL.size == N+1 {srname}'
        # assert maxchist * (maxchist - maxhist) == 0, f'CHIST.shape == 0 or MAXHIST {srname}'
        # assert conhist.shape[0] == num_constraints and maxconhist * (maxconhist - maxhist) == 0, 'CONHIST.shape[0] == num_constraints, SIZE(CONHIST, 2) == 0 or MAXHIST {srname)}'
        # assert maxfhist * (maxfhist - maxhist) == 0, f'FHIST.shape == 0 or MAXHIST {srname}'
        # assert xhist.shape[0] == num_vars and maxxhist * (maxxhist - maxhist) == 0, 'XHIST.shape[0] == N, SIZE(XHIST, 2) == 0 or MAXHIST {srname)}'
        assert all(np.isfinite(x0)), f'X0 is finite {srname}'
        assert rhobeg > 0, f'RHOBEG > 0 {srname}'

    #====================#
    # Calculation starts #
    #====================#

    # Initialize info to the default value. At return, a value different from this
    # value will indicate an abnormal return
    info = INFO_DEFAULT

    # Initialize the simplex. It will be revised during the initialization.
    sim = np.eye(num_vars, num_vars+1) * rhobeg
    sim[:, num_vars] = x0

    # Initialize the matrix simi. In most cases simi is overwritten, but not always.
    simi = np.eye(num_vars) / rhobeg

    # evaluated[j] = True iff the function/constraint of SIM[:, j] has been evaluated.
    evaluated = np.zeros(num_vars+1, dtype=bool)

    # Initialize fval
    fval = np.zeros(num_vars+1) + REALMAX
    cval = np.zeros(num_vars+1) + REALMAX
    conmat = np.zeros((num_constraints, num_vars+1)) + REALMAX


    for k in range(num_vars + 1):
        x = sim[:, num_vars].copy()
        # We will evaluate F corresponding to SIM(:, J).
        if k == 0:
            j = num_vars
            f = f0
            constr = constr0
        else:
            j = k - 1
            x[j] += rhobeg
            f, constr = evaluate(calcfc, x, m_nlcon, amat, bvec)
        cstrv = np.max(np.append(0, constr))

        # Print a message about the function/constraint evaluation according to IPRINT.
        fmsg(solver, 'Initialization', iprint, k, rhobeg, f, x, cstrv, constr)

        # Save X, F, CONSTR, CSTRV into the history.
        savehist(maxhist, x, xhist, f, fhist, cstrv, chist, constr, conhist)

        # Save F, CONSTR, and CSTRV to FVAL, CONMAT, and CVAL respectively.
        evaluated[j] = True
        fval[j] = f
        conmat[:, j] = constr
        cval[j] = cstrv

        # Check whether to exit.
        subinfo = checkbreak_con(maxfun, k, cstrv, ctol, f, ftarget, x)
        if subinfo != INFO_DEFAULT:
            info = subinfo
            break

        # Exchange the new vertex of the initial simplex with the optimal vertex if necessary.
        # This is the ONLY part that is essentially non-parallel.
        if j < num_vars and fval[j] < fval[num_vars]:
            fval[j], fval[num_vars] = fval[num_vars], fval[j]
            cval[j], cval[num_vars] = cval[num_vars], cval[j]
            conmat[:, [j, num_vars]] = conmat[:, [num_vars, j]]
            sim[:, num_vars] = x
            sim[j, :j+1] = -rhobeg  # SIM[:, :j+1] is lower triangular

    nf = np.count_nonzero(evaluated)

    if evaluated.all():
        # Initialize SIMI to the inverse of SIM[:, :num_vars]
        simi = inv(sim[:, :num_vars])

    #==================#
    # Calculation ends #
    #==================#

    # Postconditions
    if DEBUGGING:
        assert nf <= maxfun, f'NF <= MAXFUN {srname}'
        assert evaluated.size == num_vars + 1, f'EVALUATED.size == Num_vars + 1 {srname}'
        # assert chist.size == maxchist, f'CHIST.size == MAXCHIST {srname}'
        # assert conhist.shape== (num_constraints, maxconhist), f'CONHIST.shape == [M, MAXCONHIST] {srname}'
        assert conmat.shape == (num_constraints, num_vars + 1), f'CONMAT.shape = [M, N+1] {srname}'
        assert not (np.isnan(conmat).any() or np.isneginf(conmat).any()), f'CONMAT does not contain NaN/-Inf {srname}'
        assert cval.size == num_vars + 1 and not (any(cval < 0) or any(np.isnan(cval)) or any(np.isposinf(cval))), f'CVAL.shape == Num_vars+1 and CVAL does not contain negative values or NaN/+Inf {srname}'
        # assert fhist.shape == maxfhist, f'FHIST.shape == MAXFHIST {srname}'
        # assert maxfhist * (maxfhist - maxhist) == 0, f'FHIST.shape == 0 or MAXHIST {srname}'
        assert fval.size == num_vars + 1 and not (any(np.isnan(fval)) or any(np.isposinf(fval))), f'FVAL.shape == Num_vars+1 and FVAL is not NaN/+Inf {srname}'
        # assert xhist.shape == (num_vars, maxxhist), f'XHIST.shape == [N, MAXXHIST] {srname}'
        assert sim.shape == (num_vars, num_vars + 1), f'SIM.shape == [N, N+1] {srname}'
        assert np.isfinite(sim).all(), f'SIM is finite {srname}'
        assert all(np.max(abs(sim[:, :num_vars]), axis=0) > 0), f'SIM(:, 1:N) has no zero column {srname}'
        assert simi.shape == (num_vars, num_vars), f'SIMI.shape == [N, N] {srname}'
        assert np.isfinite(simi).all(), f'SIMI is finite {srname}'
        assert np.allclose(sim[:, :num_vars] @ simi, np.eye(num_vars), rtol=0.1, atol=0.1) or not all(evaluated), f'SIMI = SIM(:, 1:N)^{-1} {srname}'

    return evaluated, conmat, cval, sim, simi, fval, nf, info

################################################################################
# Method Path: scipy._external.pyprima.common.linalg.inv
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/subprojects/pyprima/pyprima/pyprima/src/pyprima/common/linalg.py
# Target:      inv
# Call Count:  56
# Total Exec:  56
# Func Size:   29
################################################################################

def inv(A):
    if not USE_NAIVE_MATH:
        return np.linalg.inv(A)
    A = A.copy()
    n = A.shape[0]
    if istril(A):
        # This case is invoked in COBYLA.
        R = A.T
        B = np.zeros((n, n))
        for i in range(n):
            B[i, i] = 1 / R[i, i]
            B[:i, i] = -matprod(B[:i, :i], R[:i, i]) / R[i, i]
        return B.T
    elif istriu(A):
        B = np.zeros((n, n))
        for i in range(n):
            B[i, i] = 1 / A[i, i]
            B[:i, i] = -matprod(B[:i, :i], A[:i, i]) / A[i, i]
    else:
        # This is NOT the best algorithm for the inverse, but since the QR subroutine is available ...
        Q, R, P = qr(A)
        R = R.T
        B = np.zeros((n, n))
        for i in range(n - 1, -1, -1):
            B[:, i] = (Q[:, i] - matprod(B[:, i + 1:n], R[i + 1:n, i])) / R[i, i]
        InvP = np.zeros(n, dtype=int)
        InvP[P] = np.linspace(0, n-1, n)
        B = B[:, InvP].T
    return B

################################################################################
# Method Path: scipy.signal._lti_conversion.ss2tf
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_lti_conversion.py
# Target:      dot
# Call Count:  54
# Total Exec:  40
# Func Size:   99
################################################################################

def ss2tf(A, B, C, D, input=0):
    r"""State-space to transfer function.

    A, B, C, D defines a linear state-space system with `p` inputs,
    `q` outputs, and `n` state variables.

    Parameters
    ----------
    A : array_like
        State (or system) matrix of shape ``(n, n)``
    B : array_like
        Input matrix of shape ``(n, p)``
    C : array_like
        Output matrix of shape ``(q, n)``
    D : array_like
        Feedthrough (or feedforward) matrix of shape ``(q, p)``
    input : int, optional
        For multiple-input systems, the index of the input to use.

    Returns
    -------
    num : 2-D ndarray
        Numerator(s) of the resulting transfer function(s). `num` has one row
        for each of the system's outputs. Each row is a sequence representation
        of the numerator polynomial.
    den : 1-D ndarray
        Denominator of the resulting transfer function(s). `den` is a sequence
        representation of the denominator polynomial.

    Notes
    -----
    Before calculating `num` and `den`, the function `abcd_normalize` is called to
    convert the parameters `A`, `B`, `C`, `D` into two-dimesional arrays of the
    same dtype. The resulting dtype will be based on NumPy's dtype promotion rules,
    except in the case where each of `A`, `B`, `C`, and `D` has integer dtype, in which
    case the resulting dtype will be the default floating point dtype of ``float64``.

    The :ref:`tutorial_signal_state_space_representation` section of the
    :ref:`user_guide` presents the corresponding definitions of continuous-time and
    disrcete time state space systems.

    Examples
    --------
    Convert the state-space representation:

    .. math::

        \dot{\textbf{x}}(t) =
        \begin{bmatrix} -2 & -1 \\ 1 & 0 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \\ 0 \end{bmatrix} \textbf{u}(t) \\

        \textbf{y}(t) = \begin{bmatrix} 1 & 2 \end{bmatrix} \textbf{x}(t) +
        \begin{bmatrix} 1 \end{bmatrix} \textbf{u}(t)

    >>> A = [[-2, -1], [1, 0]]
    >>> B = [[1], [0]]  # 2-D column vector
    >>> C = [[1, 2]]    # 2-D row vector
    >>> D = 1

    to the transfer function:

    .. math:: H(s) = \frac{s^2 + 3s + 3}{s^2 + 2s + 1}

    >>> from scipy.signal import ss2tf
    >>> ss2tf(A, B, C, D)
    (array([[1., 3., 3.]]), array([ 1.,  2.,  1.]))
    """
    # transfer function is C (sI - A)**(-1) B + D

    # Check consistency and make them all floating point 2d arrays:
    A, B, C, D = abcd_normalize(A, B, C, D)

    nout, nin = D.shape
    if input >= nin:
        raise ValueError("System does not have the input specified.")

    # make SIMO from possibly MIMO system.
    B = B[:, input:input + 1]
    D = D[:, input:input + 1]

    try:
        den = poly(A)
    except ValueError:
        den = 1

    if (B.size == 0) and (C.size == 0):
        num = np.ravel(D)
        if (D.size == 0) and (A.size == 0):
            den = []
        return num, den

    num_states = A.shape[0]
    type_test = A[:, 0] + B[:, 0] + C[0, :] + D + 0.0
    num = np.empty((nout, num_states + 1), type_test.dtype)
    for k in range(nout):
        Ck = atleast_2d(C[k, :])
        num[k] = poly(A - dot(B, Ck)) + (D[k] - 1) * den

    return num, den

################################################################################
# Method Path: scipy.stats._kde.gaussian_kde.resample
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_kde.py
# Target:      dot
# Call Count:  52
# Total Exec:  28
# Func Size:   34
################################################################################

def resample(self, size=None, seed=None):
        """Randomly sample a dataset from the estimated pdf.

        Parameters
        ----------
        size : int, optional
            The number of samples to draw.  If not provided, then the size is
            the same as the effective number of samples in the underlying
            dataset.
        seed : {None, int, `numpy.random.Generator`, `numpy.random.RandomState`}, optional
            If `seed` is None (or `np.random`), the `numpy.random.RandomState`
            singleton is used.
            If `seed` is an int, a new ``RandomState`` instance is used,
            seeded with `seed`.
            If `seed` is already a ``Generator`` or ``RandomState`` instance then
            that instance is used.

        Returns
        -------
        resample : (self.d, `size`) ndarray
            The sampled dataset.

        """ # numpy/numpydoc#87  # noqa: E501
        if size is None:
            size = int(self.neff)

        random_state = check_random_state(seed)
        norm = transpose(random_state.multivariate_normal(
            zeros((self.d,), float), self.covariance, size=size
        ))
        indices = random_state.choice(self.n, size=size, p=self.weights)
        means = self.dataset[:, indices]

        return means + norm

################################################################################
# Method Path: scipy.sparse.linalg._interface.LbfgsInvHessProduct.__mul__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_interface.py
# Target:      dot
# Call Count:  50
# Total Exec:  50
# Func Size:   6
################################################################################

def __mul__(self, x):
        """Multiplication.

        Used by the ``*`` operator. Equivalent to `dot`.
        """
        return self.dot(x)

################################################################################
# Method Path: scipy.stats._multivariate.special_ortho_group_gen.rvs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      det
# Call Count:  46
# Total Exec:  50
# Func Size:   23
################################################################################

def rvs(self, dim, size=1, random_state=None):
        """Draw random samples from SO(N).

        Parameters
        ----------
        dim : int
            Dimension of rotation space (N).
        size : int, optional
            Number of samples to draw (default 1).

        Returns
        -------
        rvs : ndarray or scalar
            Random size N-dimensional matrices, dimension (size, dim, dim)

        """
        random_state = self._get_random_state(random_state)

        q = ortho_group.rvs(dim, size, random_state)
        dets = np.linalg.det(q)
        if dim:
            q[..., 0, :] /= dets[..., np.newaxis]
        return q

################################################################################
# Method Path: scipy.stats._multivariate.multivariate_t_gen.rvs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  44
# Total Exec:  22
# Func Size:   46
################################################################################

def rvs(self, loc=None, shape=1, df=1, size=1, random_state=None):
        """Draw random samples from a multivariate t-distribution.

        Parameters
        ----------
        %(_mvt_doc_default_callparams)s
        size : int, optional
            Number of samples to draw (default 1).
        %(_doc_random_state)s

        Returns
        -------
        rvs : ndarray or scalar
            Random variates of size (`size`, `P`), where `P` is the
            dimension of the random variable.

        Examples
        --------
        >>> from scipy.stats import multivariate_t
        >>> x = [0.4, 5]
        >>> loc = [0, 1]
        >>> shape = [[1, 0.1], [0.1, 1]]
        >>> df = 7
        >>> multivariate_t.rvs(loc, shape, df)
        array([[0.93477495, 3.00408716]])

        """
        # For implementation details, see equation (3):
        #
        #    Hofert, "On Sampling from the Multivariatet Distribution", 2013
        #     http://rjournal.github.io/archive/2013-2/hofert.pdf
        #
        dim, loc, shape, df = self._process_parameters(loc, shape, df)
        if random_state is not None:
            rng = check_random_state(random_state)
        else:
            rng = self._random_state

        if np.isinf(df):
            x = np.ones(size)
        else:
            x = rng.chisquare(df, size=size) / df

        z = rng.multivariate_normal(np.zeros(dim), shape, size=size)
        samples = loc + z / np.sqrt(x)[..., None]
        return _squeeze_output(samples)

################################################################################
# Method Path: scipy.linalg._matfuncs_sqrtm._sqrtm_triu
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs_sqrtm.py
# Target:      dot
# Call Count:  40
# Total Exec:  3763
# Func Size:   87
################################################################################

def _sqrtm_triu(T, blocksize=64):
    """
    Matrix square root of an upper triangular matrix.

    This is a helper function for `sqrtm` and `logm`.

    Parameters
    ----------
    T : (N, N) array_like upper triangular
        Matrix whose square root to evaluate
    blocksize : int, optional
        If the blocksize is not degenerate with respect to the
        size of the input array, then use a blocked algorithm. (Default: 64)

    Returns
    -------
    sqrtm : (N, N) ndarray
        Value of the sqrt function at `T`

    References
    ----------
    .. [1] Edvin Deadman, Nicholas J. Higham, Rui Ralha (2013)
           "Blocked Schur Algorithms for Computing the Matrix Square Root,
           Lecture Notes in Computer Science, 7782. pp. 171-182.

    """
    T_diag = np.diag(T)
    keep_it_real = np.isrealobj(T) and np.min(T_diag, initial=0.) >= 0

    # Cast to complex as necessary + ensure double precision
    if not keep_it_real:
        T = np.asarray(T, dtype=np.complex128, order="C")
        T_diag = np.asarray(T_diag, dtype=np.complex128)
    else:
        T = np.asarray(T, dtype=np.float64, order="C")
        T_diag = np.asarray(T_diag, dtype=np.float64)

    R = np.diag(np.sqrt(T_diag))

    trsyl = get_lapack_funcs('trsyl', (R,), ilp64="preferred")

    # Compute the number of blocks to use; use at least one block.
    n, n = T.shape
    nblocks = max(n // blocksize, 1)

    # Compute the smaller of the two sizes of blocks that
    # we will actually use, and compute the number of large blocks.
    bsmall, nlarge = divmod(n, nblocks)
    blarge = bsmall + 1
    nsmall = nblocks - nlarge
    if nsmall * bsmall + nlarge * blarge != n:
        raise Exception('internal inconsistency')

    # Define the index range covered by each block.
    start_stop_pairs = []
    start = 0
    for count, size in ((nsmall, bsmall), (nlarge, blarge)):
        for i in range(count):
            start_stop_pairs.append((start, start + size))
            start += size

    # Within-block interactions (Cythonized)
    try:
        within_block_loop(R, T, start_stop_pairs, nblocks)
    except RuntimeError as e:
        raise SqrtmError(*e.args) from e

    # Between-block interactions (Cython would give no significant speedup)
    for j in range(nblocks):
        jstart, jstop = start_stop_pairs[j]
        for i in range(j-1, -1, -1):
            istart, istop = start_stop_pairs[i]
            S = T[istart:istop, jstart:jstop]
            if j - i > 1:
                S = S - R[istart:istop, istop:jstart].dot(R[istop:jstart,
                                                            jstart:jstop])

            # Invoke LAPACK.
            # For more details, see the solve_sylvester implementation
            # and the fortran dtrsyl and ztrsyl docs.
            Rii = R[istart:istop, istart:istop]
            Rjj = R[jstart:jstop, jstart:jstop]
            x, scale, info = trsyl(Rii, Rjj, S)
            R[istart:istop, jstart:jstop] = x * scale

    # Return the matrix square root.
    return R

################################################################################
# Method Path: scipy.linalg._matfuncs.signm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs.py
# Target:      dot
# Call Count:  34
# Total Exec:  187
# Func Size:   69
################################################################################

def signm(A):
    """
    Matrix sign function.

    Extension of the scalar sign(x) to matrices.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix at which to evaluate the sign function

    Returns
    -------
    signm : (N, N) ndarray
        Value of the sign function at `A`

    Examples
    --------
    >>> from scipy.linalg import signm, eigvals
    >>> a = [[1,2,3], [1,2,1], [1,1,1]]
    >>> eigvals(a)
    array([ 4.12488542+0.j, -0.76155718+0.j,  0.63667176+0.j])
    >>> eigvals(signm(a))
    array([-1.+0.j,  1.+0.j,  1.+0.j])

    """
    A = _asarray_square(A)

    def rounded_sign(x):
        rx = np.real(x)
        if rx.dtype.char == 'f':
            c = 1e3*feps*amax(x)
        else:
            c = 1e3*eps*amax(x)
        return sign((absolute(rx) > c) * rx)
    result, errest = funm(A, rounded_sign, disp=0)
    errtol = {0: 1e3*feps, 1: 1e3*eps}[_array_precision[result.dtype.char]]
    if errest < errtol:
        return result

    # Handle signm of defective matrices:

    # See "E.D.Denman and J.Leyva-Ramos, Appl.Math.Comp.,
    # 8:237-250,1981" for how to improve the following (currently a
    # rather naive) iteration process:

    # a = result # sometimes iteration converges faster but where??

    # Shifting to avoid zero eigenvalues. How to ensure that shifting does
    # not change the spectrum too much?
    vals = svd(A, compute_uv=False)
    max_sv = np.amax(vals)
    # min_nonzero_sv = vals[(vals>max_sv*errtol).tolist().count(1)-1]
    # c = 0.5/min_nonzero_sv
    c = 0.5/max_sv
    S0 = A + c*np.identity(A.shape[0])
    prev_errest = errest
    for i in range(100):
        iS0 = inv(S0)
        S0 = 0.5*(S0 + iS0)
        Pp = 0.5*(dot(S0, S0)+S0)
        errest = norm(dot(Pp, Pp)-Pp, 1)
        if errest < errtol or prev_errest == errest:
            break
        prev_errest = errest
    if not isfinite(errest) or errest >= errtol:
        print("signm result may be inaccurate, approximate err =", errest)
    return S0

################################################################################
# Method Path: scipy.signal._ltisys.StateSpaceContinuous.__mul__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  34
# Total Exec:  18
# Func Size:   59
################################################################################

def __mul__(self, other):
        """
        Post-multiply another system or a scalar

        Handles multiplication of systems in the sense of a frequency domain
        multiplication. That means, given two systems E1(s) and E2(s), their
        multiplication, H(s) = E1(s) * E2(s), means that applying H(s) to U(s)
        is equivalent to first applying E2(s), and then E1(s).

        Notes
        -----
        For SISO systems the order of system application does not matter.
        However, for MIMO systems, where the two systems are matrices, the
        order above ensures standard Matrix multiplication rules apply.
        """
        if not self._check_binop_other(other):
            return NotImplemented

        if isinstance(other, StateSpace):
            # Disallow mix of discrete and continuous systems.
            if type(other) is not type(self):
                return NotImplemented

            if self.dt != other.dt:
                raise TypeError('Cannot multiply systems with different `dt`.')

            n1 = self.A.shape[0]
            n2 = other.A.shape[0]

            # Interconnection of systems
            # x1' = A1 x1 + B1 u1
            # y1  = C1 x1 + D1 u1
            # x2' = A2 x2 + B2 y1
            # y2  = C2 x2 + D2 y1
            #
            # Plugging in with u1 = y2 yields
            # [x1']   [A1 B1*C2 ] [x1]   [B1*D2]
            # [x2'] = [0  A2    ] [x2] + [B2   ] u2
            #                    [x1]
            #  y2   = [C1 D1*C2] [x2] + D1*D2 u2
            a = np.vstack((np.hstack((self.A, np.dot(self.B, other.C))),
                           np.hstack((zeros((n2, n1)), other.A))))
            b = np.vstack((np.dot(self.B, other.D), other.B))
            c = np.hstack((self.C, np.dot(self.D, other.C)))
            d = np.dot(self.D, other.D)
        else:
            # Assume that other is a scalar / matrix
            # For post multiplication the input gets scaled
            a = self.A
            b = np.dot(self.B, other)
            c = self.C
            d = np.dot(self.D, other)

        common_dtype = np.result_type(a.dtype, b.dtype, c.dtype, d.dtype)
        return StateSpace(np.asarray(a, dtype=common_dtype),
                          np.asarray(b, dtype=common_dtype),
                          np.asarray(c, dtype=common_dtype),
                          np.asarray(d, dtype=common_dtype),
                          **self._dt_dict)

################################################################################
# Method Path: scipy.sparse.linalg._matfuncs._fragment_2_1
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_matfuncs.py
# Target:      dot
# Call Count:  30
# Total Exec:  30
# Func Size:   47
################################################################################

def _fragment_2_1(X, T, s):
    """
    A helper function for expm_2009.

    Notes
    -----
    The argument X is modified in-place, but this modification is not the same
    as the returned value of the function.
    This function also takes pains to do things in ways that are compatible
    with sparse arrays, for example by avoiding fancy indexing
    and by using methods of the matrices whenever possible instead of
    using functions of the numpy or scipy libraries themselves.

    """
    # Form X = r_m(2^-s T)
    # Replace diag(X) by exp(2^-s diag(T)).
    n = X.shape[0]
    diag_T = np.ravel(T.diagonal().copy())

    # Replace diag(X) by exp(2^-s diag(T)).
    scale = 2 ** -s
    exp_diag = np.exp(scale * diag_T)
    for k in range(n):
        X[k, k] = exp_diag[k]

    for i in range(s-1, -1, -1):
        X = X.dot(X)

        # Replace diag(X) by exp(2^-i diag(T)).
        scale = 2 ** -i
        exp_diag = np.exp(scale * diag_T)
        for k in range(n):
            X[k, k] = exp_diag[k]

        # Replace (first) superdiagonal of X by explicit formula
        # for superdiagonal of exp(2^-i T) from Eq (10.42) of
        # the author's 2008 textbook
        # Functions of Matrices: Theory and Computation.
        for k in range(n-1):
            lam_1 = scale * diag_T[k]
            lam_2 = scale * diag_T[k+1]
            t_12 = scale * T[k, k+1]
            value = _eq_10_42(lam_1, lam_2, t_12)
            X[k, k+1] = value

    # Return the updated X matrix.
    return X

################################################################################
# Method Path: scipy.linalg._solvers._solve_discrete_lyapunov_bilinear
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_solvers.py
# Target:      dot
# Call Count:  30
# Total Exec:  10
# Func Size:   13
################################################################################

def _solve_discrete_lyapunov_bilinear(a, q):
    """
    Solves the discrete Lyapunov equation using a bilinear transformation.

    This function is called by the `solve_discrete_lyapunov` function with
    `method=bilinear`. It is not supposed to be called directly.
    """
    eye = np.eye(a.shape[0])
    aH = a.conj().transpose()
    aHI_inv = inv(aH + eye)
    b = np.dot(aH - eye, aHI_inv)
    c = 2*np.dot(np.dot(inv(a + eye), q), aHI_inv)
    return solve_lyapunov(b.conj().transpose(), -c)

################################################################################
# Method Path: scipy.linalg._testutils.assert_no_overwrite
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_testutils.py
# Target:      det
# Call Count:  24
# Total Exec:  33
# Func Size:   17
################################################################################

def assert_no_overwrite(call, shapes, dtypes=None):
    """
    Test that a call does not overwrite its input arguments
    """

    if dtypes is None:
        dtypes = [np.float32, np.float64, np.complex64, np.complex128]

    for dtype in dtypes:
        for order in ["C", "F"]:
            for faker in [_id, _FakeMatrix, _FakeMatrix2]:
                orig_inputs = [_get_array(s, dtype) for s in shapes]
                inputs = [faker(x.copy(order)) for x in orig_inputs]
                call(*inputs)
                msg = f"call modified inputs [{dtype!r}, {faker!r}]"
                for a, b in zip(inputs, orig_inputs):
                    np.testing.assert_equal(a, b, err_msg=msg)

################################################################################
# Method Path: scipy.linalg._testutils.assert_no_overwrite
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_testutils.py
# Target:      inv
# Call Count:  24
# Total Exec:  33
# Func Size:   17
################################################################################

def assert_no_overwrite(call, shapes, dtypes=None):
    """
    Test that a call does not overwrite its input arguments
    """

    if dtypes is None:
        dtypes = [np.float32, np.float64, np.complex64, np.complex128]

    for dtype in dtypes:
        for order in ["C", "F"]:
            for faker in [_id, _FakeMatrix, _FakeMatrix2]:
                orig_inputs = [_get_array(s, dtype) for s in shapes]
                inputs = [faker(x.copy(order)) for x in orig_inputs]
                call(*inputs)
                msg = f"call modified inputs [{dtype!r}, {faker!r}]"
                for a, b in zip(inputs, orig_inputs):
                    np.testing.assert_equal(a, b, err_msg=msg)

################################################################################
# Method Path: scipy.linalg._expm_frechet._diff_pade3
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_expm_frechet.py
# Target:      dot
# Call Count:  24
# Total Exec:  4
# Func Size:   9
################################################################################

def _diff_pade3(A, E, ident):
    b = (120., 60., 12., 1.)
    A2 = A.dot(A)
    M2 = np.dot(A, E) + np.dot(E, A)
    U = A.dot(b[3]*A2 + b[1]*ident)
    V = b[2]*A2 + b[0]*ident
    Lu = A.dot(b[3]*M2) + E.dot(b[3]*A2 + b[1]*ident)
    Lv = b[2]*M2
    return U, V, Lu, Lv

################################################################################
# Method Path: scipy.linalg._special_matrices.pascal
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_special_matrices.py
# Target:      dot
# Call Count:  23
# Total Exec:  87
# Func Size:   84
################################################################################

def pascal(n, kind='symmetric', exact=True):
    """
    Returns the n x n Pascal matrix.

    The Pascal matrix is a matrix containing the binomial coefficients as
    its elements.

    Parameters
    ----------
    n : int
        The size of the matrix to create; that is, the result is an n x n
        matrix.
    kind : str, optional
        Must be one of 'symmetric', 'lower', or 'upper'.
        Default is 'symmetric'.
    exact : bool, optional
        If `exact` is True, the result is either an array of type
        numpy.uint64 (if n < 35) or an object array of Python long integers.
        If `exact` is False, the coefficients in the matrix are computed using
        `scipy.special.comb` with ``exact=False``. The result will be a floating
        point array, and the values in the array will not be the exact
        coefficients, but this version is much faster than ``exact=True``.

    Returns
    -------
    p : (n, n) ndarray
        The Pascal matrix.

    See Also
    --------
    invpascal

    Notes
    -----
    See https://en.wikipedia.org/wiki/Pascal_matrix for more information
    about Pascal matrices.

    .. versionadded:: 0.11.0

    Examples
    --------
    >>> from scipy.linalg import pascal
    >>> pascal(4)
    array([[ 1,  1,  1,  1],
           [ 1,  2,  3,  4],
           [ 1,  3,  6, 10],
           [ 1,  4, 10, 20]], dtype=uint64)
    >>> pascal(4, kind='lower')
    array([[1, 0, 0, 0],
           [1, 1, 0, 0],
           [1, 2, 1, 0],
           [1, 3, 3, 1]], dtype=uint64)
    >>> pascal(50)[-1, -1]
    25477612258980856902730428600
    >>> from scipy.special import comb
    >>> comb(98, 49, exact=True)
    25477612258980856902730428600

    """

    from scipy.special import comb
    if kind not in ['symmetric', 'lower', 'upper']:
        raise ValueError("kind must be 'symmetric', 'lower', or 'upper'")

    if exact:
        if n >= 35:
            L_n = np.empty((n, n), dtype=object)
            L_n.fill(0)
        else:
            L_n = np.zeros((n, n), dtype=np.uint64)
        for i in range(n):
            for j in range(i + 1):
                L_n[i, j] = comb(i, j, exact=True)
    else:
        L_n = comb(*np.ogrid[:n, :n])

    if kind == 'lower':
        p = L_n
    elif kind == 'upper':
        p = L_n.T
    else:
        p = np.dot(L_n, L_n.T)

    return p

################################################################################
# Method Path: scipy.stats._multivariate.random_correlation_gen.rvs
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_multivariate.py
# Target:      dot
# Call Count:  20
# Total Exec:  17
# Func Size:   39
################################################################################

def rvs(self, eigs, random_state=None, tol=1e-13, diag_tol=1e-7):
        """Draw random correlation matrices.

        Parameters
        ----------
        eigs : 1d ndarray
            Eigenvalues of correlation matrix
        tol : float, optional
            Tolerance for input parameter checks
        diag_tol : float, optional
            Tolerance for deviation of the diagonal of the resulting
            matrix. Default: 1e-7

        Raises
        ------
        RuntimeError
            Floating point error prevented generating a valid correlation
            matrix.

        Returns
        -------
        rvs : ndarray or scalar
            Random size N-dimensional matrices, dimension (size, dim, dim),
            each having eigenvalues eigs.

        """
        dim, eigs = self._process_parameters(eigs, tol=tol)

        random_state = self._get_random_state(random_state)

        m = ortho_group.rvs(dim, random_state=random_state)
        m = np.dot(np.dot(m, np.diag(eigs)), m.T)  # Set the trace of m
        m = self._to_corr(m)  # Carefully rotate to unit diagonal

        # Check diagonal
        if abs(m.diagonal() - 1).max() > diag_tol:
            raise RuntimeError("Failed to generate a valid correlation matrix")

        return m

################################################################################
# Method Path: scipy.linalg._solvers._solve_discrete_lyapunov_bilinear
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_solvers.py
# Target:      inv
# Call Count:  20
# Total Exec:  10
# Func Size:   13
################################################################################

def _solve_discrete_lyapunov_bilinear(a, q):
    """
    Solves the discrete Lyapunov equation using a bilinear transformation.

    This function is called by the `solve_discrete_lyapunov` function with
    `method=bilinear`. It is not supposed to be called directly.
    """
    eye = np.eye(a.shape[0])
    aH = a.conj().transpose()
    aHI_inv = inv(aH + eye)
    b = np.dot(aH - eye, aHI_inv)
    c = 2*np.dot(np.dot(inv(a + eye), q), aHI_inv)
    return solve_lyapunov(b.conj().transpose(), -c)

################################################################################
# Method Path: scipy.signal._signaltools._filtfilt_gust
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_signaltools.py
# Target:      dot
# Call Count:  19
# Total Exec:  20
# Func Size:   177
################################################################################

def _filtfilt_gust(b, a, x, axis=-1, irlen=None):
    """Forward-backward IIR filter that uses Gustafsson's method.

    Apply the IIR filter defined by ``(b,a)`` to `x` twice, first forward
    then backward, using Gustafsson's initial conditions [1]_.

    Let ``y_fb`` be the result of filtering first forward and then backward,
    and let ``y_bf`` be the result of filtering first backward then forward.
    Gustafsson's method is to compute initial conditions for the forward
    pass and the backward pass such that ``y_fb == y_bf``.

    Parameters
    ----------
    b : scalar or 1-D ndarray
        Numerator coefficients of the filter.
    a : scalar or 1-D ndarray
        Denominator coefficients of the filter.
    x : ndarray
        Data to be filtered.
    axis : int, optional
        Axis of `x` to be filtered.  Default is -1.
    irlen : int or None, optional
        The length of the nonnegligible part of the impulse response.
        If `irlen` is None, or if the length of the signal is less than
        ``2 * irlen``, then no part of the impulse response is ignored.

    Returns
    -------
    y : ndarray
        The filtered data.
    x0 : ndarray
        Initial condition for the forward filter.
    x1 : ndarray
        Initial condition for the backward filter.

    Notes
    -----
    Typically the return values `x0` and `x1` are not needed by the
    caller.  The intended use of these return values is in unit tests.

    References
    ----------
    .. [1] F. Gustaffson. Determining the initial states in forward-backward
           filtering. Transactions on Signal Processing, 46(4):988-992, 1996.

    """
    # In the comments, "Gustafsson's paper" and [1] refer to the
    # paper referenced in the docstring.

    b = np.atleast_1d(b)
    a = np.atleast_1d(a)

    order = max(len(b), len(a)) - 1
    if order == 0:
        # The filter is just scalar multiplication, with no state.
        scale = (b[0] / a[0])**2
        y = scale * x
        return y, np.array([]), np.array([])

    if axis != -1 or axis != x.ndim - 1:
        # Move the axis containing the data to the end.
        x = np.swapaxes(x, axis, x.ndim - 1)

    # n is the number of samples in the data to be filtered.
    n = x.shape[-1]

    if irlen is None or n <= 2*irlen:
        m = n
    else:
        m = irlen

    # Create Obs, the observability matrix (called O in the paper).
    # This matrix can be interpreted as the operator that propagates
    # an arbitrary initial state to the output, assuming the input is
    # zero.
    # In Gustafsson's paper, the forward and backward filters are not
    # necessarily the same, so he has both O_f and O_b.  We use the same
    # filter in both directions, so we only need O. The same comment
    # applies to S below.
    Obs = np.zeros((m, order))
    zi = np.zeros(order)
    zi[0] = 1
    Obs[:, 0] = lfilter(b, a, np.zeros(m), zi=zi)[0]
    for k in range(1, order):
        Obs[k:, k] = Obs[:-k, 0]

    # Obsr is O^R (Gustafsson's notation for row-reversed O)
    Obsr = Obs[::-1]

    # Create S.  S is the matrix that applies the filter to the reversed
    # propagated initial conditions.  That is,
    #     out = S.dot(zi)
    # is the same as
    #     tmp, _ = lfilter(b, a, zeros(), zi=zi)  # Propagate ICs.
    #     out = lfilter(b, a, tmp[::-1])          # Reverse and filter.

    # Equations (5) & (6) of [1]
    S = lfilter(b, a, Obs[::-1], axis=0)

    # Sr is S^R (row-reversed S)
    Sr = S[::-1]

    # M is [(S^R - O), (O^R - S)]
    if m == n:
        M = np.hstack((Sr - Obs, Obsr - S))
    else:
        # Matrix described in section IV of [1].
        M = np.zeros((2*m, 2*order))
        M[:m, :order] = Sr - Obs
        M[m:, order:] = Obsr - S

    # Naive forward-backward and backward-forward filters.
    # These have large transients because the filters use zero initial
    # conditions.
    y_f = lfilter(b, a, x)
    y_fb = lfilter(b, a, y_f[..., ::-1])[..., ::-1]

    y_b = lfilter(b, a, x[..., ::-1])[..., ::-1]
    y_bf = lfilter(b, a, y_b)

    delta_y_bf_fb = y_bf - y_fb
    if m == n:
        delta = delta_y_bf_fb
    else:
        start_m = delta_y_bf_fb[..., :m]
        end_m = delta_y_bf_fb[..., -m:]
        delta = np.concatenate((start_m, end_m), axis=-1)

    # ic_opt holds the "optimal" initial conditions.
    # The following code computes the result shown in the formula
    # of the paper between equations (6) and (7).
    if delta.ndim == 1:
        ic_opt = linalg.lstsq(M, delta)[0]
    else:
        # Reshape delta so it can be used as an array of multiple
        # right-hand-sides in linalg.lstsq.
        delta2d = delta.reshape(-1, delta.shape[-1]).T
        ic_opt0 = linalg.lstsq(M, delta2d)[0].T
        ic_opt = ic_opt0.reshape(delta.shape[:-1] + (M.shape[-1],))

    # Now compute the filtered signal using equation (7) of [1].
    # First, form [S^R, O^R] and call it W.
    if m == n:
        W = np.hstack((Sr, Obsr))
    else:
        W = np.zeros((2*m, 2*order))
        W[:m, :order] = Sr
        W[m:, order:] = Obsr

    # Equation (7) of [1] says
    #     Y_fb^opt = Y_fb^0 + W * [x_0^opt; x_{N-1}^opt]
    # `wic` is (almost) the product on the right.
    # W has shape (m, 2*order), and ic_opt has shape (..., 2*order),
    # so we can't use W.dot(ic_opt).  Instead, we dot ic_opt with W.T,
    # so wic has shape (..., m).
    wic = ic_opt.dot(W.T)

    # `wic` is "almost" the product of W and the optimal ICs in equation
    # (7)--if we're using a truncated impulse response (m < n), `wic`
    # contains only the adjustments required for the ends of the signal.
    # Here we form y_opt, taking this into account if necessary.
    y_opt = y_fb
    if m == n:
        y_opt += wic
    else:
        y_opt[..., :m] += wic[..., :m]
        y_opt[..., -m:] += wic[..., -m:]

    x0 = ic_opt[..., :order]
    x1 = ic_opt[..., -order:]
    if axis != -1 or axis != x.ndim - 1:
        # Restore the data axis to its original position.
        x0 = np.swapaxes(x0, axis, x.ndim - 1)
        x1 = np.swapaxes(x1, axis, x.ndim - 1)
        y_opt = np.swapaxes(y_opt, axis, x.ndim - 1)

    return y_opt, x0, x1

################################################################################
# Method Path: scipy.optimize._linesearch.line_search_armijo
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linesearch.py
# Target:      dot
# Call Count:  18
# Total Exec:  18
# Func Size:   50
################################################################################

def line_search_armijo(f, xk, pk, gfk, old_fval, args=(), c1=1e-4, alpha0=1):
    """Minimize over alpha, the function ``f(xk+alpha pk)``.

    Parameters
    ----------
    f : callable
        Function to be minimized.
    xk : array_like
        Current point.
    pk : array_like
        Search direction.
    gfk : array_like
        Gradient of `f` at point `xk`.
    old_fval : float
        Value of `f` at point `xk`.
    args : tuple, optional
        Optional arguments.
    c1 : float, optional
        Value to control stopping criterion.
    alpha0 : scalar, optional
        Value of `alpha` at start of the optimization.

    Returns
    -------
    alpha
    f_count
    f_val_at_alpha

    Notes
    -----
    Uses the interpolation algorithm (Armijo backtracking) as suggested by
    Wright and Nocedal in 'Numerical Optimization', 1999, pp. 56-57

    """
    xk = np.atleast_1d(xk)
    fc = [0]

    def phi(alpha1):
        fc[0] += 1
        return f(xk + alpha1*pk, *args)

    if old_fval is None:
        phi0 = phi(0.)
    else:
        phi0 = old_fval  # compute f(xk) -- done in past loop

    derphi0 = np.dot(gfk, pk)
    alpha, phi1 = scalar_search_armijo(phi, phi0, derphi0, c1=c1,
                                       alpha0=alpha0)
    return alpha, fc[0], phi1

################################################################################
# Method Path: scipy.linalg._matfuncs.signm
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_matfuncs.py
# Target:      inv
# Call Count:  17
# Total Exec:  187
# Func Size:   69
################################################################################

def signm(A):
    """
    Matrix sign function.

    Extension of the scalar sign(x) to matrices.

    Parameters
    ----------
    A : (N, N) array_like
        Matrix at which to evaluate the sign function

    Returns
    -------
    signm : (N, N) ndarray
        Value of the sign function at `A`

    Examples
    --------
    >>> from scipy.linalg import signm, eigvals
    >>> a = [[1,2,3], [1,2,1], [1,1,1]]
    >>> eigvals(a)
    array([ 4.12488542+0.j, -0.76155718+0.j,  0.63667176+0.j])
    >>> eigvals(signm(a))
    array([-1.+0.j,  1.+0.j,  1.+0.j])

    """
    A = _asarray_square(A)

    def rounded_sign(x):
        rx = np.real(x)
        if rx.dtype.char == 'f':
            c = 1e3*feps*amax(x)
        else:
            c = 1e3*eps*amax(x)
        return sign((absolute(rx) > c) * rx)
    result, errest = funm(A, rounded_sign, disp=0)
    errtol = {0: 1e3*feps, 1: 1e3*eps}[_array_precision[result.dtype.char]]
    if errest < errtol:
        return result

    # Handle signm of defective matrices:

    # See "E.D.Denman and J.Leyva-Ramos, Appl.Math.Comp.,
    # 8:237-250,1981" for how to improve the following (currently a
    # rather naive) iteration process:

    # a = result # sometimes iteration converges faster but where??

    # Shifting to avoid zero eigenvalues. How to ensure that shifting does
    # not change the spectrum too much?
    vals = svd(A, compute_uv=False)
    max_sv = np.amax(vals)
    # min_nonzero_sv = vals[(vals>max_sv*errtol).tolist().count(1)-1]
    # c = 0.5/min_nonzero_sv
    c = 0.5/max_sv
    S0 = A + c*np.identity(A.shape[0])
    prev_errest = errest
    for i in range(100):
        iS0 = inv(S0)
        S0 = 0.5*(S0 + iS0)
        Pp = 0.5*(dot(S0, S0)+S0)
        errest = norm(dot(Pp, Pp)-Pp, 1)
        if errest < errtol or prev_errest == errest:
            break
        prev_errest = errest
    if not isfinite(errest) or errest >= errtol:
        print("signm result may be inaccurate, approximate err =", errest)
    return S0

################################################################################
# Method Path: scipy.sparse._coo.coo_array._dense_tensordot
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/_coo.py
# Target:      dot
# Call Count:  16
# Total Exec:  16
# Func Size:   20
################################################################################

def _dense_tensordot(self, other, s_axes, o_axes):
        s_ndim = len(self.shape)
        o_ndim = len(other.shape)

        s_non_axes = [i for i in range(s_ndim) if i not in s_axes]
        s_axes_shape = [self.shape[i] for i in s_axes]
        s_non_axes_shape = [self.shape[i] for i in s_non_axes]

        o_non_axes = [i for i in range(o_ndim) if i not in o_axes]
        o_axes_shape = [other.shape[i] for i in o_axes]
        o_non_axes_shape = [other.shape[i] for i in o_non_axes]

        left = self.transpose(s_non_axes + s_axes)
        right = np.transpose(other, o_non_axes[:-1] + o_axes + o_non_axes[-1:])

        reshape_left = (*s_non_axes_shape, math.prod(s_axes_shape))
        reshape_right = (*o_non_axes_shape[:-1], math.prod(o_axes_shape),
                         *o_non_axes_shape[-1:])

        return left.reshape(reshape_left).dot(right.reshape(reshape_right))

################################################################################
# Method Path: scipy.linalg._basic.inv
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/linalg/_basic.py
# Target:      inv
# Call Count:  15
# Total Exec:  300
# Func Size:   142
################################################################################

def inv(a, overwrite_a=False, check_finite=True, *, assume_a=None, lower=False):
    r"""
    Compute the inverse of a matrix.

    If the data matrix is known to be a particular type then supplying the
    corresponding string to ``assume_a`` key chooses the dedicated solver.
    The available options are

    =============================  ================================
     general                        'general' (or 'gen')
     diagonal                       'diagonal'
     upper triangular               'upper triangular'
     lower triangular               'lower triangular'
     symmetric positive definite    'pos'
     symmetric                      'sym'
     Hermitian                      'her'
    =============================  ================================

    For the 'pos' option, only the triangle of the input matrix specified in
    the `lower` argument is used, and the other triangle is not referenced.
    Likewise, an explicit `assume_a='diagonal'` means that off-diagonal elements
    are not referenced.

    The `a` array argument may have additional "batch" dimensions prepended to the core
    shape. In this case, the array is treated as a batch of lower-dimensional slices;
    see :ref:`linalg_batch` for details.

    Parameters
    ----------
    a : array_like, shape (..., M, M)
        Square matrix (or a batch of matrices) to be inverted.
    overwrite_a : bool, optional
        Discard data in `a` (may improve performance). Default is False.
        See :ref:`tutorial_linalg_overwrite` for details.
    check_finite : bool, optional
        Whether to check that the input matrix contains only finite numbers.
        Disabling may give a performance gain, but may result in problems
        (crashes, non-termination) if the inputs do contain infinities or NaNs.
    assume_a : str, optional
        Valid entries are described above.
        If omitted or ``None``, checks are performed to identify structure so the
        appropriate solver can be called.
    lower : bool, optional
        Ignored unless `assume_a` is one of 'sym', 'her', or 'pos'. If True, the
        calculation uses only the data in the lower triangle of `a`; entries above the
        diagonal are ignored. If False (default), the calculation uses only the data in
        the upper triangle of `a`; entries below the diagonal are ignored.

    Returns
    -------
    ainv : ndarray
        Inverse of the matrix `a`.

    Raises
    ------
    LinAlgError
        If `a` is singular.
    ValueError
        If `a` is not square, or not 2D.

    Notes
    -----

    This routine checks the condition number of the `a` matrix and emits a
    `LinAlgWarning` for ill-conditioned inputs.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy import linalg
    >>> a = np.array([[1., 2.], [3., 4.]])
    >>> linalg.inv(a)
    array([[-2. ,  1. ],
           [ 1.5, -0.5]])
    >>> a @ linalg.inv(a)
    array([[ 1.,  0.],
           [ 0.,  1.]])

    The input array ``a`` may represent a single matrix or a collection (a.k.a.
    a "batch") of square matrices. For example, if ``a.shape == (4, 3, 2, 2)``, it is
    interpreted as a ``(4, 3)``-shaped batch of :math:`2\times 2` matrices.
    See :ref:`linalg_batch` for further details.
    To illustrate:

    >>> a = np.stack((np.eye(2), [[1, 2], [3, 4]]))
    >>> linalg.inv(a)
    array([[[ 1. ,  0. ],
            [ 0. ,  1. ]],
           [[-2. ,  1. ],
            [ 1.5, -0.5]]])

    Note that the structure detection runs per-slice: in the example above, each of the
    two slices will be independently discovered as being diagonal. Setting an explicit
    ``assume_a`` argument will bypass structure detection and use the provided value
    without checking:

    >>> a = np.stack((np.eye(2), [[1, 2], [3, 4]]))
    >>> linalg.inv(a, assume_a="diagonal")
    array([[[1.  , 0.  ],
            [0.  , 1.  ]],
           [[1.  , 2.  ],   # off-diagonal elements are incorrect
            [3.  , 0.25]]])
    """
    a1 = _asarray_validated(a, check_finite=check_finite)
    _deprecate_dtypes("linalg.inv", a1)

    if a1.ndim < 2:
        raise ValueError(f"Expected at least ndim=2, got {a1.ndim=}")
    if a1.shape[-1] != a1.shape[-2]:
        raise ValueError(f"Expected square matrix, got {a1.shape=}")

    # accommodate empty matrices
    if a1.size == 0:
        dt = inv(np.eye(2, dtype=a1.dtype)).dtype
        return np.empty_like(a1, dtype=dt)

    # Also check if dtype is LAPACK compatible
    a1, overwrite_a = _normalize_lapack_dtype(a1, overwrite_a)
    a1, overwrite_a = _ensure_aligned_and_native(a1, overwrite_a)

    # XXX can relax a1.ndim == 2?
    overwrite_a = overwrite_a and (a1.ndim == 2) and (a1.flags["F_CONTIGUOUS"])

    # keep the numbers in sync with C at `linalg/src/_common_array_utils.hh`
    structure = {
        None: -1,
        'general': 0, 'gen': 0,
        'diagonal': 11,
        'upper triangular': 21,
        'lower triangular': 22,
        'pos' : 101,
        'sym' : 201,
        'her' : 211,
    }[assume_a]

    # a1 is well behaved, invert it.
    inv_a, err_lst = _batched_linalg._inv(a1, structure, overwrite_a, lower)

    if err_lst:
        _format_emit_errors_warnings(err_lst)

    return inv_a

################################################################################
# Method Path: scipy.optimize._linprog_ip._ip_hsd
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_linprog_ip.py
# Target:      dot
# Call Count:  14
# Total Exec:  155
# Func Size:   265
################################################################################

def _ip_hsd(A, b, c, c0, alpha0, beta, maxiter, disp, tol, sparse, lstsq,
            sym_pos, cholesky, pc, ip, permc_spec, callback, postsolve_args):
    r"""
    Solve a linear programming problem in standard form:

    Minimize::

        c @ x

    Subject to::

        A @ x == b
            x >= 0

    using the interior point method of [4].

    Parameters
    ----------
    A : 2-D array
        2-D array such that ``A @ x``, gives the values of the equality
        constraints at ``x``.
    b : 1-D array
        1-D array of values representing the RHS of each equality constraint
        (row) in ``A`` (for standard form problem).
    c : 1-D array
        Coefficients of the linear objective function to be minimized (for
        standard form problem).
    c0 : float
        Constant term in objective function due to fixed (and eliminated)
        variables. (Purely for display.)
    alpha0 : float
        The maximal step size for Mehrota's predictor-corrector search
        direction; see :math:`\beta_3`of [4] Table 8.1
    beta : float
        The desired reduction of the path parameter :math:`\mu` (see  [6]_)
    maxiter : int
        The maximum number of iterations of the algorithm.
    disp : bool
        Set to ``True`` if indicators of optimization status are to be printed
        to the console each iteration.
    tol : float
        Termination tolerance; see [4]_ Section 4.5.
    sparse : bool
        Set to ``True`` if the problem is to be treated as sparse. However,
        the inputs ``A_eq`` and ``A_ub`` should nonetheless be provided as
        (dense) arrays rather than sparse matrices.
    lstsq : bool
        Set to ``True`` if the problem is expected to be very poorly
        conditioned. This should always be left as ``False`` unless severe
        numerical difficulties are frequently encountered, and a better option
        would be to improve the formulation of the problem.
    sym_pos : bool
        Leave ``True`` if the problem is expected to yield a well conditioned
        symmetric positive definite normal equation matrix (almost always).
    cholesky : bool
        Set to ``True`` if the normal equations are to be solved by explicit
        Cholesky decomposition followed by explicit forward/backward
        substitution. This is typically faster for moderate, dense problems
        that are numerically well-behaved.
    pc : bool
        Leave ``True`` if the predictor-corrector method of Mehrota is to be
        used. This is almost always (if not always) beneficial.
    ip : bool
        Set to ``True`` if the improved initial point suggestion due to [4]_
        Section 4.3 is desired. It's unclear whether this is beneficial.
    permc_spec : str (default = 'MMD_AT_PLUS_A')
        (Has effect only with ``sparse = True``, ``lstsq = False``, ``sym_pos =
        True``.) A matrix is factorized in each iteration of the algorithm.
        This option specifies how to permute the columns of the matrix for
        sparsity preservation. Acceptable values are:

        - ``NATURAL``: natural ordering.
        - ``MMD_ATA``: minimum degree ordering on the structure of A^T A.
        - ``MMD_AT_PLUS_A``: minimum degree ordering on the structure of A^T+A.
        - ``COLAMD``: approximate minimum degree column ordering.

        This option can impact the convergence of the
        interior point algorithm; test different values to determine which
        performs best for your problem. For more information, refer to
        ``scipy.sparse.linalg.splu``.
    callback : callable, optional
        If a callback function is provided, it will be called within each
        iteration of the algorithm. The callback function must accept a single
        `scipy.optimize.OptimizeResult` consisting of the following fields:

            x : 1-D array
                Current solution vector
            fun : float
                Current value of the objective function
            success : bool
                True only when an algorithm has completed successfully,
                so this is always False as the callback function is called
                only while the algorithm is still iterating.
            slack : 1-D array
                The values of the slack variables. Each slack variable
                corresponds to an inequality constraint. If the slack is zero,
                the corresponding constraint is active.
            con : 1-D array
                The (nominally zero) residuals of the equality constraints,
                that is, ``b - A_eq @ x``
            phase : int
                The phase of the algorithm being executed. This is always
                1 for the interior-point method because it has only one phase.
            status : int
                For revised simplex, this is always 0 because if a different
                status is detected, the algorithm terminates.
            nit : int
                The number of iterations performed.
            message : str
                A string descriptor of the exit status of the optimization.
    postsolve_args : tuple
        Data needed by _postsolve to convert the solution to the standard-form
        problem into the solution to the original problem.

    Returns
    -------
    x_hat : float
        Solution vector (for standard form problem).
    status : int
        An integer representing the exit status of the optimization::

         0 : Optimization terminated successfully
         1 : Iteration limit reached
         2 : Problem appears to be infeasible
         3 : Problem appears to be unbounded
         4 : Serious numerical difficulties encountered

    message : str
        A string descriptor of the exit status of the optimization.
    iteration : int
        The number of iterations taken to solve the problem

    References
    ----------
    .. [4] Andersen, Erling D., and Knud D. Andersen. "The MOSEK interior point
           optimizer for linear programming: an implementation of the
           homogeneous algorithm." High performance optimization. Springer US,
           2000. 197-232.
    .. [6] Freund, Robert M. "Primal-Dual Interior-Point Methods for Linear
           Programming based on Newton's Method." Unpublished Course Notes,
           March 2004. Available 2/25/2017 at:
           https://ocw.mit.edu/courses/sloan-school-of-management/15-084j-nonlinear-programming-spring-2004/lecture-notes/lec14_int_pt_mthd.pdf

    """

    iteration = 0

    # default initial point
    x, y, z, tau, kappa = _get_blind_start(A.shape)

    # first iteration is special improvement of initial point
    ip = ip if pc else False

    # [4] 4.5
    rho_p, rho_d, rho_A, rho_g, rho_mu, obj = _indicators(
        A, b, c, c0, x, y, z, tau, kappa)
    go = rho_p > tol or rho_d > tol or rho_A > tol  # we might get lucky : )

    if disp:
        _display_iter(rho_p, rho_d, rho_g, "-", rho_mu, obj, header=True)
    if callback is not None:
        x_o, fun, slack, con = _postsolve(x/tau, postsolve_args)
        res = OptimizeResult({'x': x_o, 'fun': fun, 'slack': slack,
                              'con': con, 'nit': iteration, 'phase': 1,
                              'complete': False, 'status': 0,
                              'message': "", 'success': False})
        callback(res)

    status = 0
    message = "Optimization terminated successfully."

    if sparse:
        A = sps.csc_array(A)

    while go:

        iteration += 1

        if ip:  # initial point
            # [4] Section 4.4
            gamma = 1

            def eta(g):
                return 1
        else:
            # gamma = 0 in predictor step according to [4] 4.1
            # if predictor/corrector is off, use mean of complementarity [6]
            # 5.1 / [4] Below Figure 10-4
            gamma = 0 if pc else beta * np.mean(z * x)
            # [4] Section 4.1

            def eta(g=gamma):
                return 1 - g

        try:
            # Solve [4] 8.6 and 8.7/8.13/8.23
            d_x, d_y, d_z, d_tau, d_kappa = _get_delta(
                A, b, c, x, y, z, tau, kappa, gamma, eta,
                sparse, lstsq, sym_pos, cholesky, pc, ip, permc_spec)

            if ip:  # initial point
                # [4] 4.4
                # Formula after 8.23 takes a full step regardless if this will
                # take it negative
                alpha = 1.0
                x, y, z, tau, kappa = _do_step(
                    x, y, z, tau, kappa, d_x, d_y,
                    d_z, d_tau, d_kappa, alpha)
                x[x < 1] = 1
                z[z < 1] = 1
                tau = max(1, tau)
                kappa = max(1, kappa)
                ip = False  # done with initial point
            else:
                # [4] Section 4.3
                alpha = _get_step(x, d_x, z, d_z, tau,
                                  d_tau, kappa, d_kappa, alpha0)
                # [4] Equation 8.9
                x, y, z, tau, kappa = _do_step(
                    x, y, z, tau, kappa, d_x, d_y, d_z, d_tau, d_kappa, alpha)

        except (LinAlgError, FloatingPointError,
                ValueError, ZeroDivisionError):
            # this can happen when sparse solver is used and presolve
            # is turned off. Also observed ValueError in AppVeyor Python 3.6
            # Win32 build (PR #8676). I've never seen it otherwise.
            status = 4
            message = _get_message(status)
            break

        # [4] 4.5
        rho_p, rho_d, rho_A, rho_g, rho_mu, obj = _indicators(
            A, b, c, c0, x, y, z, tau, kappa)
        go = rho_p > tol or rho_d > tol or rho_A > tol

        if disp:
            _display_iter(rho_p, rho_d, rho_g, alpha, rho_mu, obj)
        if callback is not None:
            x_o, fun, slack, con = _postsolve(x/tau, postsolve_args)
            res = OptimizeResult({'x': x_o, 'fun': fun, 'slack': slack,
                                  'con': con, 'nit': iteration, 'phase': 1,
                                  'complete': False, 'status': 0,
                                  'message': "", 'success': False})
            callback(res)

        # [4] 4.5
        inf1 = (rho_p < tol and rho_d < tol and rho_g < tol and tau < tol *
                max(1, kappa))
        inf2 = rho_mu < tol and tau < tol * min(1, kappa)
        if inf1 or inf2:
            # [4] Lemma 8.4 / Theorem 8.3
            if b.transpose().dot(y) > tol:
                status = 2
            else:  # elif c.T.dot(x) < tol: ? Probably not necessary.
                status = 3
            message = _get_message(status)
            break
        elif iteration >= maxiter:
            status = 1
            message = _get_message(status)
            break

    x_hat = x / tau
    # [4] Statement after Theorem 8.2
    return x_hat, status, message, iteration

################################################################################
# Method Path: scipy.optimize._nonlin.BroydenFirst.todense
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_nonlin.py
# Target:      inv
# Call Count:  14
# Total Exec:  14
# Func Size:   2
################################################################################

def todense(self):
        return inv(self.Gm)

################################################################################
# Method Path: scipy.signal._ltisys.StateSpaceContinuous.__rmul__
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  14
# Total Exec:  9
# Func Size:   17
################################################################################

def __rmul__(self, other):
        """Pre-multiply a scalar or matrix (but not StateSpace)"""
        if not self._check_binop_other(other) or isinstance(other, StateSpace):
            return NotImplemented

        # For pre-multiplication only the output gets scaled
        a = self.A
        b = self.B
        c = np.dot(other, self.C)
        d = np.dot(other, self.D)

        common_dtype = np.result_type(a.dtype, b.dtype, c.dtype, d.dtype)
        return StateSpace(np.asarray(a, dtype=common_dtype),
                          np.asarray(b, dtype=common_dtype),
                          np.asarray(c, dtype=common_dtype),
                          np.asarray(d, dtype=common_dtype),
                          **self._dt_dict)

################################################################################
# Method Path: scipy.signal._ltisys._KNV0_loop
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      det
# Call Count:  14
# Total Exec:  1
# Func Size:   24
################################################################################

def _KNV0_loop(ker_pole, transfer_matrix, poles, B, maxiter, rtol):
    """
    Loop over all poles one by one and apply KNV method 0 algorithm
    """
    # This method is useful only because we need to be able to call
    # _KNV0 from YT without looping over all poles, otherwise it would
    # have been fine to mix _KNV0_loop and _KNV0 in a single function
    stop = False
    nb_try = 0
    while nb_try < maxiter and not stop:
        det_transfer_matrixb = np.abs(np.linalg.det(transfer_matrix))
        for j in range(B.shape[0]):
            _KNV0(B, ker_pole, transfer_matrix, j, poles)

        det_transfer_matrix = np.max((np.sqrt(np.spacing(1)),
                                  np.abs(np.linalg.det(transfer_matrix))))
        cur_rtol = np.abs((det_transfer_matrix - det_transfer_matrixb) /
                       det_transfer_matrix)
        if cur_rtol < rtol and det_transfer_matrix > np.sqrt(np.spacing(1)):
            # Convergence test from YT page 21
            stop = True

        nb_try += 1
    return stop, cur_rtol, nb_try

################################################################################
# Method Path: scipy.signal._fir_filter_design.firls
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_fir_filter_design.py
# Target:      dot
# Call Count:  12
# Total Exec:  17
# Func Size:   218
################################################################################

def firls(numtaps, bands, desired, *, weight=None, fs=None):
    """
    FIR filter design using least-squares error minimization.

    Calculate the filter coefficients for the linear-phase finite
    impulse response (FIR) filter which has the best approximation
    to the desired frequency response described by `bands` and
    `desired` in the least squares sense (i.e., the integral of the
    weighted mean-squared error within the specified bands is
    minimized).

    Parameters
    ----------
    numtaps : int
        The number of taps in the FIR filter. `numtaps` must be odd.
    bands : array_like
        A monotonic nondecreasing sequence containing the band edges in
        Hz. All elements must be non-negative and less than or equal to
        the Nyquist frequency given by `nyq`. The bands are specified as
        frequency pairs, thus, if using a 1D array, its length must be
        even, e.g., `np.array([0, 1, 2, 3, 4, 5])`. Alternatively, the
        bands can be specified as an nx2 sized 2D array, where n is the
        number of bands, e.g, `np.array([[0, 1], [2, 3], [4, 5]])`.
    desired : array_like
        A sequence the same size as `bands` containing the desired gain
        at the start and end point of each band.
    weight : array_like, optional
        A relative weighting to give to each band region when solving
        the least squares problem. `weight` has to be half the size of
        `bands`.
    fs : float, optional
        The sampling frequency of the signal. Each frequency in `bands`
        must be between 0 and ``fs/2`` (inclusive). Default is 2.

    Returns
    -------
    coeffs : ndarray
        Coefficients of the optimal (in a least squares sense) FIR filter.

    See Also
    --------
    firwin
    firwin2
    minimum_phase
    remez

    Notes
    -----
    This implementation follows the algorithm given in [1]_.
    As noted there, least squares design has multiple advantages:

    1. Optimal in a least-squares sense.
    2. Simple, non-iterative method.
    3. The general solution can obtained by solving a linear
       system of equations.
    4. Allows the use of a frequency dependent weighting function.

    This function constructs a Type I linear phase FIR filter, which
    contains an odd number of `coeffs` satisfying for :math:`n < numtaps`:

    .. math:: coeffs(n) = coeffs(numtaps - 1 - n)

    The odd number of coefficients and filter symmetry avoid boundary
    conditions that could otherwise occur at the Nyquist and 0 frequencies
    (e.g., for Type II, III, or IV variants).

    .. versionadded:: 0.18

    References
    ----------
    .. [1] Ivan Selesnick, Linear-Phase Fir Filter Design By Least Squares.
           OpenStax CNX. Aug 9, 2005.
           https://eeweb.engineering.nyu.edu/iselesni/EL713/firls/firls.pdf

    Examples
    --------
    We want to construct a band-pass filter. Note that the behavior in the
    frequency ranges between our stop bands and pass bands is unspecified,
    and thus may overshoot depending on the parameters of our filter:

    >>> import numpy as np
    >>> from scipy import signal
    >>> import matplotlib.pyplot as plt
    >>> fig, axs = plt.subplots(2)
    >>> fs = 10.0  # Hz
    >>> desired = (0, 0, 1, 1, 0, 0)
    >>> for bi, bands in enumerate(((0, 1, 2, 3, 4, 5), (0, 1, 2, 4, 4.5, 5))):
    ...     fir_firls = signal.firls(73, bands, desired, fs=fs)
    ...     fir_remez = signal.remez(73, bands, desired[::2], fs=fs)
    ...     fir_firwin2 = signal.firwin2(73, bands, desired, fs=fs)
    ...     hs = list()
    ...     ax = axs[bi]
    ...     for fir in (fir_firls, fir_remez, fir_firwin2):
    ...         freq, response = signal.freqz(fir)
    ...         hs.append(ax.semilogy(0.5*fs*freq/np.pi, np.abs(response))[0])
    ...     for band, gains in zip(zip(bands[::2], bands[1::2]),
    ...                            zip(desired[::2], desired[1::2])):
    ...         ax.semilogy(band, np.maximum(gains, 1e-7), 'k--', linewidth=2)
    ...     if bi == 0:
    ...         ax.legend(hs, ('firls', 'remez', 'firwin2'),
    ...                   loc='lower center', frameon=False)
    ...     else:
    ...         ax.set_xlabel('Frequency (Hz)')
    ...     ax.grid(True)
    ...     ax.set(title='Band-pass %d-%d Hz' % bands[2:4], ylabel='Magnitude')
    ...
    >>> fig.tight_layout()
    >>> plt.show()

    """
    xp = array_namespace(bands, desired)
    bands = np.asarray(bands)
    desired = np.asarray(desired)

    fs = _validate_fs(fs, allow_none=True)
    fs = 2 if fs is None else fs
    nyq = 0.5 * fs

    numtaps = int(numtaps)
    if numtaps % 2 == 0 or numtaps < 1:
        raise ValueError("numtaps must be odd and >= 1")
    M = (numtaps-1) // 2

    # normalize bands 0->1 and make it 2 columns
    nyq = float(nyq)
    if nyq <= 0:
        raise ValueError(f'nyq must be positive, got {nyq} <= 0.')
    bands = np.asarray(bands).flatten() / nyq
    if len(bands) % 2 != 0:
        raise ValueError("bands must contain frequency pairs.")
    if (bands < 0).any() or (bands > 1).any():
        raise ValueError("bands must be between 0 and 1 relative to Nyquist")
    bands = bands.reshape((-1, 2))

    # check remaining params
    desired = np.asarray(desired).flatten()
    if bands.size != desired.size:
        raise ValueError(
            f"desired must have one entry per frequency, got {desired.size} "
            f"gains for {bands.size} frequencies."
        )
    desired = desired.reshape((-1, 2))
    if (np.diff(bands) <= 0).any() or (np.diff(bands[:, 0]) < 0).any():
        raise ValueError("bands must be monotonically nondecreasing and have "
                         "width > 0.")
    if (bands[:-1, 1] > bands[1:, 0]).any():
        raise ValueError("bands must not overlap.")
    if (desired < 0).any():
        raise ValueError("desired must be non-negative.")
    if weight is None:
        weight = np.ones(len(desired))
    weight = np.asarray(weight).flatten()
    if len(weight) != len(desired):
        raise ValueError("weight must be the same size as the number of "
                         f"band pairs ({len(bands)}).")
    if (weight < 0).any():
        raise ValueError("weight must be non-negative.")

    # Set up the linear matrix equation to be solved, Qa = b

    # We can express Q(k,n) = 0.5 Q1(k,n) + 0.5 Q2(k,n)
    # where Q1(k,n)=q(k-n) and Q2(k,n)=q(k+n), i.e. a Toeplitz plus Hankel.

    # We omit the factor of 0.5 above, instead adding it during coefficient
    # calculation.

    # We also omit the 1/π from both Q and b equations, as they cancel
    # during solving.

    # We have that:
    #     q(n) = 1/π ∫W(ω)cos(nω)dω (over 0->π)
    # Using our normalization ω=πf and with a constant weight W over each
    # interval f1->f2 we get:
    #     q(n) = W∫cos(πnf)df (0->1) = Wf sin(πnf)/πnf
    # integrated over each f1->f2 pair (i.e., value at f2 - value at f1).
    n = np.arange(numtaps)[:, np.newaxis, np.newaxis]
    q = np.dot(np.diff(np.sinc(bands * n) * bands, axis=2)[:, :, 0], weight)

    # Now we assemble our sum of Toeplitz and Hankel
    Q1 = toeplitz(q[:M+1])
    Q2 = hankel(q[:M+1], q[M:])
    Q = Q1 + Q2

    # Now for b(n) we have that:
    #     b(n) = 1/π ∫ W(ω)D(ω)cos(nω)dω (over 0->π)
    # Using our normalization ω=πf and with a constant weight W over each
    # interval and a linear term for D(ω) we get (over each f1->f2 interval):
    #     b(n) = W ∫ (mf+c)cos(πnf)df
    #          = f(mf+c)sin(πnf)/πnf + mf**2 cos(nπf)/(πnf)**2
    # integrated over each f1->f2 pair (i.e., value at f2 - value at f1).
    n = n[:M + 1]  # only need this many coefficients here
    # Choose m and c such that we are at the start and end weights
    m = (np.diff(desired, axis=1) / np.diff(bands, axis=1))
    c = desired[:, [0]] - bands[:, [0]] * m
    b = bands * (m*bands + c) * np.sinc(bands * n)
    # Use L'Hospital's rule here for cos(nπf)/(πnf)**2 @ n=0
    b[0] -= m * bands * bands / 2.
    b[1:] += m * np.cos(n[1:] * np.pi * bands) / (np.pi * n[1:]) ** 2
    b = np.dot(np.diff(b, axis=2)[:, :, 0], weight)

    # Now we can solve the equation
    try:  # try the fast way
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter('always')
            a = solve(Q, b, assume_a="pos", check_finite=False)
        for ww in w:
            if (ww.category == LinAlgWarning and
                    str(ww.message).startswith('Ill-conditioned matrix')):
                raise LinAlgError(str(ww.message))
    except LinAlgError:  # in case Q is rank deficient
        # This is faster than pinvh, even though we don't explicitly use
        # the symmetry here. gelsy was faster than gelsd and gelss in
        # some non-exhaustive tests.
        a = lstsq(Q, b, lapack_driver='gelsy')[0]

    # make coefficients symmetric (linear phase)
    coeffs = np.hstack((a[:0:-1], 2 * a[0], a[1:]))
    return xp.asarray(coeffs)

################################################################################
# Method Path: scipy.stats._mstats_extras._mjci_1D
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/stats/_mstats_extras.py
# Target:      dot
# Call Count:  12
# Total Exec:  2
# Func Size:   15
################################################################################

def _mjci_1D(data, p):
        data = np.sort(data.compressed())
        n = data.size
        prob = (np.array(p) * n + 0.5).astype(int)
        betacdf = beta.cdf

        mj = np.empty(len(prob), float64)
        x = np.arange(1,n+1, dtype=float64) / n
        y = x - 1./n
        for (i,m) in enumerate(prob):
            W = betacdf(x,m-1,n-m) - betacdf(y,m-1,n-m)
            C1 = np.dot(W,data)
            C2 = np.dot(W,data**2)
            mj[i] = np.sqrt(C2 - C1**2)
        return mj

################################################################################
# Method Path: scipy.signal._lti_conversion.cont2discrete
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_lti_conversion.py
# Target:      dot
# Call Count:  8
# Total Exec:  47
# Func Size:   200
################################################################################

def cont2discrete(system, dt, method="zoh", alpha=None):
    """
    Transform a continuous to a discrete state-space system.

    Parameters
    ----------
    system : a tuple describing the system or an instance of `lti`
        The following gives the number of elements in the tuple and
        the interpretation:

        * 1: (instance of `lti`)
        * 2: (num, den)
        * 3: (zeros, poles, gain)
        * 4: (A, B, C, D)

    dt : float
        The discretization time step.
    method : str, optional
        Which method to use:

        * gbt: generalized bilinear transformation
        * bilinear: Tustin's approximation ("gbt" with alpha=0.5)
        * euler: Euler (or forward differencing) method ("gbt" with alpha=0)
        * backward_diff: Backwards differencing ("gbt" with alpha=1.0)
        * zoh: zero-order hold (default)
        * foh: first-order hold (*versionadded: 1.3.0*)
        * impulse: equivalent impulse response (*versionadded: 1.3.0*)

    alpha : float within [0, 1], optional
        The generalized bilinear transformation weighting parameter, which
        should only be specified with method="gbt", and is ignored otherwise

    Returns
    -------
    sysd : tuple containing the discrete system
        Based on the input type, the output will be of the form

        * (num, den, dt)   for transfer function input
        * (zeros, poles, gain, dt)   for zeros-poles-gain input
        * (A, B, C, D, dt) for state-space system input

    Notes
    -----
    By default, the routine uses a Zero-Order Hold (zoh) method to perform
    the transformation. Alternatively, a generalized bilinear transformation
    may be used, which includes the common Tustin's bilinear approximation,
    an Euler's method technique, or a backwards differencing technique.

    The Zero-Order Hold (zoh) method is based on [1]_, the generalized bilinear
    approximation is based on [2]_ and [3]_, the First-Order Hold (foh) method
    is based on [4]_.

    References
    ----------
    .. [1] https://en.wikipedia.org/wiki/Discretization#Discretization_of_linear_state_space_models

    .. [2] http://techteach.no/publications/discretetime_signals_systems/discrete.pdf

    .. [3] G. Zhang, X. Chen, and T. Chen, Digital redesign via the generalized
        bilinear transformation, Int. J. Control, vol. 82, no. 4, pp. 741-754,
        2009.
        (https://www.mypolyuweb.hk/~magzhang/Research/ZCC09_IJC.pdf)

    .. [4] G. F. Franklin, J. D. Powell, and M. L. Workman, Digital control
        of dynamic systems, 3rd ed. Menlo Park, Calif: Addison-Wesley,
        pp. 204-206, 1998.

    Examples
    --------
    We can transform a continuous state-space system to a discrete one:

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.signal import cont2discrete, lti, dlti, dstep

    Define a continuous state-space system.

    >>> A = np.array([[0, 1],[-10., -3]])
    >>> B = np.array([[0],[10.]])
    >>> C = np.array([[1., 0]])
    >>> D = np.array([[0.]])
    >>> l_system = lti(A, B, C, D)
    >>> t, x = l_system.step(T=np.linspace(0, 5, 100))
    >>> fig, ax = plt.subplots()
    >>> ax.plot(t, x, label='Continuous', linewidth=3)

    Transform it to a discrete state-space system using several methods.

    >>> dt = 0.1
    >>> for method in ['zoh', 'bilinear', 'euler', 'backward_diff', 'foh', 'impulse']:
    ...    d_system = cont2discrete((A, B, C, D), dt, method=method)
    ...    s, x_d = dstep(d_system)
    ...    ax.step(s, np.squeeze(x_d), label=method, where='post')
    >>> ax.axis([t[0], t[-1], x[0], 1.4])
    >>> ax.legend(loc='best')
    >>> fig.tight_layout()
    >>> plt.show()

    """
    if hasattr(system, 'to_discrete') and callable(system.to_discrete):
        return system.to_discrete(dt=dt, method=method, alpha=alpha)

    if len(system) == 2:
        sysd = cont2discrete(tf2ss(system[0], system[1]), dt, method=method,
                             alpha=alpha)
        return ss2tf(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(system) == 3:
        sysd = cont2discrete(zpk2ss(system[0], system[1], system[2]), dt,
                             method=method, alpha=alpha)
        return ss2zpk(sysd[0], sysd[1], sysd[2], sysd[3]) + (dt,)
    elif len(system) == 4:
        a, b, c, d = system
    else:
        raise ValueError("First argument must either be a tuple of 2 (tf), "
                         "3 (zpk), or 4 (ss) arrays.")

    if method == 'gbt':
        if alpha is None:
            raise ValueError("Alpha parameter must be specified for the "
                             "generalized bilinear transform (gbt) method")
        elif alpha < 0 or alpha > 1:
            raise ValueError("Alpha parameter must be within the interval "
                             "[0,1] for the gbt method")

    if method == 'gbt':
        # This parameter is used repeatedly - compute once here
        ima = np.eye(a.shape[0]) - alpha*dt*a
        ad = linalg.solve(ima, np.eye(a.shape[0]) + (1.0-alpha)*dt*a)
        bd = linalg.solve(ima, dt*b)

        # Similarly solve for the output equation matrices
        cd = linalg.solve(ima.transpose(), c.transpose())
        cd = cd.transpose()
        dd = d + alpha*np.dot(c, bd)

    elif method == 'bilinear' or method == 'tustin':
        return cont2discrete(system, dt, method="gbt", alpha=0.5)

    elif method == 'euler' or method == 'forward_diff':
        return cont2discrete(system, dt, method="gbt", alpha=0.0)

    elif method == 'backward_diff':
        return cont2discrete(system, dt, method="gbt", alpha=1.0)

    elif method == 'zoh':
        # Build an exponential matrix
        em_upper = np.hstack((a, b))

        # Need to stack zeros under the a and b matrices
        em_lower = np.hstack((np.zeros((b.shape[1], a.shape[0])),
                              np.zeros((b.shape[1], b.shape[1]))))

        em = np.vstack((em_upper, em_lower))
        ms = linalg.expm(dt * em)

        # Dispose of the lower rows
        ms = ms[:a.shape[0], :]

        ad = ms[:, 0:a.shape[1]]
        bd = ms[:, a.shape[1]:]

        cd = c
        dd = d

    elif method == 'foh':
        # Size parameters for convenience
        n = a.shape[0]
        m = b.shape[1]

        # Build an exponential matrix similar to 'zoh' method
        em_upper = linalg.block_diag(np.block([a, b]) * dt, np.eye(m))
        em_lower = zeros((m, n + 2 * m))
        em = np.block([[em_upper], [em_lower]])

        ms = linalg.expm(em)

        # Get the three blocks from upper rows
        ms11 = ms[:n, 0:n]
        ms12 = ms[:n, n:n + m]
        ms13 = ms[:n, n + m:]

        ad = ms11
        bd = ms12 - ms13 + ms11 @ ms13
        cd = c
        dd = d + c @ ms13

    elif method == 'impulse':
        if not np.allclose(d, 0):
            raise ValueError("Impulse method is only applicable "
                             "to strictly proper systems")

        ad = linalg.expm(a * dt)
        bd = ad @ b * dt
        cd = c
        dd = c @ b * dt

    else:
        raise ValueError(f"Unknown transformation method '{method}'")

    return ad, bd, cd, dd, dt

################################################################################
# Method Path: scipy.integrate._bvp.solve_bvp
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_bvp.py
# Target:      dot
# Call Count:  8
# Total Exec:  33
# Func Size:   454
################################################################################

def solve_bvp(fun, bc, x, y, p=None, S=None, fun_jac=None, bc_jac=None,
              tol=1e-3, max_nodes=1000, verbose=0, bc_tol=None):
    """Solve a boundary value problem for a system of ODEs.

    This function numerically solves a first order system of ODEs subject to
    two-point boundary conditions::

        dy / dx = f(x, y, p) + S * y / (x - a), a <= x <= b
        bc(y(a), y(b), p) = 0

    Here x is a 1-D independent variable, y(x) is an n-D
    vector-valued function and p is a k-D vector of unknown
    parameters which is to be found along with y(x). For the problem to be
    determined, there must be n + k boundary conditions, i.e., bc must be an
    (n + k)-D function.

    The last singular term on the right-hand side of the system is optional.
    It is defined by an n-by-n matrix S, such that the solution must satisfy
    S y(a) = 0. This condition will be forced during iterations, so it must not
    contradict boundary conditions. See [2]_ for the explanation how this term
    is handled when solving BVPs numerically.

    Problems in a complex domain can be solved as well. In this case, y and p
    are considered to be complex, and f and bc are assumed to be complex-valued
    functions, but x stays real. Note that f and bc must be complex
    differentiable (satisfy Cauchy-Riemann equations [4]_), otherwise you
    should rewrite your problem for real and imaginary parts separately. To
    solve a problem in a complex domain, pass an initial guess for y with a
    complex data type (see below).

    Parameters
    ----------
    fun : callable
        Right-hand side of the system. The calling signature is ``fun(x, y)``,
        or ``fun(x, y, p)`` if parameters are present. All arguments are
        ndarray: ``x`` with shape (m,), ``y`` with shape (n, m), meaning that
        ``y[:, i]`` corresponds to ``x[i]``, and ``p`` with shape (k,). The
        return value must be an array with shape (n, m) and with the same
        layout as ``y``.
    bc : callable
        Function evaluating residuals of the boundary conditions. The calling
        signature is ``bc(ya, yb)``, or ``bc(ya, yb, p)`` if parameters are
        present. All arguments are ndarray: ``ya`` and ``yb`` with shape (n,),
        and ``p`` with shape (k,). The return value must be an array with
        shape (n + k,).
        The order of the returned residuals does not matter, as `solve_bvp`
        attempts to drive all residuals to zero.
    x : array_like, shape (m,)
        Initial mesh. Must be a strictly increasing sequence of real numbers
        with ``x[0]=a`` and ``x[-1]=b``.
    y : array_like, shape (n, m)
        Initial guess for the function values at the mesh nodes, ith column
        corresponds to ``x[i]``. For problems in a complex domain pass `y`
        with a complex data type (even if the initial guess is purely real).
    p : array_like with shape (k,) or None, optional
        Initial guess for the unknown parameters. If None (default), it is
        assumed that the problem doesn't depend on any parameters.
    S : array_like with shape (n, n) or None
        Matrix defining the singular term. If None (default), the problem is
        solved without the singular term.
    fun_jac : callable or None, optional
        Function computing derivatives of f with respect to y and p. The
        calling signature is ``fun_jac(x, y)``, or ``fun_jac(x, y, p)`` if
        parameters are present. The return must contain 1 or 2 elements in the
        following order:

        * df_dy : array_like with shape ``(n, n, m)``, where an element
          ``(i, j, q)`` equals to ``d f_i(x_q, y_q, p) / d (y_q)_j``.
        * df_dp : array_like with shape ``(n, k, m)``, where an element
          ``(i, j, q)`` equals to ``d f_i(x_q, y_q, p) / d p_j``.

        Here q numbers nodes at which x and y are defined, whereas i and j
        number vector components. If the problem is solved without unknown
        parameters, df_dp should not be returned.

        If `fun_jac` is None (default), the derivatives will be estimated
        by the forward finite differences.
    bc_jac : callable or None, optional
        Function computing derivatives of bc with respect to ya, yb, and p.
        The calling signature is ``bc_jac(ya, yb)``, or ``bc_jac(ya, yb, p)``
        if parameters are present. The return must contain 2 or 3 elements in
        the following order:

        * ``dbc_dya`` : array_like with shape ``(n, n)``, where an element ``(i, j)``
          equals to ``d bc_i(ya, yb, p) / d ya_j``.
        * ``dbc_dyb`` : array_like with shape ``(n, n)``, where an element ``(i, j)``
          equals to ``d bc_i(ya, yb, p) / d yb_j``.
        * ``dbc_dp`` : array_like with shape ``(n, k)``, where an element ``(i, j)``
          equals to ``d bc_i(ya, yb, p) / d p_j``.

        If the problem is solved without unknown parameters, dbc_dp should not
        be returned.

        If `bc_jac` is None (default), the derivatives will be estimated by
        the forward finite differences.
    tol : float, optional
        Desired tolerance of the solution. If we define ``r = y' - f(x, y)``,
        where y is the found solution, then the solver tries to achieve on each
        mesh interval ``norm(r / (1 + abs(f)) < tol``, where ``norm`` is
        estimated in a root mean squared sense (using a numerical quadrature
        formula). Default is 1e-3.
    max_nodes : int, optional
        Maximum allowed number of the mesh nodes. If exceeded, the algorithm
        terminates. Default is 1000.
    verbose : {0, 1, 2}, optional
        Level of algorithm's verbosity:

        * 0 (default) : work silently.
        * 1 : display a termination report.
        * 2 : display progress during iterations.
    bc_tol : float, optional
        Desired absolute tolerance for the boundary condition residuals: `bc`
        value should satisfy ``abs(bc) < bc_tol`` component-wise.
        Equals to `tol` by default. Up to 10 iterations are allowed to achieve this
        tolerance.

    Returns
    -------
    result : BVPResult
        Bunch object with the following fields defined:

        sol : PPoly
            Found solution for y as `scipy.interpolate.PPoly` instance, a C1
            continuous cubic spline.
        p : ndarray or None, shape (k,)
            Found parameters. None, if the parameters were not present in the
            problem.
        x : ndarray, shape (m,)
            Nodes of the final mesh.
        y : ndarray, shape (n, m)
            Solution values at the mesh nodes.
        yp : ndarray, shape (n, m)
            Solution derivatives at the mesh nodes.
        rms_residuals : ndarray, shape (m - 1,)
            RMS values of the relative residuals over each mesh interval (see the
            description of `tol` parameter).
        niter : int
            Number of completed iterations.
        status : int
            Reason for algorithm termination:

            * 0: The algorithm converged to the desired accuracy.
            * 1: The maximum number of mesh nodes is exceeded.
            * 2: A singular Jacobian encountered when solving the collocation system.

        message : str
            Verbal description of the termination reason.
        success : bool
            True if the algorithm converged to the desired accuracy (``status=0``).

    Notes
    -----
    This function implements a 4th order collocation algorithm with the
    control of residuals similar to [1]_. A collocation system is solved
    by a damped Newton method with an affine-invariant criterion function as
    described in [3]_.

    Note that in [1]_  integral residuals are defined without normalization
    by interval lengths. So, their definition is different by a multiplier of
    h**0.5 (h is an interval length) from the definition used here.

    .. versionadded:: 0.18.0

    References
    ----------
    .. [1] J. Kierzenka, L. F. Shampine, "A BVP Solver Based on Residual
           Control and the Maltab PSE", ACM Trans. Math. Softw., Vol. 27,
           Number 3, pp. 299-316, 2001.
    .. [2] L.F. Shampine, P. H. Muir and H. Xu, "A User-Friendly Fortran BVP
           Solver", J. Numer. Anal., Ind. Appl. Math. (JNAIAM), Vol. 1, 
           Number 2, pp. 201-217, 2006.
    .. [3] U. Ascher, R. Mattheij and R. Russell "Numerical Solution of
           Boundary Value Problems for Ordinary Differential Equations",
           Philidelphia, PA: Society for Industrial and Applied Mathematics,
           1995.
           :doi:`10.1137/1.9781611971231`
    .. [4] `Cauchy-Riemann equations
            <https://en.wikipedia.org/wiki/Cauchy-Riemann_equations>`_ on
            Wikipedia.

    Examples
    --------
    In the first example, we solve Bratu's problem::

        y'' + k * exp(y) = 0
        y(0) = y(1) = 0

    for k = 1.

    We rewrite the equation as a first-order system and implement its
    right-hand side evaluation::

        y1' = y2
        y2' = -exp(y1)

    >>> import numpy as np
    >>> def fun(x, y):
    ...     return np.vstack((y[1], -np.exp(y[0])))

    Implement evaluation of the boundary condition residuals:

    >>> def bc(ya, yb):
    ...     return np.array([ya[0], yb[0]])

    Define the initial mesh with 5 nodes:

    >>> x = np.linspace(0, 1, 5)

    This problem is known to have two solutions. To obtain both of them, we
    use two different initial guesses for y. We denote them by subscripts
    a and b.

    >>> y_a = np.zeros((2, x.size))
    >>> y_b = np.zeros((2, x.size))
    >>> y_b[0] = 3

    Now we are ready to run the solver.

    >>> from scipy.integrate import solve_bvp
    >>> res_a = solve_bvp(fun, bc, x, y_a)
    >>> res_b = solve_bvp(fun, bc, x, y_b)

    Let's plot the two found solutions. We take an advantage of having the
    solution in a spline form to produce a smooth plot.

    >>> x_plot = np.linspace(0, 1, 100)
    >>> y_plot_a = res_a.sol(x_plot)[0]
    >>> y_plot_b = res_b.sol(x_plot)[0]
    >>> import matplotlib.pyplot as plt
    >>> plt.plot(x_plot, y_plot_a, label='y_a')
    >>> plt.plot(x_plot, y_plot_b, label='y_b')
    >>> plt.legend()
    >>> plt.xlabel("x")
    >>> plt.ylabel("y")
    >>> plt.show()

    We see that the two solutions have similar shape, but differ in scale
    significantly.

    In the second example, we solve a simple Sturm-Liouville problem::

        y'' + k**2 * y = 0
        y(0) = y(1) = 0

    It is known that a non-trivial solution y = A * sin(k * x) is possible for
    k = pi * n, where n is an integer. To establish the normalization constant
    A = 1 we add a boundary condition::

        y'(0) = k

    Again, we rewrite our equation as a first-order system and implement its
    right-hand side evaluation::

        y1' = y2
        y2' = -k**2 * y1

    >>> def fun(x, y, p):
    ...     k = p[0]
    ...     return np.vstack((y[1], -k**2 * y[0]))

    Note that parameters p are passed as a vector (with one element in our
    case).

    Implement the boundary conditions:

    >>> def bc(ya, yb, p):
    ...     k = p[0]
    ...     return np.array([ya[0], yb[0], ya[1] - k])

    Set up the initial mesh and guess for y. We aim to find the solution for
    k = 2 * pi, to achieve that we set values of y to approximately follow
    sin(2 * pi * x):

    >>> x = np.linspace(0, 1, 5)
    >>> y = np.zeros((2, x.size))
    >>> y[0, 1] = 1
    >>> y[0, 3] = -1

    Run the solver with 6 as an initial guess for k.

    >>> sol = solve_bvp(fun, bc, x, y, p=[6])

    We see that the found k is approximately correct:

    >>> sol.p[0]
    6.28329460046

    And, finally, plot the solution to see the anticipated sinusoid:

    >>> x_plot = np.linspace(0, 1, 100)
    >>> y_plot = sol.sol(x_plot)[0]
    >>> plt.plot(x_plot, y_plot)
    >>> plt.xlabel("x")
    >>> plt.ylabel("y")
    >>> plt.show()
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("`x` must be 1 dimensional.")
    h = np.diff(x)
    if np.any(h <= 0):
        raise ValueError("`x` must be strictly increasing.")
    a = x[0]

    y = np.asarray(y)
    if np.issubdtype(y.dtype, np.complexfloating):
        dtype = complex
    else:
        dtype = float
    y = y.astype(dtype, copy=False)

    if y.ndim != 2:
        raise ValueError("`y` must be 2 dimensional.")
    if y.shape[1] != x.shape[0]:
        raise ValueError(f"`y` is expected to have {x.shape[0]} columns, but actually "
                         f"has {y.shape[1]}.")

    if p is None:
        p = np.array([])
    else:
        p = np.asarray(p, dtype=dtype)
    if p.ndim != 1:
        raise ValueError("`p` must be 1 dimensional.")

    if tol < 100 * EPS:
        warn(f"`tol` is too low, setting to {100 * EPS:.2e}", stacklevel=2)
        tol = 100 * EPS

    if verbose not in [0, 1, 2]:
        raise ValueError("`verbose` must be in [0, 1, 2].")

    n = y.shape[0]
    k = p.shape[0]

    if S is not None:
        S = np.asarray(S, dtype=dtype)
        if S.shape != (n, n):
            raise ValueError(f"`S` is expected to have shape {(n, n)}, "
                             f"but actually has {S.shape}")

        # Compute I - S^+ S to impose necessary boundary conditions.
        B = np.identity(n) - np.dot(pinv(S), S)

        y[:, 0] = np.dot(B, y[:, 0])

        # Compute (I - S)^+ to correct derivatives at x=a.
        D = pinv(np.identity(n) - S)
    else:
        B = None
        D = None

    if bc_tol is None:
        bc_tol = tol

    # Maximum number of iterations
    max_iteration = 10

    fun_wrapped, bc_wrapped, fun_jac_wrapped, bc_jac_wrapped = wrap_functions(
        fun, bc, fun_jac, bc_jac, k, a, S, D, dtype)

    f = fun_wrapped(x, y, p)
    if f.shape != y.shape:
        raise ValueError(f"`fun` return is expected to have shape {y.shape}, "
                         f"but actually has {f.shape}.")

    bc_res = bc_wrapped(y[:, 0], y[:, -1], p)
    if bc_res.shape != (n + k,):
        raise ValueError(f"`bc` return is expected to have shape {(n + k,)}, "
                         f"but actually has {bc_res.shape}.")

    status = 0
    iteration = 0
    if verbose == 2:
        print_iteration_header()

    while True:
        m = x.shape[0]

        col_fun, jac_sys = prepare_sys(n, m, k, fun_wrapped, bc_wrapped,
                                       fun_jac_wrapped, bc_jac_wrapped, x, h)
        y, p, singular = solve_newton(n, m, h, col_fun, bc_wrapped, jac_sys,
                                      y, p, B, tol, bc_tol)
        iteration += 1

        col_res, y_middle, f, f_middle = collocation_fun(fun_wrapped, y,
                                                         p, x, h)
        bc_res = bc_wrapped(y[:, 0], y[:, -1], p)
        max_bc_res = np.max(abs(bc_res))

        # This relation is not trivial, but can be verified.
        r_middle = 1.5 * col_res / h
        sol = create_spline(y, f, x, h)
        rms_res = estimate_rms_residuals(fun_wrapped, sol, x, h, p,
                                         r_middle, f_middle)
        max_rms_res = np.max(rms_res)

        if singular:
            status = 2
            break

        insert_1, = np.nonzero((rms_res > tol) & (rms_res < 100 * tol))
        insert_2, = np.nonzero(rms_res >= 100 * tol)
        nodes_added = insert_1.shape[0] + 2 * insert_2.shape[0]

        if m + nodes_added > max_nodes:
            status = 1
            if verbose == 2:
                nodes_added = f"({nodes_added})"
                print_iteration_progress(iteration, max_rms_res, max_bc_res,
                                         m, nodes_added)
            break

        if verbose == 2:
            print_iteration_progress(iteration, max_rms_res, max_bc_res, m,
                                     nodes_added)

        if nodes_added > 0:
            x = modify_mesh(x, insert_1, insert_2)
            h = np.diff(x)
            y = sol(x)
        elif max_bc_res <= bc_tol:
            status = 0
            break
        elif iteration >= max_iteration:
            status = 3
            break

    if verbose > 0:
        if status == 0:
            print(f"Solved in {iteration} iterations, number of nodes {x.shape[0]}. \n"
                  f"Maximum relative residual: {max_rms_res:.2e} \n"
                  f"Maximum boundary residual: {max_bc_res:.2e}")
        elif status == 1:
            print(f"Number of nodes is exceeded after iteration {iteration}. \n"
                  f"Maximum relative residual: {max_rms_res:.2e} \n"
                  f"Maximum boundary residual: {max_bc_res:.2e}")
        elif status == 2:
            print("Singular Jacobian encountered when solving the collocation "
                  f"system on iteration {iteration}. \n"
                  f"Maximum relative residual: {max_rms_res:.2e} \n"
                  f"Maximum boundary residual: {max_bc_res:.2e}")
        elif status == 3:
            print("The solver was unable to satisfy boundary conditions "
                  f"tolerance on iteration {iteration}. \n"
                  f"Maximum relative residual: {max_rms_res:.2e} \n"
                  f"Maximum boundary residual: {max_bc_res:.2e}")

    if p.size == 0:
        p = None

    return BVPResult(sol=sol, p=p, x=x, y=y, yp=f, rms_residuals=rms_res,
                     niter=iteration, status=status,
                     message=TERMINATION_MESSAGES[status], success=status == 0)

################################################################################
# Method Path: scipy.integrate._cubature.cubature
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_cubature.py
# Target:      inv
# Call Count:  5
# Total Exec:  246
# Func Size:   417
################################################################################

# Error reading or parsing /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_cubature.py: invalid syntax (<unknown>, line 27)

################################################################################
# Method Path: scipy.signal._ltisys.lsim
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/signal/_ltisys.py
# Target:      dot
# Call Count:  5
# Total Exec:  89
# Func Size:   213
################################################################################

def lsim(system, U, T, X0=None, interp=True):
    """
    Simulate output of a continuous-time linear system.

    Parameters
    ----------
    system : an instance of the LTI class or a tuple describing the system
        The following gives the number of elements in the tuple and
        the interpretation:

        * 1: (instance of `lti`)
        * 2: (num, den)
        * 3: (zeros, poles, gain)
        * 4: (A, B, C, D)

    U : array_like
        An input array describing the input at each time `T`
        (interpolation is assumed between given times).  If there are
        multiple inputs, then each column of the rank-2 array
        represents an input.  If U = 0 or None, a zero input is used.
    T : array_like
        The time steps at which the input is defined and at which the
        output is desired.  Must be nonnegative, increasing, and equally spaced.
    X0 : array_like, optional
        The initial conditions on the state vector (zero by default).
    interp : bool, optional
        Whether to use linear (True, the default) or zero-order-hold (False)
        interpolation for the input array.

    Returns
    -------
    T : 1D ndarray
        Time values for the output.
    yout : 1D ndarray
        System response.
    xout : ndarray
        Time evolution of the state vector.

    Notes
    -----
    If (num, den) is passed in for ``system``, coefficients for both the
    numerator and denominator should be specified in descending exponent
    order (e.g. ``s^2 + 3s + 5`` would be represented as ``[1, 3, 5]``).

    Examples
    --------
    We'll use `lsim` to simulate an analog Bessel filter applied to
    a signal.

    >>> import numpy as np
    >>> from scipy.signal import bessel, lsim
    >>> import matplotlib.pyplot as plt

    Create a low-pass Bessel filter with a cutoff of 12 Hz.

    >>> b, a = bessel(N=5, Wn=2*np.pi*12, btype='lowpass', analog=True)

    Generate data to which the filter is applied.

    >>> t = np.linspace(0, 1.25, 500, endpoint=False)

    The input signal is the sum of three sinusoidal curves, with
    frequencies 4 Hz, 40 Hz, and 80 Hz.  The filter should mostly
    eliminate the 40 Hz and 80 Hz components, leaving just the 4 Hz signal.

    >>> u = (np.cos(2*np.pi*4*t) + 0.6*np.sin(2*np.pi*40*t) +
    ...      0.5*np.cos(2*np.pi*80*t))

    Simulate the filter with `lsim`.

    >>> tout, yout, xout = lsim((b, a), U=u, T=t)

    Plot the result.

    >>> plt.plot(t, u, 'r', alpha=0.5, linewidth=1, label='input')
    >>> plt.plot(tout, yout, 'k', linewidth=1.5, label='output')
    >>> plt.legend(loc='best', shadow=True, framealpha=1)
    >>> plt.grid(alpha=0.3)
    >>> plt.xlabel('t')
    >>> plt.show()

    In a second example, we simulate a double integrator ``y'' = u``, with
    a constant input ``u = 1``.  We'll use the state space representation
    of the integrator.

    >>> from scipy.signal import lti
    >>> A = np.array([[0.0, 1.0], [0.0, 0.0]])
    >>> B = np.array([[0.0], [1.0]])
    >>> C = np.array([[1.0, 0.0]])
    >>> D = 0.0
    >>> system = lti(A, B, C, D)

    `t` and `u` define the time and input signal for the system to
    be simulated.

    >>> t = np.linspace(0, 5, num=50)
    >>> u = np.ones_like(t)

    Compute the simulation, and then plot `y`.  As expected, the plot shows
    the curve ``y = 0.5*t**2``.

    >>> tout, y, x = lsim(system, u, t)
    >>> plt.plot(t, y)
    >>> plt.grid(alpha=0.3)
    >>> plt.xlabel('t')
    >>> plt.show()

    """
    if isinstance(system, lti):
        sys = system._as_ss()
    elif isinstance(system, dlti):
        raise AttributeError('lsim can only be used with continuous-time '
                             'systems.')
    else:
        sys = lti(*system)._as_ss()
    T = atleast_1d(T)
    if len(T.shape) != 1:
        raise ValueError("T must be a rank-1 array.")

    A, B, C, D = map(np.asarray, (sys.A, sys.B, sys.C, sys.D))
    n_states = A.shape[0]
    n_inputs = B.shape[1]

    n_steps = T.size
    if X0 is None:
        X0 = zeros(n_states, sys.A.dtype)
    xout = np.empty((n_steps, n_states), sys.A.dtype)

    if T[0] == 0:
        xout[0] = X0
    elif T[0] > 0:
        # step forward to initial time, with zero input
        xout[0] = dot(X0, linalg.expm(transpose(A) * T[0]))
    else:
        raise ValueError("Initial time must be nonnegative")

    no_input = (U is None or
                (isinstance(U, int | float) and U == 0.) or
                not np.any(U))

    if n_steps == 1:
        yout = squeeze(xout @ C.T)
        if not no_input:
            yout += squeeze(U @ D.T)
        return T, yout, squeeze(xout)

    dt = T[1] - T[0]
    if not np.allclose(np.diff(T), dt):
        raise ValueError("Time steps are not equally spaced.")

    if no_input:
        # Zero input: just use matrix exponential
        # take transpose because state is a row vector
        expAT_dt = linalg.expm(A.T * dt)
        for i in range(1, n_steps):
            xout[i] = xout[i-1] @ expAT_dt
        yout = squeeze(xout @ C.T)
        return T, yout, squeeze(xout)

    # Nonzero input
    U = atleast_1d(U)
    if U.ndim == 1:
        U = U[:, np.newaxis]

    if U.shape[0] != n_steps:
        raise ValueError("U must have the same number of rows "
                         "as elements in T.")

    if U.shape[1] != n_inputs:
        raise ValueError("System does not define that many inputs.")

    if not interp:
        # Zero-order hold
        # Algorithm: to integrate from time 0 to time dt, we solve
        #   xdot = A x + B u,  x(0) = x0
        #   udot = 0,          u(0) = u0.
        #
        # Solution is
        #   [ x(dt) ]       [ A*dt   B*dt ] [ x0 ]
        #   [ u(dt) ] = exp [  0     0    ] [ u0 ]
        M = np.vstack([np.hstack([A * dt, B * dt]),
                       np.zeros((n_inputs, n_states + n_inputs))])
        # transpose everything because the state and input are row vectors
        expMT = linalg.expm(M.T)
        Ad = expMT[:n_states, :n_states]
        Bd = expMT[n_states:, :n_states]
        for i in range(1, n_steps):
            xout[i] = xout[i-1] @ Ad + U[i-1] @ Bd
    else:
        # Linear interpolation between steps
        # Algorithm: to integrate from time 0 to time dt, with linear
        # interpolation between inputs u(0) = u0 and u(dt) = u1, we solve
        #   xdot = A x + B u,        x(0) = x0
        #   udot = (u1 - u0) / dt,   u(0) = u0.
        #
        # Solution is
        #   [ x(dt) ]       [ A*dt  B*dt  0 ] [  x0   ]
        #   [ u(dt) ] = exp [  0     0    I ] [  u0   ]
        #   [u1 - u0]       [  0     0    0 ] [u1 - u0]
        M = np.vstack([np.hstack([A * dt, B * dt,
                                  np.zeros((n_states, n_inputs))]),
                       np.hstack([np.zeros((n_inputs, n_states + n_inputs)),
                                  np.identity(n_inputs)]),
                       np.zeros((n_inputs, n_states + 2 * n_inputs))])
        expMT = linalg.expm(M.T)
        Ad = expMT[:n_states, :n_states]
        Bd1 = expMT[n_states+n_inputs:, :n_states]
        Bd0 = expMT[n_states:n_states + n_inputs, :n_states] - Bd1
        for i in range(1, n_steps):
            xout[i] = xout[i-1] @ Ad + U[i-1] @ Bd0 + U[i] @ Bd1

    yout = squeeze(xout @ C.T) + squeeze(U @ D.T)
    return T, yout, squeeze(xout)

################################################################################
# Method Path: scipy.optimize._optimize.check_grad
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_optimize.py
# Target:      dot
# Call Count:  4
# Total Exec:  9
# Func Size:   91
################################################################################

def check_grad(func, grad, x0, *args, epsilon=_epsilon, direction='all', rng=None):
    r"""Check the correctness of a gradient function by comparing it against a
    (forward) finite-difference approximation of the gradient.

    Parameters
    ----------
    func : callable ``func(x0, *args)``
        Function whose derivative is to be checked.
    grad : callable ``grad(x0, *args)``
        Jacobian of `func`.
    x0 : ndarray
        Points to check `grad` against forward difference approximation of grad
        using `func`.
    *args : optional
        Extra arguments passed to `func` and `grad`.
    epsilon : float, optional
        Step size used for the finite difference approximation. It defaults to
        ``sqrt(np.finfo(float).eps)``, which is approximately 1.49e-08.
    direction : str, optional
        If set to ``'random'``, then gradients along a random vector
        are used to check `grad` against forward difference approximation
        using `func`. By default it is ``'all'``, in which case, all
        the one hot direction vectors are considered to check `grad`.
        If `func` is a vector valued function then only ``'all'`` can be used.
    rng : `numpy.random.Generator`, optional
        Pseudorandom number generator state. When `rng` is None, a new
        `numpy.random.Generator` is created using entropy from the
        operating system. Types other than `numpy.random.Generator` are
        passed to `numpy.random.default_rng` to instantiate a ``Generator``.

        The random numbers generated affect the random vector along which gradients
        are computed to check ``grad``. Note that `rng` is only used when `direction`
        argument is set to `'random'`.

    Returns
    -------
    err : float
        The square root of the sum of squares (i.e., the 2-norm) of the
        difference between ``grad(x0, *args)`` and the finite difference
        approximation of `grad` using func at the points `x0`.

    See Also
    --------
    approx_fprime

    Examples
    --------
    >>> import numpy as np
    >>> def func(x):
    ...     return x[0]**2 - 0.5 * x[1]**3
    >>> def grad(x):
    ...     return [2 * x[0], -1.5 * x[1]**2]
    >>> from scipy.optimize import check_grad
    >>> check_grad(func, grad, [1.5, -1.5])
    2.9802322387695312e-08  # may vary
    >>> rng = np.random.default_rng()
    >>> check_grad(func, grad, [1.5, -1.5],
    ...             direction='random', seed=rng)
    2.9802322387695312e-08

    """
    step = epsilon
    x0 = np.asarray(x0)

    def g(w, func, x0, v, *args):
        return func(x0 + w*v, *args)

    if direction == 'random':
        _grad = np.asanyarray(grad(x0, *args))
        if _grad.ndim > 1:
            raise ValueError("'random' can only be used with scalar valued"
                             " func")
        rng_gen = check_random_state(rng)
        v = rng_gen.standard_normal(size=(x0.shape))
        _args = (func, x0, v) + args
        _func = g
        vars = np.zeros((1,))
        analytical_grad = np.dot(_grad, v)
    elif direction == 'all':
        _args = args
        _func = func
        vars = x0
        analytical_grad = grad(x0, *args)
    else:
        raise ValueError(f"{direction} is not a valid string for "
                         "``direction`` argument")

    return np.sqrt(np.sum(np.abs(
        (analytical_grad - approx_fprime(vars, _func, step, *_args))**2
    )))

################################################################################
# Method Path: scipy.sparse.linalg._matfuncs.ProductOperator._rmatvec
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/sparse/linalg/_matfuncs.py
# Target:      dot
# Call Count:  4
# Total Exec:  2
# Func Size:   5
################################################################################

def _rmatvec(self, x):
        x = x.ravel()
        for A in self._operator_sequence:
            x = A.T.dot(x)
        return x

################################################################################
# Method Path: scipy.optimize._direct_py._func_wrap
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/optimize/_direct_py.py
# Target:      inv
# Call Count:  2
# Total Exec:  51785
# Func Size:   8
################################################################################

def _func_wrap(x, args=None):
        x = np.asarray(x)
        if args is None:
            f = func(x)
        else:
            f = func(x, *args)
        # always return a float
        return np.asarray(f).item()

################################################################################
# Method Path: scipy.integrate._bvp.fun_jac_wrapped
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_bvp.py
# Target:      dot
# Call Count:  2
# Total Exec:  4
# Func Size:   9
################################################################################

def fun_jac_wrapped(x, y, p):
                df_dy, df_dp = fun_jac_p(x, y, p)
                if x[0] == a:
                    df_dy[:, :, 0] = np.dot(D, df_dy[:, :, 0])
                    df_dy[:, :, 1:] += Sr / (x[1:] - a)
                else:
                    df_dy += Sr / (x - a)

                return df_dy, df_dp

################################################################################
# Method Path: scipy.integrate._ivp.rk.DOP853._estimate_error
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/scipy/integrate/_ivp/rk.py
# Target:      dot
# Call Count:  2
# Total Exec:  1
# Func Size:   8
################################################################################

def _estimate_error(self, K, h):  # Left for testing purposes.
        err5 = np.dot(K.T, self.E5)
        err3 = np.dot(K.T, self.E3)
        denom = np.hypot(np.abs(err5), 0.1 * np.abs(err3))
        correction_factor = np.ones_like(err5)
        mask = denom > 0
        correction_factor[mask] = np.abs(err5[mask]) / denom[mask]
        return h * err5 * correction_factor

################################################################################
# Method Path: scipy._external.pyprima.common._project._project
# File Path:   /Users/muyeedahmed/Desktop/Gitcode/ForkedLibrary/scipy/subprojects/pyprima/pyprima/pyprima/src/pyprima/common/_project.py
# Target:      dot
# Call Count:  1
# Total Exec:  33
# Func Size:   157
################################################################################

def _project(x0, lb, ub, constraints):
    """Projection of the initial guess onto the feasible set.

    Parameters
    ----------
    x0: ndarray, shape (n,)
        The same as in prepdfo.
    lb: ndarray, shape (n,)
        The same as in prepdfo.
    ub: ndarray, shape (n,)
        The same as in prepdfo.
    constraints: dict
        The general constraints of the problem, defined as a dictionary with
        fields:
            linear: LinearConstraint
                The linear constraints of the problem.
            nonlinear: dict
                The nonlinear constraints of the problem. When ``_project`` is called, the nonlinear constraints are
                None.

    Returns
    -------
    result: OptimizeResult
        The result of the projection.

    Authors
    -------
    Tom M. RAGONNEAU (ragonneau.github.io)
    and Zaikun ZHANG (www.zhangzk.net)

    Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    """
    invoker = 'prima'

    # Validate x0.
    if isinstance(x0, scalar_types):
        x0_c = [x0]
    elif hasattr(x0, '__len__'):
        x0_c = x0
    else:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should be a vector.'.format(invoker))
    try:
        x0_c = np.asarray(x0_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should contain only scalars.'.format(invoker))
    if len(x0_c.shape) != 1:
        raise ValueError('{}: UNEXPECTED ERROR: x0 should be a vector.'.format(invoker))
    lenx0 = x0_c.size

    # Validate lb.
    if isinstance(lb, scalar_types):
        lb_c = [lb]
    elif hasattr(lb, '__len__'):
        lb_c = lb
    else:
        raise ValueError('{}: UNEXPECTED ERROR: lb should be a vector.'.format(invoker))
    try:
        lb_c = np.asarray(lb_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: lb should contain only scalars.'.format(invoker))
    if len(lb_c.shape) != 1 or lb_c.size != lenx0:
        raise ValueError('{}: UNEXPECTED ERROR: the size of lb is inconsistent with x0.'.format(invoker))

    # Validate ub.
    if isinstance(ub, scalar_types):
        ub_c = [ub]
    elif hasattr(ub, '__len__'):
        ub_c = ub
    else:
        raise ValueError('{}: UNEXPECTED ERROR: ub should be a vector.'.format(invoker))
    try:
        ub_c = np.asarray(ub_c, dtype=np.float64)
    except ValueError:
        raise ValueError('{}: UNEXPECTED ERROR: ub should contain only scalars.'.format(invoker))
    if len(ub_c.shape) != 1 or ub_c.size != lenx0:
        raise ValueError('{}: UNEXPECTED ERROR: the size of ub is inconsistent with x0.'.format(invoker))

    # Validate constraints.
    if not isinstance(constraints, dict) or not ({'linear', 'nonlinear'} <= set(constraints.keys())) or \
            not (isinstance(constraints['linear'], LinearConstraint) or constraints['linear'] is None):
        # the nonlinear constraints will not be taken into account in this function and are, therefore, not validated
        raise ValueError('{}: UNEXPECTED ERROR: The constraints are ill-defined.'.format(invoker))

    max_con = 1e20  # Decide whether an inequality constraint can be ignored

    # Project onto the feasible set.
    if constraints['linear'] is None:
        # Direct projection onto the bound constraints
        x_proj = np.nanmin((np.nanmax((x0_c, lb_c), axis=0), ub_c), axis=0)
        return OptimizeResult(x=x_proj)
    elif all(np.less_equal(np.abs(constraints['linear'].ub - constraints['linear'].lb), eps)) and \
            np.max(lb_c) <= -max_con and np.min(ub_c) >= max_con:
        # The linear constraints are all equality constraints. The projection can therefore be done by solving the
        # least-squares problem: min ||A*x - (b - A*x_0)||.
        a = constraints['linear'].A
        b = (constraints['linear'].lb + constraints['linear'].ub) / 2
        xi, _, _, _ = np.linalg.lstsq(a, b - np.dot(a, x0_c), rcond=None)

        # The problem is not bounded. However, if the least-square solver returned values bigger in absolute value
        # than max_con, they will be reduced to this bound.
        x_proj = np.nanmin((np.nanmax((x0_c + xi, lb_c), axis=0), ub_c), axis=0)

        return OptimizeResult(x=x_proj)

    if constraints['linear'] is not None:
        try:
            # Project the initial guess onto the linear constraints via SciPy.
            from scipy.optimize import minimize
            from scipy.optimize import Bounds as ScipyBounds
            from scipy.optimize import LinearConstraint as ScipyLinearConstraint

            linear = constraints['linear']

            # To be more efficient, SciPy asks to separate the equality and the inequality constraints into two
            # different LinearConstraint structures
            pc_args_ineq, pc_args_eq = dict(), dict()
            pc_args_ineq['A'], pc_args_eq['A'] = np.asarray([[]]), np.asarray([[]])
            pc_args_ineq['A'] = pc_args_ineq['A'].reshape(0, linear.A.shape[1])
            pc_args_eq['A'] = pc_args_eq['A'].reshape(0, linear.A.shape[1])
            pc_args_ineq['lb'], pc_args_eq['lb'] = np.asarray([]), np.asarray([])
            pc_args_ineq['ub'], pc_args_eq['ub'] = np.asarray([]), np.asarray([])

            for i in range(linear.lb.size):
                if linear.lb[i] != linear.ub[i]:
                    pc_args_ineq['A'] = np.concatenate((pc_args_ineq['A'], linear.A[i:i+1, :]), axis=0)
                    pc_args_ineq['lb'] = np.r_[pc_args_ineq['lb'], linear.lb[i]]
                    pc_args_ineq['ub'] = np.r_[pc_args_ineq['ub'], linear.ub[i]]
                else:
                    pc_args_eq['A'] = np.concatenate((pc_args_eq['A'], linear.A[i:i+1, :]), axis=0)
                    pc_args_eq['lb'] = np.r_[pc_args_eq['lb'], linear.lb[i]]
                    pc_args_eq['ub'] = np.r_[pc_args_eq['ub'], linear.ub[i]]

            if pc_args_ineq['A'].size > 0 and pc_args_ineq['lb'].size > 0 and pc_args_eq['lb'].size > 0:
                project_constraints = [ScipyLinearConstraint(**pc_args_ineq), ScipyLinearConstraint(**pc_args_eq)]
            elif pc_args_ineq['A'].size > 0 and pc_args_ineq['lb'].size > 0:
                project_constraints = ScipyLinearConstraint(**pc_args_ineq)
            elif pc_args_eq['A'].size > 0:
                project_constraints = ScipyLinearConstraint(**pc_args_eq)
            else:
                project_constraints = ()

            # Perform the actual projection.
            ax_ineq = np.dot(pc_args_ineq['A'], x0_c)
            ax_eq = np.dot(pc_args_eq['A'], x0_c)
            if np.greater(ax_ineq, pc_args_ineq['ub']).any() or np.greater(pc_args_ineq['lb'], ax_ineq).any() or \
                    np.not_equal(ax_eq, pc_args_eq['lb']).any() or \
                    np.greater(x0_c, ub_c).any() or np.greater(lb_c, x0_c).any():
                return minimize(lambda x: np.dot(x - x0_c, x - x0_c) / 2, x0_c, jac=lambda x: (x - x0_c),
                                bounds=ScipyBounds(lb_c, ub_c), constraints=project_constraints)
            else:
                # Do not perform any projection if the initial guess is feasible.
                return OptimizeResult(x=x0_c)

        except ImportError:
            return OptimizeResult(x=x0_c)

    return OptimizeResult(x=x0_c)
