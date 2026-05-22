def _line_search_wolfe1(
    f,
    fprime,
    xk,
    pk,
    gfk=None,
    old_fval=None,
    old_old_fval=None,
    args=(),
    c1=1e-4,
    c2=0.9,
    amax=50,
    amin=1e-8,
    xtol=1e-14,
):
    """
    Same as `scalar_search_wolfe1` but do a line search to direction `pk`
    """
    if gfk is None:
        gfk = fprime(xk, *args)

    gval = [gfk]
    gc = [0]
    fc = [0]

    def phi(s):
        fc[0] += 1
        return f(xk + s * pk, *args)

    def derphi(s):
        gval[0] = fprime(xk + s * pk, *args)
        gc[0] += 1
        return gval[0] @ pk

    derphi0 = gfk @ pk

    stp, fval, old_fval = scalar_search_wolfe1(
        phi,
        derphi,
        old_fval,
        old_old_fval,
        derphi0,
        c1=c1,
        c2=c2,
        amax=amax,
        amin=amin,
        xtol=xtol,
    )

    return stp, fc[0], gc[0], fval, old_fval, gval[0]