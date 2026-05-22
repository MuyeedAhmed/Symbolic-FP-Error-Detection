    def derphi(s):
        gval[0] = fprime(xk + s * pk, *args)
        gc[0] += 1
        return gval[0] @ pk