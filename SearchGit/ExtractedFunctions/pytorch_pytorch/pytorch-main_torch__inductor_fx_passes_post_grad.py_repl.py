    def repl(inp, x1, x2, alpha, beta):
        mm_result = x1 @ x2
        if alpha != 1:
            mm_result = alpha * mm_result
        if beta != 1:
            inp = beta * inp
        return inp + mm_result