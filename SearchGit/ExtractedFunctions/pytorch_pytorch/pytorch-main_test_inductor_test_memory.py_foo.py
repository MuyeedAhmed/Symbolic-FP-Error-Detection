        def foo(inp, inp2):
            inp = inp @ inp
            inp = inp.view(2, -1, 256)
            x = inp[0]
            y = inp[1]
            x, y = torch._foreach_add([x, y], 1.0)
            out = x.sum()
            out2 = y.sum(dim=-1)

            return out, out2, inp2 @ inp2