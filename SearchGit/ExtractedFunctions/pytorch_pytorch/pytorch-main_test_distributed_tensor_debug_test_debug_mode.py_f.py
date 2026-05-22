        def f(dA, dB):
            dy = dA @ dB
            loss = dy.sum()
            loss.backward()
            return dA.grad, dB.grad