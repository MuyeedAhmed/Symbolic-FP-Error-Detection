            def f(x, y):
                x1 = x + 1
                y1 = y + 1
                y_cpu = y1.cpu() + 1
                z = x1 + y1 + x @ y
                u = (y_cpu.to("cuda") + 2) @ y + 3
                u_cpu = u.cpu() + 2
                return z + u_cpu.to("cuda")