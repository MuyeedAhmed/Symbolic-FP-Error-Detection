            def foo(x):
                # partition 1
                x1 = x @ x
                y1 = x1 + 1
                z_cpu = y1.cpu() + 1
                # partition 2
                # partition 2 should reuse the fused triton kernel generated
                # in partition 1
                x2 = z_cpu.to("cuda") @ z_cpu.to("cuda")
                y2 = x2 + 1
                return y1, y2