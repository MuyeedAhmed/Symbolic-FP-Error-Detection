        def fn(x, y):
            z = x @ y
            with torch.profiler.record_function("this_is_ignored"):
                z = z + 1
            return z.sum()