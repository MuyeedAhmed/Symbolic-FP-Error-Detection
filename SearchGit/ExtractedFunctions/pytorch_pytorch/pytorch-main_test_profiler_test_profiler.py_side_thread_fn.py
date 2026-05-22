        def side_thread_fn():
            x, y = prep_inputs()
            profiling_started.wait()
            for _ in range(n_rep):
                _ = x @ y
            profiling_ended.set()