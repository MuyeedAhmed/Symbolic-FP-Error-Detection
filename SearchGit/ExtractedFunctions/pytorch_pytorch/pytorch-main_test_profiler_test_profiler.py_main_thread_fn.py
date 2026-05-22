        def main_thread_fn(profile_all_threads, returned_events):
            x, y = prep_inputs()
            experimental_config = torch._C._profiler._ExperimentalConfig(
                profile_all_threads=profile_all_threads
            )
            with torch.profiler.profile(
                experimental_config=experimental_config, record_shapes=True
            ) as p:
                profiling_started.set()
                for _ in range(n_rep):
                    _ = x @ y
                profiling_ended.wait()
            returned_events.append(p.events())