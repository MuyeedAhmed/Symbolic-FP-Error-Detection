        def main_with_thread_fn(profile_all_threads):
            x, y = prep_inputs()
            experimental_config = torch._C._profiler._ExperimentalConfig(
                profile_all_threads=profile_all_threads
            )
            with torch.profiler.profile(
                experimental_config=experimental_config, record_shapes=True
            ) as p:
                side_thread = threading.Thread(target=side_thread_fn)
                side_thread.start()
                for _ in range(n_rep):
                    _ = x @ y
                side_thread.join()
            return p.events()