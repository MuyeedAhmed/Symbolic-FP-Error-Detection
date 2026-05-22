        def shape_tracker_op(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            # This runs during benchmarking with REAL values
            shapes_seen.append(x.shape[0])
            return x @ weight