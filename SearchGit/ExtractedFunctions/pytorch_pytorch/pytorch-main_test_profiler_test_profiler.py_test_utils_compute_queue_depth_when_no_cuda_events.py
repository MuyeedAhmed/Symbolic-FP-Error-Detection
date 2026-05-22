    def test_utils_compute_queue_depth_when_no_cuda_events(self):
        # For traces with only cpu events, we expect empty queue depth list
        x = torch.ones((100, 100))
        with profile() as prof:
            for _ in range(5):
                x = x @ x
        basic_evaluation = _utils.BasicEvaluation(prof.profiler)
        self.assertFalse(basic_evaluation.compute_queue_depth())