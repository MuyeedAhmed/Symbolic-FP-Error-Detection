    def _nested_rpc_call_backward_error(t1, t2, dst):
        t1 = t1 * t2
        t2 = t1 + t2
        res = rpc.rpc_sync(
            worker_name(dst),
            DistAutogradTest._python_udf_with_backward_error,
            args=(t1, t2),
        )
        return torch.linalg.multi_dot([t1, t2, res])