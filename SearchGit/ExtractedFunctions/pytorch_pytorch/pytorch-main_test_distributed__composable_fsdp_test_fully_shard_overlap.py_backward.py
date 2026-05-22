    def backward(ctx, grad_output: torch.Tensor):
        (input, weight) = ctx.saved_tensors
        torch.get_device_module(device_type)._sleep(
            int(2 * ctx.sleep_ms * get_cycles_per_ms())
        )
        grad_input = grad_output @ weight.T
        grad_weight = input.T @ grad_output
        return grad_input, grad_weight, None