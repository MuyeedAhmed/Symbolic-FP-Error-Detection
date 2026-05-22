    def calculate_weight(
        self,
        weight,
        key,
        strength,
        strength_model,
        offset,
        function,
        intermediate_dtype=torch.float32,
        original_weight=None,
    ):
        v = self.weights
        blocks = v[0]
        rescale = v[1]
        alpha = v[2]
        if alpha is None:
            alpha = 0
        dora_scale = v[3]

        blocks = comfy.model_management.cast_to_device(
            blocks, weight.device, intermediate_dtype
        )
        if rescale is not None:
            rescale = comfy.model_management.cast_to_device(
                rescale, weight.device, intermediate_dtype
            )

        block_num, block_size, *_ = blocks.shape

        try:
            # Get r
            I = torch.eye(block_size, device=blocks.device, dtype=blocks.dtype)
            # for Q = -Q^T
            q = blocks - blocks.transpose(1, 2)
            normed_q = q
            if alpha > 0:  # alpha in oft/boft is for constraint
                q_norm = torch.norm(q) + 1e-8
                if q_norm > alpha:
                    normed_q = q * alpha / q_norm
            # use float() to prevent unsupported type in .inverse()
            r = (I + normed_q) @ (I - normed_q).float().inverse()
            r = r.to(weight)
            # Create I in weight's dtype for the einsum
            I_w = torch.eye(block_size, device=weight.device, dtype=weight.dtype)
            _, *shape = weight.shape
            lora_diff = torch.einsum(
                "k n m, k n ... -> k m ...",
                (r * strength) - strength * I_w,
                weight.view(block_num, block_size, *shape),
            ).view(-1, *shape)
            if dora_scale is not None:
                weight = weight_decompose(
                    dora_scale,
                    weight,
                    lora_diff,
                    alpha,
                    strength,
                    intermediate_dtype,
                    function,
                )
            else:
                weight += function((strength * lora_diff).type(weight.dtype))
        except Exception as e:
            logging.error("ERROR {} {} {}".format(self.name, key, e))
        return weight