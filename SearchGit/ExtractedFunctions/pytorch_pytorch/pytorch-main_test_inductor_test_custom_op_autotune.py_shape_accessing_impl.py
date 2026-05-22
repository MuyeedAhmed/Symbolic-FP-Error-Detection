        def shape_accessing_impl(
            mat1: torch.Tensor, mat2: torch.Tensor
        ) -> torch.Tensor:
            m, k = mat1.shape  # Shape access that would break naive make_fx
            n = mat2.shape[1]
            k_splits = 4
            if k % k_splits == 0:
                k_parts = k // k_splits
                a = torch.permute(mat1.reshape(m, k_splits, k_parts), (1, 0, 2))
                b = mat2.reshape(k_splits, k_parts, n)
                return torch.sum(torch.bmm(a, b), dim=0)
            return mat1 @ mat2