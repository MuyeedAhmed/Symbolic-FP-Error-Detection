        def unsafe_impl(mat1: torch.Tensor, mat2: torch.Tensor) -> torch.Tensor:
            m = mat1.shape[0]
            if m % 4 == 0:
                n = mat2.shape[1]
                k = mat1.shape[1]
                m_parts = m // 4
                a = mat1.reshape(4, m_parts, k)
                result = torch.bmm(a, mat2.unsqueeze(0).expand(4, -1, -1))
                return result.reshape(m, n)
            return mat1 @ mat2