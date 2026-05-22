        def fn(x, weight):
            qkv = x @ weight
            qkv = qkv.view(B, M, 3, H, D)
            q, k, v = qkv.unbind(2)
            q, k, v = (t.transpose(1, 2) for t in (q, k, v))
            return flex_attention(q, k, v, kernel_options={"BACKEND": "FLASH"})