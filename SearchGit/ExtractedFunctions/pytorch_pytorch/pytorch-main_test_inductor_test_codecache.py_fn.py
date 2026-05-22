        def fn():
            inp = torch.randn(32, 50, 768, device=GPU_TYPE)
            weight = torch.randn(768, 768, device=GPU_TYPE)
            layer = torch.nn.LayerNorm(768, device=GPU_TYPE)
            return layer(inp @ weight)