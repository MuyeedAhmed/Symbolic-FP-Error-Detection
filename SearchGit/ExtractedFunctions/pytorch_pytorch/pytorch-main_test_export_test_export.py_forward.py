            def forward(self, x: torch.Tensor, y: torch.Tensor, zs: list[torch.Tensor]):
                return x[:3] + y @ torch.cat(zs)