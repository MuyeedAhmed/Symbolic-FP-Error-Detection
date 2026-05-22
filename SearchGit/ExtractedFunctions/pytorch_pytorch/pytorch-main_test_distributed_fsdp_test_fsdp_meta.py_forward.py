            def forward(self, x: torch.Tensor) -> torch.Tensor:
                return x @ self.weight