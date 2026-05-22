    def forward(self, xs: tuple[torch.Tensor, torch.Tensor]) -> torch.Tensor:
        if len(xs) != 2:
            raise AssertionError(f"Expects 2 args but got {len(xs)}")
        x, y = xs
        z = x + y
        z = z @ self.w
        return nn.functional.relu(z)