            def forward(self, x):
                return checkpoint(
                    lambda x: torch.sin(x @ self.w), x, use_reentrant=False
                )