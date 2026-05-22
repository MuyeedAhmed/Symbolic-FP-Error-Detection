            def forward(self, x):
                # Inner DDP forward runs first (like DP embedding lookup)
                x = self.inner_ddp(x)
                # After inner DDP exits, _active_ddp_module must still be set
                # for DDPOptimizer to split this region
                x = x @ self.weight
                return x