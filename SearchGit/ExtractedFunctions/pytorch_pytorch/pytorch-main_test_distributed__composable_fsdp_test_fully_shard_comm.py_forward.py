            def forward(self, x: torch.Tensor):
                y = F.relu(x @ self.weight)
                # NOTE: This all-reduce is not differentiable and is included
                # to exercise the overlap.
                work = dist.all_reduce(y, group=self.mesh.get_group(), async_op=True)
                return y, work