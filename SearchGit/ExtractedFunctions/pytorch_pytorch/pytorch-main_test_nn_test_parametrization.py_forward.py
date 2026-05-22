            def forward(self, x, y):
                # Form a rank-1 matrix from a pair of vectors
                return x.unsqueeze(-1) @ y.unsqueeze(-2)