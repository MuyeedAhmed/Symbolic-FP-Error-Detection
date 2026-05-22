            def forward(self, x, y, a, b):
                u0, s0 = a.item(), b.item()
                u_max = max(u0, 15)
                # construct the equality rule Max(15, u0) == s0 * Max(15, u0)
                torch._check(u_max == s0 * u_max)
                # size x - [Max(u0, 15), 64]
                x = x.expand(u_max, *x.shape).clone()
                return x @ y