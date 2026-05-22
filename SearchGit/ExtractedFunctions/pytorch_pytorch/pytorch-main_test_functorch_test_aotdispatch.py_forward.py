            def forward(self, x):
                self.buf.add_(1)
                return (self.w @ x).sum() + self.buf.sum()