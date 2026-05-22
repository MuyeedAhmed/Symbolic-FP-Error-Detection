            def forward(self, x):
                return self.layer(x @ self.weight + self.buf1)