        def forward(self, x):
            return self.base(x @ self.weight.T)