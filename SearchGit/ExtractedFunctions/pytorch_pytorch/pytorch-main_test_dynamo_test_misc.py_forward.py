            def forward(self, x):
                return x @ self.other_attr + self.queue[-1]