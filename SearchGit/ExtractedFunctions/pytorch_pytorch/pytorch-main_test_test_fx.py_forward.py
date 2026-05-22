            def forward(self, x):
                return torch.dot(self.W, x)