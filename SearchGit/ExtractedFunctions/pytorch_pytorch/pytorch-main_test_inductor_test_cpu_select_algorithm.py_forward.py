            def forward(self, x, weight):
                return x @ weight.permute(2, 0, 1)