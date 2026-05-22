    def forward(self, x):
        x = x @ self.weight
        x = self.bias_module0(x)
        x = self.relu(x)  # ensure biases have different gradients
        x = self.bias_module1(x)
        return x