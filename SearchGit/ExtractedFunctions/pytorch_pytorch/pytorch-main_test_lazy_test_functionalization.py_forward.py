                def forward(self, x):
                    return x @ self.fc1.weight.transpose(0, 1)