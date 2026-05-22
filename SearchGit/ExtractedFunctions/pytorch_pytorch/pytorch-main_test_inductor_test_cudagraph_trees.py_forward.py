                def forward(self, x):
                    x1 = x + 1
                    y1 = x + 2
                    y_cpu = y1 + 1
                    z = x @ y1
                    inp = x1 + y1 + z + y_cpu
                    return self.linear(inp)