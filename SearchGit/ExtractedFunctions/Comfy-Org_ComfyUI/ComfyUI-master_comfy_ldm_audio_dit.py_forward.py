    def forward(self, input):
        f = 2 * math.pi * input @ comfy.ops.cast_to_input(self.weight.T, input)
        return torch.cat([f.cos(), f.sin()], dim=-1)