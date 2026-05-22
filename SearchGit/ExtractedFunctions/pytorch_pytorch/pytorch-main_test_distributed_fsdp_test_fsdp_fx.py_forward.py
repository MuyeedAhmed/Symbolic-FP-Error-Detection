    def forward(self, x: torch.Tensor, run_all_layers: bool) -> torch.Tensor:
        z = self.relu(self.layer0(x))
        z = self.relu(self.layer2(z))
        z = z @ self.weight1
        if run_all_layers:
            z = self.relu(self.layer1(z))
            z = z @ self.weight2
            # Use `layer0` twice to check the handling of multiplicity in the
            # saved data structures
            z = self.relu(self.layer0(x))
        return z