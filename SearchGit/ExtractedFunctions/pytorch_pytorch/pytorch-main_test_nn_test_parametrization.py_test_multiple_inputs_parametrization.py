    def test_multiple_inputs_parametrization(self):
        # A parametrization with several outputs
        class RankOne(nn.Module):
            def forward(self, x, y):
                # Form a rank-1 matrix from a pair of vectors
                return x.unsqueeze(-1) @ y.unsqueeze(-2)

            def right_inverse(self, Y):
                # We project the given matrix onto the rank 1 matrices
                U, S, Vh = torch.linalg.svd(Y, full_matrices=False)
                # S is ordered in a decreasing way.
                s0_sqrt = S[0].sqrt().unsqueeze(-1)
                return U[..., :, 0] * s0_sqrt, Vh[..., 0, :] * s0_sqrt

        # Simple parametrisation
        class Double(nn.Module):
            def forward(self, x):
                return 2.0 * x

            def right_inverse(self, w):
                return 0.5 * w

        model = nn.Linear(3, 3)
        # Test one parametrization
        parametrize.register_parametrization(model, "weight", RankOne())
        self.assertTrue(hasattr(model, "parametrizations"))
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertTrue(hasattr(model.parametrizations.weight, "original0"))
        self.assertIn("original0", model.parametrizations.weight._parameters)
        self.assertTrue(hasattr(model.parametrizations.weight, "original1"))
        self.assertIn("original1", model.parametrizations.weight._parameters)
        self.assertFalse(parametrize.is_parametrized(model, "bias"))
        self.assertNotIn("weight", model._parameters)
        # Result should be rank 1
        self.assertEqual(torch.linalg.matrix_rank(model.weight).item(), 1)

        with self.assertRaisesRegex(ValueError, "leave_parametrized=False"):
            # Cannot remove a parametrization with multiple inputs and not leave it parametrized
            parametrize.remove_parametrizations(
                model, "weight", leave_parametrized=False
            )
        # Remove parametrization and check consistency
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=True)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.__class__, nn.Linear)
        self.assertFalse(parametrize.is_parametrized(model))
        self.assertEqual(torch.linalg.matrix_rank(model.weight).item(), 1)
        self.assertIn("weight", model._parameters)

        # Registering parametrizations with one input on top of one with multiple inputs should work
        init_weight = model.weight.clone()
        parametrize.register_parametrization(model, "weight", RankOne())
        # Projecting a rank 1 matrix onto the matrices of rank one does not change the matrix
        self.assertEqual(init_weight, model.weight)
        parametrize.register_parametrization(model, "weight", Double())
        # The matrix now is twice the initial matrix
        self.assertEqual(2.0 * init_weight, model.weight)
        # Multiplying by a scalar does not change the rank
        self.assertEqual(torch.linalg.matrix_rank(model.weight).item(), 1)

        # The model has now three parameters
        self.assertEqual(len(list(model.parameters())), 3)

        sgd = torch.optim.SGD(model.parameters(), lr=0.1)

        # Test backward. Should not throw
        for _ in range(2):
            sgd.zero_grad()
            loss = (model.weight.T @ model.bias).sum()
            loss.backward()
            sgd.step()

        # Same drill as before, removing should work as expected
        with self.assertRaisesRegex(ValueError, "leave_parametrized=False"):
            # Cannot remove a parametrization with multiple inputs and not leave it parametrized
            parametrize.remove_parametrizations(
                model, "weight", leave_parametrized=False
            )
        # Remove parametrization and check consistency
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=True)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.__class__, nn.Linear)
        self.assertFalse(parametrize.is_parametrized(model))
        self.assertEqual(torch.linalg.matrix_rank(model.weight).item(), 1)
        self.assertIn("weight", model._parameters)

        # The model has now two parameters
        self.assertEqual(len(list(model.parameters())), 2)

        # Test backward. Should not throw
        sgd = torch.optim.SGD(model.parameters(), lr=0.1)
        for _ in range(2):
            sgd.zero_grad()
            loss = (model.weight.T @ model.bias).sum()
            loss.backward()
            sgd.step()