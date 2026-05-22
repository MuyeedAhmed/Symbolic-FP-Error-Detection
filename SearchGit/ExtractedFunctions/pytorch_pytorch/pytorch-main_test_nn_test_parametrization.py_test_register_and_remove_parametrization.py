    def test_register_and_remove_parametrization(self):
        r"""Test that it is possible to add a few parametrizations
        on a parameter or a buffer and that removing them restores the initial state
        It also tests that backpropagating through them works as expected
        """

        # Define a couple matrix parametrizations
        class Skew(nn.Module):
            def forward(self, X):
                X = X.tril(-1)
                return X - X.T

        class Orthogonal(nn.Module):
            def forward(self, X):
                # Cayley map
                # If X is skew-symmetric it returns an orthogonal matrix
                Id = torch.eye(X.size(0), device=X.device)
                # We call contiguous because solve returns a tensor with strides that are Fortran-contiguous
                # and autograd raises a performance warning.
                # This happens when we remove the parametrization with leave_parametrized=True,
                # which does a set_ with a non-contiguous tensor while the gradient is contiguous
                return torch.linalg.solve(Id + X, Id - X).contiguous()

        class Resize(nn.Module):
            def forward(self, X):
                return X[[0]]

        class NoResize(nn.Module):
            def forward(self, X):
                return X

        # Define a couple vector parametrizations
        class FirstZero(nn.Module):
            def forward(self, x):
                return torch.cat([x.new_zeros(1), x[1:]])

        class LastZero(nn.Module):
            def forward(self, x):
                return torch.cat([x[:-1], x.new_zeros(1)])

        model = nn.Linear(8, 8)
        initial_weight_id = id(model.weight)
        initial_bias_id = id(model.bias)
        initial_model = deepcopy(model)

        # Test unsafe flag
        with self.assertRaisesRegex(
            ValueError,
            "Registering a parametrization may not change the shape of the tensor",
        ):
            parametrize.register_parametrization(
                model, "weight", Resize()
            )  # default unsafe = False
            model(torch.ones(8, 8))

        # One parametrization with unsafe=True
        parametrize.register_parametrization(model, "weight", Resize(), unsafe=True)
        self.assertTrue(hasattr(model, "parametrizations"))
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertFalse(parametrize.is_parametrized(model, "bias"))
        self.assertNotIn("weight", model._parameters)
        self.assertTrue(model.weight.shape[0] == 1)
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=False)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.weight, initial_model.weight)
        self.assertEqual(id(model.weight), initial_weight_id)
        self.assertEqual(model.__class__, nn.Linear)

        # Two parametrizations with unsafe=True
        parametrize.register_parametrization(model, "weight", Resize(), unsafe=True)
        parametrize.register_parametrization(model, "weight", NoResize(), unsafe=False)
        self.assertTrue(hasattr(model, "parametrizations"))
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertFalse(parametrize.is_parametrized(model, "bias"))
        self.assertNotIn("weight", model._parameters)
        self.assertTrue(model.weight.shape[0] == 1)
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=False)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.weight, initial_model.weight)
        self.assertEqual(id(model.weight), initial_weight_id)
        self.assertEqual(model.__class__, nn.Linear)

        # Test unsafe flag doesn't change expected behavior
        parametrize.register_parametrization(model, "weight", Skew(), unsafe=True)
        self.assertTrue(hasattr(model, "parametrizations"))
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertFalse(parametrize.is_parametrized(model, "bias"))
        self.assertNotIn("weight", model._parameters)
        # Result should be skew-symmetric
        A = model.weight
        self.assertEqual(A, -A.T)
        if get_swap_module_params_on_conversion():
            # When using the swap_tensors path, this is needed so that the autograd
            # graph is not alive anymore.
            del A
        # Remove and check consistency
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=False)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.weight, initial_model.weight)
        self.assertEqual(id(model.weight), initial_weight_id)
        self.assertEqual(model.__class__, nn.Linear)

        # Test one parametrization
        parametrize.register_parametrization(model, "weight", Skew())
        self.assertTrue(hasattr(model, "parametrizations"))
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertFalse(parametrize.is_parametrized(model, "bias"))
        self.assertNotIn("weight", model._parameters)
        # Result should be skew-symmetric
        A = model.weight
        self.assertEqual(A, -A.T)
        if get_swap_module_params_on_conversion():
            # When using the swap_tensors path, this is needed so that the autograd
            # graph is not alive anymore.
            del A
        # Remove and check consistency
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=False)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.weight, initial_model.weight)
        self.assertEqual(id(model.weight), initial_weight_id)
        self.assertEqual(model.__class__, nn.Linear)

        # Test two parametrizations at the same time and removing them
        parametrize.register_parametrization(model, "weight", Skew())
        parametrize.register_parametrization(model, "weight", Orthogonal())
        # Result should be orthogonal
        X = model.weight
        Id = torch.eye(X.size(0), device=X.device)
        self.assertEqual(X.T @ X, Id)
        if get_swap_module_params_on_conversion():
            # When using the swap_tensors path, this is needed so that the autograd
            # graph is not alive anymore.
            del X
        # Structure tests
        self.assertTrue(hasattr(model, "parametrizations"))
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertFalse(parametrize.is_parametrized(model, "bias"))
        self.assertIn("weight", model.parametrizations)
        self.assertNotIn("weight", model._parameters)
        # Remove
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=False)
        self.assertEqual(model.weight, initial_model.weight)
        self.assertEqual(id(model.weight), initial_weight_id)
        self.assertFalse(hasattr(model, "parametrizations"))
        self.assertEqual(model.__class__, nn.Linear)

        # Add everything
        parametrize.register_parametrization(model, "weight", Skew())
        parametrize.register_parametrization(model, "weight", Orthogonal())
        parametrize.register_parametrization(model, "bias", FirstZero())
        parametrize.register_parametrization(model, "bias", LastZero())

        # Basic tests
        self.assertTrue(parametrize.is_parametrized(model))
        self.assertTrue(parametrize.is_parametrized(model, "weight"))
        self.assertTrue(parametrize.is_parametrized(model, "bias"))
        self.assertEqual(model.bias[0].item(), 0.0)
        self.assertEqual(model.bias[-1].item(), 0.0)
        self.assertEqual(len(list(model.parameters())), 2)  # Nothing weird has happened
        # Should not throw

        sgd = torch.optim.SGD(model.parameters(), lr=0.01)

        weight_copy = model.weight.clone()
        bias_copy = model.bias.clone()
        sgd.zero_grad()
        (model.weight.T @ model.bias).sum().backward()
        sgd.step()
        self.assertNotEqual(model.weight, weight_copy)
        self.assertNotEqual(model.bias, bias_copy)

        # Remove first parametrization.
        # Check that the model is still parametrized and so is the second parameter
        parametrize.remove_parametrizations(model, "weight", leave_parametrized=False)
        self.assertTrue(parametrize.is_parametrized(model))  # Still parametrized
        self.assertFalse(
            parametrize.is_parametrized(model, "weight")
        )  # Parametrization removed
        self.assertTrue(
            parametrize.is_parametrized(model, "bias")
        )  # Still parametrized
        self.assertEqual(model.bias[0].item(), 0.0)  # Still parametrized
        self.assertEqual(model.bias[-1].item(), 0.0)  # Still parametrized
        self.assertNotEqual(model.weight, initial_model.weight)  # Has been updated
        self.assertEqual(id(model.weight), initial_weight_id)  # Keeps the same id
        self.assertEqual(len(list(model.parameters())), 2)  # Nothing weird has happened
        # Should not throw
        weight_copy = model.weight.clone()
        bias_copy = model.bias.clone()
        sgd.zero_grad()
        (model.weight.T @ model.bias).sum().backward()
        sgd.step()
        self.assertNotEqual(model.weight, weight_copy)
        self.assertNotEqual(model.bias, bias_copy)

        # Remove the second parametrization.
        # Check that the module is not parametrized
        parametrize.remove_parametrizations(model, "bias", leave_parametrized=False)
        self.assertFalse(parametrize.is_parametrized(model))  # Not parametrized
        self.assertNotEqual(model.bias, initial_model.bias)  # Has been updated
        self.assertNotEqual(model.bias[0].item(), 0.0)  # Not parametrized
        self.assertNotEqual(model.bias[-1].item(), 0.0)  # Not parametrized
        self.assertEqual(id(model.bias), initial_bias_id)  # Keeps the same id
        self.assertFalse(
            hasattr(model, "parametrizations")
        )  # Not parametrized the module
        self.assertEqual(model.__class__, nn.Linear)  # Resores the previous class
        self.assertEqual(len(list(model.parameters())), 2)  # Nothing weird has happeed

        # Should not throw things are updated
        weight_copy = model.weight.clone()
        bias_copy = model.bias.clone()
        sgd.zero_grad()
        (model.weight.T @ model.bias).sum().backward()
        sgd.step()
        self.assertNotEqual(model.weight, weight_copy)
        self.assertNotEqual(model.bias, bias_copy)
        if get_swap_module_params_on_conversion():
            # When using the swap_tensors path, this is needed so that the autograd
            # graph is not alive anymore.
            del weight_copy, bias_copy

        # Test leave_parametrized=True
        for _ in range(2):
            parametrize.register_parametrization(model, "weight", Skew())
            parametrize.register_parametrization(model, "weight", Orthogonal())
            parametrize.remove_parametrizations(
                model, "weight", leave_parametrized=True
            )
            # We didn't change the dtype nor had multiple inputs, so the id should be the same
            self.assertEqual(id(model.weight), initial_weight_id)
            self.assertEqual(id(model.bias), initial_bias_id)

            # Should not throw. Things are updated
            weight_copy = model.weight.clone()
            bias_copy = model.bias.clone()
            sgd.zero_grad()
            (model.weight.T @ model.bias).sum().backward()
            sgd.step()
            self.assertNotEqual(model.weight, weight_copy)
            self.assertNotEqual(model.bias, bias_copy)
            if get_swap_module_params_on_conversion():
                # When using the swap_tensors path, this is needed so that the autograd
                # graph is not alive anymore.
                del weight_copy, bias_copy