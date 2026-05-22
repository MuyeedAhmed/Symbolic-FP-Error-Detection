    def test_leaf_function_constant_tensor_closure_error(self, backend):
        constant_weight = torch.randn(3, 3)

        @leaf_function
        def constant_closure_forward(x):
            return (x @ constant_weight,)

        @constant_closure_forward.register_fake
        def constant_closure_forward_fake(x):
            return (x @ constant_weight,)

        class ConstantClosureModule(torch.nn.Module):
            def __init__(self):
                super().__init__()

            def forward(self, x):
                return constant_closure_forward(x)

        mod = ConstantClosureModule()
        x = torch.randn(3, 3, requires_grad=True)

        result = mod(x)
        expected = x @ constant_weight
        self.assertEqual(result[0], expected)

        compiled_mod = torch.compile(mod, backend=backend, fullgraph=True)
        with self.assertRaisesRegex(
            Exception, "Please convert all Tensors to FakeTensors"
        ):
            compiled_mod(x)