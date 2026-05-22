    def test_matmul(self):
        operator = self.module
        self.assertRaises(TypeError, operator.matmul)
        self.assertRaises(TypeError, operator.matmul, 42, 42)
        class M:
            def __matmul__(self, other):
                return other - 1
        self.assertEqual(M() @ 42, 41)