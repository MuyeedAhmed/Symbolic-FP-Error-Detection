    def test_func_slots(self):
        cls = _testlimitedcapi.type_from_slots("matmul_123")
        self.assertEqual(cls() @ None, 123)