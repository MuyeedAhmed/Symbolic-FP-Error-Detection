    def test_schema_check_mode_functionality_list_input(self):
        a = torch.rand((3, 3))
        b = torch.rand((3, 3))
        c = torch.rand((3, 3))
        expected = torch.linalg.multi_dot([a, b, c])
        with SchemaCheckMode():
            actual = torch.linalg.multi_dot([a, b, c])
        self.assertEqual(expected, actual)