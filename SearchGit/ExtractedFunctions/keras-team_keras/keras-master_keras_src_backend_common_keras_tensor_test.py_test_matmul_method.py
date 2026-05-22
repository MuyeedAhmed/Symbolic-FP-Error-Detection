    def test_matmul_method(self, mock_method):
        y = Mock()
        self._test_binary_op_method(mock_method, y, lambda x, y: x @ y)