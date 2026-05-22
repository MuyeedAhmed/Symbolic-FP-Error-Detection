    def test_rmatmul_method(self, mock_symbolic_call):
        """Test __rmatmul__ method"""
        mock_tensor = Mock()
        mock_symbolic_call.return_value = mock_tensor
        x = keras_tensor.KerasTensor(shape=(3, 4), dtype="float32")
        y = Mock()
        result = y @ x
        mock_symbolic_call.assert_called_once_with(y, x)
        self.assertEqual(result, mock_tensor)