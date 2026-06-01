            def _expert_linear_impl(self, input, weight, bias, i):
                if isinstance(weight, QuantizedTensor):
                    qw = self._expert_qt_from(weight, i)
                else:
                    qw = weight[i]
                b = cast_to_input(bias[i], input, copy=False) if bias is not None else None

                if isinstance(qw, QuantizedTensor):
                    use_fast = (
                        not self._full_precision_mm
                        and qw.layout_cls.supports_fast_matmul()
                        and input.dim() == 2
                    )
                    if use_fast:
                        qin = QuantizedTensor.from_float(input, self.layout_type)
                        return torch.nn.functional.linear(qin, qw, b)
                    out = input @ qw.dequantize().t()
                    return out + b if b is not None else out
                return torch.nn.functional.linear(input, qw, b)