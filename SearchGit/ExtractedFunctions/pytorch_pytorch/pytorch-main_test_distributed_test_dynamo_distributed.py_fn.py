        def fn(output, input, w):
            result = dist._all_gather_base(output, input, async_op=False)
            assert result is None  # noqa: S101
            return output @ w