        def test_multiple_fusions(x):
            y = x.to(torch.float)
            return (y @ y).relu()