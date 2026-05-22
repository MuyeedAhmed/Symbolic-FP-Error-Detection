        def func_3(x):
            x = x @ x
            x = torch.ops._c10d_functional.all_reduce(x, "sum", "0")
            return torch.ops._c10d_functional.wait_tensor(x)