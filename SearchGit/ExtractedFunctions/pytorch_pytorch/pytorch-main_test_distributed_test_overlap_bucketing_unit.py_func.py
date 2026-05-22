        def func(a, b, c):
            group_name = "0"
            group_size = 1

            ag1 = torch.ops._c10d_functional.all_gather_into_tensor(
                a, group_size, group_name
            )
            b = torch.ops.prims.convert_element_type.default(b, torch.float16)
            ag1_out = torch.ops._c10d_functional.wait_tensor(ag1)

            ag2 = torch.ops._c10d_functional.all_gather_into_tensor(
                b, group_size, group_name
            )
            ag3 = torch.ops._c10d_functional.all_gather_into_tensor(
                c, group_size, group_name
            )

            mm = ag1_out @ ag1_out

            ag2_out = torch.ops._c10d_functional.wait_tensor(ag2)
            ag3_out = torch.ops._c10d_functional.wait_tensor(ag3)

            return ag1_out, ag2_out, ag3_out, mm