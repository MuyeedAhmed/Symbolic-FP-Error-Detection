    def _test_fsdp_not_all_outputs_used_in_loss(
        self, sharding_strategy: ShardingStrategy
    ):
        class MyModule(nn.Module):
            def __init__(self) -> None:
                super().__init__()
                self.lin1 = nn.Linear(4, 4)
                self.lin2 = nn.Linear(4, 4)

            def forward(self, x):
                a = self.lin1(x)
                b = self.lin2(x)
                return (a, b)

        def _check_resharded(fsdp_module):
            handle = fsdp_module._handle
            if not handle:
                return
            param = handle.flat_param
            if handle.uses_sharded_strategy:
                full_param = param._full_param_padded
                self.assertEqual(full_param.storage().size(), 0)

            self.assertEqual(param.data_ptr(), param._local_shard.data_ptr())

        def _check_equal(local, fsdp):
            with FSDP.summon_full_params(fsdp):
                for p1, p2 in zip(fsdp.parameters(), local.parameters()):
                    torch.testing.assert_close(p1, p2)

        fsdp_ctor = functools.partial(FSDP, sharding_strategy=sharding_strategy)
        m = MyModule().to(device=device_type)
        m_local = deepcopy(m)
        local_m = m_local
        prev_params = [p.clone() for p in m_local.parameters()]

        m.lin1 = fsdp_ctor(m.lin1)
        m = fsdp_ctor(m)
        _check_equal(m_local, m)

        opt = torch.optim.SGD(m.parameters(), lr=1e-3)
        opt_local = torch.optim.SGD(local_m.parameters(), lr=1e-3)

        for i in range(6):
            t = torch.ones(4, device=device_type)
            a, b = m(t)
            local_a, local_b = local_m(t)
            if i < 2:
                # use both params in loss computation. Later,
                # b will go unused and we check grads are the
                # same as local training.
                loss = (a @ b).sum()
                loss_local = (local_a @ local_b).sum()
            else:
                loss = a.sum()
                loss_local = local_a.sum()

            loss.backward()
            loss_local.backward()
            _check_resharded(m)
            opt.step()
            opt_local.step()
            _check_equal(m_local, m)
            # Ensure at least some change from previous params, otherwise
            # above check would be vacuously true.
            self.assertTrue(
                any(
                    not torch.equal(p1, p2)
                    for p1, p2 in zip(prev_params, m_local.parameters())
                )
            )
            prev_params = [p.clone() for p in local_m.parameters()]
            opt.zero_grad()
            opt_local.zero_grad()

        dist.barrier()