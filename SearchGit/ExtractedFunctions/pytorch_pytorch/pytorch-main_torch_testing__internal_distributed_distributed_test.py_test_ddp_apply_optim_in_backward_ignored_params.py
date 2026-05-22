        def test_ddp_apply_optim_in_backward_ignored_params(self):
            torch.cuda.set_device(self.rank)
            for init_before in [True, False]:
                with self.subTest(init_before=init_before):
                    torch.manual_seed(self.rank)
                    torch.cuda.manual_seed(self.rank)
                    model = TwoLinLayerNet()
                    # Parameters to ignore are in the format {module_name}.{param_name}
                    params_to_ignore = ["a.weight"]
                    torch.nn.parallel.DistributedDataParallel._set_params_and_buffers_to_ignore_for_model(
                        model, params_to_ignore
                    )
                    if init_before:
                        _apply_optimizer_in_backward(
                            optimizer_class=torch.optim.SGD,
                            params=model.parameters(),
                            optimizer_kwargs={"lr": 0.03},
                        )
                    net = torch.nn.parallel.DistributedDataParallel(
                        model.cuda(self.rank),
                        device_ids=[self.rank],
                    )
                    if not init_before:
                        _apply_optimizer_in_backward(
                            optimizer_class=torch.optim.SGD,
                            params=model.parameters(),
                            optimizer_kwargs={"lr": 0.03},
                        )
                    inp = torch.randn(1, 10)
                    a, b = net(inp)
                    (a.transpose(0, 1) @ b).sum().backward()
                    # a.weight did not go through allreduce, so optimizer acted on local
                    # gradient, which should be different across ranks. Remaining params
                    # should be equal.
                    models = [None for _ in range(dist.get_world_size())]
                    dist.all_gather_object(models, model)
                    rank0_model, remainder = models[0], models[1:]
                    for m in remainder:
                        self.assertNotEqual(rank0_model.a.weight, m.a.weight)
                        self.assertEqual(
                            list(rank0_model.b.parameters()), list(m.b.parameters())
                        )
                        self.assertEqual(rank0_model.a.bias, m.a.bias)