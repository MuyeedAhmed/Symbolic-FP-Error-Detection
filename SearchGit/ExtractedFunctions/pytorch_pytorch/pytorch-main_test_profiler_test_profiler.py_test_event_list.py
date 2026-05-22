    def test_event_list(self):
        # AFAIK event list is part of legacy profiler and/or used when kineto is not available.
        # This test has basic sanity checks to test against obvious regressions.
        x, y = (torch.rand((4, 4), requires_grad=True, device="cuda") for _ in range(2))
        with profile(with_stack=True) as p:
            z = (x @ y).relu().sum()
            z.backward()

        event_list = torch.autograd.profiler_util.EventList(p.events())
        # event_list._build_tree()

        with TemporaryFileName(mode="w+") as fname:
            event_list.export_chrome_trace(fname)
            with open(fname) as f:
                json.load(f)

        event_list.table()