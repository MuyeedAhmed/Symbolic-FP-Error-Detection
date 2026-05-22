    def test_profiler_name_pattern(self):
        x = torch.ones((100, 100))
        with profile() as prof:
            for _ in range(5):
                x = x @ x
                x = x + x
        matched_events = NamePattern(prof, "aten::mm").matched_events()
        output = "\n".join([f"{event.name}" for event in matched_events])
        self.assertExpectedInline(
            output,
            """\
aten::mm
aten::mm
aten::mm
aten::mm
aten::mm""",
        )