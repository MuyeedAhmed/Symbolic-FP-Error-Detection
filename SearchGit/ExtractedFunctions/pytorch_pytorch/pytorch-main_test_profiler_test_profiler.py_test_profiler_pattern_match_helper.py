    def test_profiler_pattern_match_helper(self):
        x = torch.ones((100, 100))
        with profile() as prof:
            for _ in range(5):
                x = x @ x
                x = x + x
        event_tree = prof.profiler.kineto_results.experimental_event_tree()
        pattern = Pattern(prof)
        self.assertEqual([], pattern.siblings_of(event_tree[0])[0])
        self.assertEqual(event_tree[1:], pattern.siblings_of(event_tree[0])[1])
        child_nodes = event_tree[0].children
        self.assertEqual([], pattern.siblings_of(child_nodes[0])[0])
        self.assertEqual(child_nodes[1:], pattern.siblings_of(child_nodes[0])[1])
        self.assertEqual(
            event_tree[0], pattern.root_of(event_tree[0].children[0].children[0])
        )
        self.assertEqual(None, pattern.next_of(event_tree[-1]))
        self.assertEqual(event_tree[1], pattern.next_of(event_tree[0]))
        self.assertEqual(event_tree[0], pattern.prev_of(event_tree[1]))