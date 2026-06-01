        def check(tensors, out, msg):
            with self.assertRaisesRegex(RuntimeError, msg):
                torch.linalg.multi_dot(tensors, out=out)