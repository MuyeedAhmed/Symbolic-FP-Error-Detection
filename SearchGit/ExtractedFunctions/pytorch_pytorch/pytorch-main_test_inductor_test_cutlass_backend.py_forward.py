            def forward(self, a, b, extra_args):
                acc = a @ b
                return acc, op(acc.relu(), *extra_args)