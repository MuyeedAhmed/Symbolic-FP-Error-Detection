            def forward(self, inp, *, kwarg=None):
                if kwarg is not None:
                    inp = inp @ kwarg
                return self.a(inp)