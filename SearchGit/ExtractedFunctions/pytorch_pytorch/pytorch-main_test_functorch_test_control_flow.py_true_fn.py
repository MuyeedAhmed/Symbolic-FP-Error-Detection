                def true_fn(x):
                    x.add_(1)
                    return x.sin() @ self.buf