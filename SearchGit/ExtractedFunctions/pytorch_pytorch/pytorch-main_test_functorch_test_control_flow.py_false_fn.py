                def false_fn(x):
                    self.buf.add_(1)
                    return x.sin() @ self.buf