                def fn(x):
                    return torch.nn.functional.dropout(torch.sin(x @ self.w), p=0.5)