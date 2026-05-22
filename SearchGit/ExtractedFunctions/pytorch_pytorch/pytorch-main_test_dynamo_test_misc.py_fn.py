        def fn(x, targets):
            # Forward computation before detach (e.g. transformer layers)
            h = x * 2 + 1
            x_detached = h.detach().requires_grad_()
            chunksz = x_detached.shape[0] // 2
            total_loss = torch.tensor(0.0)
            for start in range(0, x_detached.shape[0], chunksz):
                chunk = x_detached[start : start + chunksz]
                chunk_targets = targets[start : start + chunksz]
                logits = chunk @ torch.eye(chunk.shape[-1])
                loss = torch.nn.functional.cross_entropy(logits, chunk_targets)
                loss.backward()
                total_loss = total_loss + loss.detach()
            # Propagate chunked grad back through the forward computation
            h.backward(x_detached.grad)
            return x.grad, total_loss