        def fn(x, w):
            s = torch.cuda.Stream()
            e = torch.cuda.Event()

            h = x
            for _ in range(4):
                h = h @ w
            e.record()

            with torch.cuda.stream(s):
                e.wait()
                out = torch.relu(h) + 1.0

            s.synchronize()
            return out