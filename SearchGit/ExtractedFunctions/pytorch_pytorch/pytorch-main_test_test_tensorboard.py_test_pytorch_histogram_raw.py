    def test_pytorch_histogram_raw(self):
        with self.createSummaryWriter() as w:
            num = 50
            floats = make_np(torch.rand((num,)))
            bins = [0.0, 0.25, 0.5, 0.75, 1.0]
            counts, limits = np.histogram(floats, bins)
            sum_sq = floats.dot(floats).item()
            w.add_histogram_raw(
                "float histogram raw",
                min=floats.min().item(),
                max=floats.max().item(),
                num=num,
                sum=floats.sum().item(),
                sum_squares=sum_sq,
                bucket_limits=limits[1:].tolist(),
                bucket_counts=counts.tolist(),
            )

            ints = make_np(torch.randint(0, 100, (num,)))
            bins = [0, 25, 50, 75, 100]
            counts, limits = np.histogram(ints, bins)
            sum_sq = ints.dot(ints).item()
            w.add_histogram_raw(
                "int histogram raw",
                min=ints.min().item(),
                max=ints.max().item(),
                num=num,
                sum=ints.sum().item(),
                sum_squares=sum_sq,
                bucket_limits=limits[1:].tolist(),
                bucket_counts=counts.tolist(),
            )

            ints = torch.tensor(range(100)).float()
            nbins = 100
            counts = torch.histc(ints, bins=nbins, min=0, max=99)
            limits = torch.tensor(range(nbins))
            sum_sq = ints.dot(ints).item()
            w.add_histogram_raw(
                "int histogram raw",
                min=ints.min().item(),
                max=ints.max().item(),
                num=num,
                sum=ints.sum().item(),
                sum_squares=sum_sq,
                bucket_limits=limits.tolist(),
                bucket_counts=counts.tolist(),
            )