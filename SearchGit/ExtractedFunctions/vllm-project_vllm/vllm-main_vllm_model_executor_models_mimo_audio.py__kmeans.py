def _kmeans(samples, num_clusters: int, num_iters: int = 10):
    dim, dtype = samples.shape[-1], samples.dtype

    means = _sample_vectors(samples, num_clusters)

    for _ in range(num_iters):
        dists = -(
            samples.pow(2).sum(1, keepdim=True)
            - 2 * samples @ means.t()
            + means.t().pow(2).sum(0, keepdim=True)
        )

        buckets = dists.max(dim=-1).indices
        bins = torch.bincount(buckets, minlength=num_clusters)

        new_means = buckets.new_zeros(num_clusters, dim, dtype=dtype)
        new_means = new_means.scatter_add_(
            0, repeat(buckets, "n -> n d", d=dim), samples
        )

        if dist.is_initialized():
            dist.all_reduce(bins, op=dist.ReduceOp.SUM)
            dist.all_reduce(new_means, op=dist.ReduceOp.SUM)

        zero_mask = bins == 0
        bins_min_clamped = bins.masked_fill(zero_mask, 1)

        new_means = new_means / bins_min_clamped[..., None]

        means = torch.where(zero_mask[..., None], means, new_means)

    return means, bins