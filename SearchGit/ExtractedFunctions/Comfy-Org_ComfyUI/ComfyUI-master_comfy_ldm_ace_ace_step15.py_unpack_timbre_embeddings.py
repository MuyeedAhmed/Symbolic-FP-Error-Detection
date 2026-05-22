    def unpack_timbre_embeddings(self, timbre_embs_packed, refer_audio_order_mask):
        N, d = timbre_embs_packed.shape
        device = timbre_embs_packed.device
        B = N
        counts = torch.bincount(refer_audio_order_mask, minlength=B)
        max_count = counts.max().item()

        sorted_indices = torch.argsort(
            refer_audio_order_mask * N + torch.arange(N, device=device),
            stable=True
        )
        sorted_batch_ids = refer_audio_order_mask[sorted_indices]

        positions = torch.arange(N, device=device)
        batch_starts = torch.cat([torch.tensor([0], device=device), torch.cumsum(counts, dim=0)[:-1]])
        positions_in_sorted = positions - batch_starts[sorted_batch_ids]

        inverse_indices = torch.empty_like(sorted_indices)
        inverse_indices[sorted_indices] = torch.arange(N, device=device)
        positions_in_batch = positions_in_sorted[inverse_indices]

        indices_2d = refer_audio_order_mask * max_count + positions_in_batch
        one_hot = F.one_hot(indices_2d, num_classes=B * max_count).to(timbre_embs_packed.dtype)

        timbre_embs_flat = one_hot.t() @ timbre_embs_packed
        timbre_embs_unpack = timbre_embs_flat.view(B, max_count, d)

        mask_flat = (one_hot.sum(dim=0) > 0).long()
        new_mask = mask_flat.view(B, max_count)
        return timbre_embs_unpack, new_mask