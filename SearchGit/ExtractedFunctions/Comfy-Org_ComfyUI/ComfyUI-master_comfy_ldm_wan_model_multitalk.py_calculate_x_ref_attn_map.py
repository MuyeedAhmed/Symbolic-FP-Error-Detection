def calculate_x_ref_attn_map(visual_q, ref_k, ref_target_masks, split_num=8):
    scale = 1.0 / visual_q.shape[-1] ** 0.5
    visual_q = visual_q.transpose(1, 2) * scale

    B, H, x_seqlens, K = visual_q.shape

    x_ref_attn_maps = []
    for class_idx, ref_target_mask in enumerate(ref_target_masks):
        ref_target_mask = ref_target_mask.view(1, 1, 1, -1)

        x_ref_attnmap = torch.zeros(B, H, x_seqlens, device=visual_q.device, dtype=visual_q.dtype)
        chunk_size = min(max(x_seqlens // split_num, 1), x_seqlens)

        for i in range(0, x_seqlens, chunk_size):
            end_i = min(i + chunk_size, x_seqlens)

            attn_chunk = visual_q[:, :, i:end_i] @ ref_k.permute(0, 2, 3, 1)  # B, H, chunk, ref_seqlens

            # Apply softmax
            attn_max = attn_chunk.max(dim=-1, keepdim=True).values
            attn_chunk = (attn_chunk - attn_max).exp()
            attn_sum = attn_chunk.sum(dim=-1, keepdim=True)
            attn_chunk = attn_chunk / (attn_sum + 1e-8)

            # Apply mask and sum
            masked_attn = attn_chunk * ref_target_mask
            x_ref_attnmap[:, :, i:end_i] = masked_attn.sum(-1) / (ref_target_mask.sum() + 1e-8)

            del attn_chunk, masked_attn

        # Average across heads
        x_ref_attnmap = x_ref_attnmap.mean(dim=1)  # B, x_seqlens
        x_ref_attn_maps.append(x_ref_attnmap)

    del visual_q, ref_k

    return torch.cat(x_ref_attn_maps, dim=0)