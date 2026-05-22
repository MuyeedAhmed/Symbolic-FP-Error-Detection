def use_torch(
    hidden_states,
    token_lora_mapping,
    topk_ids,
    lora_a_stacked,
    lora_b_stacked,
    top_k_num,
    num_slices=1,
):
    outputs = []
    for i in range(hidden_states.shape[0]):
        slice_tensors = []
        for slice_id in range(num_slices):
            lora_idx = token_lora_mapping[i]
            expert_ids = topk_ids[i]
            lora_a = lora_a_stacked[slice_id][lora_idx][expert_ids]
            lora_b = lora_b_stacked[slice_id][lora_idx][expert_ids]
            tensors = [
                hidden_states[i] @ lora_a[x].T @ lora_b[x].T for x in range(top_k_num)
            ]
            slice_tensors.append(torch.stack(tensors, dim=0))

        outputs.append(torch.concat(slice_tensors, dim=-1))
    return torch.stack(outputs, dim=0)