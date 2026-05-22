def is_torch_fp16_available_on_device(device: str) -> bool:
    if not is_torch_available():
        return False

    if is_torch_hpu_available():
        if is_habana_gaudi1():
            return False
        else:
            return True

    import torch

    try:
        x = torch.zeros(2, 2, dtype=torch.float16, device=device)
        _ = x @ x
        # At this moment, let's be strict of the check: check if `LayerNorm` is also supported on device, because many
        # models use this layer.
        batch, sentence_length, embedding_dim = 3, 4, 5
        embedding = torch.randn(batch, sentence_length, embedding_dim, dtype=torch.float16, device=device)
        layer_norm = torch.nn.LayerNorm(embedding_dim, dtype=torch.float16, device=device)
        _ = layer_norm(embedding)
        return True
    except Exception:
        return False