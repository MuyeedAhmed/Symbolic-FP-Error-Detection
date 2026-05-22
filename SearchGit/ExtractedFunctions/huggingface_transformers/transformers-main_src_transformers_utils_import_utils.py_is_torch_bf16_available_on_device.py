def is_torch_bf16_available_on_device(device: str) -> bool:
    if not is_torch_available():
        return False

    import torch

    if device == "cuda":
        return is_torch_bf16_gpu_available()

    if device == "hpu":
        return True

    try:
        x = torch.zeros(2, 2, dtype=torch.bfloat16, device=device)
        _ = x @ x
        return True
    except Exception:
        return False