def make_symmetric_pd_matrices(*shape, device, dtype):
    if shape[-1] != shape[-2]:
        raise AssertionError(f"expected square matrix, got shape[-1]={shape[-1]} != shape[-2]={shape[-2]}")
    t = make_tensor(shape, device=device, dtype=dtype)
    i = torch.eye(shape[-1], device=device, dtype=dtype) * 1e-5
    return t @ t.mT + i