def lstsq(a, b, rcond=None):
    a = convert_to_tensor(a)
    b = convert_to_tensor(b)
    # `torch.linalg.lstsq` defaults to the QR-based `gelsy` driver on CPU,
    # which interprets `rcond` differently from numpy and can return
    # large errors on near-rank-deficient inputs. Pin the SVD-based
    # `gelsd` driver when an explicit `rcond` is supplied so the result
    # matches numpy's `linalg.lstsq` semantics. CUDA only supports the
    # `gels` driver, which has no rcond cutoff, so fall back to SVD via
    # `pinv` to match numpy semantics there.
    if rcond is not None:
        if a.device.type == "cuda":
            return torch.linalg.pinv(a, rtol=rcond) @ b
        return torch.linalg.lstsq(a, b, rcond=rcond, driver="gelsd")[0]
    return torch.linalg.lstsq(a, b, rcond=rcond)[0]