def math_attention(
    query: torch.Tensor,
    key: torch.Tensor,
    value: torch.Tensor,
    score_mod: Callable,
    block_mask: tuple,
    scale: float,
    kernel_options: dict[str, Any],
    score_mod_other_buffers: tuple = (),
    mask_mod_other_buffers: tuple = (),
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Eager implementation

    This implementation uses vmap to vectorize the score_mod function over the batch, head, m, and n dimensions.
    We then apply the vectorized score_mod function to the scores matrix. Each wrap of vmap applies one of the
    batch, head, m, or n dimensions. We need to apply vmap 4 times to vectorized over all 4 dimensions.

    Args:
        query: The query tensor
        key: The key tensor
        value: The value tensor
        score_mod: The score_mod function
        other_buffers: Other buffers that are passed to the score_mod function

    Notes:
        Query and Keys are dtype cast up to float64 (if query.dtype is float64) and float32 otherwise.
        Scores and Values are dtype cast to input query.dtype at the end.
    """
    # broadcast query & key along head dim for GQA
    G = query.size(1) // key.size(1)
    value = torch.repeat_interleave(value, G, dim=1)
    key = torch.repeat_interleave(key, G, dim=1)

    Bq, Bkv = query.size(0), key.size(0)
    if not ((Bq == Bkv) or (Bq > 1 and Bkv == 1)):
        raise RuntimeError(f"Bq and Bkv must broadcast. Got Bq={Bq} and Bkv={Bkv}")

    key = key.expand((Bq, *key.size()[1:]))
    value = value.expand((Bq, *value.size()[1:]))

    _, post_mod_scores = _math_attention_inner(
        query,
        key,
        value,
        score_mod,
        block_mask,
        scale,
        kernel_options,
        score_mod_other_buffers,
        mask_mod_other_buffers,
    )

    # Set fully masked rows' sumexp to 0.0
    logsumexp = post_mod_scores.logsumexp(dim=-1)
    masked_rows = torch.all(post_mod_scores == -float("inf"), dim=-1)
    logsumexp = torch.where(masked_rows, -float("inf"), logsumexp)

    # working precision will be used so no need to cast to fp32
    max_scores = torch.max(post_mod_scores, dim=-1)[0]

    post_mod_scores = torch._safe_softmax(post_mod_scores, dim=-1)

    # NB: kernel computes in ln2 space, we always convert back at the top level op, so
    # for math impl we divide by log(2) because we will multiply by log(2)

    return (
        post_mod_scores.to(query.dtype) @ value.to(query.dtype),
        logsumexp / math.log(2),
        max_scores / math.log(2),
    )