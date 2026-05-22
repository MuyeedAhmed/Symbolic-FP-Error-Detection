def convert_paligemma_checkpoint(
    checkpoint_path,
    tokenizer_model_file,
    pytorch_dump_folder_path,
    variant: str,
    precision: str,
    do_convert_weights=False,
):
    """
    Read checkpoints from flax npz files, rename/reshape, send result to state dict and verify logits if needed.
    """
    config = get_paligemma_config(variant, precision=precision)
    if do_convert_weights:
        if variant == "2b-test":
            # for the test model, the vocabulary was smaller
            tokenizer_id = "google/gemma-2b"
            tokenizer = AutoTokenizer.from_pretrained(tokenizer_id)
        else:
            tokenizer_class = GemmaTokenizer if GemmaTokenizerFast is None else GemmaTokenizerFast
            tokenizer = tokenizer_class(tokenizer_model_file)
        image_token = AddedToken("<image>", normalized=False, special=True)
        tokens_to_add = {"additional_special_tokens": [image_token]}
        tokenizer.add_special_tokens(tokens_to_add)

        # tokenizer.padding_side = 'right' # uncomment for testing purposes only.

        image_processor = SiglipImageProcessor.from_pretrained("google/siglip-so400m-patch14-384")
        image_processor.size = {"width": config.vision_config.image_size, "height": config.vision_config.image_size}
        image_processor.image_seq_length = config.vision_config.num_image_tokens

        processor = PaliGemmaProcessor(image_processor=image_processor, tokenizer=tokenizer)
        data = load(checkpoint_path)
        state_dict = flatten_nested_dict(data)
        del data
        state_dict_transformers = slice_state_dict(state_dict, config)
        del state_dict

        model = PaliGemmaForConditionalGeneration(config).to(device).eval()
        model.load_state_dict(state_dict_transformers)
        del state_dict_transformers

    else:
        processor = PaliGemmaProcessor.from_pretrained(pytorch_dump_folder_path)
        model = (
            PaliGemmaForConditionalGeneration.from_pretrained(pytorch_dump_folder_path, attn_implementation="sdpa")
            .to(device)
            .eval()
        )
    model.config.text_config._attn_implementation = "sdpa"

    # model expansion to get random embeds of image tokens
    pad_shape = 64  # for performance reasons
    pre_expansion_embeddings = model.language_model.model.embed_tokens.weight.data
    mu = torch.mean(pre_expansion_embeddings, dim=0).float()
    n = pre_expansion_embeddings.size()[0]
    sigma = ((pre_expansion_embeddings - mu).T @ (pre_expansion_embeddings - mu)) / n
    dist = torch.distributions.multivariate_normal.MultivariateNormal(mu, covariance_matrix=1e-5 * sigma)

    # We add an image token so we resize the model
    model.resize_token_embeddings(config.text_config.vocab_size + 2, pad_shape)
    model.language_model.model.embed_tokens.weight.data[257152:] = torch.stack(
        tuple(dist.sample() for _ in range(model.language_model.model.embed_tokens.weight.data[257152:].shape[0])),
        dim=0,
    )
    model.language_model.lm_head.weight.data[257152:] = torch.stack(
        tuple(dist.sample() for _ in range(model.language_model.lm_head.weight.data[257152:].shape[0])),
        dim=0,
    )

    model.save_pretrained(pytorch_dump_folder_path, max_shard_size="2GB")
    processor.save_pretrained(pytorch_dump_folder_path)