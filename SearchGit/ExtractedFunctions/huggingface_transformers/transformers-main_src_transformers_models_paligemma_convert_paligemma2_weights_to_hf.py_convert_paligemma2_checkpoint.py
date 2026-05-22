def convert_paligemma2_checkpoint(
    checkpoint_path,
    pytorch_dump_folder_path,
    variant: str,
    precision: str,
    do_convert_weights=False,
):
    """
    Read checkpoints from flax npz files, rename/reshape, send result to state dict and verify logits if needed.
    """
    config = get_paligemma2_config(variant, precision=precision)
    if do_convert_weights:
        tokenizer_id = "google/paligemma-3b-pt-224"  # same tokenizer as paligemma 1
        tokenizer = AutoTokenizer.from_pretrained(tokenizer_id)
        image_token = AddedToken("<image>", normalized=False, special=True)
        tokens_to_add = {"additional_special_tokens": [image_token]}
        tokenizer.add_special_tokens(tokens_to_add)

        # tokenizer.padding_side = 'right' # uncomment for testing purposes only.

        image_processor = SiglipImageProcessor.from_pretrained("google/paligemma-3b-pt-224")
        image_processor.size = {"width": config.vision_config.image_size, "height": config.vision_config.image_size}
        image_processor.image_seq_length = config.vision_config.num_image_tokens

        processor = PaliGemmaProcessor(image_processor=image_processor, tokenizer=tokenizer)
        data = jnp.load(checkpoint_path)
        state_dict = flatten_nested_dict(data, precision=precision)
        del data
        state_dict_transformers = slice_state_dict(state_dict, config)
        del state_dict
        del config.hidden_size  # this key is unused
        model = PaliGemmaForConditionalGeneration(config).to(device).eval()
        model.load_state_dict(state_dict_transformers)
        del state_dict_transformers
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
        # convert to needed precision

        model.to(DTYPES[precision])
        model.save_pretrained(pytorch_dump_folder_path)
        processor.save_pretrained(pytorch_dump_folder_path)

    else:
        processor = PaliGemmaProcessor.from_pretrained(pytorch_dump_folder_path, do_rescale=False)
        model = (
            PaliGemmaForConditionalGeneration.from_pretrained(pytorch_dump_folder_path, attn_implementation="sdpa")
            .to(device)
            .eval()
        )