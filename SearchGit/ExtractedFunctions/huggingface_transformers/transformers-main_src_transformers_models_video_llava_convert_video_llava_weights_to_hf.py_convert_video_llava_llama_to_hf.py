def convert_video_llava_llama_to_hf(text_model_id, vision_model_id, output_hub_path, old_state_dict_id):
    torch.set_default_dtype(torch.float16)
    text_config = AutoConfig.from_pretrained(text_model_id)

    tokenizer = AutoTokenizer.from_pretrained(text_model_id)
    tokenizer.add_tokens(AddedToken("<image>", special=True, normalized=False), special_tokens=True)
    tokenizer.add_tokens(AddedToken("<video>", special=True, normalized=False), special_tokens=True)
    tokenizer.add_special_tokens({"pad_token": "<pad>"})
    tokenizer.padding_side = "left"

    image_processor = VideoLlavaImageProcessor.from_pretrained(vision_model_id)

    processor = VideoLlavaProcessor(tokenizer=tokenizer, image_processor=image_processor)

    config = VideoLlavaConfig(text_config=text_config)
    config.pad_token_id = 32002

    with torch.device("meta"):
        model = VideoLlavaForConditionalGeneration(config)

    model_state_dict = set(model.state_dict().keys())

    # Pad to 64 for performance reasons
    pad_shape = 64
    state_dict_temp = "pytorch_model-0000{i}-of-00002.bin"
    for shard in range(1, 3):
        state_dict_path = hf_hub_download(old_state_dict_id, state_dict_temp.format(i=shard))
        state_dict = torch.load(state_dict_path, map_location="cpu", weights_only=True)
        state_dict = convert_state_dict_to_hf(state_dict)
        model.load_state_dict(state_dict, strict=False, assign=True)
        model_state_dict -= set(state_dict.keys())

    if len(model_state_dict) > 0:
        raise RuntimeError(f"Missing keys in state dict: {model_state_dict}")

    pre_expansion_embeddings = model.language_model.model.embed_tokens.weight.data
    mu = torch.mean(pre_expansion_embeddings, dim=0).float()
    n = pre_expansion_embeddings.size()[0]
    sigma = ((pre_expansion_embeddings - mu).T @ (pre_expansion_embeddings - mu)) / n
    dist = torch.distributions.multivariate_normal.MultivariateNormal(mu, covariance_matrix=1e-5 * sigma)

    # We add an image and video token so we resize the model
    model.resize_token_embeddings(config.text_config.vocab_size + 3, pad_shape)
    model.language_model.model.embed_tokens.weight.data[32000:] = torch.stack(
        tuple(dist.sample() for _ in range(model.language_model.model.embed_tokens.weight.data[32000:].shape[0])),
        dim=0,
    )
    model.language_model.lm_head.weight.data[32000:] = torch.stack(
        tuple(dist.sample() for _ in range(model.language_model.lm_head.weight.data[32000:].shape[0])),
        dim=0,
    )

    model.push_to_hub(output_hub_path)
    processor.push_to_hub(output_hub_path)