def _convert_lora_to_linear(module: LoraLayer, adapter_name: str = "default"):
    base_layer = module.get_base_layer()
    weight = base_layer.weight

    assert isinstance(weight, bnb.nn.Params4bit)
    quant_state = weight.quant_state
    original_dtype = quant_state.dtype

    w_dq = dequantize_4bit(weight.data, quant_state).float()
    lora_delta = (
        module.lora_B[adapter_name].weight
        @ module.lora_A[adapter_name].weight
        * module.scaling[adapter_name]
    )
    w_dq += lora_delta.float()
    w_dq = w_dq.to(original_dtype)

    new_module = torch.nn.Linear(
        w_dq.shape[1], w_dq.shape[0], bias = module.base_layer.bias is not None
    )
    new_module.weight.data = torch.nn.Parameter(w_dq, requires_grad = False)
    if module.lora_bias[adapter_name]:
        bias_data = module.base_layer.bias.data + module.lora_B[adapter_name].bias
        new_module.bias.data = torch.nn.Parameter(bias_data, requires_grad = False)
    return new_module