        def _merge_weights(
            weights: Iterable[tuple[str, torch.Tensor]],
        ) -> Iterable[tuple[str, torch.Tensor]]:
            for name, tensor in weights:
                clean_name = name
                if clean_name.startswith("model."):
                    clean_name = clean_name[len("model.") :]

                if clean_name in lora_pairs:
                    pair = lora_pairs[clean_name]
                    if "A" in pair and "B" in pair:
                        lora_A = pair["A"].to(device=tensor.device, dtype=tensor.dtype)
                        lora_B = pair["B"].to(device=tensor.device, dtype=tensor.dtype)
                        tensor = tensor + (lora_B @ lora_A) * scaling
                yield name, tensor