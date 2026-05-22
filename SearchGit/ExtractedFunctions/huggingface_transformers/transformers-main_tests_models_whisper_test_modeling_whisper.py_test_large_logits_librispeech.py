    def test_large_logits_librispeech(self):
        set_seed(42)

        torch_device = "cpu"
        model = WhisperModel.from_pretrained("openai/whisper-large")
        model.to(torch_device)

        input_speech = self._load_datasamples(1)

        processor = WhisperProcessor.from_pretrained("openai/whisper-large")
        processed_inputs = processor(
            audio=input_speech,
            text="This part of the speech",
            add_special_tokens=False,
            return_tensors="pt",
            sampling_rate=16_000,
        )
        input_features = processed_inputs.input_features.to(torch_device)
        decoder_input_ids = processed_inputs.labels.to(torch_device)

        logits = model(
            input_features,
            decoder_input_ids=decoder_input_ids,
            output_hidden_states=False,
            output_attentions=False,
            use_cache=False,
        )

        logits = logits.last_hidden_state @ model.decoder.embed_tokens.weight.T

        # fmt: off
        EXPECTED_LOGITS = torch.tensor(
            [
                2.1382, 0.9381, 4.4671, 3.5589, 2.4022, 3.8576, -0.6521, 2.5472,
                1.8301, 1.9957, 2.3432, 1.4678, 0.5459, 2.2597, 1.5179, 2.5357,
                1.1624, 0.6194, 1.0757, 1.8259, 2.4076, 1.6601, 2.3503, 1.3376,
                1.9891, 1.8635, 3.8931, 5.3699, 4.4772, 3.9184
            ]
        )
        # fmt: on

        torch.testing.assert_close(logits[0, 0, :30].cpu(), EXPECTED_LOGITS, rtol=1e-4, atol=1e-4)