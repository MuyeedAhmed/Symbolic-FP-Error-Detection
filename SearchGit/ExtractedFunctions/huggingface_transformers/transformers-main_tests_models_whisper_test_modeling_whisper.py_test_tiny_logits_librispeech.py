    def test_tiny_logits_librispeech(self):
        torch_device = "cpu"
        set_seed(42)
        model = WhisperModel.from_pretrained("openai/whisper-tiny")
        model.to(torch_device)
        input_speech = self._load_datasamples(1)
        feature_extractor = WhisperFeatureExtractor()
        input_features = feature_extractor(input_speech, return_tensors="pt", sampling_rate=16_000).input_features

        with torch.no_grad():
            logits = model(
                input_features,
                decoder_input_ids=torch.tensor([[50258, 50259, 50359]]),
                output_hidden_states=False,
                output_attentions=False,
                return_dict=False,
                use_cache=False,
            )

        # fmt: off
        EXPECTED_LOGITS = torch.tensor(
            [
                2.9892, -6.7607, 5.7348, 3.6096, 0.2152, -5.7321, 4.8855, -1.6407,
                0.2823, -1.5718, 10.4269, 3.4427, 0.0219, -8.0612, 3.4784, 8.4246,
                4.0575, -2.2864, 11.1084, 0.9963, 0.9884, -8.5154, -3.5469, -9.3713,
                0.9786, 3.5435, 7.4850, -5.2579, -1.4366, 10.4841
            ]
        )
        # fmt: on
        torch.testing.assert_close(logits[0][0, 0, :30].cpu(), EXPECTED_LOGITS, rtol=1e-4, atol=1e-4)

        # fmt: off
        EXPECTED_GENERATION = torch.tensor(
            [
                -1.4651, -2.6944, 2.7821, 2.3793, 4.0738, 0.0188, -3.3203, 1.9836,
                0.0520, 0.7095, 1.1063, 0.2952, -3.6786, -0.5249, 0.3105, 4.7691,
                1.1562, 1.3046, 0.5810, -0.3624, 1.7006, 1.3424, 0.9817, 2.1958,
                1.8775, -5.7046, -0.7679, 4.0113, 2.6848, 2.8609
            ]
        )
        # fmt: on

        head_logits = logits[0] @ model.decoder.embed_tokens.weight.T
        torch.testing.assert_close(head_logits[0, 0, :30].cpu(), EXPECTED_GENERATION, rtol=1e-4, atol=1e-4)