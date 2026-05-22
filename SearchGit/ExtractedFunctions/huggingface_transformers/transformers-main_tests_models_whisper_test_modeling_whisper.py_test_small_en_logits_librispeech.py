    def test_small_en_logits_librispeech(self):
        set_seed(42)
        torch_device = "cpu"
        model = WhisperModel.from_pretrained("openai/whisper-small.en")
        model.to(torch_device)

        input_speech = self._load_datasamples(1)

        feature_extractor = WhisperFeatureExtractor()
        input_features = feature_extractor(input_speech, return_tensors="pt").input_features.to(torch_device)

        logits = model(
            input_features,
            decoder_input_ids=torch.tensor([[model.config.decoder_start_token_id]]),
            output_hidden_states=False,
            output_attentions=False,
            use_cache=False,
        )

        logits = logits.last_hidden_state @ model.decoder.embed_tokens.weight.T

        # fmt: off
        EXPECTED_LOGITS = torch.tensor(
            [
                -3.6784, -7.7211, -9.5070, -11.9286, -7.6489, -9.7026, -5.6188,
                -8.0104, -4.6238, -5.1833, -9.0485, -3.4079, -5.4874, -2.6935,
                -6.3479, -7.3398, -6.9558, -7.6867, -7.4748, -8.3463, -9.9781,
                -10.8389, -10.3105, -11.7201, -9.7261, -7.1590, -5.9272, -12.4509,
                -11.1146, -8.1918
            ]
        )
        # fmt: on
        torch.testing.assert_close(logits[0, 0, :30].cpu(), EXPECTED_LOGITS, rtol=1e-4, atol=1e-4)