    def _clip_inference(self, images, texts):
        """Perform CLIP inference to calculate similarity between images and text prompts.

        Args:
            images (list[PIL.Image]): List of source images, each should be PIL.Image with RGB channel order.
            texts (list[str]): List of prompt texts, each should be a string object.

        Returns:
            (torch.Tensor): Similarity matrix between given images and texts with shape (M, N).
        """
        from ultralytics.nn.text_model import CLIP

        if not hasattr(self, "clip"):
            self.clip = CLIP("ViT-B/32", device=self.device)
        images = torch.stack([self.clip.image_preprocess(image).to(self.device) for image in images])
        image_features = self.clip.encode_image(images)
        text_features = self.clip.encode_text(self.clip.tokenize(texts))
        return text_features @ image_features.T  # (M, N)