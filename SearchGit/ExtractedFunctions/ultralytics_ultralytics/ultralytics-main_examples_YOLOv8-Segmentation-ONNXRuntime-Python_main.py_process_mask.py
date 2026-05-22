    def process_mask(
        self, protos: torch.Tensor, masks_in: torch.Tensor, bboxes: torch.Tensor, shape: tuple[int, int]
    ) -> torch.Tensor:
        """Process prototype masks with predicted mask coefficients to generate instance segmentation masks.

        Args:
            protos (torch.Tensor): Prototype masks with shape (mask_dim, mask_h, mask_w).
            masks_in (torch.Tensor): Predicted mask coefficients with shape (N, mask_dim), where N is number of
                detections.
            bboxes (torch.Tensor): Bounding boxes with shape (N, 4), where N is number of detections.
            shape (tuple[int, int]): The size of the input image as (height, width).

        Returns:
            (torch.Tensor): Binary segmentation masks with shape (N, height, width).
        """
        c, mh, mw = protos.shape  # CHW
        masks = (masks_in @ protos.float().view(c, -1)).view(-1, mh, mw)  # Matrix multiplication
        masks = ops.scale_masks(masks[None], shape)[0]  # Scale masks to original image size
        masks = ops.crop_mask(masks, bboxes)  # Crop masks to bounding boxes
        return masks.gt_(0.0)  # Convert to binary masks