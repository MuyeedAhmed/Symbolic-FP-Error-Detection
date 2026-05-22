    def single_mask_loss(self, gt_mask, pred, proto, xyxy, area):
        """Calculates and normalizes single mask loss for YOLOv5 between predicted and ground truth masks."""
        pred_mask = (pred @ proto.view(self.nm, -1)).view(-1, *proto.shape[1:])  # (n,32) @ (32,80,80) -> (n,80,80)
        loss = F.binary_cross_entropy_with_logits(pred_mask, gt_mask, reduction="none")
        return (crop_mask(loss, xyxy).mean(dim=(1, 2)) / area).mean()