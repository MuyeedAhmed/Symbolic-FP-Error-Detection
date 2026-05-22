    def _inverse(self, y):
        if not (-y.dim() <= self.dim < y.dim()):
            raise AssertionError(
                f"dim {self.dim} out of range for tensor with {y.dim()} dimensions"
            )
        if y.size(self.dim) != len(self.transforms):
            raise AssertionError(
                f"y.size({self.dim}) = {y.size(self.dim)} must equal len(transforms) {len(self.transforms)}"
            )
        xslices = []
        for yslice, trans in zip(self._slice(y), self.transforms):
            xslices.append(trans.inv(yslice))
        return torch.stack(xslices, dim=self.dim)