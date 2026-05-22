def dot(self, other):
    if self.is_complex():
        if self.is_conj():
            if other.is_conj():
                return torch.dot(self.conj(), other.conj()).conj()
            else:
                return torch.vdot(self.conj(), other)
        elif other.is_conj():
            return torch.vdot(other.conj(), self)

    return (self * other).sum()