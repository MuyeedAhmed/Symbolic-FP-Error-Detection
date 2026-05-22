def vdot(self, other):
    if not self.is_complex():
        return torch.dot(self, other)

    if self.is_conj():
        if other.is_conj():
            return torch.vdot(other.conj(), self.conj())
        else:
            return torch.dot(self.conj(), other)
    elif other.is_conj():
        return torch.dot(self, other.conj()).conj()

    # The decomposition fails if you do self.conj()... not sure why
    return (self.conj_physical() * other).sum()