    def _build_pe(w, h, dim=256, temp=10000.):
        assert dim % 4 == 0
        gw = torch.arange(w, dtype=torch.float32)
        gh = torch.arange(h, dtype=torch.float32)
        gw, gh = torch.meshgrid(gw, gh, indexing='ij')
        pdim  = dim // 4
        omega = 1. / (temp ** (torch.arange(pdim, dtype=torch.float32) / pdim))
        ow = gw.flatten()[:, None] @ omega[None]
        oh = gh.flatten()[:, None] @ omega[None]
        return torch.cat([ow.sin(), ow.cos(), oh.sin(), oh.cos()], 1)[None]