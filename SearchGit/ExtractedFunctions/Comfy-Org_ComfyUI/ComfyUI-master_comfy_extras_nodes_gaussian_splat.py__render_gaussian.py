def _render_gaussian(xyz, rgb, opacity, scale, rot, width, height, splat_scale, bg, camera_info,
                     sharpen=1.0, headlight_shading=0.0, render_style="color"):
    # Perspective-correct anisotropic gaussian splat rasterizer. Each splat is weighted by its 3D Gaussian's
    # peak along each pixel's ray (AAA / Hahlbohm), composited front-to-back across depth slabs. `render_style`
    # selects the image: color / clay / depth / normal. Returns (image HxWx3, coverage mask HxW) on CPU.
    dev = comfy.model_management.get_torch_device()
    t = lambda a: torch.as_tensor(a, dtype=torch.float32, device=dev)
    idev, idtype = comfy.model_management.intermediate_device(), comfy.model_management.intermediate_dtype()
    xyz, rgb, opacity = t(xyz), t(rgb).clamp(0, 1), t(opacity).reshape(-1)
    scale, rot = t(scale) * float(splat_scale), t(rot)
    do_linear = render_style == "color"  # colour blends in linear light, re-encoded at the end
    if do_linear:
        rgb = _srgb_to_linear(rgb)
    flat = width * height
    bg_t = t(bg)
    bg_comp = _srgb_to_linear(bg_t) if do_linear else bg_t  # background blended in the same space as the splats
    need_depth = render_style == "depth"
    need_normal = render_style in ("normal", "clay") or headlight_shading > 0

    def background_only():  # no splats to rasterize -> just the background + empty mask
        img = bg_t.expand(height, width, 3) if render_style == "color" else torch.zeros(height, width, 3, device=dev)
        return img.to(idev, idtype), torch.zeros(height, width, device=idev, dtype=idtype)

    if xyz.shape[0] == 0:  # empty input (e.g. all culled by opacity_threshold)
        return background_only()

    eye, target, right, up, fwd = _camera_basis(camera_info, dev)  # all camera state comes from camera_info
    W = torch.stack([right, up, fwd], 0)                           # rows = camera axes (world -> camera)
    cam = (xyz - eye) @ W.T
    fov = float(camera_info.get("fov", 0) or 0) or 35.0
    zoom = float(camera_info.get("zoom", 1.0) or 1.0)              # three.js digital zoom: scales the focal length
    is_ortho = str(camera_info.get("cameraType", "")).lower().startswith("ortho")
    xc, yc, zc = cam.unbind(-1)

    keep = zc > 1e-2
    xc, yc, zc, rgb, opacity, scale, rot = (a[keep] for a in (xc, yc, zc, rgb, opacity, scale, rot))
    if xc.shape[0] == 0:  # nothing in front of the camera -> background only
        return background_only()
    if render_style == "clay":
        rgb = torch.full_like(rgb, 0.75)  # neutral albedo -> shading shows pure geometry

    f = (min(width, height) / 2) / math.tan(math.radians(fov) / 2) * zoom  # fov over the smaller axis, x camera zoom
    cx0, cy0 = width / 2, height / 2

    # Camera-space 3D covariance per splat: Sigma = (W Rq) diag(scale^2) (W Rq)^T, plus a tiny relative
    # regularizer for a stable inverse (a pixel-size Mip low-pass would over-thicken flat surfels and blur).
    Mw = W[None] @ _quat_to_mat(rot)  # (N,3,3) world -> camera
    cam_cov = (Mw * scale.square()[:, None, :]) @ Mw.transpose(1, 2)
    cam_cov = cam_cov + (cam_cov.diagonal(dim1=-2, dim2=-1).mean(-1) * 1e-3)[:, None, None] * torch.eye(3, device=dev)

    # Perspective-correct weighting: peak of the 3D Gaussian along each pixel ray. Precompute Si, Si@mu, mu^T Si mu.
    mu = torch.stack([xc, yc, zc], -1)
    si = torch.linalg.inv(cam_cov)
    simu = (si @ mu[:, :, None])[:, :, 0]  # (N,3)
    musimu = (mu * simu).sum(-1)           # (N,)
    s00, s01, s02 = si[:, 0, 0], si[:, 0, 1], si[:, 0, 2]
    s11, s12, s22 = si[:, 1, 1], si[:, 1, 2], si[:, 2, 2]
    simu0, simu1, simu2 = simu.unbind(-1)
    if need_normal:  # surfel normal = thinnest axis, oriented toward camera
        nrm = Mw[torch.arange(Mw.shape[0], device=dev), :, scale.argmin(-1)]  # (N,3) camera-space normal
        nrm = nrm * torch.where(nrm[:, 2:3] > 0, -1.0, 1.0)                   # flip so nz <= 0 (faces camera)

    # Screen centre (exact) + footprint radius from the affine 2D projection (used only to size the kernel).
    # The image is +y-down, so the projection's y row is unflipped - it matches the splat frame's +Y.
    jm = torch.zeros(xc.shape[0], 2, 3, device=dev)
    if is_ortho:                                              # parallel projection: screen = s * (xc, yc)
        s = f / float((target - eye).norm().clamp_min(1e-6))  # pixels per world unit at the target plane
        cx, cy = cx0 + s * xc, cy0 + s * yc
        jm[:, 0, 0] = s
        jm[:, 1, 1] = s
    else:  # perspective: screen = f * (xc, yc) / zc
        invz = 1.0 / zc
        cx, cy = cx0 + f * xc * invz, cy0 + f * yc * invz
        jm[:, 0, 0], jm[:, 0, 2] = f * invz, -f * xc * invz.square()
        jm[:, 1, 1], jm[:, 1, 2] = f * invz, -f * yc * invz.square()
    cov2 = jm @ cam_cov @ jm.transpose(1, 2)
    a, b, c = cov2[:, 0, 0], cov2[:, 0, 1], cov2[:, 1, 1]
    max_eig = (a + c) * 0.5 + (((a - c) * 0.5).square() + b * b).clamp_min(0).sqrt()
    radius = 3.0 * max_eig.clamp_min(1e-8).sqrt()
    K = int(min(max(24, min(width, height) // 16), max(2, math.ceil(_quantile(radius, 0.995).item()))))

    # Per-splat kernel size: bucket splats by radius into a coarse ladder of window sizes (global K stays the cap) so
    # small splats (the bulk of it) use a small window.
    levels = [L for L in (16, 64, 256) if L < K] + [K]
    levels_t = torch.tensor(levels, device=dev, dtype=torch.float32)
    grids = []
    for L in levels:
        rng = torch.arange(-L, L + 1, device=dev, dtype=torch.float32)
        gy, gx = torch.meshgrid(rng, rng, indexing="ij")
        grids.append((gx.reshape(-1), gy.reshape(-1)))
    blevel = torch.bucketize(radius * (4.0 / 3.0), levels_t).clamp_(max=len(levels) - 1)  # window >= ~4 sigma

    n = zc.shape[0]
    ns = int(min(256, max(1, n // 1000)))                      # depth slabs: 1 per ~1000 splats, capped
    nl = len(levels)
    order = torch.argsort(zc)                                  # front (small zc) -> back -> defines the slabs
    bounds = torch.linspace(0, n, ns + 1, device=dev).round().long()
    rank = torch.empty(n, dtype=torch.long, device=dev)
    rank[order] = torch.arange(n, device=dev)                  # depth rank of each splat
    slab_id = (torch.searchsorted(bounds, rank, right=True) - 1).clamp_(0, ns - 1)
    key = slab_id * nl + blevel                                # group by slab, then kernel level (order-free within)
    order = torch.argsort(key)
    key = key[order]

    cxr, cyr = cx[order].round(), cy[order].round()
    s00, s01, s02 = s00[order], s01[order], s02[order]
    s11, s12, s22 = s11[order], s12[order], s22[order]
    s01b, s02b, s12b = s01 * 2, s02 * 2, s12 * 2               # doubled cross terms for the fused quadratic forms
    simu0, simu1, simu2, musimu = simu0[order], simu1[order], simu2[order], musimu[order]
    opacity, rgb = opacity[order], rgb[order]
    zc_o = zc[order] if need_depth else None
    nrm_o = nrm[order] if need_normal else None
    mux_o, muy_o, muz_o = (xc[order], yc[order], zc[order]) if is_ortho else (None, None, None)

    # Pack the per-splat scalars into one tensor so each chunk slices once
    common = [cxr, cyr, s00, s11, s22, s01b, s02b, s12b, opacity]
    pstack = torch.stack(common + ([s02, s12, mux_o, muy_o, muz_o] if is_ortho else [simu0, simu1, simu2, musimu]))

    # Precompute the (slab, level) run table on-GPU and pull it to the CPU once
    starts = torch.cat([torch.zeros(1, dtype=torch.long, device=dev), (key[1:] != key[:-1]).nonzero().flatten() + 1])
    ks = key[starts]
    run_lo = starts.tolist() + [n]
    run_lev = (ks % nl).tolist()
    run_slab = torch.div(ks, nl, rounding_mode="floor").tolist()
    slab_runs = [[] for _ in range(ns)]
    for r in range(len(run_lev)):
        slab_runs[run_slab[r]].append((run_lo[r], run_lo[r + 1], run_lev[r]))

    def splat(lo, hi, ox, oy):  # -> pixel idx (m,M), alpha (m,M); weight = 3D Gaussian peak along each pixel's ray
        cols = pstack[:, lo:hi, None].unbind(0)
        cxr_, cyr_, a00, a11, a22, b01, b02, b12, opa = cols[:9]   # a* = Si components; b* = 2 * cross terms
        px = cxr_ + ox[None, :]
        py = cyr_ + oy[None, :]
        valid = (px >= 0) & (px < width) & (py >= 0) & (py < height)
        if is_ortho:  # parallel ray (0,0,1) from screen point (X, Y, 0); rz constant per splat
            c02, c12, mx, my, mz = cols[9:]
            rx = (px - cx0) / s - mx
            ry = (py - cy0) / s - my
            rz = -mz
            a22rz = a22 * rz
            inx = torch.addcmul(b02 * rz, a00, rx).addcmul_(b01, ry)          # a00 rx + b01 ry + b02 rz
            rSr = torch.addcmul(a22rz * rz, rx, inx).addcmul_(ry, torch.addcmul(b12 * rz, a11, ry))
            dsr = torch.addcmul(a22rz, c02, rx).addcmul_(c12, ry)
            q = torch.addcdiv(rSr, dsr * dsr, a22.clamp_min(1e-12), value=-1).clamp_min_(0)
        else:  # perspective ray (dx,dy,1) through the camera origin
            su0, su1, su2, mus = cols[9:]
            dx, dy = (px - cx0) / f, (py - cy0) / f
            dsid = torch.addcmul(a22, dx, torch.addcmul(b02, a00, dx))        # a22 + dx*(a00 dx + b02)
            dsid = dsid.addcmul_(dy, torch.addcmul(b12, a11, dy))             # + dy*(a11 dy + b12)
            dsid = dsid.addcmul_(b01 * dx, dy)                               # + (2 s01) dx dy
            dsimu = torch.addcmul(su2, dx, su0).addcmul_(dy, su1)
            q = torch.addcdiv(mus, dsimu * dsimu, dsid.clamp_min(1e-12), value=-1).clamp_min_(0)
        alpha = (opa * torch.exp(-0.5 * q) * valid).clamp_(0, 0.999)
        idx = py.long().clamp(0, height - 1) * width + px.long().clamp(0, width - 1)
        return idx, alpha

    # Front-to-back compositing over the depth slabs set up above. Within a slab the accumulation is a pure
    # sum (order-independent), so splats are grouped by kernel level and each level uses its own tight window.
    sharp = sharpen != 1.0  # winner-take-more colour blend: dominant splat shows more
    cacc = torch.zeros((flat, 3), device=dev)
    trans = torch.ones((flat,), device=dev)
    a_buf = torch.zeros((flat,), device=dev)                            # sum alpha -> colour/depth/normal weight (alpha-weighted mean)
    tau_buf = torch.zeros((flat,), device=dev)                          # sum -ln(1-alpha) -> slab opacity = 1-prod(1-alpha)
    crgb = torch.zeros((flat, 3), device=dev)                           # sum alpha^p * rgb -> slab colour
    wbuf = torch.zeros((flat,), device=dev) if sharp else None          # sum alpha^p -> colour normalizer (sharp only)
    dacc = torch.zeros((flat,), device=dev) if need_depth else None     # front-weighted depth
    nacc = torch.zeros((flat, 3), device=dev) if need_normal else None  # front-weighted camera-space normal
    zslab = torch.zeros((flat,), device=dev) if need_depth else None
    nslab = torch.zeros((flat, 3), device=dev) if need_normal else None
    stale = 0  # consecutive fully-occluded slabs -> early-out
    for si in range(ns):
        runs = slab_runs[si]
        if not runs:
            continue
        a_buf.zero_()
        tau_buf.zero_()
        crgb.zero_()
        if sharp:
            wbuf.zero_()
        if need_depth:
            zslab.zero_()
        if need_normal:
            nslab.zero_()
        for r_lo, r_hi, li in runs:            # contiguous same-kernel-level runs in this slab
            ox, oy = grids[li]
            ch = max(2048, 10_000_000 // ox.shape[0])          # splats/chunk, bounded by this level's kernel size
            for lo in range(r_lo, r_hi, ch):
                hi = min(lo + ch, r_hi)
                idx, alpha = splat(lo, hi, ox, oy)
                idx, af = idx.reshape(-1), alpha.reshape(-1)
                a_buf.index_add_(0, idx, af)
                tau_buf.index_add_(0, idx, (-torch.log1p(-alpha)).reshape(-1))      # -ln(1-alpha), correct opacity merge
                apw = alpha.pow(sharpen) if sharp else alpha                        # bias colour toward the highest-alpha splat
                crgb.index_add_(0, idx, (apw[:, :, None] * rgb[lo:hi, None, :]).reshape(-1, 3))
                if sharp:
                    wbuf.index_add_(0, idx, apw.reshape(-1))
                if need_depth:
                    zslab.index_add_(0, idx, (alpha * zc_o[lo:hi, None]).reshape(-1))
                if need_normal:
                    nslab.index_add_(0, idx, (alpha[:, :, None] * nrm_o[lo:hi, None, :]).reshape(-1, 3))
        slab_a = 1 - torch.exp(-tau_buf)   # 1 - prod(1-alpha): true opacity of the slab's splats
        front = trans * slab_a
        denom = wbuf if sharp else a_buf
        cacc.addcmul_(front[:, None], crgb / denom.clamp_min(1e-8)[:, None])  # cacc += front * (crgb/denom)
        if need_depth or need_normal:
            ainv = a_buf.clamp_min(1e-8)   # alpha-weighted-mean normalizer (depth/normal only)
            if need_depth:
                dacc.addcmul_(front, zslab / ainv)
            if need_normal:
                nacc.addcmul_(front[:, None], nslab / ainv[:, None])
        trans.mul_(1 - slab_a)
        if si % 8 == 7:                    # checkpoint every 8 slabs (a per-slab GPU sync would cost more)
            if float(front.max()) < 1e-3:  # this checkpoint slab is fully occluded by what is in front
                stale += 1
                if stale >= 2:             # two occluded checkpoints running -> the rest are too -> stop
                    break
            else:
                stale = 0

    cov = 1 - trans
    covg = cov.reshape(height, width)
    covm = covg > 0.5 if render_style in ("depth", "normal") else None  # silhouette mask (depth/normal styles only)
    depth_map = (dacc / cov.clamp_min(1e-6)).reshape(height, width) if need_depth else None
    nrm_map = None
    if need_normal:
        # Per-splat surfel normals are jittery, so do a masked blur
        nb = nacc.reshape(height, width, 3).permute(2, 0, 1)[None]
        cb = cov.reshape(1, 1, height, width)
        nb, cb = _gauss_blur(nb, 1.2, dev), _gauss_blur(cb, 1.2, dev)
        normal = (nb / cb.clamp_min(1e-6))[0].permute(1, 2, 0)
        nrm_map = normal / normal.norm(dim=-1, keepdim=True).clamp_min(1e-6)

    if render_style == "depth":  # near = bright, far = dark, 0 off-object
        d = torch.zeros(height, width, device=dev)
        if bool(covm.any()):
            lo, hi = depth_map[covm].min(), depth_map[covm].max()
            d = torch.where(covm, ((hi - depth_map) / (hi - lo).clamp_min(1e-6)).clamp(0, 1), d)
        img = d[:, :, None].expand(height, width, 3)
    elif render_style == "normal":  # OpenGL normal map: +X right, +Y up, +Z to viewer
        enc = (nrm_map * t([1.0, -1.0, -1.0]) * 0.5 + 0.5).clamp(0, 1)
        img = enc * covm[:, :, None]
    else:  # color / clay
        img = cacc.reshape(height, width, 3)
        if render_style == "clay":  # studio key light + ambient -> sculpted matte look
            kl = t([-0.4, -0.7, -0.6])  # key from screen upper-left, angled toward the viewer
            kl = kl / kl.norm()
            hl = (0.5 * (nrm_map * kl).sum(-1) + 0.5).clamp(0, 1)  # half-Lambert: soft terminator, no harsh dark side
            img = img * (0.35 + 0.65 * hl * hl)[:, :, None]        # ambient floor + diffuse key
        elif headlight_shading > 0:  # camera headlight: darken faces turned from view
            k = float(headlight_shading)
            ndotl = (-nrm_map[:, :, 2]).clamp(0, 1)
            img = img * (1 - 0.6 * k + 0.6 * k * ndotl)[:, :, None]
        img = img.addcmul_(trans.reshape(height, width, 1), bg_comp)
        if do_linear:  # back to display space after linear compositing
            img = _linear_to_srgb(img)
    return img.clamp(0, 1).to(idev, idtype), covg.clamp(0, 1).to(idev, idtype)