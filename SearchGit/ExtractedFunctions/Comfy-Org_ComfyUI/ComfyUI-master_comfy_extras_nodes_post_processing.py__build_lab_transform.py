    def _build_lab_transform(cls, image_target, image_ref, device, stats_mode, target_index, is_reinhard):
        """Build transform parameters for Lab-based methods. Returns a transform function."""
        eps = 1e-6
        B, H, W, C = image_target.shape
        B_ref = image_ref.shape[0]
        single_ref = B_ref == 1
        HW = H * W
        HW_ref = image_ref.shape[1] * image_ref.shape[2]

        # Precompute ref stats
        if single_ref or stats_mode in ('uniform', 'target_frame'):
            ref_mean, ref_sc = cls._pool_stats(image_ref, device, is_reinhard, eps)

        # Uniform/target_frame: precompute single affine transform
        if stats_mode in ('uniform', 'target_frame'):
            if stats_mode == 'target_frame':
                ti = min(target_index, B - 1)
                s_lab = cls._to_lab(image_target, ti, device).view(C, -1)
                s_mean, s_sc = cls._frame_stats(s_lab, HW, is_reinhard, eps)
            else:
                s_mean, s_sc = cls._pool_stats(image_target, device, is_reinhard, eps)

            if is_reinhard:
                scale = ref_sc / s_sc
                offset = ref_mean - scale * s_mean
                return lambda src_flat, **_: src_flat * scale + offset
            T = cls._mkl_matrix(s_sc, ref_sc, eps)
            offset = ref_mean - T @ s_mean
            return lambda src_flat, **_: T @ src_flat + offset

        # per_frame
        def per_frame_transform(src_flat, frame_idx):
            s_mean, s_sc = cls._frame_stats(src_flat, HW, is_reinhard, eps)

            if single_ref:
                r_mean, r_sc = ref_mean, ref_sc
            else:
                ri = min(frame_idx, B_ref - 1)
                r_mean, r_sc = cls._frame_stats(cls._to_lab(image_ref, ri, device).view(C, -1), HW_ref, is_reinhard, eps)

            centered = src_flat - s_mean
            if is_reinhard:
                return centered * (r_sc / s_sc) + r_mean
            T = cls._mkl_matrix(centered @ centered.T / HW, r_sc, eps)
            return T @ centered + r_mean

        return per_frame_transform