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