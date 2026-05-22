        def _run_test(x_hp, x_recipe, x_fp8, x_scales, x_scales_original,
                      y_hp, y_recipe, y_fp8, y_scales, y_scales_original):

            # Calculate actual F8 mm
            out_scaled_mm = scaled_mm_wrap(
                x_fp8,
                y_fp8.t(),
                scale_a=x_scales.reciprocal(),
                scale_recipe_a=x_recipe,
                # Note: No more .t() on scale_b, not necessary.
                scale_b=y_scales.reciprocal(),
                scale_recipe_b=y_recipe,
                out_dtype=output_dtype,
            )

            # Calculate emulated F8 mm
            out_emulated = mm_float8_emulated_block(
                x_fp8,
                x_scales_original,
                y_fp8.t(),
                y_scales_original.t(),
                output_dtype
            )

            cosine_sim = torch.nn.functional.cosine_similarity(
                out_emulated.flatten().float(), (x @ y.t()).flatten().float(), dim=0
            )
            self.assertGreaterEqual(float(cosine_sim), 0.999)

            cosine_sim = torch.nn.functional.cosine_similarity(
                out_scaled_mm.flatten().float(), out_emulated.flatten().float(), dim=0
            )
            self.assertGreaterEqual(float(cosine_sim), 0.999)

            if output_dtype in {torch.bfloat16, torch.float16}:
                atol, rtol = 6e-1, 7e-2
            else:
                atol, rtol = 7e-1, 2e-3

            self.assertEqual(out_scaled_mm, out_emulated.to(output_dtype), atol=atol, rtol=rtol)

            # One last check against the full-precision reference, to ensure we
            # didn't mess up the scaling itself and made the test trivial.
            cosine_sim = torch.nn.functional.cosine_similarity(
                out_scaled_mm.flatten().float(), (x @ y.t()).flatten().float(), dim=0
            )
            self.assertGreaterEqual(float(cosine_sim), 0.999)