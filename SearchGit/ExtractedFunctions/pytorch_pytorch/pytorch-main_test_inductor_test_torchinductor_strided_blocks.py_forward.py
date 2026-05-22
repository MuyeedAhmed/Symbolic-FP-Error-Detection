            def forward(
                self,
                sub_dense0_w,  # f16[64, 80]
                sub_dense0_b,  # f16[64]
                s25,  # Sym(s25) - batch size
                s70,  # Sym(s70) - sequence length
                input_feat,  # f16[s25, s70, 80]
                sub_conv0_w,  # f16[64, 64, 5]
                sub_conv0_b,  # f16[64]
                sub_conv1_w,  # f16[32, 64, 5]
                sub_conv1_b,  # f16[32]
                sub_dense1_w,  # f16[64, 32]
                sub_dense1_b,  # f16[64]
                rot_inv_freq,  # f32[8]
                rot_attn_scale,  # f64[] cpu
                l0_norm_ff1_w,  # f16[64]
                l0_ff1_lin1_w,  # f16[256, 64]
                l0_ff1_lin2_w,  # f16[64, 256]
                l0_ff_res0,  # f64[] cpu
                l0_ff_res1,  # f64[] cpu
                l0_norm_attn_w,  # f16[64]
                l0_q_w,  # f16[64, 64]
                l0_k_w,  # f16[64, 64]
                l0_v_w,  # f16[64, 64]
                l0_o_w,  # f16[64, 64]
                l0_norm_conv_w,  # f16[64]
                l0_pw_conv1_w,  # f16[128, 64, 1]
                l0_dw_conv_w,  # f16[64, 1, 8]
                l0_bn_mean,  # f16[64]
                l0_bn_var,  # f16[64]
                l0_bn_w,  # f16[64]
                l0_bn_b,  # f16[64]
                l0_pw_conv2_w,  # f16[64, 64, 1]
                l0_conv_res0,  # f64[] cpu
                l0_conv_res1,  # f64[] cpu
                l0_norm_ff2_w,  # f16[64]
                l0_ff2_lin1_w,  # f16[256, 64]
                l0_ff2_lin2_w,  # f16[64, 256]
                l0_norm_out_w,  # f16[64]
                l1_norm_ff1_w,  # f16[64]
                l1_ff1_lin1_w,  # f16[256, 64]
                l1_ff1_lin2_w,  # f16[64, 256]
                l1_norm_attn_w,  # f16[64]
                l1_q_w,  # f16[64, 64]
                l1_k_w,  # f16[64, 64]
                l1_v_w,  # f16[64, 64]
                l1_o_w,  # f16[64, 64]
                l1_norm_conv_w,  # f16[64]
                l1_pw_conv1_w,  # f16[128, 64, 1]
                l1_dw_conv_w,  # f16[64, 1, 8]
                l1_bn_mean,  # f16[64]
                l1_bn_var,  # f16[64]
                l1_bn_w,  # f16[64]
                l1_bn_b,  # f16[64]
                l1_pw_conv2_w,  # f16[64, 64, 1]
                l1_norm_ff2_w,  # f16[64]
                l1_ff2_lin1_w,  # f16[256, 64]
                l1_ff2_lin2_w,  # f16[64, 256]
                l1_norm_out_w,  # f16[64]
                out_norm_w,  # f16[64]
            ):
                # Subsampler: dense_0 + relu
                linear = torch._C._nn.linear(input_feat, sub_dense0_w, sub_dense0_b)
                hidden_states = F.relu(linear, inplace=False)

                # Transpose for conv
                hidden_states_1 = hidden_states.transpose(1, 2)

                # Conv_0 with stride=2
                conv1d = torch.conv1d(
                    hidden_states_1, sub_conv0_w, sub_conv0_b, (2,), (0,), (1,), 1
                )
                hidden_states_2 = F.relu(conv1d, inplace=False)

                # Conv_1 with stride=2
                conv1d_1 = torch.conv1d(
                    hidden_states_2, sub_conv1_w, sub_conv1_b, (2,), (0,), (1,), 1
                )
                hidden_states_3 = F.relu(conv1d_1, inplace=False)

                # Transpose back
                hidden_states_4 = hidden_states_3.transpose(1, 2)

                # Dense_1
                hidden_states_5 = torch._C._nn.linear(
                    hidden_states_4, sub_dense1_w, sub_dense1_b
                )

                # Compute output sequence length: ((s70 - 5) // 4) - 1
                # Note: In the dynamo graph, torch.sym_sum is used, but we use direct arithmetic here
                sym_sum = s70 - 5
                floordiv = sym_sum // 4
                sym_sum_1 = floordiv - 1

                # Create position ids
                arange = torch.arange(sym_sum_1, device=GPU_TYPE)
                unsqueeze = arange.unsqueeze(0)

                # Rotary embedding computation
                getitem_3 = rot_inv_freq[(None, slice(None, None, None), None)]
                float_1 = getitem_3.float()
                expand = float_1.expand(1, -1, 1)
                inv_freq_expanded = expand.to(torch.device(GPU_TYPE, index=0))

                getitem_6 = unsqueeze[
                    (slice(None, None, None), None, slice(None, None, None))
                ]
                position_ids_expanded = getitem_6.float()

                # Rotary embedding frequency computation (autocast removed for tracing)
                float_3 = inv_freq_expanded.float()
                float_4 = position_ids_expanded.float()
                matmul = float_3 @ float_4
                freqs = matmul.transpose(1, 2)

                emb = torch.cat((freqs, freqs), dim=-1)

                cos = emb.cos()
                item = rot_attn_scale.item()
                cos_1 = cos * item

                sin = emb.sin()
                sin_1 = sin * item

                cos_2 = cos_1.to(dtype=torch.float16)
                sin_2 = sin_1.to(dtype=torch.float16)

                # Dropout (no-op in eval mode)
                hidden_states_6 = F.dropout(hidden_states_5, p=0.1, training=False)
                cos_3 = F.dropout(cos_2, p=0.0, training=False)
                sin_3 = F.dropout(sin_2, p=0.0, training=False)

                # Create attention mask
                cache_position = torch.arange(
                    sym_sum_1, device=GPU_TYPE, dtype=torch.int64
                )
                arange_4 = torch.arange(sym_sum_1, device=GPU_TYPE)

                q_indices = cache_position[(None, None, slice(None, None, None), None)]
                attention_mask = q_indices >= 0
                attention_mask_1 = attention_mask.expand(s25, -1, sym_sum_1, sym_sum_1)

                # ============ LAYER 0 ============
                # Feed forward 1
                layer_norm = F.layer_norm(
                    hidden_states_6, (64,), l0_norm_ff1_w, None, 1e-06
                )
                linear_2 = torch._C._nn.linear(layer_norm, l0_ff1_lin1_w, None)
                hidden_states_7 = F.silu(linear_2)
                hidden_states_8 = F.dropout(hidden_states_7, p=0.1, training=False)
                hidden_states_9 = torch._C._nn.linear(
                    hidden_states_8, l0_ff1_lin2_w, None
                )

                # Residual connection with weights
                item_5 = l0_ff_res0.item()
                mul_2 = item_5 * hidden_states_6
                item_6 = l0_ff_res1.item()
                mul_3 = item_6 * hidden_states_9
                hidden_states_10 = mul_2 + mul_3

                # Self attention
                normalized_hidden_states = F.layer_norm(
                    hidden_states_10, (64,), l0_norm_attn_w, None, 1e-06
                )

                linear_4 = torch._C._nn.linear(normalized_hidden_states, l0_q_w, None)
                view = linear_4.view((s25, sym_sum_1, -1, 16))
                query_states = view.transpose(1, 2)

                linear_5 = torch._C._nn.linear(normalized_hidden_states, l0_k_w, None)
                view_1 = linear_5.view((s25, sym_sum_1, -1, 16))
                key_states = view_1.transpose(1, 2)

                linear_6 = torch._C._nn.linear(normalized_hidden_states, l0_v_w, None)
                view_2 = linear_6.view((s25, sym_sum_1, -1, 16))
                value_states = view_2.transpose(1, 2)

                # Apply rotary embeddings
                cos_4 = cos_3.unsqueeze(1)
                sin_4 = sin_3.unsqueeze(1)

                mul_4 = query_states * cos_4
                x1 = query_states[(Ellipsis, slice(None, 8, None))]
                x2 = query_states[(Ellipsis, slice(8, None, None))]
                neg = -x2
                cat_1 = torch.cat((neg, x1), dim=-1)
                mul_5 = cat_1 * sin_4
                q_embed = mul_4 + mul_5

                mul_6 = key_states * cos_4
                x1_1 = key_states[(Ellipsis, slice(None, 8, None))]
                x2_1 = key_states[(Ellipsis, slice(8, None, None))]
                neg_1 = -x2_1
                cat_2 = torch.cat((neg_1, x1_1), dim=-1)
                mul_7 = cat_2 * sin_4
                k_embed = mul_6 + mul_7

                # SDPA
                attn_output = torch._C._nn.scaled_dot_product_attention(
                    q_embed,
                    k_embed,
                    value_states,
                    attn_mask=attention_mask_1,
                    dropout_p=0.0,
                    scale=0.25,
                    is_causal=False,
                )

                transpose_6 = attn_output.transpose(1, 2)
                attn_output_1 = transpose_6.contiguous()
                reshape = attn_output_1.reshape(s25, sym_sum_1, -1)
                attn_output_2 = reshape.contiguous()
                attn_output_3 = torch._C._nn.linear(attn_output_2, l0_o_w, None)

                hidden_states_11 = hidden_states_10 + attn_output_3

                # Convolution module
                layer_norm_2 = F.layer_norm(
                    hidden_states_11, (64,), l0_norm_conv_w, None, 1e-06
                )
                hidden_states_12 = layer_norm_2.transpose(1, 2)
                hidden_states_13 = torch.conv1d(
                    hidden_states_12, l0_pw_conv1_w, None, (1,), (0,), (1,), 1
                )
                hidden_states_14 = F.glu(hidden_states_13, dim=1)

                invert = ~attention_mask_1
                all_masked_rows = torch.all(invert, dim=2)
                hidden_states_15 = hidden_states_14.masked_fill(all_masked_rows, 0.0)

                hidden_states_16 = torch.conv1d(
                    hidden_states_15, l0_dw_conv_w, None, (1,), "same", (1,), 64
                )
                hidden_states_17 = F.batch_norm(
                    hidden_states_16,
                    l0_bn_mean,
                    l0_bn_var,
                    l0_bn_w,
                    l0_bn_b,
                    False,
                    0.01,
                    1e-05,
                )
                hidden_states_18 = F.silu(hidden_states_17)
                hidden_states_19 = torch.conv1d(
                    hidden_states_18, l0_pw_conv2_w, None, (1,), (0,), (1,), 1
                )
                conv_output = hidden_states_19.transpose(1, 2)

                # Conv residual
                item_12 = l0_conv_res0.item()
                item_13 = l0_conv_res1.item()
                mul_8 = item_12 * hidden_states_11
                mul_9 = item_13 * conv_output
                hidden_states_20 = mul_8 + mul_9

                # Feed forward 2
                layer_norm_3 = F.layer_norm(
                    hidden_states_20, (64,), l0_norm_ff2_w, None, 1e-06
                )
                linear_8 = torch._C._nn.linear(layer_norm_3, l0_ff2_lin1_w, None)
                hidden_states_21 = F.silu(linear_8)
                hidden_states_22 = F.dropout(hidden_states_21, p=0.1, training=False)
                hidden_states_23 = torch._C._nn.linear(
                    hidden_states_22, l0_ff2_lin2_w, None
                )

                mul_10 = item_5 * hidden_states_20
                mul_11 = item_6 * hidden_states_23
                hidden_states_24 = mul_10 + mul_11

                hidden_states_25 = F.layer_norm(
                    hidden_states_24, (64,), l0_norm_out_w, None, 1e-06
                )

                # ============ LAYER 1 ============
                # Feed forward 1
                layer_norm_5 = F.layer_norm(
                    hidden_states_25, (64,), l1_norm_ff1_w, None, 1e-06
                )
                linear_10 = torch._C._nn.linear(layer_norm_5, l1_ff1_lin1_w, None)
                hidden_states_26 = F.silu(linear_10)
                hidden_states_27 = F.dropout(hidden_states_26, p=0.1, training=False)
                hidden_states_28 = torch._C._nn.linear(
                    hidden_states_27, l1_ff1_lin2_w, None
                )

                mul_12 = item_5 * hidden_states_25
                mul_13 = item_6 * hidden_states_28
                hidden_states_29 = mul_12 + mul_13

                # Self attention
                normalized_hidden_states_1 = F.layer_norm(
                    hidden_states_29, (64,), l1_norm_attn_w, None, 1e-06
                )

                linear_12 = torch._C._nn.linear(
                    normalized_hidden_states_1, l1_q_w, None
                )
                view_3 = linear_12.view((s25, sym_sum_1, -1, 16))
                query_states_1 = view_3.transpose(1, 2)

                linear_13 = torch._C._nn.linear(
                    normalized_hidden_states_1, l1_k_w, None
                )
                view_4 = linear_13.view((s25, sym_sum_1, -1, 16))
                key_states_1 = view_4.transpose(1, 2)

                linear_14 = torch._C._nn.linear(
                    normalized_hidden_states_1, l1_v_w, None
                )
                view_5 = linear_14.view((s25, sym_sum_1, -1, 16))
                value_states_1 = view_5.transpose(1, 2)

                # Apply rotary embeddings
                cos_5 = cos_3.unsqueeze(1)
                sin_5 = sin_3.unsqueeze(1)

                mul_14 = query_states_1 * cos_5
                x1_2 = query_states_1[(Ellipsis, slice(None, 8, None))]
                x2_2 = query_states_1[(Ellipsis, slice(8, None, None))]
                neg_2 = -x2_2
                cat_3 = torch.cat((neg_2, x1_2), dim=-1)
                mul_15 = cat_3 * sin_5
                q_embed_1 = mul_14 + mul_15

                mul_16 = key_states_1 * cos_5
                x1_3 = key_states_1[(Ellipsis, slice(None, 8, None))]
                x2_3 = key_states_1[(Ellipsis, slice(8, None, None))]
                neg_3 = -x2_3
                cat_4 = torch.cat((neg_3, x1_3), dim=-1)
                mul_17 = cat_4 * sin_5
                k_embed_1 = mul_16 + mul_17

                # SDPA
                attn_output_4 = torch._C._nn.scaled_dot_product_attention(
                    q_embed_1,
                    k_embed_1,
                    value_states_1,
                    attn_mask=attention_mask_1,
                    dropout_p=0.0,
                    scale=0.25,
                    is_causal=False,
                )

                transpose_12 = attn_output_4.transpose(1, 2)
                attn_output_5 = transpose_12.contiguous()
                reshape_1 = attn_output_5.reshape(s25, sym_sum_1, -1)
                attn_output_6 = reshape_1.contiguous()
                attn_output_7 = torch._C._nn.linear(attn_output_6, l1_o_w, None)

                hidden_states_30 = hidden_states_29 + attn_output_7

                # Convolution module
                layer_norm_7 = F.layer_norm(
                    hidden_states_30, (64,), l1_norm_conv_w, None, 1e-06
                )
                hidden_states_31 = layer_norm_7.transpose(1, 2)
                hidden_states_32 = torch.conv1d(
                    hidden_states_31, l1_pw_conv1_w, None, (1,), (0,), (1,), 1
                )
                hidden_states_33 = F.glu(hidden_states_32, dim=1)

                invert_1 = ~attention_mask_1
                all_masked_rows_1 = torch.all(invert_1, dim=2)
                hidden_states_34 = hidden_states_33.masked_fill(all_masked_rows_1, 0.0)

                hidden_states_35 = torch.conv1d(
                    hidden_states_34, l1_dw_conv_w, None, (1,), "same", (1,), 64
                )
                hidden_states_36 = F.batch_norm(
                    hidden_states_35,
                    l1_bn_mean,
                    l1_bn_var,
                    l1_bn_w,
                    l1_bn_b,
                    False,
                    0.01,
                    1e-05,
                )
                hidden_states_37 = F.silu(hidden_states_36)
                hidden_states_38 = torch.conv1d(
                    hidden_states_37, l1_pw_conv2_w, None, (1,), (0,), (1,), 1
                )
                conv_output_1 = hidden_states_38.transpose(1, 2)

                # Conv residual
                mul_18 = item_12 * hidden_states_30
                mul_19 = item_13 * conv_output_1
                hidden_states_39 = mul_18 + mul_19

                # Feed forward 2
                layer_norm_8 = F.layer_norm(
                    hidden_states_39, (64,), l1_norm_ff2_w, None, 1e-06
                )
                linear_16 = torch._C._nn.linear(layer_norm_8, l1_ff2_lin1_w, None)
                hidden_states_40 = F.silu(linear_16)
                hidden_states_41 = F.dropout(hidden_states_40, p=0.1, training=False)
                hidden_states_42 = torch._C._nn.linear(
                    hidden_states_41, l1_ff2_lin2_w, None
                )

                mul_20 = item_5 * hidden_states_39
                mul_21 = item_6 * hidden_states_42
                hidden_states_43 = mul_20 + mul_21

                hidden_states_44 = F.layer_norm(
                    hidden_states_43, (64,), l1_norm_out_w, None, 1e-06
                )

                # Final output norm
                hidden_states_45 = F.layer_norm(
                    hidden_states_44, (64,), out_norm_w, None, 1e-06
                )

                return (hidden_states_45,)