    def mlstm_chunkwise_parallel_fw_H(
        matQ: torch.Tensor,
        matK: torch.Tensor,
        matV: torch.Tensor,
        # these states must be all states up to the last chunk, i.e. :-1
        matC_states: torch.Tensor,
        vecN_states: torch.Tensor,
        scaMinter_states: torch.Tensor,
        vecI: torch.Tensor,
        vecB: torch.Tensor,
        qk_scale: float,
        chunk_size: int = 64,
        num_chunks: int = 1,
        eps: float = 1e-6,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        _device = matQ.device
        nc = num_chunks
        batch_size, nh, dqk, dhv = matC_states.shape
        dhqk = dqk // nc
        matC_k_states = matC_states.view(batch_size, nh, nc, dhqk, dhv)
        vecN_k_states = vecN_states.view(batch_size, nh, nc, dhqk)
        scaMinter_k_states = scaMinter_states

        matQ = matQ.view(batch_size, nh, nc, chunk_size, dhqk)
        matK = matK.view(batch_size, nh, nc, chunk_size, dhqk)
        matV = matV.view(batch_size, nh, nc, chunk_size, dhv)

        ltr = torch.tril(
            torch.ones(
                (chunk_size, chunk_size),
                dtype=torch.bool,
                device=_device,
            )
        )

        # Compute intra chunk contribution: H_intra
        matF_logsig_chunk = vecB[:, :, :, :, None] - vecB[:, :, :, None, :]

        matF_logsig_mask_chunk = torch.where(ltr, matF_logsig_chunk, -float("inf"))

        matLogD_chunk = matF_logsig_mask_chunk + vecI[:, :, :, None, :]

        # max_state intra
        vecMintra_k = torch.max(matLogD_chunk, dim=-1, keepdim=False).values

        # max_state combined
        vecM_b_inter = vecB + scaMinter_k_states[:, :, :, None]
        vecM_k_combine = torch.maximum(vecM_b_inter, vecMintra_k)

        vecM_k_combine = vecM_k_combine[:, :, :, :, None]
        vecM_b_inter = vecM_b_inter[:, :, :, :, None]

        matLogD_stabilized_chunk = matLogD_chunk - vecM_k_combine
        matD_chunk = torch.exp(matLogD_stabilized_chunk)

        matS_chunk = (matQ @ matK.transpose(-2, -1)) * qk_scale

        matM_chunk = matS_chunk * matD_chunk

        # ? Combine H_intra with H_inter
        vecBbar = torch.exp(vecM_b_inter - vecM_k_combine)
        matQ_chunk_gated = matQ * vecBbar * qk_scale

        matNumerator_common = matQ_chunk_gated @ matC_k_states + matM_chunk @ matV

        vecDenom_l_common = matQ_chunk_gated @ vecN_k_states.unsqueeze(-1) + matM_chunk.sum(dim=-1, keepdim=True)

        vecDenom_max_common = torch.maximum(torch.abs(vecDenom_l_common), torch.exp(-vecM_k_combine))

        matH_k_chunk = matNumerator_common / (vecDenom_max_common + eps)

        matH_out = matH_k_chunk.view(batch_size, nh, nc * chunk_size, dhv)

        # we need the denominator and the overall max state for the backward pass
        vecN_out = vecDenom_max_common.reshape(batch_size, nh, nc * chunk_size)
        vecM_out = vecM_k_combine.reshape(batch_size, nh, nc * chunk_size)
        return matH_out, vecN_out, vecM_out