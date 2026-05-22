    def mlstm_recurrent_step_native(
        query: torch.Tensor,
        key: torch.Tensor,
        value: torch.Tensor,
        igate: torch.Tensor,
        fgate: torch.Tensor,
        cstate: torch.Tensor,
        nstate: torch.Tensor,
        mstate: torch.Tensor,
        eps: float = 1e-6,
        dtype_state: torch.dtype = torch.float32,
        **kwargs,
    ) -> tuple[torch.Tensor, tuple[torch.Tensor, torch.Tensor, torch.Tensor]]:
        """This is a single step of the mLSTM operation in recurrent form."""
        dtype_qkv = query.dtype
        matC_old = cstate.to(dtype=dtype_state)
        vecN_old = nstate.to(dtype=dtype_state)
        scaM_old = mstate.to(dtype=dtype_state)

        batch_size, nh, dhqk = query.shape
        _, _, dhhv = value.shape
        if query.shape != key.shape:
            raise ValueError("query and key must have the same shape")
        if matC_old.shape != (batch_size, nh, dhqk, dhhv):
            raise ValueError(f"matC_old has wrong shape, got {matC_old.shape}")
        if vecN_old.shape != (batch_size, nh, dhqk):
            raise ValueError(f"vecN_old has wrong shape, got {vecN_old.shape}")
        if scaM_old.shape != (batch_size, nh, 1):
            raise ValueError(f"scaM_old has wrong shape, got {scaM_old.shape}")
        if igate.shape != (batch_size, nh, 1):
            raise ValueError(f"scaI has wrong shape, got {igate.shape}")
        if fgate.shape != (batch_size, nh, 1):
            raise ValueError(f"scaF has wrong shape, got {fgate.shape}")

        # gates
        scaF_log = torch.nn.functional.logsigmoid(fgate)

        # update rule
        scaM_state_new = torch.max(scaF_log + scaM_old, igate)

        scaF_act = torch.exp(scaF_log + scaM_old - scaM_state_new)
        scaI_act = torch.exp(igate - scaM_state_new)

        vecQ_scaled = query * (dhqk ** (-0.5))
        matC_state_new = scaF_act[:, :, :, None] * matC_old.clone() + scaI_act[:, :, :, None] * (
            key[:, :, :, None] @ value[:, :, None, :]
        )
        vecN_state_new = scaF_act * vecN_old.clone() + scaI_act * key
        h_num = vecQ_scaled[:, :, None, :] @ matC_state_new.to(dtype=dtype_qkv)
        h_num = h_num.squeeze(2).to(dtype=dtype_state)

        qn_dotproduct = vecQ_scaled[:, :, None, :] @ vecN_state_new[:, :, :, None].to(dtype=dtype_qkv)
        qn_dotproduct = qn_dotproduct.squeeze(2)
        max_val = torch.exp(-scaM_state_new)
        h_denom = (torch.maximum(qn_dotproduct.abs(), max_val) + eps).to(dtype=dtype_state)
        h = h_num / h_denom

        h = h.to(dtype=dtype_qkv)
        matC_state_new = matC_state_new.to(dtype=dtype_state)
        vecN_state_new = vecN_state_new.to(dtype=dtype_state)
        scaM_state_new = scaM_state_new.to(dtype=dtype_state)
        return h, (matC_state_new, vecN_state_new, scaM_state_new)