        def f(x, w, targets):
            logits = (x * 2) @ w
            max_logits = logits.amax(dim=-1)
            shifted = logits - max_logits.unsqueeze(-1)
            exp_shifted = shifted.exp()
            sum_exp = exp_shifted.sum(dim=-1)
            log_probs = shifted - sum_exp.log().unsqueeze(-1)
            target_log_probs = log_probs.gather(1, targets.unsqueeze(1))
            loss = -target_log_probs.squeeze(1).sum() / M
            loss.backward()
            return loss