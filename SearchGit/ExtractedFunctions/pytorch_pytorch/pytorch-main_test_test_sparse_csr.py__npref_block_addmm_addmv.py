def _npref_block_addmm_addmv(c, a, b, alpha, beta):
    return alpha * (a @ b) + beta * c