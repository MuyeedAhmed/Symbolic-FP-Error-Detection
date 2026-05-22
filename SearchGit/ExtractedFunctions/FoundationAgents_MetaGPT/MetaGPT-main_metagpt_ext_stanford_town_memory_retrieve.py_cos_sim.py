def cos_sim(a, b):
    """
    计算余弦相似度
    """
    return dot(a, b) / (norm(a) * norm(b))