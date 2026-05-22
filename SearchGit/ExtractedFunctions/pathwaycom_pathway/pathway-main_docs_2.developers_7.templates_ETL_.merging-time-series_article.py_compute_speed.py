def compute_speed(t_prev, position_prev, t_next, position_next):
    try:
        _, _, distance_2d = g.inv(
            position_prev[1], position_prev[0], position_next[1], position_next[0]
        )
    except:
        return 0.0
    return float(distance_2d / (t_next - t_prev))