  def run_rnn(rnn, x0, ntime_steps, input_time=None):
    rs = np.zeros([N,ntime_steps])
    x_tm1 = x0
    r_tm1 = np.tanh(x0)
    tau = rnn['tau']
    dt = rnn['dt']
    alpha = (1.0-dt/tau)
    W = dt/tau*rnn['W']*rnn['g']
    Bin = dt/tau*rnn['Bin']
    Bin2 = dt/tau*rnn['Bin2']
    b = dt/tau*rnn['b']

    us = np.zeros([1, ntime_steps])
    for t in range(ntime_steps):
      x_t = alpha*x_tm1 + np.dot(W,r_tm1) + b
      if input_time is not None and t == input_time:
        us[0,t] = input_magnitude
        x_t += Bin * us[0,t] # DCS is this what was used?
      r_t = np.tanh(x_t)
      x_tm1 = x_t
      r_tm1 = r_t
      rs[:,t] = r_t
    return rs, us