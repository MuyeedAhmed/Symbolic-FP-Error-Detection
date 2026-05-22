def generate_data(rnn, T, E, x0s=None, P_sxn=None, input_magnitude=0.0,
                  input_times=None):
  """ Generates data from an randomly initialized RNN.
  Args:
    rnn: the rnn
    T: Time in seconds to run (divided by rnn['dt'] to get steps, rounded down.
    E: total number of examples
    S: number of samples (subsampling N)
  Returns:
    A list of length E of NxT tensors of the network being run.
  """
  N = rnn['N']
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

  if P_sxn is None:
    P_sxn = np.eye(N)
  ntime_steps = int(T / rnn['dt'])
  data_e = []
  inputs_e = []
  for e in range(E):
    input_time = input_times[e] if input_times is not None else None
    r_nxt, u_uxt = run_rnn(rnn, x0s[:,e], ntime_steps, input_time)
    r_sxt = np.dot(P_sxn, r_nxt)
    inputs_e.append(u_uxt)
    data_e.append(r_sxt)

  S = P_sxn.shape[0]
  data_e = normalize_rates(data_e, E, S)

  return data_e, x0s, inputs_e