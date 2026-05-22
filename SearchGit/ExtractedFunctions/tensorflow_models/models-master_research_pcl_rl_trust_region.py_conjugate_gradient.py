def conjugate_gradient(f_Ax, b, cg_iters=10, residual_tol=1e-10):
  p = b.copy()
  r = b.copy()
  x = np.zeros_like(b)
  rdotr = r.dot(r)
  for i in xrange(cg_iters):
    z = f_Ax(p)
    v = rdotr / p.dot(z)
    x += v * p
    r -= v * z
    newrdotr = r.dot(r)
    mu = newrdotr / rdotr
    p = r + mu * p
    rdotr = newrdotr
    if rdotr < residual_tol:
      break
  return x