  def car_space_to_ff(self, x, y, z):
    car_space_projective = np.column_stack((x, y, z)).T

    ep = self.extrinsics_matrix.dot(car_space_projective)
    kep = self.intrinsic.dot(ep)
    return (kep[:-1, :] / kep[-1, :]).T