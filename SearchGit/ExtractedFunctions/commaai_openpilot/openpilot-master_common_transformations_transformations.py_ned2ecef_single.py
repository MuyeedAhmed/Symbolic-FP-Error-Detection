  def ned2ecef_single(self, ned):
    """
    Convert a single NED point to ECEF coordinates.
    """
    return self.ned2ecef_matrix @ ned + self.init_ecef