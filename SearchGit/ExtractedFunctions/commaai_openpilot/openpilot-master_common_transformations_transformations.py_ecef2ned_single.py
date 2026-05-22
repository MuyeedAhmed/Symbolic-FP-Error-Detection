  def ecef2ned_single(self, ecef):
    """
    Convert a single ECEF point to NED coordinates relative to the origin.
    """
    return self.ecef2ned_matrix @ (ecef - self.init_ecef)