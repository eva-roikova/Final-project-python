import numpy as np
class GaussianPulse:
    """Function for the construction of a Gaussian pulse"""
    # Parameters:
    #  t - time range
    #  t0 - typically 0 because peak is placed intentionaly into 0
    #  tp - pulse duration in I(t)

    # Returns:
    #  E_t_in - pulse electric field
    #  I_t_in - pulse intensity


    def __init__(self,t,t0,tp):
        self.t = t
        self.t0 = t0
        self.tp = tp

    def efield(self):
        self.sig = self.tp / (2 * np.sqrt(2 * np.log(2)))  # standard deviation
        self.A = 1 / np.sqrt(np.dot(self.sig, np.sqrt(2 * np.pi))) # amplitude
        self.E_t_in = np.dot(self.A, np.exp(-1 / 2 * (self.t - self.t0) ** 2 / (2 * self.sig ** 2)))  # gaussian pulse (E field (t))
        return self.E_t_in

    def intensity(self):
        self.I_t_in = np.abs(self.E_t_in) ** 2
        return self.I_t_in