import numpy as np

class Sellmeier:
    """Caculation of fiber refractive index n based on Sellmeier eqn. n^2-1 = sum_k B_k* wavelength^2/( wavelength^2-C_k)"""
    # Parameters:
    #   B,C - coefficients empirically evaluated for individual materials (dispersion formula)
    #   wavelength - range of the wavelengths for which the dispersion formula is applicable

    # Returns:
    #   n - refractive index
    

    def __init__(self,B,C,wavelength):
        self.B = B
        self.C = C
        self.wavelength = wavelength

    def refractiveindex(self):
        n = np.sqrt(np.maximum(0, 1 + (self.B[0] * self.wavelength ** 2) / (self.wavelength ** 2 - self.C[0]) +
                                     (self.B[1] * self.wavelength ** 2) / (self.wavelength ** 2 - self.C[1]) +
                                     (self.B[2] * self.wavelength ** 2) / (self.wavelength ** 2 - self.C[2])))
        return n
