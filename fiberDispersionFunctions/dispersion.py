import numpy as np

class Dispersion:
    """Compute electric field, spectrum and spectral phase after passing through dispersive material"""
    # Parameters:
    #   S_w - input spectrum
    #   phi_w - input spectral phase
    #   alfa - absorption coefficient (generally function)
    #   n - refractive index (generally function - could be calculated by Sellmeier eqn.)
    #   L - thickness of the material (length of the fiber)
    #   w - angular frequency
    #   w0 - carrier angular frequency (n is calculated in "original" wavelengths/freq. it's necessary to "shift" w by w0 in electric field calculation)
    #   c - speed of light

    # Returns:
    # Ew_out - electric field
    # Sw_out - spectrum
    # phiw_out - spectral phase
    def __init__(self,S_w,phi_w,alfa,n,L,w,w0,c):
        self.S_w = S_w
        self.phi_w = phi_w
        self.alfa = alfa
        self.L = L
        self.n = n
        self.w = w
        self.c = c
        self.w0 = w0

    def compute(self):
        Ew_out = np.sqrt(self.S_w) * np.exp(-1j * self.phi_w) * np.exp(-self.alfa * self.L / 2) * np.exp(1j * self.n * (self.w + self.w0) / self.c * self.L)
        Ew_out = np.nan_to_num(Ew_out) # electric field
        Sw_out = abs(Ew_out)** 2  # spectrum
        phiw_out = np.angle(Ew_out)  # spectral phase
        return Ew_out, Sw_out, phiw_out