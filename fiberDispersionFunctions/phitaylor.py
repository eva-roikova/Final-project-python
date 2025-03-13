import numpy as np
import math

class PhiwTaylor:
    """Phase Taylor expansion in Frequency domain"""
    # Parameters:
    #   w - angular frequency range
    #   w0 - carrier freq.
    #   phase coefficients:
    #       phi_0 - absolute phase (same for frequency and time domains)
    #       phi_1... GD = dphi/dw;
    #       phi_2... GDD = dGD/dw = d2phi/dw2
    #       phi_3... TOD = dGDD/dw = d3phi/dw3
    #       etc...
    # phiw_taylor = phi_0 + phi_1.*(w-w0)./factorial(1) +  phi_2.*(w-w0).^2./factorial(2) + phi_3.*(w-w0).^3./factorial(3) +...;

    # Returns:
    #   phiw_taylor - reconstructed phase

    def __init__(self,w,w0,phi_tot):
        self.w = w
        self.w0 = w0
        self.phi_tot = phi_tot

    def calculate(self):
        length_phi = np.shape(self.phi_tot)[0]
        length_w = np.shape(self.w)[0]
        phiw_taylor_table = np.zeros([length_phi,length_w])
        for i in range(length_phi):
            phiw_taylor_table[i,:] = (self.phi_tot[i] * (self.w-self.w0) ** i / math.factorial(i))

        phiw_taylor = np.sum(phiw_taylor_table, axis=0)
        return phiw_taylor
