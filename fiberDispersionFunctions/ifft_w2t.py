import numpy as np

class InverseFourierTransform:
    """iFFT of the electric field in the frequency domain to the time domain"""
    # Parameters:
    #   S_w - spectrum
    #   phi_w - spectral phase

    # Returns:
    #   E_t - electric field in the time domain
    #   I_t - intensity in the time domain
    #   phi_t - phase

    def __init__(self,S_w,phi_w):
        self.S_w = S_w
        self.phi_w = phi_w

    def computeIFFT(self):
        E_w = np.sqrt(np.abs(self.S_w)) * np.exp(-1j * self.phi_w) # Electric field in the spectral domain
        E_t = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(E_w))) # Electric field in the time domain (inverse fourier transform from spectral domain)
        I_t = np.abs(E_t)**2 # Intensity of the electric field
        phi_t = np.angle(E_t) # phase of the field
        return E_t,I_t,phi_t
