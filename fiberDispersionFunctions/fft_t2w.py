import numpy as np

class FourierTransform:
    """FFT of the electric field in the time domain to the frequency domain"""
    # Parameters:
    #   E_t - Electric field envelope
    #   Phi_t - Phase in the time domain
    #   t - time

    # Returns:
    #   E_w - Electric field in the frequency domain
    #   S_w - Spectrum (Power Spectrum)
    #   phi_w - Spectral phase
    #   w - Angular frequency range

    def __init__(self,E_t,phi_t,t):
        self.t = t
        self.Fs = 1 / (t[1] - t[0])  # Sampling frequency
        self.N = len(t)              # Signal length
        self.f = np.linspace(-self.Fs / 2, self.Fs / 2, self.N)  # Frequency range
        self.w = self.f * 2 * np.pi  # Angular frequency range
        self.E_t = E_t
        self.phi_t = phi_t

    def FFT_t2w(self):
        Phi_t = np.exp(-1j * self.phi_t)  # Convert phase to complex form
        # Perform FFT with fftshift and ifftshift for proper centering
        E_w = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(np.abs(self.E_t) * Phi_t))) # Electric field
        S_w = np.abs(E_w) ** 2  # Spectrum
        phi_w = np.angle(E_w)   # Spectral phase

        return E_w, S_w, phi_w, self.w  # Returns angular frequency range as well
