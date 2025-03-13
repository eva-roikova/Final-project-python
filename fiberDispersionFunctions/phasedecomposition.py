import numpy as np
from scipy.integrate import cumtrapz, trapz
from math import factorial

class PhaseDecomposition:
    """Computes the phase dispersion coefficients using polynomial fitting."""
    # Parameters:
    #   w - Angular frequency range
    #   S_w - Spectrum (power spectrum)
    #   phi_w - Spectral phase
    #   degree - Degree of polynomial fit

    # Returns:
    #    Phase dispersion coefficients (phi_0, GD, GDD, TOD, etc.)

    def __init__(self, w, S_w, phi_w, degree):
        self.w = w
        self.S_w = S_w
        self.phi_w = phi_w
        self.degree = degree
        self.phase_Di = self.compute_phase_dispersion()

    def compute_phase_dispersion(self):
        # Compute 2.5% threshold of integral of S(w)
        d = trapz(self.S_w) * 0.025

        # Compute cumulative integral of S(w)
        S_wint = cumtrapz(self.S_w, initial=0)

        # Find indices where S_wint is within the 2.5% - 97.5% range
        iw = np.where((d < S_wint) & (S_wint < trapz(self.S_w) - d))[0]

        # Polynomial fitting of the phase in the selected range
        p = np.polyfit(self.w[iw], np.unwrap(self.phi_w[iw]), self.degree)

        # Convert polynomial coefficients to Taylor series coefficients
        phase_Di = np.zeros(len(p))
        for i in range(len(p)):
            phase_Di[i] = p[-(i + 1)] * factorial(i)  # Convert polyfit coefficients

        return phase_Di

    def get_phase_coefficients(self):
        """Returns the computed phase dispersion coefficients."""
        return self.phase_Di
