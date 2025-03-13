"""
==========================================================================
Dispersion calculation for the optical fiber
==========================================================================

 Author: Eva Roikova
 Organisation: CERN
 Date: 10/03/25

Properties of the fiber
    a_core - core radius [m]
    a_clad - cladding radius [m]
    alpha - fiber absorption 
    L - length of the fiber [m]
    Parameters for refractive index calculation (Sellmeier eqn): 
        B,C - core
        Bcl,Ccl - cladding
    
Laser parameters
    carrier_wavelength - central wavelength [m]
    tp - pulse duration
    
Another parameters  
    c - speed of light [m/s]
    
- data are saved in "fiber_parameters.json"
----------------------------------------------------------------------------
"""
import numpy as np
import matplotlib.pyplot as plt
import json
from fiberDispersionFunctions import GaussianPulse as GP
from fiberDispersionFunctions import FourierTransform as FT
from fiberDispersionFunctions import PhaseDecomposition as PDE
from fiberDispersionFunctions import Dispersion
from fiberDispersionFunctions import PhiwTaylor
from fiberDispersionFunctions import InverseFourierTransform as IFT
from fiberDispersionFunctions import Sellmeier


"""
--------------------------------------------------------------------------
Fiber parameters
--------------------------------------------------------------------------
"""
with open("fiber_parameters.json") as f:
    fiber_data = json.load(f)
    
a_core = fiber_data["a_core"]    #core radius
a_clad = fiber_data["a_clad"]    # cladding radius
c = fiber_data["c"]     #speed of light [m/s]
carrier_wavelength = fiber_data["lambda0"]         # [m] central wavelength of carrier 
f0 = c/ carrier_wavelength                # [Hz] carrier freq 
w0 = 2*np.pi*f0                  # [rad/s] angular frequency

T = 1e-16
t = np.arange(-1000e-12,1000e-12,T) # time vector

#range for plotting in time domain
dt = 10e-13
it = np.where((-dt < t) & (t < dt))[0]
del dt

"""
--------------------------------------------------------------------------
Pulse definition
--------------------------------------------------------------------------
"""
# Gaussian envelope (pulse)
tp = fiber_data["tp"]    #time of the pulse [s] in INTENSITY!!
pulse = GP(t, 0, tp)
E_t_in = pulse.efield() # input laser electric field
I_t_in = pulse.intensity() # input laser intensity
Phi_t_in = np.zeros(np.size(t)) # input phase

# Plot of input laser pulse intensity in time domain
fig1, ax1 = plt.subplots()
ax1.plot(t[it], I_t_in[it])
ax1.set_xlabel("Time [s]")
ax1.set_ylabel("Intensity [a.u.]")
ax1.set_title("Input laser pulse")
plt.show()

# From time domain to frequency domain
fft_transform = FT(E_t_in,Phi_t_in,t)
E_w, S_w, phi_w, w = fft_transform.FFT_t2w()

dw = 5e13
iw = np.where((w > -dw) & (w < dw))[0]  # wavelength range applicable for silica

"""
--------------------------------------------------------------------------
Refractive index calculation
--------------------------------------------------------------------------
"""
L = fiber_data["L"]         # length of the fiber [m]
wavelength = 2*np.pi*c/(w+w0)*1e6      # wavelength in um!! due to Sellmeier
wavelength_min = fiber_data["wavelength_min"]
wavelength_max = fiber_data["wavelength_max"]

il = np.where((wavelength_min< wavelength) & (wavelength<wavelength_max))[0]    # wavelength range applicable for silica
il1 = np.where((0.630< wavelength) & (wavelength<0.800))[0]    # wavelength range for plotting

blankn = np.zeros(np.shape(wavelength))   # blanking for the refractive index - zero outside of the range (il)
blankn[il] = 1

# Parameters for Sellmeier eqn. n^2-1 = sum_k B_k* wavelength^2/( wavelength^2-C_k)
# core
B = np.array(fiber_data["B"])
C = np.array(fiber_data["C"])**2
# cladding
Bcl = np.array(fiber_data["Bcl"])
Ccl = np.array(fiber_data["Ccl"])**2

# refractive index (material of the core and the cladding)
core = Sellmeier(B,C,wavelength)
n0_core = core.refractiveindex()
n_core = n0_core*blankn
cladding = Sellmeier(Bcl,Ccl,wavelength)
n0_clad = cladding.refractiveindex()
n_clad = n0_clad*blankn

# effective index (waveguide)
# Marcuse approximation
V = 2.*np.pi /(wavelength*1e-6) * a_core * np.sqrt(n_core**2-n_clad**2)
V = np.nan_to_num(V)
epsilon = 1e-10 # to prevent division by zero
u = 2.4048 # empirical value
bV = (1+(u / (V + epsilon))**2) / (1+(u / (V + epsilon))**2 + 2 / (V + epsilon) * np.sqrt(np.maximum(0,1+1 / (V+ epsilon)**2)))
n_eff = np.sqrt(np.maximum(0,n_clad**2 + bV * (n_core**2-n_clad**2)))

# refractive index plot - wavelength
fig2, ax2 = plt.subplots(2,1)
ax2[0].plot(wavelength[il1], n_core[il1], label='Core')
ax2[0].plot(wavelength[il1], n_clad[il1], label='Cladding')
ax2[0].plot(wavelength[il1], n_eff[il1], label='Effective')
ax2[0].set_xlabel(r'$\lambda$ [$\mu m$]')
ax2[0].set_ylabel('Refractive index n')
ax2[0].set_title('Refractive index distribution')
ax2[0].legend()

# refractive index plot - frequency
ax2[1].plot((w+w0)[il1], n_core[il1], label='Core')
ax2[1].plot((w+w0)[il1], n_clad[il1], label='Cladding')
ax2[1].plot((w+w0)[il1], n_eff[il1], label='Effective')
ax2[1].set_xlabel(r'$\omega$ [Hz]')
ax2[1].set_ylabel('Refractive index n')
ax2[1].legend()
plt.tight_layout()
plt.show()

"""
--------------------------------------------------------------------------
Absorption
--------------------------------------------------------------------------
"""
alpha = fiber_data["alpha"] #theoreticaly could be an array, for now just zeros


"""
--------------------------------------------------------------------------
Chromatic dispersion
--------------------------------------------------------------------------
"""
fibre_dispersion = Dispersion(S_w,phi_w,alpha,n_eff,L,w,w0,c)
E_w_out_tot, S_w_out_tot, phi_w_out_tot = fibre_dispersion.compute() #calculation of fiber dispersion effect on laser pulse

# plots of spectrum and spectral phase after propagation through the fiber
fig3, ax3 = plt.subplots()
ax3.plot(w[iw],np.unwrap(phi_w_out_tot)[iw], color='tab:blue')
ax3.set_xlabel(r'$\omega$ [Hz]')
ax3.set_ylabel('Total spectral phase [rad]', color='tab:blue')
ax3.set_title('Chromatic dispersion')
ax3.tick_params(axis='y', labelcolor='tab:blue')
ax3.set_title('Spectral domain')

ax31 = ax3.twinx()
ax31.plot(w[iw], S_w_out_tot[iw], color='tab:orange')
ax31.set_ylabel(r'Spectrum S($\omega$) [a.u.]', color='tab:orange')
ax31.tick_params(axis='y', labelcolor='tab:orange')
plt.tight_layout()
plt.show()

"""
--------------------------------------------------------------------------
Phase decomposition 
--------------------------------------------------------------------------
"""
order = 4

# dispersion coefficient
phase_decomp_out = PDE(w,S_w_out_tot,phi_w_out_tot,order)
phase_coeff_out_tot = phase_decomp_out.get_phase_coefficients() #decomposition of phase into the polynomial coefficients


GD = phase_coeff_out_tot[1] * 1e15
GDD = phase_coeff_out_tot[2] * 1e30
TOD = phase_coeff_out_tot[3] * 1e45
FOD = phase_coeff_out_tot[4] * 1e60


print('phase coefficients: \n',
      'Group delay GD = ', f"{GD:.2e}", '[fs]', '\n',
      'Group delay dispersion GDD = ',f"{GDD:.2e}", '[fs^2]', '\n',
      'Third order dispersion TOD = ', f"{TOD:.2e}", '[fs^3]', '\n',
      'Fourth order dispersion FOD = ', f"{FOD:.2e}", '[fs^4]', '\n')

phi_tot_compressed = np.array([0, phase_coeff_out_tot[1]*0,  phase_coeff_out_tot[2]*0,  phase_coeff_out_tot[3]]) # correction of the dispersion (Fourier limited pulse)
phi_out_taylor = PhiwTaylor(w,0,phi_tot_compressed)
phi_w_out_centr= phi_out_taylor.calculate() # reconstruction of the phase from corrected coefficients (as a Taylor series)


"""
--------------------------------------------------------------------------
From spectral domain back to time domain
--------------------------------------------------------------------------
"""
ifft_transform = IFT(S_w_out_tot,phi_w_out_centr)
E_t_out,I_t_out,phi_t_out = ifft_transform.computeIFFT()    #Transition from the spectral domain to the time domain is achieved by inverse Fourier transform

# plot of the input pulse and Fourier limited output pulse after propagation through optical fiber
fig4, ax4 = plt.subplots()
ax4.plot(t[it],I_t_in[it], label='input pulse')
ax4.plot(t[it],I_t_out[it], label='output pulse')
ax4.set_xlabel('Time [s]')
ax4.set_ylabel('Intensity [a.u.]')
ax4.set_title('Reconstruction of the pulse')
ax4.legend()
plt.show()



