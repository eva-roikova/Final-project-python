# Final-project-python
Calculation of dispersion for step-index fibers 

- Implementation of existing Matlab code into Python
- main.py is the main script
- fiberDispersionFunctions is a package with used functions

The script is working with the parameters:
  Properties of the fiber:
      a_core - core radius [m]
      a_clad - cladding radius [m]
      alpha - fiber absorption 
      L - length of the fiber [m]
      Parameters for refractive index calculation (Sellmeier eqn): 
          B, C - core
          Bcl, Ccl - cladding
      
  Laser parameters:
      carrier_wavelength - central wavelength [m]
      tp - pulse duration
      
  Another parameter:  
      c - speed of light [m/s]
      
  - parameters are saved in "fiber_parameters.json"


The Input Gaussian laser pulse is defined in the time domain by its electric field, intensity, and phase. It is then transformed into the frequency domain by a Fourier transform, where it is characterized by the spectrum and spectral phase. Then, the dispersive medium (optical fiber) is applied, which changes the spectral phase. The dispersed pulse is then transformed from the frequency domain into the time domain by inverse Fourier transform.

The knowledge of phase dispersion orders (phase coefficients) as group delay, group delay dispersion, third order dispersion, and higher orders is crucial for dispersion compensation, which can be done by different materials before the propagation through the fiber after that. Which can lead to the production of narrow pulses.
