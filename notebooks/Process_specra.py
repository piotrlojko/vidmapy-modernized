"""
This script plots the synthetic spectra and converts them into CROCODILE-friendly format:

wavelength_Angstrom    flux
3000.000    1.234567e-10
...    ...

Note: no errors are provided for synthetic specra!
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u

plt.rcParams['figure.dpi'] = 150
plt.rcParams['axes.linewidth'] = 1
plt.rcParams['axes.titlepad'] = 15
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['axes.formatter.use_mathtext'] = True
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['axes.labelsize'] = 16
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'
plt.rcParams['savefig.dpi'] = 300
plt.rcParams["savefig.directory"] = os.getcwd()


class spectrum:
    # arguments are passed via the constructor method
    def __init__(self, path):  # this is a constructor method that initializes automatically the attributes of the class when a `spectrum` object is created
        self.wave, self.Edd_flux, self.Edd_continuum, self.normalized_flux = np.loadtxt(path, unpack=True, dtype=float, comments="#")

save = False  
if not save:
    print("Spectra will not be saved to text files. Set `save = True` inside the script to save them.")
    
try:
    hot_star_spectrum = spectrum('hot_star_spectrum.txt')
except FileNotFoundError:
    print("File `hot_star_spectrum.txt` not found. Generate it first by running the proper python script.")

try:
    blap_spectrum = spectrum('blap_spectrum.txt')
except FileNotFoundError:
    print("File `blap_spectrum.txt` not found. Generate it first by running the proper python script.")

try:
    vega_spectrum = spectrum('vega_spectrum.txt')
except FileNotFoundError:
    print("File `vega_spectrum.txt` not found. Generate it first by running the proper python script.")
    
try:
    sun_spectrum = spectrum('sun_spectrum.txt')
except FileNotFoundError:
    print("File `sun_spectrum.txt` not found. Generate it first by running the proper python script.")

#---------------------------------------------
fig, ax = plt.subplots(figsize=(13, 6))
ax.plot(hot_star_spectrum.wave, hot_star_spectrum.Edd_flux, label="Hot Star (Teff=29000K, logg=4.5) [R=50k]", color='blue', alpha=0.7, lw=0.7)
ax.plot(hot_star_spectrum.wave, hot_star_spectrum.Edd_continuum, color='blue', alpha=0.7, lw=0.7, ls='--')

ax.plot(blap_spectrum.wave, blap_spectrum.Edd_flux, label="BLAP (Teff=29000K, logg=4.5, He-rich, N-rich, C/O/Si/Fe-poor) [R=50k]", color='orange', alpha=0.7, lw=0.7)
ax.plot(blap_spectrum.wave, blap_spectrum.Edd_continuum, color='orange', alpha=0.7, lw=0.7, ls='--')

ax.plot(vega_spectrum.wave, vega_spectrum.Edd_flux, label="Vega (Teff=9550K, logg=4.1) [R=30k]", color='green', alpha=0.7, lw=0.7)
ax.plot(vega_spectrum.wave, vega_spectrum.Edd_continuum, color='green', alpha=0.7, lw=0.7, ls='--')

ax.plot(sun_spectrum.wave, sun_spectrum.Edd_flux, label="Sun (Teff=5772K, logg=4.438) [R=30k]", color='red', alpha=0.7, lw=0.7)
ax.plot(sun_spectrum.wave, sun_spectrum.Edd_continuum, color='red', alpha=0.7, lw=0.7, ls='--')
ax.set_yscale('log')
ax.set_ylim(1e-3, None)

ax.set_xlabel("Wavelength (Å)")
ax.set_ylabel(rf"Eddington flux, $H_\lambda$ ($\rm erg\,cm^{{-2}}\,s^{{-1}}\,\AA^{{-1}}$)")
ax.legend()
plt.show()


#---------------------------------------------
offsets = np.linspace(0, 1, 4)  # offsets for the normalized Edd_fluxes
fig, ax = plt.subplots(figsize=(13, 6))
ax.plot(hot_star_spectrum.wave, hot_star_spectrum.normalized_flux + offsets[0], label="Hot Star (Teff=29000K, logg=4.5) [R=50k]", color='blue', alpha=0.7, lw=0.7)
ax.plot(blap_spectrum.wave, blap_spectrum.normalized_flux + offsets[1], label="BLAP (Teff=29000K, logg=4.5, He-rich, N-rich, C/O/Si/Fe-poor) [R=50k]", color='orange', alpha=0.7, lw=0.7)
ax.plot(vega_spectrum.wave, vega_spectrum.normalized_flux + offsets[2], label="Vega (Teff=9550K, logg=4.1) [R=30k]", color='green', alpha=0.7, lw=0.7)
ax.plot(sun_spectrum.wave, sun_spectrum.normalized_flux + offsets[3], label="Sun (Teff=5772K, logg=4.438) [R=30k]", color='red', alpha=0.7, lw=0.7)

ax.set_xlabel("Wavelength (Å)")
ax.set_ylabel(rf"Normalized flux + offset [-]")
ax.legend()
plt.show()


if save:
    np.savetxt('hot_star_Teff+29000_logg+4.5_m+0.0_vmic+10_vsini+20_R+50000_air.txt', np.column_stack((hot_star_spectrum.wave, 4*np.pi*hot_star_spectrum.Edd_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)
    np.savetxt('norm_hot_star_Teff+29000_logg+4.5_m+0.0_vmic+10_vsini+20_R+50000_air.txt', np.column_stack((hot_star_spectrum.wave, hot_star_spectrum.normalized_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)

    np.savetxt('blap_Teff+29000_logg+4.5_m+0.0_vmic+10_vsini+20_R+50000_air.txt', np.column_stack((blap_spectrum.wave, 4*np.pi*blap_spectrum.Edd_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)
    np.savetxt('norm_blap_Teff+29000_logg+4.5_m+0.0_vmic+10_vsini+20_R+50000_air.txt', np.column_stack((blap_spectrum.wave, blap_spectrum.normalized_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)

    np.savetxt('vega_Teff+9550_logg+4.1_m-0.5_vmic+3_vsini+21_R+30000_air.txt', np.column_stack((vega_spectrum.wave, 4*np.pi*vega_spectrum.Edd_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)
    np.savetxt('norm_vega_Teff+9550_logg+4.1_m-0.5_vmic+3_vsini+21_R+30000_air.txt', np.column_stack((vega_spectrum.wave, vega_spectrum.normalized_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)

    np.savetxt('sun_Teff+5772_logg+4.438_m+0.0_vmic+1_vsini+2_R+30000_air.txt', np.column_stack((sun_spectrum.wave, 4*np.pi*sun_spectrum.Edd_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)
    np.savetxt('norm_sun_Teff+5772_logg+4.438_m+0.0_vmic+1_vsini+2_R+30000_air.txt', np.column_stack((sun_spectrum.wave, sun_spectrum.normalized_flux)), fmt='%.6e', header='wavelength_Angstrom    flux', comments=None)

    print("Spectra saved to text files.")

#---------------------------------------------
sys.exit()

flux_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA**-1

hot_star_F_tot = 4*np.pi*np.trapezoid(hot_star_spectrum.Edd_flux * flux_unit, hot_star_spectrum.wave * u.AA) # Edd_flux in the SYNTHE output files is the Eddington Edd_flux, H
blap_F_tot = 4*np.pi*np.trapezoid(blap_spectrum.Edd_flux * flux_unit, blap_spectrum.wave * u.AA)
vega_F_tot = 4*np.pi*np.trapezoid(vega_spectrum.Edd_flux * flux_unit, vega_spectrum.wave * u.AA)
sun_F_tot = 4*np.pi*np.trapezoid(sun_spectrum.Edd_flux * flux_unit, sun_spectrum.wave * u.AA)


hot_star_and_blap_theoretical_F_tot = (c.sigma_sb.cgs * (29e3*u.K)**4).to(flux_unit*u.AA)
vega_theoretical_F_tot = (c.sigma_sb.cgs * (9550*u.K)**4).to(flux_unit*u.AA)
sun_theoretical_F_tot = (c.sigma_sb.cgs * (5772*u.K)**4).to(flux_unit*u.AA)

print(f"hot star calculated (in {hot_star_spectrum.wave[0]:.0f} - {hot_star_spectrum.wave[-1]:.0f} Å) vs theoretical F_tot: {hot_star_F_tot:.6e} / {hot_star_and_blap_theoretical_F_tot:.6e}")
print(f"Vega calculated (in {vega_spectrum.wave[0]:.0f} - {vega_spectrum.wave[-1]:.0f} Å) vs theoretical F_tot: {vega_F_tot:.6e} / {vega_theoretical_F_tot:.6e}")
print(f"Sun calculated (in {sun_spectrum.wave[0]:.0f} - {sun_spectrum.wave[-1]:.0f} Å) vs theoretical F_tot: {sun_F_tot:.6e} / {sun_theoretical_F_tot:.6e}")