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
        self.wave, self.flux, self.continuum, self.normed_flux = np.loadtxt(path, unpack=True, dtype=float, comments="#")

hot_star_spectrum = spectrum("hot_star_spectrum.txt")
blap_spectrum = spectrum("blap_spectrum.txt")

### # Temp:
### py_kurucz_hot_star_spectrum_wave, py_kurucz_hot_star_spectrum_flux = np.loadtxt("t29000g4.75_mh+0.00_am+0.00_300_1300.spec", unpack=True, usecols=[0,1], dtype=float, comments="#")
### py_kurucz_hot_star_spectrum_wave *= 10  # convert from nm to Å
### py_kurucz_hot_star_spectrum_flux *= 10  # convert from erg/s/cm^2/nm to erg/s/cm^2/Å


fig, ax = plt.subplots(figsize=(13, 6))
ax.plot(hot_star_spectrum.wave, hot_star_spectrum.flux, label="Hot Star (Teff=29000K, logg=4.5)", color='blue', alpha=0.7)
ax.plot(blap_spectrum.wave, blap_spectrum.flux, label="BLAP (Teff=29000K, logg=4.5, He-rich, N-rich, C/O/Si/Fe-poor)", color='orange', alpha=0.7)
# Temp:
# ax.plot(py_kurucz_hot_star_spectrum_wave, py_kurucz_hot_star_spectrum_flux, label="Kurucz (Teff=29000K, logg=4.75)", color='darkgreen', alpha=0.7)

ax.set_xlabel("Wavelength (Å)")
ax.set_ylabel(rf"Eddington Flux, $H_\lambda$ ($\rm erg\,cm^{{-2}}\,s^{{-1}}\,\AA^{{-1}}$)")
ax.legend()
ax.set_title("Spectral Comparison")
plt.show()

flux_unit = u.erg * u.s**-1 * u.cm**-2 * u.AA**-1

hot_star_F_tot = 4*np.pi*np.trapezoid(hot_star_spectrum.flux * flux_unit, hot_star_spectrum.wave * u.AA) # flux in the SYNTHE output files is the Eddington flux, H
blap_F_tot = 4*np.pi*np.trapezoid(blap_spectrum.flux * flux_unit, blap_spectrum.wave * u.AA)

theoretical_F_tot = (c.sigma_sb.cgs * (29e3*u.K)**4).to(flux_unit*u.AA)

print(hot_star_F_tot, blap_F_tot)
print(f"Theoretical F_tot: {theoretical_F_tot:.6e}")

sys.exit()


ratio_hot_star = hot_star_F_tot / theoretical_F_tot.value
ratio_blap = blap_F_tot / theoretical_F_tot.value
print(f"Hot star: F_tot / (sigma_sb * Teff^4) = {ratio_hot_star:.6f}")
print(f"BLAP: F_tot / (sigma_sb * Teff^4) = {ratio_blap:.6f}")