import os, sys
import numpy as np
import matplotlib.pyplot as plt
import copy

from vidmapy.kurucz.atlas import Atlas
from vidmapy.kurucz.synthe import Synthe
from vidmapy.kurucz.parameters import Parameters


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



vega_params = Parameters(teff=9550,
                        logg=4.1,
                        metallicity=-0.5,
                        microturbulence=3,
                        vsini=21,
                        resolution=30000,
                        wave_min=500,
                        wave_max=13000)


atlas_worker= Atlas()
vega_atmosphere_model = atlas_worker.get_model(vega_params)
print(vega_atmosphere_model) # controlled by the __str__() special method
synthe_worker = Synthe()
vega_spectrum = synthe_worker.get_spectrum(vega_atmosphere_model, quiet=False)



fig, ax = plt.subplots(figsize=(13, 6))
ax.plot(vega_spectrum.wave, vega_spectrum.flux, label="vega (Teff=9550K, logg=4.1)", color='darkorange', alpha=0.7)
ax.set_xlabel("Wavelength (Å)")
ax.set_ylabel(rf"Eddington Flux, $H_\lambda$ ($\rm erg\,cm^{{-2}}\,s^{{-1}}\,\AA^{{-1}}$)")
ax.legend()
ax.set_title("Spectral Comparison")
plt.show()


# Save:
#----------------------------------------------------------------------------------------------------------------------------------------------
vega_data = np.column_stack((vega_spectrum.wave, vega_spectrum.flux, vega_spectrum.continuum, vega_spectrum.normed_flux))
vega_metadata_and_header = f"""# Spectrum generated from the Kurucz ATLAS and SYNTHE codes via the VidmaPy package.
# Parameters: Teff={vega_params.teff} [K]
# logg={vega_params.logg} [dex]
# [M/H]={vega_params.metallicity} [dex]
# microturbulence={vega_params.microturbulence} [km/s]
# vsini={vega_params.vsini} [km/s]
# spectral resolution={vega_params.resolution} [km/s]
# wavelength range=({vega_params.wave_min}-{vega_params.wave_max}) [Å]
#
# wavelength [Å]    Eddington flux [erg/cm^2/s/Å]    continuum [erg/cm^2/s/Å]    normalized_flux [-]
"""
np.savetxt("vega_spectrum.txt", vega_data, delimiter="    ", header=vega_metadata_and_header, comments='')


print("Spectra saved to 'vega_spectrum.txt'.")
sys.exit(0)



"""
_reference_composition = [
    [  1, 0.9204 , 'H' , 'Hydrogen'],
    [  2, 0.07834, 'He', 'Helium'],
    [  3, -10.94 , 'Li', 'Lithium'],
    [  4, -10.64 , 'Be', 'Beryllium'],
    [  5, -9.49  , 'B' , 'Boron'],
    [  6, -3.52  , 'C' , 'Carbon'],
    [  7, -4.12  , 'N' , 'Nitrogen'],
    [  8, -3.21  , 'O' , 'Oxygen'],
    [  9, -7.48  , 'F' , 'Fluorine'],
    [ 10, -3.96  , 'Ne', 'Neon'],
    [ 11, -5.71  , 'Na', 'Sodium'],
    [ 12, -4.46  , 'Mg', 'Magnesium'],
    [ 13, -5.57  , 'Al', 'Aluminum'],
    [ 14, -4.49  , 'Si', 'Silicon'],
    [ 15, -6.59  , 'P' , 'Phosphorus'],
    [ 16, -4.71  , 'S' , 'Sulfur'],
    [ 17, -6.54  , 'Cl', 'Chlorine'],
    [ 18, -5.64  , 'Ar', 'Argon'],
    [ 19, -6.92  , 'K' , 'Potassium'],
    [ 20, -5.68  , 'Ca', 'Calcium'],
    [ 21, -8.87  , 'Sc', 'Scandium'],
    [ 22, -7.02  , 'Ti', 'Titanium'],
    [ 23, -8.04  , 'V' , 'Vanadium'],
    [ 24, -6.37  , 'Cr', 'Chromium'],
    [ 25, -6.65  , 'Mn', 'Manganese'],
    [ 26, -4.54  , 'Fe', 'Iron'],
    [ 27, -7.12  , 'Co', 'Cobalt'],
    [ 28, -5.79  , 'Ni', 'Nickel'],
    [ 29, -7.83  , 'Cu', 'Copper'],
    [ 30, -7.44  , 'Zn', 'Zinc'],
    [ 31, -9.16  , 'Ga', 'Gallium'],
    [ 32, -8.63  , 'Ge', 'Germanium'],
    [ 33, -9.67  , 'As', 'Arsenic'],
    [ 34, -8.63  , 'Se', 'Selenium'],
    [ 35, -9.41  , 'Br', 'Bromine'],
    [ 36, -8.73  , 'Kr', 'Krypton'],
    [ 37, -9.44  , 'Rb', 'Rubidium'],
    [ 38, -9.07  , 'Sr', 'Strontium'],
    [ 39, -9.8   , 'Y' , 'Yttrium'],
    [ 40, -9.44  , 'Zr', 'Zirconium'],
    [ 41, -10.62 , 'Nb', 'Niobium'],
    [ 42, -10.12 , 'Mo', 'Molybdenum'],
    [ 43, -20.0  , 'Tc', 'Technetium'],
    [ 44, -10.2  , 'Ru', 'Ruthenium'],
    [ 45, -10.92 , 'Rh', 'Rhodium'],
    [ 46, -10.35 , 'Pd', 'Palladium'],
    [ 47, -11.1  , 'Ag', 'Silver'],
    [ 48, -10.27 , 'Cd', 'Cadmium'],
    [ 49, -10.38 , 'In', 'Indium'],
    [ 50, -10.04 , 'Sn', 'Tin'],
    [ 51, -11.04 , 'Sb', 'Antimony'],
    [ 52, -9.8   , 'Te', 'Tellurium'],
    [ 53, -10.53 , 'I' , 'Iodine'],
    [ 54, -9.87  , 'Xe', 'Xenon'],
    [ 55, -10.91 , 'Cs', 'Cesium'],
    [ 56, -9.91  , 'Ba', 'Barium'],
    [ 57, -10.87 , 'La', 'Lanthanum'],
    [ 58, -10.46 , 'Ce', 'Cerium'],
    [ 59, -11.33 , 'Pr', 'Praseodymium'],
    [ 60, -10.54 , 'Nd', 'Neodymium'],
    [ 61, -20.0  , 'Pm', 'Promethium'],
    [ 62, -11.03 , 'Sm', 'Samarium'],
    [ 63, -11.53 , 'Eu', 'Europium'],
    [ 64, -10.92 , 'Gd', 'Gadolinium'],
    [ 65, -11.69 , 'Tb', 'Terbium'],
    [ 66, -10.9  , 'Dy', 'Dysprosium'],
    [ 67, -11.78 , 'Ho', 'Holmium'],
    [ 68, -11.11 , 'Er', 'Erbium'],
    [ 69, -12.04 , 'Tm', 'Thulium'],
    [ 70, -10.96 , 'Yb', 'Ytterbium'],
    [ 71, -11.98 , 'Lu', 'Lutetium'],
    [ 72, -11.16 , 'Hf', 'Hafnium'],
    [ 73, -12.17 , 'Ta', 'Tantalum'],
    [ 74, -10.93 , 'W' , 'Tungsten'],
    [ 75, -11.76 , 'Re', 'Rhenium'],
    [ 76, -10.59 , 'Os', 'Osmium'],
    [ 77, -10.69 , 'Ir', 'Iridium'],
    [ 78, -10.24 , 'Pt', 'Platinum'],
    [ 79, -11.03 , 'Au', 'Gold'],
    [ 80, -10.91 , 'Hg', 'Mercury'],
    [ 81, -11.14 , 'Tl', 'Thallium'],
    [ 82, -10.09 , 'Pb', 'Lead'],
    [ 83, -11.33 , 'Bi', 'Bismuth'],
    [ 84, -20.0  , 'Po', 'Polonium'],
    [ 85, -20.0  , 'At', 'Astatine'],
    [ 86, -20.0  , 'Rn', 'Radon'],
    [ 87, -20.0  , 'Fr', 'Francium'],
    [ 88, -20.0  , 'Ra', 'Radium'],
    [ 89, -20.0  , 'Ac', 'Actinium'],
    [ 90, -11.95 , 'Th', 'Thorium'],
    [ 91, -20.0  , 'Pa', 'Protactinium'],
    [ 92, -12.54 , 'U' , 'Uranium'],
    [ 93, -20.0  , 'Np', 'Neptunium'],
    [ 94, -20.0  , 'Pu', 'Plutonium'],
    [ 95, -20.0  , 'Am', 'Americium'],
    [ 96, -20.0  , 'Cm', 'Curium'],
    [ 97, -20.0  , 'Bk', 'Berkelium'],
    [ 98, -20.0  , 'Cf', 'Californium'],
    [ 99, -20.0  , 'Es', 'Einsteinium'],
    [100, None   , 'Fm', 'Fermium'],
    [101, None   , 'Md', 'Mendelevium'],
    [102, None   , 'No', 'Nobelium'],
    [103, None   , 'Lr', 'Lawrencium'],
    [104, None   , 'Rf', 'Rutherfordium'],
    [105, None   , 'Db', 'Dubnium'],
    [106, None   , 'Sg', 'Seaborgium'],
    [107, None   , 'Bh', 'Bohrium'],
    [108, None   , 'Hs', 'Hassium'],
    [109, None   , 'Mt', 'Meitnerium'],
    [110, None   , 'Ds', 'Darmstadtium'],
    [111, None   , 'Rg', 'Roentgenium'],
    [112, None   , 'Cn', 'Copernicium'],
    [113, None   , 'Nh', 'Nihonium'],
    [114, None   , 'Fl', 'Flerovium'],
    [115, None   , 'Mc', 'Moscovium'],
    [116, None   , 'Lv', 'Livermorium'],
    [117, None   , 'Ts', 'Tennessine'],
    [118, None   , 'Og', 'Oganesson']
    ]
"""
