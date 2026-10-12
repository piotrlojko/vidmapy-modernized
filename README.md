> [!IMPORTANT]
> ### This repository contains a modern fork of the original `VidmaPy` software by ***Tomasz Różański*** and ***Collaborators***. The code is functional and tested.

> [!CAUTION]
> ### This code is co-created by AI. If you really don't like AI-generated software, you should consider using the original one.

# VidmaPy - wrapper for Kurucz's ATLAS/SYNTHE codes

VidmaPy is a python module that enable user to use the power of ATLAS/SYNTHE codes with just few lines of python code.

## Key features
* Easy access to ATLAS and SYNTHE codes through Python
* Currently works only under LINUX based systems (because of compiled codes included)

### Example code
```python
from vidmapy.kurucz.atlas import Atlas
from vidmapy.kurucz.synthe import Synthe
from vidmapy.kurucz.parameters import Parameters
import matplotlib.pyplot as plt

atlas_worker = Atlas()
model = atlas_worker.get_model(Parameters(teff=8000., logg=4.0, metallicity=0.0, microturbulence=2.0))

synthe_worker = Synthe()
spectrum = synthe_worker.get_spectrum(model, Parameters(wave_min=4500., wave_max=5000., vsini=50.))

# Plot results
f, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(model.structure["RHOX"],model.structure["T"])
ax1.set_title('Atmosphere model')
ax1.set_xlabel("RHOX")
ax1.set_ylabel("T [K]")

ax2.plot(spectrum.wave, spectrum.normed_flux)
ax2.set_title('Spectrum')
ax2.set_xlabel("Wavelength [$\AA$]")
ax2.set_ylabel("Normalized Flux")

plt.show()
```
![Output plot of example code](docs/img/example_code_output.png)

### Units of the computed spectrum
* `spectrum.wave`: wavelength [Å]
* `spectrum.flux`, `spectrum.continuum`: Eddington flux per unit wavelength, $H_\lambda$ [erg cm<sup>-2</sup> s<sup>-1</sup> Å<sup>-1</sup> sr<sup>-1</sup>]; the surface flux is $F_\lambda = 4\pi H_\lambda$
* `spectrum.flux_nu`, `spectrum.continuum_nu`: Eddington flux per unit frequency, $H_\nu$ [erg cm<sup>-2</sup> s<sup>-1</sup> Hz<sup>-1</sup> sr<sup>-1</sup>], as written by SYNTHE, with $H_\lambda = H_\nu c/\lambda^2$
* `spectrum.normed_flux`: `flux / continuum`

### Notes on opacities, microturbulence and molecular lines
* ATLAS9 **cannot** rescale the line opacity of its opacity distribution functions (ODFs), so the `Atlas` code uses the scaled-solar ODF (`{p|m}XXbigN.bdf`) and Rosseland table (`kap{p|m}XX.ros`) from `atomic_data/ODF/NEW` that is closest to the requested `metallicity` and `microturbulence`, but warns when no exact match exists. **Changing individual abundances alters the continuum opacity and the equation of state but does not modify the ODF line opacity.**
* A model is taken directly from the grid only when all its parameters match the requested ones. This includes the `microturbulence` (given in km/s). The default grid was computed with the microturbulent velocity of **2 km/s**. `Synthe` computes the spectrum with exactly the requested `microturbulence`, whether the model comes from the grid or from ATLAS.
* By default `Synthe` uses only the atomic lines of `gfall08oct17.dat`. Specifying the `molecules` keyword parameter in `synthe_worker.get_spectrum(model, molecules=True)` adds every molecular line list in `atomic_data/molecules` (C<sub>2</sub>, CH, CN, CO, H<sub>2</sub>, MgH, NH, OH, SiH, SiO), which is important for cool stars, e.g. for the solar flux below about 4500 Å. **Be aware that setting `molecules=True` comes at the cost of a much longer computation.**

> [!CAUTION]
> ### Known issues
> * TiO and H<sub>2</sub>O lines are not included.
> * With molecular lines, `synthe.exe` crashes on large wavelength ranges. To avoid this, the code now computes the spectrum in segments spanning **at most a factor of 1.5 in wavelength**  and joins the resulting segments.
> * `synthe.exe` fails when the computed wavelength range starts in the far wing of Lyman α, between about 1220 and 1980 Å. With `molecules=True` the segments avoid this interval by default. **Otherwise `wave_min` should be set above this interval.**

## Getting Started

### Prerequisites

* Python3
* Conda - usefull but not necessary

### Download

Two steps:
* Clone the repository or download it as the .zip file:
  - Clone by:
    ```
    git clone https://github.com/piotrlojko/vidmapy-modernized.git
    ```
* Download and untar the `atom_data_kurucz.tar.xz` with necessary atomic data in `/vidmapy/kurucz/atomic_data/`:
  - Can be downloaded from : [atomic data + example spectra](https://drive.google.com/drive/folders/1arYoJCUaMMkw-8Fd1f93Pan1DwQLbdON?usp=sharing)
  - And extracted via:
    ```
    tar -xf atom_data_kurucz.tar.xz -C ~/vidmapy-modernized/vidmapy/kurucz/atomic_data/
    ```
  
### Install

* Move to the top-level vidmapy directory:
```
cd vidmapy-modernized
```
* Install VidmaPy in currently activated environment:
```
pip install .
```
or if you want to develop VidmaPy:
```
pip install -e .
```
  

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

## Acknowledgments
* Tomasz Różański
* Ewa Niemczura
