#!/usr/bin/env python3

from vidmapy.kurucz.spectrum import Spectrum, SPEED_OF_LIGHT
# from vidmapy.kurucz.parameters import Parameters
import numpy as np
import pytest

PLANCK = 6.62607015e-27            # [erg s]
BOLTZMANN = 1.380649e-16           # [erg K^-1]
STEFAN_BOLTZMANN = 5.670374419e-5  # [erg cm^-2 s^-1 K^-4]
WIEN = 2.897771955e7               # [A K]
X_PEAK_NU = 2.821439372            # h nu_max / k T of B_nu

def _write_blackbody_synthe_spectrum(path, teff, wave_min=1.e3, wave_max=9.e5, resolution=2000.):
    # Same layout as br_spec.dat written by syntoascanga, FORMAT(F11.4,2E20.8,F13.6):
    # wavelength [A], H_nu and continuum H_nu [erg cm^-2 s^-1 Hz^-1 sr^-1], normalized flux:
    ratio = 1. + 1./resolution
    wave = wave_min * ratio**np.arange(int(np.log(wave_max/wave_min)/np.log(ratio)) + 1)
    nu = SPEED_OF_LIGHT / wave
    planck_nu = 2.*PLANCK*nu**3/(SPEED_OF_LIGHT*1.e-8)**2 / np.expm1(PLANCK*nu/(BOLTZMANN*teff))
    continuum_nu = planck_nu / 4. # Eddington flux of a blackbody
    normed_flux = 1. - 0.5*np.exp(-0.5*((wave - 6562.8)/2.)**2)
    with open(path, 'w') as f:
        for row in zip(wave, normed_flux*continuum_nu, continuum_nu, normed_flux):
            f.write("{:11.4f}{:20.8E}{:20.8E}{:13.6f}\n".format(*row))

@pytest.fixture
def blackbody_spectrum(tmp_path):
    teff = 5772.
    path = tmp_path / "br_spec.dat"
    _write_blackbody_synthe_spectrum(path, teff)
    spectrum = Spectrum()
    spectrum.read_spectrum_file(path)
    return spectrum, teff, path

def test_create_spectrum_object():
    spectrum = Spectrum()
    assert True

def test_flux_peak_follows_wien_law(blackbody_spectrum):
    spectrum, teff, _ = blackbody_spectrum
    assert spectrum.wave[np.argmax(spectrum.continuum)] == pytest.approx(WIEN/teff, rel=1e-3)

def test_flux_integrates_to_bolometric_flux(blackbody_spectrum):
    spectrum, teff, _ = blackbody_spectrum
    bolometric_flux = 4.*np.pi*np.trapezoid(spectrum.continuum, spectrum.wave)
    assert bolometric_flux == pytest.approx(STEFAN_BOLTZMANN*teff**4, rel=1e-3)

def test_flux_per_unit_frequency_is_kept(blackbody_spectrum):
    spectrum, teff, path = blackbody_spectrum
    data = np.loadtxt(path)
    np.testing.assert_array_equal(spectrum.flux_nu, data[:,1])
    np.testing.assert_array_equal(spectrum.continuum_nu, data[:,2])
    # H_nu peaks ~1.76 times redward of H_lambda:
    wave_peak_nu = PLANCK*SPEED_OF_LIGHT/(X_PEAK_NU*BOLTZMANN*teff)
    assert spectrum.wave[np.argmax(spectrum.continuum_nu)] == pytest.approx(wave_peak_nu, rel=1e-3)

def test_normed_flux_is_unchanged(blackbody_spectrum):
    spectrum, _, path = blackbody_spectrum
    np.testing.assert_array_equal(spectrum.normed_flux, np.loadtxt(path)[:,3])
    np.testing.assert_allclose(spectrum.flux/spectrum.continuum, spectrum.normed_flux, atol=1e-6)

def main():
    pass

if __name__ == '__main__':
    main()
