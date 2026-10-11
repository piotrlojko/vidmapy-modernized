#!/usr/bin/env python3

"""
Implements class that contains all spectrum data:
- wave, flux, contrinuum, normed_flux
- parameters of spectrum         if availible
- lines identification dataframe if availible

Units:
- wave                  : wavelength [A]
- flux, continuum       : Eddington flux per unit wavelength, H_lambda [erg cm^-2 s^-1 A^-1 sr^-1]
- flux_nu, continuum_nu : Eddington flux per unit frequency, H_nu [erg cm^-2 s^-1 Hz^-1 sr^-1],
                          as written by SYNTHE
- normed_flux           : flux / continuum [-]
The physical surface flux is F = 4 * pi * H.
"""


from vidmapy.kurucz.parameters import Parameters
from vidmapy.kurucz.utility_functions import _read_fixed_width_data
import numpy as np
import pandas as pd
import os

SPEED_OF_LIGHT = 2.99792458e18 # [A s^-1]


class Spectrum:
    def __init__(self):
        self.wave = None
        self.flux = None 
        self.continuum = None 
        self.normed_flux = None

        self.flux_nu = None
        self.continuum_nu = None

        self.parameters = None

        self.lines_identification = None

    @classmethod
    def from_synthe_spectrum(cls, parameters, directory, spectrum_file_name="br_spec.dat", linelist_file_name="line_list.dat"):
        sp = cls()
        sp.read_lines_identification_file(os.path.join(directory,linelist_file_name))
        sp.read_spectrum_file(os.path.join(directory,spectrum_file_name))
        sp.parameters = parameters
        return sp

    @classmethod
    def from_segments(cls, segments, parameters):
        # Join spectra of consecutive wavelength segments, [(spectrum, wave_min, wave_max), ...],
        # each contributing wave_min <= wave < wave_max (wave <= wave_max for the last one):
        def in_segment(wave, i, wave_min, wave_max):
            return (wave >= wave_min) & ((wave <= wave_max) if i == len(segments) - 1 else (wave < wave_max))
        sp = cls()
        for attribute in ["wave", "flux", "continuum", "normed_flux", "flux_nu", "continuum_nu"]:
            setattr(sp, attribute, np.concatenate([getattr(segment, attribute)[in_segment(segment.wave, i, wave_min, wave_max)]
                                                   for i, (segment, wave_min, wave_max) in enumerate(segments)]))
        sp.lines_identification = pd.concat([segment.lines_identification[in_segment(segment.lines_identification["wave"], i, wave_min, wave_max)]
                                             for i, (segment, wave_min, wave_max) in enumerate(segments)], ignore_index=True)
        sp.parameters = parameters
        return sp

    def read_lines_identification_file(self, path):
        format_string = "++++++++++^^^^^^^+++++^^^^^^^^^^^^+++++^^^^^^^^^^^^+++++++++^^^^^^^^+++++++"
        column_names = ["wave","loggf?","J1?","E1?","J2?","E2?","atom_symbol","literature","strength"]
        self.lines_identification = _read_fixed_width_data(path, format_string, column_names=column_names)
        self.lines_identification["wave"] *= 10.
    
    def read_spectrum_file(self, path):
        df = pd.read_csv(path, header=None, index_col=False, sep=r'\s+')
        self.wave = df[0].values
        # SYNTHE writes the Eddington flux per unit frequency, H_nu, so convert it
        # to the flux per unit wavelength, H_lambda = H_nu * c / lambda^2:
        self.flux_nu = df[1].values
        self.continuum_nu = df[2].values
        self.flux = self._per_unit_wavelength(self.flux_nu)
        self.continuum = self._per_unit_wavelength(self.continuum_nu)
        self.normed_flux = df[3].values

    def _per_unit_wavelength(self, flux_nu):
        return flux_nu * SPEED_OF_LIGHT / self.wave**2

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"({self.parameters.teff}, {self.parameters.logg}, {self.parameters.metallicity}, ({min(self.wave)}, {max(self.wave):.2f}))"

def main():
    pass

if __name__ == '__main__':
    main()
