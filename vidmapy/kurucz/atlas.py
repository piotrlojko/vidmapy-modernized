#!/usr/bin/env python3

"""
Implements easy-to-use interface to Kurucz's ATLAS code

Needs to be use together with class Parameters.
eg.
    p = Parameters() # All default

    wa = Atlas()
    m = wa.get_model(p)
"""

from vidmapy.kurucz.model import Model
from vidmapy.kurucz.models_grid import Grid
from vidmapy.kurucz.model_definition import ModelDefinition
from vidmapy.kurucz.parameters import Parameters
import os
import re
import subprocess
import tempfile
import copy
import stat
import warnings

class Atlas:
    def __init__(self):
        self.parameters = Parameters()
        self._kurucz_directory = os.path.dirname(os.path.abspath(__file__))
        self.default_grid = os.path.join(self._kurucz_directory,"grids","ap00k2odfnew.txt")
        self.set_grid(path=self.default_grid)
        
        self._kurucz_bin_path = self._kurucz_directory
        self._atomic_data_path = os.path.join(self._kurucz_directory,"atomic_data")
        # self._kurucz_bin_path = "/usr/local/kurucz/"
        self._atlas_name = "atlas9mem_newodf.exe"

    def get_report(self):
        return f"TEFF  = {self.parameters.teff:5.0f}\n"\
               f"LOGG  = {self.parameters.logg:4.2f}\n"\
               f"[M/H] = {self.parameters.metallicity:4.2f}"

    @property
    def teff(self):
        return self.parameters.teff
        
    @property
    def logg(self):
        return self.parameters.logg

    @property
    def metallicity(self):
        return self.parameters.metallicity  

    @teff.setter
    def teff(self, teff):
        self.parameters.teff = teff

    @logg.setter
    def logg(self, logg):
        self.parameters.logg = logg

    @metallicity.setter
    def metallicity(self, metallicity):
        self.parameters.metallicity = metallicity

    def set_grid(self, path):
        self._grid = Grid(path)

    def get_model(self, parameters):
        self.parameters = copy.deepcopy(parameters)
        initial_model = self._get_initial_model()
        if initial_model.is_model_has_this(self.parameters):
            return initial_model
        else:
            return self._create_temporary_directory_and_run_ATLAS(initial_model)

    def _get_initial_model(self):
        model_as_string = self._grid.get_best_model(self.parameters)
        return Model.from_string(model_as_string, self.parameters)

    def _create_temporary_directory_and_run_ATLAS(self, initial_model):
        with tempfile.TemporaryDirectory(prefix="atlas_") as tmpdirname:
            model = self._compute_model(tmpdirname, initial_model)
        return model

    def _compute_model(self, tmpdirname, initial_model):
        md = ModelDefinition()
        odf_path, rosseland_path = self._get_opacity_files(self.parameters)
        # Link necessary data:
        os.symlink(rosseland_path, os.path.join(tmpdirname,"fort.1"))
        os.symlink(odf_path, os.path.join(tmpdirname,"fort.9"))
        os.symlink(os.path.join(self._atomic_data_path,"lines","molecules.dat"), os.path.join(tmpdirname,"fort.2"))
        # Create intial model:
        initial_model.save_model(os.path.join(tmpdirname, "fort.3"))
        
        # Run ATLAS
        process = subprocess.Popen([os.path.join(self._kurucz_bin_path, "bin", self._atlas_name)],
                                    stdin=subprocess.PIPE,
                                    stdout=subprocess.PIPE,
                                    encoding='utf8',
                                    cwd=tmpdirname)
        outs, errs = process.communicate(md(self.parameters))
        process.wait()
        
        return Model.with_extended_set_of_parameters(os.path.join(tmpdirname,"fort.7"), self.parameters)

    def _get_opacity_files(self, parameters):
        # ATLAS9 cannot rescale the line opacity of an ODF, so use the scaled-solar ODF
        # ({p|m}<10*|[M/H]|>big<vturb>.bdf) closest in metallicity, then in microturbulence,
        # and the Rosseland opacity table of the same metallicity:
        odf_directory = os.path.join(self._atomic_data_path, "ODF", "NEW")
        odfs = []
        for file_name in sorted(os.listdir(odf_directory)):
            match = re.fullmatch(r"([pm])(\d\d)big(\d+)\.bdf", file_name)
            if match:
                metallicity = (1 if match.group(1) == "p" else -1) * int(match.group(2)) / 10.
                odfs.append((metallicity, float(match.group(3)), file_name))
        if not odfs:
            raise FileNotFoundError(f"No opacity distribution functions (*big*.bdf) in {odf_directory}")

        metallicity = min(sorted({odf[0] for odf in odfs}), key=lambda m: abs(m - parameters.metallicity))
        _, microturbulence, odf_name = min([odf for odf in odfs if odf[0] == metallicity],
                                           key=lambda odf: abs(odf[1] - parameters.microturbulence))
        if abs(metallicity - parameters.metallicity) > 1e-3 or abs(microturbulence - parameters.microturbulence) > 1e-3:
            warnings.warn(f"No ODF for [M/H] = {parameters.metallicity:+.2f}, vturb = {parameters.microturbulence:g} km/s "
                          f"in {odf_directory}; ATLAS uses {odf_name} ([M/H] = {metallicity:+.1f}, vturb = {microturbulence:g} km/s)")

        rosseland_path = os.path.join(odf_directory, f"kap{odf_name[:3]}.ros")
        if not os.path.exists(rosseland_path):
            raise FileNotFoundError(f"No Rosseland opacity table {rosseland_path} for {odf_name}")
        return os.path.join(odf_directory, odf_name), rosseland_path

def main():
    pass

if __name__ == '__main__':
    main()