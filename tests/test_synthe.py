#!/usr/bin/env python3

from vidmapy.kurucz.synthe import Synthe
from vidmapy.kurucz.parameters import Parameters
"""
"""

# TEST INTERFACE
def test_create_synthe_object():
    worker = Synthe()
    assert True 

# TEST IMPLEMENTATION
import os
from unittest.mock import MagicMock
from vidmapy.kurucz import synthe
from vidmapy.kurucz.model import Model
from vidmapy.kurucz.models_grid import Grid

def grid_model():
    # vturb = 2 km/s model from the default grid:
    path = os.path.join(os.path.dirname(synthe.__file__), "grids", "ap00k2odfnew.txt")
    return Model.from_string(Grid(path).get_best_model(Parameters(teff=5750, logg=4.5)), Parameters())

def test_synbeg_input_has_no_microturbulence():
    worker = Synthe()
    worker.model = grid_model()
    # AIRorVAC WLBEG WLEND RESOLU TURBV ...:
    assert worker._get_synbeg_input(worker.model).split()[4] == "0.00"

def test_microturbulence_passed_through_model():
    worker = Synthe()
    worker._create_temp_direcotry_and_run_SYNTHE = MagicMock(return_value="spectrum")
    model = grid_model()

    worker.get_spectrum(model, Parameters(microturbulence=3.0, wave_min=5000, wave_max=5100), quiet=True)
    synthe_model, molecules = worker._create_temp_direcotry_and_run_SYNTHE.call_args.args
    assert synthe_model.parameters.microturbulence == 3.0
    assert all(synthe_model.structure["VTURB"] == 3.0e5)
    assert molecules is False
    assert all(model.structure["VTURB"] == 2.0e5)

def test_molecular_lines_read_on_request(monkeypatch):
    worker = Synthe()
    for step in ["_run_xnfpelsyn", "_run_synbeg", "_run_rline2", "_run_synthe", "_run_spectrv", "_run_rotate", "_run_broaden", "_run_syntoascanga"]:
        monkeypatch.setattr(worker, step, lambda *args, **kwargs: None)
    read_molecules = []
    monkeypatch.setattr(worker, "_run_rmolecasc", lambda tmpdirname, model, molecule_file: read_molecules.append(molecule_file))
    monkeypatch.setattr(worker, "_get_molecule_files", lambda: ["chax.dat", "cnax.dat"])
    monkeypatch.setattr(synthe.Spectrum, "from_synthe_spectrum", lambda *args, **kwargs: "spectrum")
    worker.model = MagicMock()

    worker._compute_spectrum("tmp", worker.model, molecules=False)
    assert read_molecules == []
    worker._compute_spectrum("tmp", worker.model, molecules=True)
    assert read_molecules == ["chax.dat", "cnax.dat"]

def test_molecule_files(tmp_path):
    (tmp_path / "molecules").mkdir()
    for name in ["cnax.dat", "chax.dat", "readme.txt"]:
        (tmp_path / "molecules" / name).touch()
    worker = Synthe()
    worker._atomic_data_path = str(tmp_path)
    assert worker._get_molecule_files() == ["chax.dat", "cnax.dat"]


def main():
    pass

if __name__ == '__main__':
    main()