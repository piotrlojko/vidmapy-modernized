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
import numpy as np
import pandas as pd
import pytest
from unittest.mock import MagicMock
from vidmapy.kurucz.spectrum import Spectrum
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

def fake_synthe_run(calls):
    # Spectrum on the logarithmic grid of SYNTHE, (1 + 1/600000)**i A, over the range of the model:
    def run(model, molecules=False):
        calls.append((model.parameters.wave_min, model.parameters.wave_max, molecules))
        step = np.log1p(1/600000)
        i = np.arange(np.ceil(np.log(model.parameters.wave_min)/step), np.floor(np.log(model.parameters.wave_max)/step) + 1)
        spectrum = Spectrum()
        spectrum.wave = np.exp(i*step)
        for attribute in ["flux", "continuum", "normed_flux", "flux_nu", "continuum_nu"]:
            setattr(spectrum, attribute, np.sqrt(spectrum.wave))
        spectrum.lines_identification = pd.DataFrame({"wave": spectrum.wave[::5000]})
        return spectrum
    return run

def test_spectrum_with_molecules_computed_in_segments():
    worker = Synthe()
    calls = []
    worker._create_temp_direcotry_and_run_SYNTHE = fake_synthe_run(calls)
    model = grid_model()
    model.parameters.wave_min, model.parameters.wave_max = 4000., 9000.

    spectrum = worker.get_spectrum(model, molecules=True)
    # 9000/4000 = 1.5**2, i.e. two segments, each extended by the margin:
    assert [c[2] for c in calls] == [True, True]
    assert calls[0][0] == pytest.approx(4000.*(1 - synthe.SEGMENT_MARGIN))
    assert calls[0][1] == pytest.approx(6000.*(1 + synthe.SEGMENT_MARGIN))
    assert calls[1][0] == pytest.approx(6000.*(1 - synthe.SEGMENT_MARGIN))
    assert calls[1][1] == pytest.approx(9000.*(1 + synthe.SEGMENT_MARGIN))
    # Joined spectrum covers the requested range without gaps or repeated points:
    step = np.diff(np.log(spectrum.wave))
    assert spectrum.wave[0] >= 4000. and spectrum.wave[-1] <= 9000.
    assert spectrum.wave[0]/4000. < 1 + 2e-6 and 9000./spectrum.wave[-1] < 1 + 2e-6
    assert np.allclose(step, np.log1p(1/600000), rtol=1e-6)
    assert np.array_equal(spectrum.flux, np.sqrt(spectrum.wave))
    assert spectrum.lines_identification["wave"].is_monotonic_increasing
    assert spectrum.lines_identification["wave"].between(4000., 9000.).all()
    assert spectrum.parameters.wave_min == 4000. and spectrum.parameters.wave_max == 9000.

@pytest.mark.parametrize("wave_min,wave_max", [(500., 13000.), (1150., 10000.), (1100., 1500.), (1500., 1800.),
                                               (1990., 3000.), (2000., 10000.), (4000., 9000.)])
def test_segments_avoid_lyman_alpha_wing(wave_min, wave_max):
    worker = Synthe()
    segments = worker._get_segments(wave_min, wave_max)
    assert segments[0][0] == wave_min
    assert segments[-1][1] == pytest.approx(wave_max)
    assert all(a[1] == pytest.approx(b[0]) and a[0] < a[1] for a, b in zip(segments[:-1], segments[1:]))
    for segment_min, segment_max in segments:
        computed_min, computed_max = worker._get_computed_range(segment_min, segment_max)
        assert not synthe.LYMAN_ALPHA_WING[0] < computed_min < synthe.LYMAN_ALPHA_WING[1]
        assert computed_min < segment_min and computed_max > segment_max
        # Longest range that synthe.exe completed with all molecular lines was a factor 1.85:
        assert computed_max/computed_min < 1.8

def test_failing_program_raises():
    with pytest.raises(RuntimeError, match="false failed"):
        Synthe()._call_external_code(".", "/bin/false")

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