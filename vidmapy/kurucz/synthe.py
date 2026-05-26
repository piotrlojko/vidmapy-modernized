#!/usr/bin/env python3

"""
Implements easy-to-use interface to Kurucz's SYNTHE code

Needs to be use together with class Atlas and Parameters.
eg.
    p = Parameters() # All default

    wa = Atlas()
    m = wa.get_model(p)
    
    ws = Synthe()
    spectrum = ws.get_spectrum(m)
"""

from vidmapy.kurucz.spectrum import Spectrum
from vidmapy.kurucz.parameters import Parameters 

import glob
import hashlib
import os
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
import copy

class Synthe:
    def __init__(self):

        self._kurucz_directory = os.path.dirname(os.path.abspath(__file__))
        # self._kurucz_bin_path = "/usr/local/kurucz/"
        self._kurucz_bin_path = os.path.join(self._kurucz_directory, 'bin')
        self._atomic_data_path = os.path.join(self._kurucz_directory, "atomic_data")
        self._line_list_filename = "gfall08oct17.dat"
        self._line_list_path = os.path.join(self._atomic_data_path, "lines", self._line_list_filename)
        self._pre_synthe_cache_version = "v1"

    def get_spectrum(self, model, parameters=None, quiet=False, work_dir=None, reuse_workdir=False, cache_dir=None):
        self.model = copy.deepcopy(model)
        if parameters is not None:
            self.model.parameters.get_synthe_parameters(parameters,quiet=quiet)
        spectrum = self._create_temp_direcotry_and_run_SYNTHE(
            self.model,
            work_dir=work_dir,
            reuse_workdir=reuse_workdir,
            cache_dir=cache_dir
        )
        return spectrum
    
    def prepare_work_dir(self, work_dir):
        os.makedirs(work_dir, exist_ok=True)
        self._ensure_atomic_symlinks(work_dir)

    def _create_temp_direcotry_and_run_SYNTHE(self, model, work_dir=None, reuse_workdir=False, cache_dir=None):
        if work_dir is None or not reuse_workdir:
            temp_kwargs = {"prefix": "synthe_"}
            if work_dir is not None:
                temp_kwargs["dir"] = work_dir
            with tempfile.TemporaryDirectory(**temp_kwargs) as tmpdirname:
                spectrum = self._compute_spectrum(tmpdirname, model, cache_dir=cache_dir, reuse_workdir=False)
            return spectrum
        os.makedirs(work_dir, exist_ok=True)
        return self._compute_spectrum(work_dir, model, cache_dir=cache_dir, reuse_workdir=True)

    def _compute_spectrum(self, tmpdirname, model, cache_dir=None, reuse_workdir=False):
        if reuse_workdir:
            self._cleanup_workdir(tmpdirname)
        self._ensure_atomic_symlinks(tmpdirname)
        cache_key = None
        cache_hit = False
        if cache_dir is not None:
            cache_key = self._get_pre_synthe_cache_key(model)
            cache_hit = self._restore_pre_synthe_cache(tmpdirname, cache_dir, cache_key)
        if not cache_hit:
            self._run_xnfpelsyn(tmpdirname, model)
            self._run_synbeg(tmpdirname, model)
            self._run_rline2(tmpdirname, model) #  rgfalllinesnew.for  ?
            if cache_dir is not None:
                self._store_pre_synthe_cache(tmpdirname, cache_dir, cache_key)
        # self._run_rmolecasc(tmpdirname, model) # Optionally! Include molecular lines
        self._run_synthe(tmpdirname, model)
        self._run_spectrv(tmpdirname, model)
        self._run_rotate(tmpdirname, model)
        self._run_broaden(tmpdirname, model)
        self._run_syntoascanga(tmpdirname, model) #converfsynnmtoa.exe ?
        if reuse_workdir:
            self._ensure_atomic_symlinks(tmpdirname)

        return Spectrum.from_synthe_spectrum(self.model.parameters, tmpdirname)

    def _run_xnfpelsyn(self, tmpdirname, model):
        # Prepare lines data
        self._ensure_symlink(
            os.path.join(self._atomic_data_path,"lines","he1tables.dat"),
            os.path.join(tmpdirname,"fort.18")
        )
        self._ensure_symlink(
            os.path.join(self._atomic_data_path,"lines","molecules.dat"),
            os.path.join(tmpdirname,"fort.2")
        )
        self._ensure_symlink(
            os.path.join(self._atomic_data_path,"lines","continua.dat"),
            os.path.join(tmpdirname,"fort.17")
        )

        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "xnfpelsyn.exe"), 
                                self._extend_model_for_SYNTHE(model)
                                )

    def _run_synbeg(self, tmpdirname, model):
        self._call_external_code(tmpdirname, 
                                os.path.join(self._kurucz_bin_path, "synbeg.exe"), 
                                self._get_synbeg_input(model)
                                )

    def _run_rline2(self, tmpdirname, model):
        # os.symlink(os.path.join(self._atomic_data_path,"lines","gf0600.100"), os.path.join(tmpdirname,"fort.11"))
        os.symlink(self._line_list_path, os.path.join(tmpdirname,"fort.11"))
        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "rline2.exe")
                                )
        # self._call_external_code(tmpdirname,
        #                         os.path.join(self._kurucz_bin_path, "rgfalllinesnew.exe")
        #                         )                               
        os.remove(os.path.join(tmpdirname,"fort.11"))

    def _run_rmolecasc(self, tmpdirname, model, molecule_file="coax.dat"):
        os.symlink(os.path.join(self._atomic_data_path, "molecules", molecule_file), os.path.join(tmpdirname, "fort.11"))
        self._call_external_code(tmpdirname,
                        os.path.join(self._kurucz_bin_path, "rmolecasc.exe")
                        )
        os.remove(os.path.join(tmpdirname,"fort.11"))
    

    def _run_synthe(self, tmpdirname, model):
        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "synthe.exe")
                                )

    def _run_spectrv(self, tmpdirname, model):
        self._save_to_file(os.path.join(tmpdirname,"fort.5"), self._extend_model_for_SYNTHE(model))
        self._save_to_file(os.path.join(tmpdirname,"fort.25"), self._get_spectrv_string())
        
        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "spectrv.exe")
                                )      
    
        os.rename(os.path.join(tmpdirname,"fort.7"), os.path.join(tmpdirname,"spec.bin"))
        os.symlink(os.path.join(tmpdirname,"spec.bin"), os.path.join(tmpdirname,"fort.1"))        

    def _run_rotate(self, tmpdirname, model):
        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "rotate.exe"),
                                self._get_rotate_string()
                                )
        
        os.rename(os.path.join(tmpdirname,"ROT1"), os.path.join(tmpdirname,"spec.bin"))
        os.symlink(os.path.join(tmpdirname,"spec.bin"), os.path.join(tmpdirname,"fort.21"))
        os.remove(os.path.join(tmpdirname,"fort.1"))
        os.remove(os.path.join(tmpdirname,"fort.5"))

    def _run_broaden(self, tmpdirname, model):
        os.symlink(os.path.join(tmpdirname,"br_spec.bin"), os.path.join(tmpdirname,"fort.22"))
        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "broaden.exe"),
                                self._get_broaden_string()
                                )

    def _run_syntoascanga(self, tmpdirname, model):
        for f in glob.glob(os.path.join(tmpdirname,"fort.*")):
            if os.path.basename(f) in {"fort.17", "fort.18"}:
                continue
            os.remove(f)
        self._ensure_symlink(os.path.join(tmpdirname,"br_spec.bin"), os.path.join(tmpdirname,"fort.1"))
        self._ensure_symlink(os.path.join(tmpdirname,"br_spec.dat"), os.path.join(tmpdirname,"fort.2"))
        self._ensure_symlink(os.path.join(tmpdirname,"line_list.dat"), os.path.join(tmpdirname,"fort.3"))

        self._call_external_code(tmpdirname,
                                os.path.join(self._kurucz_bin_path, "syntoascanga.exe")
                                )
        # self._call_external_code(tmpdirname,
        #                         os.path.join(self._kurucz_bin_path, "converfsynnmtoa.exe")
        #                         )                            
    def _call_external_code(self, cwd, program_path, input_data=None):
        process = subprocess.Popen([program_path],
                            stdin=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            encoding='utf8',
                            cwd=cwd)
        outs, errs = process.communicate(input_data)
        process.wait()

    def _ensure_atomic_symlinks(self, tmpdirname):
        self._ensure_symlink(
            os.path.join(self._atomic_data_path,"lines","he1tables.dat"),
            os.path.join(tmpdirname,"fort.18")
        )
        self._ensure_symlink(
            os.path.join(self._atomic_data_path,"lines","molecules.dat"),
            os.path.join(tmpdirname,"fort.2")
        )
        self._ensure_symlink(
            os.path.join(self._atomic_data_path,"lines","continua.dat"),
            os.path.join(tmpdirname,"fort.17")
        )

    def _ensure_symlink(self, source, destination):
        if os.path.islink(destination):
            if os.readlink(destination) == source:
                return
            os.remove(destination)
        elif os.path.exists(destination):
            os.remove(destination)
        os.symlink(source, destination)

    def _cleanup_workdir(self, tmpdirname):
        preserve = {
            os.path.join(tmpdirname, "fort.2"): os.path.join(self._atomic_data_path,"lines","molecules.dat"),
            os.path.join(tmpdirname, "fort.17"): os.path.join(self._atomic_data_path,"lines","continua.dat"),
            os.path.join(tmpdirname, "fort.18"): os.path.join(self._atomic_data_path,"lines","he1tables.dat"),
        }
        for f in glob.glob(os.path.join(tmpdirname,"fort.*")):
            if f in preserve and self._is_atomic_symlink(f, preserve[f]):
                continue
            try:
                os.remove(f)
            except FileNotFoundError:
                pass
        for name in ["spec.bin", "br_spec.bin", "br_spec.dat", "line_list.dat", "ROT1"]:
            path = os.path.join(tmpdirname, name)
            if os.path.exists(path):
                os.remove(path)

    def _is_atomic_symlink(self, path, target):
        return os.path.islink(path) and os.readlink(path) == target

    def _extend_model_for_SYNTHE(self, model):
        s = [
            "SURFACE INTENSI 17 1.,.9,.8,.7,.6,.5,.4,.3,.25,.2,.15,.125,.1,.075,.05,.025,.01\n",
            "ITERATIONS 1 PRINT 2 PUNCH 2\n",
            "CORRECTION OFF\n",
            "PRESSURE OFF\n",
            "MOLECULES ON\n",
            "READ MOLECULES\n",
            model.get_model_string()
        ]
        return "".join(s)

    def _get_synbeg_input(self, model):
        s = [
            f"AIR        {self.model.parameters.wave_min/10.:6.1f}    {self.model.parameters.wave_max/10.:6.1f}    600000.    {self.model.parameters.microturbulence:5.2f}    0     30    .0001     1    0\n",
            "AIRorVAC  WLBEG     WLEND     RESOLU    TURBV  IFNLTE LINOUT CUTOFF        NREAD\n"
        ]
        return "".join(s)

    def _save_to_file(self, path, string):
        with open(path,'w') as f:
            f.write(string)

    def _get_spectrv_string(self):
        s = [
        "0.0       0.        1.        0.        0.        0.        0.        0.\n",
        "0.\n",
        "RHOXJ     R1        R101      PH1       PC1       PSI1      PRDDOP    PRDPOW\n"
        ]
        return "".join(s)

    def _get_rotate_string(self):
        s = [
            "    1\n",
            f"{self.model.parameters.vsini:3.0f}."
        ]
        return "".join(s)

    def _get_broaden_string(self):
        return f"GAUSSIAN  {self.model.parameters.resolution:^8.1f}  RESOLUTION"

    def _get_pre_synthe_cache_key(self, model):
        hasher = hashlib.sha256()
        hasher.update(self._pre_synthe_cache_version.encode("utf-8"))
        hasher.update(self._extend_model_for_SYNTHE(model).encode("utf-8"))
        hasher.update(self._get_synbeg_input(model).encode("utf-8"))
        for path in [
            self._line_list_path,
            os.path.join(self._atomic_data_path,"lines","he1tables.dat"),
            os.path.join(self._atomic_data_path,"lines","molecules.dat"),
            os.path.join(self._atomic_data_path,"lines","continua.dat"),
        ]:
            hasher.update(self._get_file_cache_token(path).encode("utf-8"))
        return hasher.hexdigest()

    def _get_file_cache_token(self, path):
        try:
            stat = os.stat(path)
            return f"{path}:{stat.st_size}:{int(stat.st_mtime)}"
        except FileNotFoundError:
            return f"{path}:missing"

    def _restore_pre_synthe_cache(self, tmpdirname, cache_dir, cache_key):
        if cache_key is None:
            return False
        cache_path = os.path.join(cache_dir, cache_key)
        if not os.path.isdir(cache_path):
            return False
        for name in os.listdir(cache_path):
            src = os.path.join(cache_path, name)
            if not os.path.isfile(src):
                continue
            dest = os.path.join(tmpdirname, name)
            if os.path.exists(dest):
                os.remove(dest)
            shutil.copy2(src, dest)
        return True

    def _store_pre_synthe_cache(self, tmpdirname, cache_dir, cache_key):
        if cache_key is None:
            return
        os.makedirs(cache_dir, exist_ok=True)
        cache_path = os.path.join(cache_dir, cache_key)
        if os.path.isdir(cache_path):
            return
        os.makedirs(cache_path, exist_ok=True)
        for name in os.listdir(tmpdirname):
            src = os.path.join(tmpdirname, name)
            if os.path.isfile(src) and not os.path.islink(src):
                shutil.copy2(src, os.path.join(cache_path, name))

def _spectrum_worker(args):
    model, parameters, quiet, cache_dir = args
    worker = Synthe()
    return worker.get_spectrum(model, parameters=parameters, quiet=quiet, cache_dir=cache_dir)

def get_spectra_parallel(models, parameters_list=None, max_workers=None, quiet=True, cache_dir=None):
    """
    Compute multiple spectra in parallel using separate processes.
    """
    if parameters_list is None:
        parameters_list = [None] * len(models)
    if len(models) != len(parameters_list):
        raise ValueError("models and parameters_list must have the same length")
    if cache_dir is not None:
        os.makedirs(cache_dir, exist_ok=True)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(_spectrum_worker, (model, parameters, quiet, cache_dir))
            for model, parameters in zip(models, parameters_list)
        ]
        return [future.result() for future in futures]

def main():
    # http://wwwuser.oats.inaf.it/castelli/sources/linuxcodes.html
    # http://wwwuser.oats.inaf.it/castelli/sources/synthe/examples/synthenop.html
    # http://wwwuser.oats.inaf.it/castelli/sources/atlas9/ap00t10000g40k2odfnew.com
    # http://wwwuser.oats.inaf.it/castelli/sources/atlas9codes.html
    # http://wwwuser.oats.inaf.it/castelli/grids/gridp00k2odfnew/ap00k2tab.html
    pass

if __name__ == '__main__':
    main()