import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)

args = parser.parse_args()

# Dictionary to define starting APL (from lipid parametrization)
eq_apl_dict = {
    'DMPC': 0.57742,
    'DPPC': 0.519,
    'DOPC': 0.655,
    'POPC': 0.62,
}

cwd = os.getcwd()
N_save = int(5e4)
N_frames = 5200
Lx = 12
Ly = Lx
area_per_lipid = eq_apl_dict['POPC']
N_lipids = int(np.ceil(Lx*Ly/area_per_lipid)*2)

sysname = args.name
charge_termini = 'None' if sysname[2] in ['R','K'] else 'C'

residues_file = f'{cwd}/input/residues.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [Lx, Ly, 22.], # nm
  temp = 297.15,
  ionic = 0.15, # molar
  pH = 7,
  topol = 'shift_ref_bead',
  report_potential_energy = True,
  slab_outer = 15,
  bilayer_eq = True,
  friction = 0.01,
  pressure_coupling = True,
  pressure = [0,0,0],

  # RUNTIME SETTINGS
  #gpu_id = args.gpu_id,
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  logfreq = N_save,
  steps = N_frames*N_save, # number of simulation steps
  steps_eq = 20*N_save,
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # 'CUDA' or 'OpenCL'
  restart = 'checkpoint',
  frestart = 'restart.chk',
  verbose = True,
)

# PATH
path = f'{cwd}/{sysname}'
output_path = 'data'
subprocess.run(f'mkdir -p {path}',shell=True)
subprocess.run(f'mkdir -p {output_path}',shell=True)

analyses = f"""
from calvados.analysis import SlabAnalysis, calc_bilayer_prop, calc_com_profiles

slab = SlabAnalysis(name="{sysname:s}", input_path="{path:s}",
                    output_path="{output_path:s}", ref_name="bilayer",
                    ref_chains = (1,{int(N_lipids):d}),
                    client_names = ["peptides"],
                    client_chain_list = [(0,0)],
                    verbose=True)

slab.center(start=200, center_target='ref')
slab.calc_profiles()
slab.calc_concentrations()
calc_bilayer_prop(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}")
chainid_dict = dict(peptides=(0, 0), bilayer=(1, {int(N_lipids):d}))
calc_com_profiles(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}",residues_file="{residues_file:s}",
                  chainid_dict=chainid_dict)
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = False, # apply restraints
  charge_termini = 'None', # charge N or C or both

  # INPUT
  ffasta = f'{cwd}/input/fastalib.fasta', # input fasta file
  fresidues = residues_file, # residue definitions
)
components.add(name='POPC', molecule_type='lipid', nmol=int(N_lipids))
components.add(name=args.name, molecule_type='protein', nmol=1, ref_bead=2, pos_bead=[0,0,2.5], fixed_bead=True)
components.write(path,name='components.yaml')
