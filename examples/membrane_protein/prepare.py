import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
parser.add_argument('--replica',nargs='?',required=True,type=int)
#parser.add_argument('--ref_bead',nargs='?',required=True,type=int)
#parser.add_argument('--tmd_sel',required=True,type=str)
#parser.add_argument('--gpu_id',nargs='?',required=True,type=int)
args = parser.parse_args()

prot_data = {
    "AEROLYSIN":  {"ref_bead": 262, "tmd_sel": "resid 216 to 282 or resid 640 to 706 or resid 1064 to 1130 or resid 1488 to 1554 or resid 1912 to 1978 or resid 2336 to 2402 or resid 2760 to 2826"},
    "KALP25_IDR": {"ref_bead": 12,  "tmd_sel": "resid 3 to 22"},
    "RHD_IDR":    {"ref_bead": 125, "tmd_sel": "resid 15 to 29 or resid 35 to 49 or resid 118 to 131 or resid 139 to 150"},
    "RHO":        {"ref_bead": 52,  "tmd_sel": "resid 34 to 321"},
    "ADRA2A":     {"ref_bead": 134, "tmd_sel": "resid 46 to 185 or resid 208 to 245 or resid 379 to 458"},
    "hGHR":       {"ref_bead": 270, "tmd_sel": "resid 265 to 287"},
    "EGFR":       {"ref_bead": 633, "tmd_sel": "resid 622 to 643"},
    "NHE1_CHP1":  {"ref_bead": 390, "tmd_sel": "resid 99 to 591 or resid 914 to 1406"}
}

ref_bead = prot_data[args.name]['ref_bead']
tmd_sel = prot_data[args.name]['tmd_sel']

cwd = os.getcwd()
N_save = int(5e4)
N_frames = 1010
Lx = 30
Ly = Lx
area_per_lipid = .63
N_lipids = int(np.ceil(Lx*Ly/area_per_lipid)*2)

sysname = f'{args.name:s}_{args.replica:d}'
residues_file = f'{cwd}/input/residues.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [Lx, Ly, 100.], # nm
  temp = 293.15,
  ionic = 0.15, # molar
  pH = 7.4,
  topol = 'shift_ref_bead',
  bilayer_eq = True,
  friction = 0.01,
  pressure_coupling = True,
  pressure = [0,0,0],

  # RUNTIME SETTINGS
  #gpu_id = args.gpu_id,
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
  logfreq = 100*N_save,
  steps = N_frames*N_save, # number of simulation steps
  steps_eq = 20*N_save,
  runtime = 0, # overwrites 'steps' keyword if > 0
  platform = 'CUDA', # 'CUDA'
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
from calvados.analysis import calc_membrane_profiles, calc_domain_angles, calc_domain_rgs

calc_membrane_profiles("{path}","{sysname}","{output_path}","{residues_file}","{tmd_sel}",10)
angle_sels = dict(TMD="{tmd_sel}", D1="resid 52 to 145", D2="resid 146 to 252")
calc_domain_angles("{path}","{sysname}","{output_path}",angle_sels,0.5)
rg_sels = dict(ECD="resid 1 to 252", ICD="resid 313 to 639", ICD_GFP="resid 313 to 863", FL="resid 1 to 863")
calc_domain_rgs("{path}","{sysname}","{output_path}","{residues_file}",rg_sels,0.5)
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
  fdomains = f'{cwd}/input/domains.yaml', # domain definitions (harmonic restraints)
  pdb_folder = f'{cwd}/input', # directory for pdb and PAE files
)
components.add(name='POPC', molecule_type='lipid', nmol=N_lipids)
components.add(name=args.name, restraint=True, charge_termini='both', ref_bead=ref_bead)
components.write(path,name='components.yaml')

