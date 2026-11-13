import os
from calvados.cfg import Config, Job, Components
import subprocess
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from Bio import SeqIO

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',required=True,type=str)
parser.add_argument('--charge_termini', nargs='?',default='None',type=str, choices=['None','N','C','both'])
#parser.add_argument('--frac',nargs='?',required=True,type=int)
#parser.add_argument('--secondary',nargs='?',default='DPPC',type=str)
#parser.add_argument('--sigma',nargs='?',default=None,type=float)
#parser.add_argument('--omega',nargs='?',default=None,type=float)
#parser.add_argument('--gpu_id',nargs='?',required=True,type=int)

args = parser.parse_args()

# Dictionary to define starting APL (from lipid parametrization)
eq_apl_dict = {
    'DMPC': 0.57742,
    'DPPC': 0.519,
    'DOPC': 0.655,
    'POPC': 0.64,
}

cwd = os.getcwd()
N_save = int(1e5)
N_frames = 500
Lx = 25
Ly = Lx
area_per_lipid = eq_apl_dict['POPC']
N_lipids = int(np.ceil(Lx*Ly/area_per_lipid)*2)

if args.charge_termini == 'None':
    sysname = f'{args.name}'
else:
    sysname = f'{args.name}_charge_{args.charge_termini}'

#if args.frac < 100:
#    if args.sigma is not None:
#        if args.omega is not None:
#            sysname = f'{args.name}_{args.frac}_{args.secondary}_{100-args.frac}_sigma_{args.sigma}_omega_{args.omega}'
#        else:
#            sysname = f'{args.name}_{args.frac}_{args.secondary}_{100-args.frac}_sigma_{args.sigma}'
#    else:
#        if args.omega is not None:
#            sysname = f'{args.name}_{args.frac}_{args.secondary}_{100-args.frac}_omega_{args.omega}'
#        else:
#            sysname = f'{args.name}_{args.frac}_{args.secondary}_{100-args.frac}'
#else:
#    if args.sigma is not None:
#        if args.omega is not None:
#            sysname = f'{args.name}_{args.frac}_sigma_{args.sigma}_omega_{args.omega}'
#        else:
#            sysname = f'{args.name}_{args.frac}_sigma_{args.sigma}'
#    else:
#        if args.omega is not None:
#            sysname = f'{args.name}_{args.frac}_omega_{args.omega}'
#        else:
#            sysname = f'{args.name}_{args.frac}'

residues_file = f'{cwd}/input/residues.csv'

config = Config(
  # GENERAL
  sysname = sysname, # name of simulation system
  box = [Lx, Ly, 60.], # nm
  temp = 297.15,
  ionic = 0.15, # molar
  pH = 7,
  topol = 'slab',
  slab_outer = 15,
  bilayer_eq = True,
  friction = 0.01,
  pressure_coupling = True,
  pressure = [0,0,0],

  # RUNTIME SETTINGS
  #gpu_id = args.gpu_id,
  wfreq = N_save, # dcd writing frequency, 1 = 10fs
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
from calvados.analysis import SlabAnalysis, calc_bilayer_prop

slab = SlabAnalysis(name="{sysname:s}", input_path="{path:s}",
                    output_path="{output_path:s}", ref_name="bilayer", 
                    ref_chains = (100,{int(N_lipids+100-1):d}),
                    client_names = ["protein"],
                    client_chain_list = [(0,99)],
                    verbose=True)

slab.center(start=400, center_target='ref')
slab.calc_profiles()
slab.calc_concentrations()
calc_bilayer_prop(path="{path:s}",sysname="{sysname:s}",output_path="{output_path:s}")
"""

config.write(path,name='config.yaml',analyses=analyses)

components = Components(
  # Defaults
  molecule_type = 'protein',
  nmol = 1, # number of molecules
  restraint = False, # apply restraints
  charge_termini = args.charge_termini, # charge N or C or both

  # INPUT
  ffasta = f'{cwd}/input/fastalib.fasta', # input fasta file
  fresidues = residues_file, # residue definitions
)
components.add(name='POPC', molecule_type='lipid', nmol=int(N_lipids))
components.add(name=args.name, molecule_type='protein', nmol=100)
#if args.frac < 100:
#    components.add(name=args.secondary, molecule_type='lipid', nmol=int((100-args.frac)*N_lipids/100))
components.write(path,name='components.yaml')
