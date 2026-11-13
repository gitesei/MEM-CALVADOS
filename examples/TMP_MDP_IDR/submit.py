import os
import time
import shutil
import subprocess
import mdtraj as md
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH --gres=gpu:1
#SBATCH --partition=gpu
#SBATCH -t 7:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -A naiss2026-3-355-gpu
#SBATCH -J {{name_1}}_{{name_2}}_{{replica}}
#SBATCH -o {{name_1}}_{{name_2}}_Sos1_{{replica}}/out
#SBATCH -e {{name_1}}_{{name_2}}_Sos1_{{replica}}/err

module load GPU/Miniforge/26.3.2-2-eb

mamba activate memcal

python prepare.py --name_1 {{name_1}} --name_2 {{name_2}} --replica {{replica}}

python {{name_1}}_{{name_2}}_Sos1_{{replica}}/run.py --path {{name_1}}_{{name_2}}_Sos1_{{replica}} &

wait
""")

frames_per_batch = 99
for name_1, name_2 in [['LAT', 'Grb2'],['pLAT', 'Grb2']]:
    for replica in [0,1,2]:
        if not os.path.isdir(f'{name_1}_{name_2}_Sos1_{replica}'):
            os.mkdir(f'{name_1}_{name_2}_Sos1_{replica}')
        t = md.load(f'{name_1}_{name_2}_Sos1_{replica}/{name_1}_{name_2}_Sos1_{replica}.dcd', top=f'{name_1}_{name_2}_Sos1_{replica}/top.pdb')
        n_keep = (t.n_frames // frames_per_batch) * frames_per_batch
        print(replica,t.n_frames,t.n_frames // frames_per_batch)
        t[:n_keep].save_dcd(f'{name_1}_{name_2}_Sos1_{replica}/{name_1}_{name_2}_Sos1_{replica}.dcd')
        print(replica,md.load(f'{name_1}_{name_2}_Sos1_{replica}/{name_1}_{name_2}_Sos1_{replica}.dcd', top=f'{name_1}_{name_2}_Sos1_{replica}/top.pdb').n_frames)
        with open(f'{name_1}_{name_2}_{replica}.sh', 'w') as submit:
            submit.write(submission.render(name_1=name_1,name_2=name_2,replica=replica))
        subprocess.run(['sbatch',f'{name_1}_{name_2}_{replica}.sh'])
        time.sleep(0.5)

"""
for name_1, name_2 in [['LAT', 'Grb2'],['pLAT', 'Grb2']][:]:
    for replica in [0,1,2]:
        if not os.path.isdir(f'{name_1}_{name_2}_Sos1_{replica}'):
            os.mkdir(f'{name_1}_{name_2}_Sos1_{replica}')
        with open(f'{name_1}_{name_2}_{replica}.sh', 'w') as submit:
            submit.write(submission.render(name_1=name_1,name_2=name_2,replica=replica))
        subprocess.run(['sbatch',f'{name_1}_{name_2}_{replica}.sh'])
        time.sleep(0.5)
"""
