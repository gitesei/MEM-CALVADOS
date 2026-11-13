import os
import time
import shutil
import subprocess
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p gpua100i
#SBATCH --gres=gpu:1
#SBATCH -A lu2025-2-64
#SBATCH -t 24:0:0
#SBATCH -J {{sysname}}
#SBATCH -o {{sysname}}/out
#SBATCH -e {{sysname}}/err

source /home/riccardosaltutti/.bashrc
conda activate calvados
module load CUDA/12.0.0

python prepare.py --name {{name}} --charge_termini {{charge}} 
python {{sysname}}/run.py --path {{sysname}}""")

peptides = ['WLALL', 'WLFLL', 'WLRLL']

#--charge_termini {{charge}}
charges = ['None']

for name in peptides:
   for charge in charges:
        if charge == 'None':
            sysname = f'{name}'
        else:
            sysname = f'{name}_charge_{charge}'
        if not os.path.isdir(sysname):
            os.mkdir(sysname)
        with open(f'{sysname}.sh', 'w') as submit:
            submit.write(submission.render(name=name, charge=charge, sysname=sysname))
        subprocess.run(['sbatch',f'{sysname}.sh'])
        time.sleep(0.5)
