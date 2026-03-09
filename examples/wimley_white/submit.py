import os
import time
import shutil
import subprocess
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -p gpua100
#SBATCH --gres=gpu:1
#SBATCH -A lu2025-2-64
#SBATCH -t 24:0:0
#SBATCH -J {{name}}
#SBATCH -o {{name}}/out
#SBATCH -e {{name}}/err

source /home/gtesei00/.bashrc
conda activate calipids
module load CUDA/12.0.0

python prepare.py --name {{name}}
python {{name}}/run.py --path {{name}}""")

peptides = ['WLKLL']

for name in peptides:
    if not os.path.isdir(name):
        os.mkdir(name)
    with open(f'{name}.sh', 'w') as submit:
        submit.write(submission.render(name=name))
    subprocess.run(['sbatch',f'{name}.sh'])
    time.sleep(0.5)
