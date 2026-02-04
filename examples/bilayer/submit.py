import subprocess
import os
import pandas as pd
import numpy as np
import mdtraj as md
import time
import shutil
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH --job-name={{path}}
#SBATCH --nodes=1
#SBATCH --partition=qgpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=18
#SBATCH -t 4:00:00
#SBATCH -o {{path}}/out
#SBATCH -e {{path}}/err

source /home/gitesei/.bashrc
module load gcc/11.2.0 openmpi/4.0.3 cuda/11.2.0
conda activate calvados

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python prepare.py --name {{name}} --temp {{temp}} --alpha {{alpha}}
python {{path}}/run.py --path {{path}}""")

for name in ['POPC']:
    for alpha in [0.4]:
        for temp in [293,303,323]:
            path = f'{name:s}_{temp:d}_{alpha:.2f}'
            if not os.path.isdir(f'{path:s}'):
                os.mkdir(f'{path:s}')
            with open(f'{path:s}/submit.sh', 'w') as submit:
                submit.write(submission.render(name=name,temp=temp,alpha=alpha,path=path))
            subprocess.run(['sbatch',f'{path:s}/submit.sh'])
            print(f'{path:s}')
            time.sleep(.6)
