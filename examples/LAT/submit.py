import os
import time
import shutil
import subprocess
from jinja2 import Template

submission = Template("""#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -p gpua100i
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH -A lu2025-2-64
#SBATCH -t 30:0:0
#SBATCH -J {{name}}_{{replica}}
#SBATCH -o {{name}}_{{replica}}.out
#SBATCH -e {{name}}_{{replica}}.err

source /home/gtesei00/.bashrc
conda activate calipids
module load CUDA/12.0.0

python prepare.py --name {{name}} --replica {{replica}}

python {{name}}_{{replica}}/run.py --path {{name}}_{{replica}}
""")

for name in ['LAT']:
    for replica in [0]:
        if not os.path.isdir(f'{name:s}_{replica:d}'):
            os.mkdir(f'{name:s}_{replica:d}')
        with open(f'{name:s}_{replica:d}.sh', 'w') as submit:
            submit.write(submission.render(name=name,replica=replica))
        subprocess.run(['sbatch',f'{name:s}_{replica:d}.sh'])
        time.sleep(0.5)
