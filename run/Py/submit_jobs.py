#!/usr/bin/env python3
# load functions -------------------------------------------------------------------------
import sys
all_args = sys.argv[1:]
source_code_at = all_args[0]
sys.path.append(source_code_at)
from Py.libs import *
from Py.func import *

# Define paths ---------------------------------------------------------------------------
input_paths = read_json(all_args[1])
save_at = all_args[2]
if not os.path.exists(save_at):
    os.makedirs(save_at, exist_ok = True)
task_name = all_args[3]
out_paths = input_paths

# Define logs ----------------------------------------------------------------------------
logs_at = f'{save_at}/{task_name}.log'
logging.basicConfig(filename=logs_at, level=logging.DEBUG, filemode='w')
logger = logging.getLogger(__name__)
logging.info(f'Attempt at {datetime.now().strftime("%d/%m/%Y %H:%M:%S")}')

# get the paths for scripts --------------------------------------------------------------
master_script_paths = out_paths['master_script_path']
logging.info("master script paths read")

#import subprocess
# set up the SSH command to execute on the "qg" node
#ssh_command = f"ssh gogna@slurm-login 'screen -S run -p 0 -X stuff \"{echo "hi"}\\n\"'"
# execute the SSH command using subprocess
#subprocess.call(ssh_command, shell=True)

# Submit jobs ----------------------------------------------------------------------------
check = pd.read_csv(master_script_paths, sep = "\t", names = ['path', 'jobs'])
logging.info("Script data was read")

# Finish off the script ------------------------------------------------------------------
write_json(out_paths, f'{all_args[4]}')
print(f'{task_name} completed successfully. Paths written at {all_args[4]}')