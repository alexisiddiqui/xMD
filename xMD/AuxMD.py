# auxiliary functions for xMD
import os
from copy import deepcopy
import glob
import shutil
import subprocess
import pandas as pd
import argparse

import os
import subprocess

import os
import subprocess
import os
import subprocess

def run_MD(md_mdp: str,
           input_path: str,
           topo_path: str,
           tpr_path: str,
           gmx: str,
           restraints: str = None,
           gpu: bool = False):
    if restraints is None:
        restraints = input_path
    
    grompp_command = [gmx, "grompp",
                      "-f", md_mdp,
                      "-c", input_path,
                      "-p", topo_path,
                      "-o", tpr_path,
                      "-r", restraints,
                      "-maxwarn", "1",
                      "-v"]
    subprocess.run(grompp_command, check=True)

    # Construct the mdrun command
    mdrun_command = [gmx, "mdrun",
                     "-v",
                     "-deffnm", tpr_path.replace(".tpr", ""),
                     "-pin", "on"]

    if gpu:
        # Get the number of GPUs
        ngpus = int(os.environ.get('SLURM_GPUS_ON_NODE', '0'))
        gpu_ids = "".join(str(i) for i in range(ngpus))
        if ngpus > 1:
            mdrun_command.extend(["-gpu_id", gpu_ids])
            print(f"Using {ngpus} GPUs with IDs: {gpu_ids}")
        elif ngpus == 1:
            mdrun_command.extend(["-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu"])
            print(f"Using {ngpus} GPUs with IDs: {gpu_ids}")
        else:
            mdrun_command.extend(["-nb", "gpu", "-pme", "gpu", "-bonded", "gpu", "-update", "gpu"])
            print("GPU flag set, but no GPUs detected in the SLURM environment.")

    print(f"Running command: {' '.join(mdrun_command)}")
    subprocess.run(mdrun_command, check=True)
    
    input_path = tpr_path.replace(".tpr", ".gro")
    return input_path

def traj_to_pdb(traj_file: str,
                top_path: str,
                pdb_path: str,
                gmx: str="gmx"):
    pdbout_command = [gmx, "trjconv", 
                        "-f", traj_file,
                        "-s", top_path,
                        "-o", pdb_path]

    subprocess.run(pdbout_command, input=b"1\n", check=True)
    print("PDB file written to: ", pdb_path)  
    
      

def PBC_conversion():

    pass
