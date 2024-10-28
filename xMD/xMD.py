import os
from copy import deepcopy
import glob
import shutil
import subprocess
from abc import ABC, abstractmethod
import pandas as pd
import argparse
from xMD.MD_Experiment import MD_Experiment
from xMD.MD_Settings import GROMACS_Settings
from xMD.AuxMD import run_MD, traj_to_pdb

class xMD(MD_Experiment):
    def __init__(self, settings: GROMACS_Settings, name=None, pdbcode: str = None, rep=None):
        print("Initializing xMD...")
        super().__init__(settings, name, pdbcode, rep)
        print(f"xMD initialized with name: {name}, pdbcode: {pdbcode}, rep: {rep}")

    def run_experiment(self, search=None, config_files=None, topology_files=None, restraints=None, rep=None, md_steps:int=None):
        print("Starting run_experiment...")
        self.set_replicate(rep)
        print(f"Replicate set to: {self.rep_no}")
        
        print("Preparing simulation...")
        self.prepare_simulation(search, config_files=config_files, topology_files=topology_files, restraints=restraints)
        
        if md_steps is None:
            md_steps = len(self.config_files)
        print(f"MD steps: {md_steps}")
        
        if len(self.config_files) == 1:
            self.config_files = self.config_files * md_steps
        print(f"Config files: {self.config_files}")
        
        assert len(self.config_files) == md_steps, "Number of config files must match number of steps"
        
        print("Running MD step...")
        tpr_path = self.run_MD_step()
        print(f"TPR path: {tpr_path}")
        
        print("Preparing analysis...")
        traj_file, pdb_top_file = self.prepare_analysis(tpr_path=tpr_path)
        print(f"Trajectory file: {traj_file}, PDB topology file: {pdb_top_file}")
        
        print("Running analysis...")
        self.run_analysis(traj_file=traj_file, tpr_path=tpr_path, pdb_top=pdb_top_file)
        print("Analysis complete.")

    def run_MD_step(self):
        print("Starting run_MD_step...")
        md_mdp, input_path, topo_path, tpr_path = super().run_MD_step()
        print(f"MD MDP: {md_mdp}, Input path: {input_path}, Topology path: {topo_path}, TPR path: {tpr_path}")
        
        assert isinstance(md_mdp, list), "md_mdp must be a list of mdp files"
        
        self.set_trajectory_number()
        print(f"Trajectory number set to: {self.traj_no}")
        
        for idx, mdp in enumerate(md_mdp):
            print(f"Running MD for MDP {idx + 1}/{len(md_mdp)}: {mdp}")
            input_path = run_MD(mdp,
                                input_path,
                                topo_path,
                                tpr_path,
                                self.gmx[0],
                                self.restraints[idx],
                                self.settings.gpu)
            print(f"New input path: {input_path}")
            
            if idx < len(md_mdp) - 1:
                self.traj_no += 1
                print(f"Incrementing trajectory number to: {self.traj_no}")
                _, _, _, tpr_path = super().run_MD_step()
                print(f"Updated TPR path: {tpr_path}")
        
        print(f"Returning final TPR path: {tpr_path}")
        return tpr_path

    def run_analysis(self, traj_file=None, tpr_path=None, pdb_top=None):
        print("Starting run_analysis...")
        traj_file, tpr_path, pdb_path = super().run_analysis(traj_file, tpr_path, pdb_top)
        print(f"Trajectory file: {traj_file}, TPR path: {tpr_path}, PDB path: {pdb_path}")
        
        print("Converting trajectory to PDB...")
        traj_to_pdb(traj_file,
                    pdb_top,
                    pdb_path,
                    self.gmx[0])
        print("Trajectory conversion complete.")