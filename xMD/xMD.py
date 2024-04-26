# xMD class for managing MD experiments
# MD_Experiment class handles all the functions for running MD experiments
# This class desribes the routines for running MD experiments

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
        super().__init__(settings, name, pdbcode, rep)


    def run_experiment(self, 
                       search=None, 
                       config_files=None, 
                       topology_files=None, 
                       restraints=None,
                       rep=None, 
                       md_steps:int=None):
        """
        This will run the experiment for the trial.
        suffix is the suffix for the initial topology files. 
        Should revert back to the suffix set in the settings.
        """
        ### TODO more flexibile setup of experiment
        # how do we make sure settings are not overwritten by this method?
        # TODO setup tensorboard logging
        # write parser for gromacs outputs
        # write tensorboard logger
        # use threading to run logger concurrently
        self.set_replicate(rep)
        self.prepare_simulation(search,
                                config_files=config_files,
                                topology_files=topology_files,
                                restraints=restraints)
        if md_steps is None:
            md_steps = len(self.config_files)
        if len(self.config_files) == 1:
            self.config_files = self.config_files * md_steps
        assert len(self.config_files) == md_steps, "Number of config files must match number of steps"

        tpr_path = self.run_MD_step()
        traj_file, pdb_top_file = self.prepare_analysis(tpr_path=tpr_path)
        self.run_analysis(traj_file=traj_file, tpr_path=tpr_path, pdb_top=pdb_top_file)

    ## TODO add repeat steps - run for as many mdp files are provided.
    def run_MD_step(self): # add MIG compatibility, MIG can only be used on the final config file provided
        """
        This will run the steps of MD for the trial.
        Retruns the tpr file name.
        """

        md_mdp, input_path, topo_path, tpr_path = super().run_MD_step()

        assert isinstance(md_mdp, list), "md_mdp must be a list of mdp files"

        self.set_trajectory_number() 

        for idx, mdp in enumerate(md_mdp):
            input_path = run_MD(mdp, 
                                input_path, 
                                topo_path, 
                                tpr_path, 
                                self.gmx[0],
                                self.restraints[idx],
                                self.settings.gpu)
            # the last step should not be repeated
            if idx < len(md_mdp)-1:
                self.traj_no += 1
                _,_,_, tpr_path = super().run_MD_step() 
    
            # add log file to tensorboard as text
        return tpr_path

### TODO create analysis commands - add data to the dataframe
    def run_analysis(self, traj_file=None, tpr_path=None, pdb_top=None):

        
        traj_file, tpr_path, pdb_path = super().run_analysis(traj_file, tpr_path, pdb_top)

        # change this to - concatenate and align the xtc files to the first frame of the trajectory
        traj_to_pdb(traj_file,
                    pdb_top,
                    pdb_path,
                    self.gmx[0])
