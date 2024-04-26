import os
from copy import deepcopy
import glob
import shutil
import subprocess
from abc import ABC, abstractmethod
import pandas as pd
import argparse
import time
import pickle
from .Experiment_ABC import Experiment
from .MD_Settings import GROMACS_Settings

class MD_Experiment(Experiment):
    def __init__(self,settings: GROMACS_Settings, name=None, pdbcode=None, rep=None):
        super().__init__(settings, name, pdbcode, rep)
        self.set_mdrun_gmx()
        if __name__ == "__main__":
            self.args = self.check_args()
        self.set_environs()

    def set_mdrun_gmx(self, mpi_on: bool = None):
        if mpi_on is None:
            mpi_on = self.settings.gmx_mpi_on

        if mpi_on:
            self.gmx = ("gmx_mpi","gmx")
        else:
            self.gmx = ("gmx","gmx_mpi")


    def check_args(self):
        """
        Checks the arguments for the trial.
        """
        parser = argparse.ArgumentParser()

        parser.add_argument("-R", "--replicate", 
                            dest="replicate",
                            help="Replicate number", type=int,
                            default=0)
        
        ### TODO setup traj continuation
        # parser.add_argument("-T", "--trajectory", 
        #                     dest="traj_no", 
        #                     help="Trajectory number", type=int)
                            
        
        parser.add_argument("-P", "--pdbcode", 
                            dest="pdbcode", 
                            help="PDB code", type=str)
        
        parser.add_argument("-N", "--name",
                            dest="name",
                            help="Name of the trial", type=str)
        
        parser.add_argument("-S", "--suffix",
                            dest="suffix",
                            help="string for the trajectory files", type=str)
        
        parser.add_argument("-s", "--search",
                            dest="search",
                            help="search string for the top files", type=str)


        args = parser.parse_args()

        print("Arguments: ", args)

        if args.replicate is not None:
            self.set_replicate(args.replicate)

        # if args.traj_no is not None:
        #     self.traj_no = args.traj_no

        if args.pdbcode is not None:
            self.settings.pdbcode = args.pdbcode

        if args.name is not None:
            self.name = args.name
        
        if args.suffix is not None:
            self.settings.suffix = args.suffix

        if args.search is not None:
            self.settings.search = args.search

        return args
    
    def pbc_conversion(self, tpr_path):
        """
        Converts the trajectory file to correct for pbc.
        Returns the corrected trajectory file name.
        """
        gmx = self.gmx[0]
        # change this to find all xtc files in the directory - concatenate them and then correct PBC

        traj_file = tpr_path.replace(".tpr", ".xtc")
        traj_file1 = traj_file.split(".")[-2] + self.settings.pbc_extensions[0] + ".xtc"
        traj_file2 = traj_file.split(".")[-2] + self.settings.pbc_extensions[1] + ".xtc"
        
        trjconv_command1 = [gmx, "trjconv",
                             "-f", traj_file, 
                             "-s", tpr_path, 
                             *self.settings.pbc_commands[0], 
                             "-o", traj_file1]
        
        trjconv_command2 = [gmx, "trjconv", 
                             "-f", traj_file1, 
                             "-s", tpr_path, 
                             *self.settings.pbc_commands[1], 
                             "-o", traj_file2]
        
        print("Running trjconv command 1: ", trjconv_command1)
        subprocess.run(trjconv_command1, input=b"1\n0\n", check=True)

        print("Running trjconv command 2: ", trjconv_command2)
        subprocess.run(trjconv_command2, input=b"1\n0\n", check=True)

        pdb_file = traj_file2.replace(".xtc", ".pdb")
        
        trjconv_command3 = [gmx, "trjconv", 
                             "-f", traj_file1, 
                             "-s", tpr_path, 
                             *self.settings.pbc_commands[1], 
                             "-o", pdb_file,
                                "-dump", "0"]

        print("Running trjconv command 3: ", trjconv_command3)
        subprocess.run(trjconv_command3, input=b"1\n0\n", check=True)  

        return traj_file2, pdb_file

### TODO sort out the trajfile naming
    def prepare_analysis(self, tpr_path):
        """
        This will prepare the analysis for the trial.
        """
        traj_file2, pdb_file = self.pbc_conversion(tpr_path)
        # TODO: concatenate trajecotry files 

        return traj_file2, pdb_file

    def prepare_simulation(self, search=None, 
                           config_files: list = None, 
                           topology_files: list = None, 
                           restraints: list = None):
        """
        This will prepare the simulation for the trial.
        """

        if restraints is not None:

            if isinstance(restraints, str):
                restraints = [restraints]*len(config_files)
        else:
            restraints = [None]*len(config_files)

        assert len(restraints) == len(config_files), "Number of restraints must match number of config files"
        self.restraints = restraints

        self.prepare_config(config_files)
        self.prepare_input_files(search, topology_files)
        self.load_input_files()

    def save_experiment(self, save_name=None):
        super().save_experiment(save_name=save_name)

        unix_time = int(time.time())
        if save_name is None:
            save_name = [self.settings.parent,
                        self.settings.pdbcode, 
                        self.name,
                        str(self.rep_no),
                        str(unix_time)]
            save_name = "_".join(save_name)+".pkl"
            
            log_dir = self.dirs[self.settings.logs_directory]
            
            save_path = os.path.join(log_dir, save_name)

            with open(save_path, 'wb') as f:
                pickle.dump(self, f)
        print("Saved experiment to: ", save_path)
        return save_path


    def set_environs(self):
        """
        Sets the environment variables for the trial.
        """
        os.environ[self.settings.environ] = self.settings.environ_path
        print("Environment variables set: ", self.settings.environ, self.settings.environ_path)

        
    def prepare_TB_writer(self):
        """
        This will prepare the TB writer for the trial.
        """
        pass

    @abstractmethod
    def run_MD_step(self):
        """
        Abstract method for running a step of MD for the trial.
        Handles the preparation of input files.
        Returns:
            None
        """
        
        topo_files = [f for f in self.topology_files if f.endswith(".top")]
        gro_files = [f for f in self.topology_files if f.endswith(".gro")]
        tpr_name = "_".join([self.settings.suffix, 
                            self.settings.pdbcode, 
                            str(self.traj_no)]) + ".tpr"
        
        tpr_path = os.path.join(self.dirs[self.settings.data_directory], 
                                self.settings.rep_directory + str(self.rep_no), 
                                tpr_name)
        
        md_mdps = [os.path.join(self.settings.config, config) for config in self.config_files]
        topo_name = topo_files[0]
        topo_path = os.path.join(self.dirs[self.settings.data_directory], 
                                self.settings.rep_directory + str(self.rep_no), 
                                topo_name)
        
        gro_name = gro_files[0]
        input_path = os.path.join(self.dirs[self.settings.data_directory], 
                                self.settings.rep_directory + str(self.rep_no), 
                                gro_name)

        return md_mdps, input_path, topo_path, tpr_path
    
    @abstractmethod
    def run_analysis(self, traj_file=None, tpr_path=None, pdb_top=None):
        ### This will take args from settings for what analyses to run
        # for now we just convert to pdb
        #if traj file is not given try to rebuild name
        if tpr_path is None:
            tpr_name = "_".join([self.settings.suffix, 
                                            self.settings.pdbcode, 
                                            str(self.traj_no)]) + ".tpr"
            tpr_path = os.path.join(self.dirs[self.settings.data_directory], 
                                    self.settings.rep_directory + str(self.rep_no), 
                                    tpr_name)

        if traj_file is None:
            traj_file = tpr_path.replace(".tpr", ".xtc")
            traj_file = traj_file.split(".")[-2] + self.settings.pbc_extensions[1] + ".xtc"
            pdb_name = str(self.rep_no) + "_" + tpr_name.replace(".tpr", ".pdb")
            pdb_name = pdb_name.split(os.sep)[-1]
        else:
            pdb_name = str(self.rep_no) + "_" + traj_file.replace(".xtc", ".pdb")
            pdb_name = pdb_name.split(os.sep)[-1]


        pdb_path = os.path.join(self.dirs[self.settings.viz], pdb_name)

        return traj_file, tpr_path, pdb_path
        

