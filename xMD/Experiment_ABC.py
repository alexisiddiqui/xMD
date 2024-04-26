
import os
from copy import deepcopy
import glob
import pandas as pd
import time
import shutil
import pickle
from abc import ABC, abstractmethod
from .MD_Settings import Settings
### Abstract method for the MD and Docking experiment classes

class Experiment(ABC):
    """
    This class represents an experimental trial.
    WHat do we need in this class?
    This class will handle some of the folder checking as well as the settings used for the run.
    """    
    def __init__(self,settings: Settings, name=None, pdbcode: str = None, rep=None):
        super().__init__()
        self.settings = settings
        if name is not None:
            self.name = name
        else:
            self.name = self.settings.trial_name
        self.dataframe  = pd.DataFrame()
        self.dirs = {dir : None for dir in self.settings.dirs_to_create}
        self.trajectories = {}
        self.rep_no : int = None
        self.traj_no : int = 0
        self.set_replicate(rep)

        self.config_files = []
        self.topology_files = []
        self.restraints = []
        if pdbcode is not None:
            self.settings.pdbcode = pdbcode
        self.generate_path_structure(self.name)


    def generate_path_structure(self, trial_name=None):
        """
        Generates the path structure for the trial.
        Sets self.dirs.
        Returns the data directory.
        """
        if trial_name is None:
            trial_name = self.name

        for dir in self.dirs.keys():
            dirs_to_join = [dir, 
                            self.settings.parent,
                            self.settings.pdbcode, 
                            trial_name]
            trial_dir = os.path.join(*dirs_to_join)
            print(f"Trial directory {dir}: ", trial_dir)
            self.dirs[dir] = trial_dir

        return self.dirs[self.settings.data_directory] 

    def create_directory_structure(self, overwrite=False):
        """
        Creates the directory structure for the trial.
        Checks if the trial directory exists. If not, creates it.
        Otherwise, adds a number to the end of the trial name and tries again.
        """
        if not overwrite:
            trial_name = self.name
            trial_dir = self.generate_path_structure(trial_name)
            trial_exists = os.path.exists(trial_dir)
            trial_number = 1

            while trial_exists:
                trial_name = self.name + str(trial_number)
                trial_dir = self.generate_path_structure(trial_name)
                trial_exists = os.path.exists(trial_dir)
                trial_number += 1

            self.name = trial_name

        self.create_directories()
        print("Created directories for trial: ", self.name)

    @abstractmethod
    def save_experiment(self, save_name=None):
        """
        Writes the settings of the experiment (contents of class) to a pickle file. Logs???
        """
        unix_time = int(time.time())
        if save_name is not None:
            save_name = save_name+"_"+str(unix_time)+".pkl"
            save_path = os.path.join(self.settings.logs_directory, save_name)

            with open(save_path, 'wb') as f:
                pickle.dump(self, f)
                print("Saving experiment to: ", save_path)
                return save_path
            

    ### TODO changed to abstract method - will this save settings? need to also be able to save and load a human readable file
    # @abstractmethod
    def load_experiment(self, latest=False, idx=None, load_path=None):
        """
        Loads the settings of the experiment from a pickle file.
        """
        if load_path is not None:
            print("Attempting to load experiment from: ", load_path)
            with open(load_path, 'rb') as f:
                print("Loading experiment from: ", load_path)
                return pickle.load(f)
            # print("Loaded object type:", type(loaded_obj))
            # return deepcopy(loaded_obj)
        # If no explicit path is provided
        search_dir = self.dirs[self.settings.logs_directory]
        print("Searching for experiment files in: ", search_dir)

        pkl_files = glob.glob(os.path.join(search_dir, "*.pkl"))
        print("Found files: ", pkl_files)

        if not pkl_files:
            print("No experiment files found.")
            raise FileNotFoundError
        pkl_files = sorted(pkl_files, key=os.path.getctime)
        if latest is True:
            print("Loading latest experiment.")
            file = pkl_files[-1]
        elif idx is not None:
            print(f"Loading {idx} experiment.")
            file = pkl_files[idx]
        else:
            print("Loading first experiment.")
            file = pkl_files[0]

        load_path = file    
        print("Loading experiment from: ", load_path)
        with open(load_path, 'rb') as f:
            return pickle.load(f)
        # return deepcopy(loaded_obj)


    def create_directories(self):
        # we are using the ABC method here
        data_dir = self.generate_path_structure()
        for dir in self.settings.dirs_to_create:
            #split path back into list
            dirs_to_join = data_dir.split(os.sep)
            dirs_to_join[0] = dir
            trial_dir = os.path.join(*dirs_to_join)
            
            os.makedirs(trial_dir, exist_ok=True)
        # then we can create the MD specific directories
        rep_dirs = [self.settings.rep_directory + 
                          str(i) for i in range(1,self.settings.replicates+1)]

        data_dir = self.generate_path_structure(self.name)

        for directory in rep_dirs:
            print("Creating directory: ", directory, end=" ")
            path = os.path.join(data_dir, directory)
            os.makedirs(path, exist_ok=True)

    def prepare_config(self, file_names=None):
        """
        Prepares the configuration files for the trial based on file names. 
        Otherwise grabs all files in the config folder.
        """
        config_files = os.listdir(os.path.join(self.settings.config))
        print("Config files: :", config_files)
        if file_names is not None:
            config_files = file_names

        self.config_files = config_files
        print("Loading config files: ", config_files)

    ### TODO add index files
    def prepare_input_files(self, search=None, file_names=None):
        """
        Prepares the topology files for the trial.
        """
        ### TODO refactor this to use self.dirs
        topology_files = os.listdir(os.path.join(self.settings.topology))
        topology_files = [file 
                        for file in topology_files 
                        if self.settings.pdbcode in file.split(".")[-2]]  
     
        print(f"Topology files for {self.settings.pdbcode}: ", topology_files)
        if search is None:
            search = self.settings.search

        if search is not None:
            topology_files = [file 
                              for file in topology_files 
                              if search in file.split(".")[-2]]
            
        if file_names is not None:
            topology_files = [file for file in topology_files if file in file_names]
        
        self.topology_files = topology_files
        print("Loading topology files: ", self.topology_files)

### This is mostly to check remote directories: return to this later.
    def check_all_trajectory_files(self, data_dir=None, traj_extension=None):
        """
        Checks if the trajectory files exist. 
        For a given trial goes into each replicate and checks for file with the correct extension.
        Adds the name to the trajectories dictionary. If a data dir is given, it will check it.
        """
        if data_dir is None:
            data_dir = self.dirs[self.settings.data_directory]

        if traj_extension is None:
            traj_extension = self.settings.search_traj
        print("Checking trajectory files in: ", data_dir)

        for _dir in os.listdir(data_dir):
            if os.path.isdir(os.path.join(data_dir, _dir)):
                rep = _dir
                self.trajectories[rep] = []

                for file in os.listdir(os.path.join(data_dir, rep)):
                    if traj_extension in file and "#" not in file:
                        self.trajectories[rep].append(file)
        print("Trajectory files: ", self.trajectories)
        return self.trajectories
    

### TODO later: This is for loading remote trajectories into local folder
    def load_trajectory_files(self, data_dir=None):
        """
        This is for syncing the trajectory files to the current experiment.
        This is meant to be for when files need to be transferred from the scratch disk to the home disk.
        This will copy all files from the data_dir to the data dict.
        """


        ### One function for local copy

        ### One function for remote copy


    def load_input_files(self, rep=None):
        """
        Copies the topology files for the trial.
        Files must be local.
        """
        if rep is None:
            rep = self.rep_no

        if rep is None:
            raise ValueError("Replicate number not set.")

        for file in self.topology_files:
            file_path = os.path.join(self.settings.topology, file)
            ### TODO refactor this to use self.dirs
            destination = os.path.join(self.dirs[self.settings.data_directory],
                                       self.settings.rep_directory + str(rep),
                                       file)
            shutil.copyfile(file_path, destination)
  
    def set_replicate(self, rep=None):
        """
        Sets the replicate for the trial.
        """
        if rep is not None:
            self.rep_no = int(rep)
        
        # check self.trajectories[rep] exists
        if self.trajectories.get(self.settings.rep_directory + str(self.rep_no)) is None:
            self.trajectories[self.settings.rep_directory + str(self.rep_no)] = []

        print("Replicate number: ", self.rep_no)

    def set_trajectory_number(self, trajectory_number=None, suffix=None):
        """
        Sets the trajectory number for the trial.
        """
        if trajectory_number is None:
            self.traj_no = self.find_latest_trajectory(suffix=suffix)
        else:
            self.traj_no = int(trajectory_number)
        print("Trajectory number: ", self.traj_no)

    def find_latest_trajectory(self, suffix=None, rep=None):
        """
        Finds the latest trajectory for the current trial. Uses the trajectory settings.
        Does not check converted trajectories.
        """
        if suffix is None:
            suffix = self.settings.suffix

        if rep is None:
            rep = self.rep_no

        if rep is None:
            raise ValueError("Replicate number not set.")

        self.check_all_trajectory_files()

        trajectory_files = self.trajectories[self.settings.rep_directory + str(rep)]
        #remove entries which dont contain the suffix
        trajectory_files = [file for file in trajectory_files if suffix in file]

        #sort the files by the number at the end
        trajectory_files = sorted(trajectory_files, 
                                  key=lambda x: int(x.split("_")[-1].split(".")[0]))

        # return traj_no
        try:
            return int(trajectory_files[-1].split("_")[-1].split(".")[0])
        except:
            return 0
