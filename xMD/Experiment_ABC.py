import os
from copy import deepcopy
import glob
import pandas as pd
import time
import shutil
import pickle
from abc import ABC, abstractmethod
from .MD_Settings import Settings

class Experiment(ABC):
    def __init__(self, settings: Settings, name=None, pdbcode: str = None, rep=None):
        super().__init__()
        print("Initializing Experiment...")
        self.settings = settings
        if name is not None:
            self.name = name
        else:
            self.name = self.settings.trial_name
        print(f"Experiment name: {self.name}")
        self.dataframe = pd.DataFrame()
        self.dirs = {dir: None for dir in self.settings.dirs_to_create}
        self.trajectories = {}
        self.rep_no: int = None
        self.traj_no: int = 0
        self.set_replicate(rep)
        print(f"Replicate set to: {self.rep_no}")

        self.config_files = []
        self.topology_files = []
        self.restraints = []
        if pdbcode is not None:
            self.settings.pdbcode = pdbcode
        print(f"PDB code: {self.settings.pdbcode}")
        self.generate_path_structure(self.name)

    def generate_path_structure(self, trial_name=None):
        print("Generating path structure...")
        if trial_name is None:
            trial_name = self.name

        for dir in self.dirs.keys():
            dirs_to_join = [dir, 
                            self.settings.parent,
                            self.settings.pdbcode, 
                            trial_name]
            trial_dir = os.path.join(*dirs_to_join)
            print(f"Trial directory {dir}: {trial_dir}")
            self.dirs[dir] = trial_dir

        return self.dirs[self.settings.data_directory]

    def create_directory_structure(self, overwrite=False):
        print("Creating directory structure...")
        if not overwrite:
            trial_name = self.name
            trial_dir = self.generate_path_structure(trial_name)
            trial_exists = os.path.exists(trial_dir)
            trial_number = 1

            while trial_exists:
                print(f"Trial {trial_name} already exists. Trying new name...")
                trial_name = self.name + str(trial_number)
                trial_dir = self.generate_path_structure(trial_name)
                trial_exists = os.path.exists(trial_dir)
                trial_number += 1

            self.name = trial_name

        self.create_directories()
        print(f"Created directories for trial: {self.name}")

    @abstractmethod
    def save_experiment(self, save_name=None):
        print("Saving experiment...")
        unix_time = int(time.time())
        if save_name is not None:
            save_name = save_name+"_"+str(unix_time)+".pkl"
            save_path = os.path.join(self.settings.logs_directory, save_name)

            with open(save_path, 'wb') as f:
                pickle.dump(self, f)
                print(f"Saved experiment to: {save_path}")
                return save_path

    def load_experiment(self, latest=False, idx=None, load_path=None):
        print("Loading experiment...")
        if load_path is not None:
            print(f"Attempting to load experiment from: {load_path}")
            with open(load_path, 'rb') as f:
                print(f"Loading experiment from: {load_path}")
                return pickle.load(f)

        search_dir = self.dirs[self.settings.logs_directory]
        print(f"Searching for experiment files in: {search_dir}")

        pkl_files = glob.glob(os.path.join(search_dir, "*.pkl"))
        print(f"Found files: {pkl_files}")

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
        print(f"Loading experiment from: {load_path}")
        with open(load_path, 'rb') as f:
            return pickle.load(f)

    def create_directories(self):
        print("Creating directories...")
        data_dir = self.generate_path_structure()
        for dir in self.settings.dirs_to_create:
            dirs_to_join = data_dir.split(os.sep)
            dirs_to_join[0] = dir
            trial_dir = os.path.join(*dirs_to_join)
            
            os.makedirs(trial_dir, exist_ok=True)
            print(f"Created directory: {trial_dir}")

        rep_dirs = [self.settings.rep_directory + 
                    str(i) for i in range(1, self.settings.replicates+1)]

        data_dir = self.generate_path_structure(self.name)

        for directory in rep_dirs:
            path = os.path.join(data_dir, directory)
            os.makedirs(path, exist_ok=True)
            print(f"Created directory: {path}")

    def prepare_config(self, file_names=None):
        print("Preparing config files...")
        config_files = os.listdir(os.path.join(self.settings.config))
        print(f"Config files: {config_files}")
        if file_names is not None:
            config_files = file_names

        self.config_files = config_files
        print(f"Loading config files: {config_files}")

    def prepare_input_files(self, search=None, file_names=None):
        print("Preparing input files...")
        topology_files = os.listdir(os.path.join(self.settings.topology))
        topology_files = [file for file in topology_files]  
     
        print(f"Topology files for {self.settings.pdbcode}: {topology_files}")
        if search is None:
            search = self.settings.search

        if search is not None:
            _topology_files = [file 
                              for file in topology_files 
                              if search in file.split(".")[-2]]
        if len(_topology_files) == 0:
            _topology_files = [file 
                               for file in topology_files 
                               if self.settings.pdbcode in file.split(".")[-2]]

        if file_names is not None:
            _topology_files = [file for file in topology_files if file in file_names]

        if len(_topology_files) > 0:
            topology_files = _topology_files
        else:
            print("No topology files found.")
            raise FileNotFoundError

        self.topology_files = topology_files
        print(f"Loading topology files: {self.topology_files}")

    def check_all_trajectory_files(self, data_dir=None, traj_extension=None):
        print("Checking trajectory files...")
        if data_dir is None:
            data_dir = self.dirs[self.settings.data_directory]

        if traj_extension is None:
            traj_extension = self.settings.search_traj
        print(f"Checking trajectory files in: {data_dir}")

        for _dir in os.listdir(data_dir):
            if os.path.isdir(os.path.join(data_dir, _dir)):
                rep = _dir
                self.trajectories[rep] = []

                for file in os.listdir(os.path.join(data_dir, rep)):
                    if traj_extension in file and "#" not in file:
                        self.trajectories[rep].append(file)
        print(f"Trajectory files: {self.trajectories}")
        return self.trajectories

    def load_input_files(self, rep=None):
        print("Loading input files...")
        if rep is None:
            rep = self.rep_no

        if rep is None:
            raise ValueError("Replicate number not set.")

        for file in self.topology_files:
            file_path = os.path.join(self.settings.topology, file)
            destination = os.path.join(self.dirs[self.settings.data_directory],
                                       self.settings.rep_directory + str(rep),
                                       file)
            shutil.copyfile(file_path, destination)
            print(f"Copied {file} to {destination}")
  
    def set_replicate(self, rep=None):
        print("Setting replicate...")
        if rep is not None:
            self.rep_no = int(rep)
        
        if self.trajectories.get(self.settings.rep_directory + str(self.rep_no)) is None:
            self.trajectories[self.settings.rep_directory + str(self.rep_no)] = []

        print(f"Replicate number: {self.rep_no}")

    def set_trajectory_number(self, trajectory_number=None, suffix=None):
        print("Setting trajectory number...")
        if trajectory_number is None:
            self.traj_no = self.find_latest_trajectory(suffix=suffix)
        else:
            self.traj_no = int(trajectory_number)
        print(f"Trajectory number: {self.traj_no}")

    def find_latest_trajectory(self, suffix=None, rep=None):
        print("Finding latest trajectory...")
        if suffix is None:
            suffix = self.settings.suffix

        if rep is None:
            rep = self.rep_no

        if rep is None:
            raise ValueError("Replicate number not set.")

        self.check_all_trajectory_files()

        trajectory_files = self.trajectories[self.settings.rep_directory + str(rep)]
        trajectory_files = [file for file in trajectory_files if suffix in file]

        trajectory_files = sorted(trajectory_files, 
                                  key=lambda x: int(x.split("_")[-1].split(".")[0]))

        try:
            latest_traj = int(trajectory_files[-1].split("_")[-1].split(".")[0])
            print(f"Latest trajectory: {latest_traj}")
            return latest_traj
        except:
            print("No trajectories found.")
            return 0