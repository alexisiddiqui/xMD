import os
import MDAnalysis as mda

def concatenate_trajs(MD_dir, stride=10, output_dir=None, traj_index=6, traj_suffix='-mol.xtc', replicate_dir_pref="R_", top_suffix="_solv_ions.gro"):
    """
    Concatenates multiple molecular dynamics trajectories from different replicates into a single trajectory file.
    
    Parameters:
    -----------
    MD_dir : str
        Path to the main directory containing replicate folders
    stride : int, optional
        Step size for reading frames (default: 10)
    output_dir : str, optional
        Directory for saving output (default: parent directory of MD_dir)
    traj_index : int, optional
        Index identifying trajectory files (default: 6)
    traj_suffix : str, optional
        Suffix identifying trajectory files (default: '-mol')
    replicate_dir_pref : str, optional
        Prefix for replicate directories (default: "R_")
    top_prefix : str, optional
        Prefix for topology file (default: "_solv_ions.gro")
    """
    print(f"\nInitiating trajectory concatenation from directory: {MD_dir}")
    print(f"Using stride: {stride}")
    
    # Set up output directory and filename
    if output_dir is None:
        output_dir = os.path.dirname(MD_dir)
    output_name = os.path.basename(MD_dir) + "_concatenated.xtc"
    output_path = os.path.join(output_dir, output_name)
    print(f"Output will be saved to: {output_path}")
    
    # Find all replicate directories
    replicate_dirs = []
    for dir in os.listdir(MD_dir):
        if os.path.isdir(os.path.join(MD_dir, dir)) and replicate_dir_pref in dir:
            replicate_dirs.append(os.path.join(MD_dir, dir))
    
    print(f"\nFound {len(replicate_dirs)} replicate directories:")
    for dir in replicate_dirs:
        print(f"  - {dir}")
    
    if len(replicate_dirs) == 0:
        raise ValueError(f"No replicate directories found in {MD_dir}")
    replicate_dirs.sort()
    # Find topology file in first replicate directory
    top_dir = replicate_dirs[0]
    top_file = None

    for file in os.listdir(top_dir):
        if top_suffix in file:
            top_file = file
            top_path = os.path.join(top_dir, top_file)
            break
    
    if top_file is None:
        raise ValueError(f"No topology file found with prefix '{top_suffix}' in {top_dir}")
    print(f"\nUsing topology file: {top_path}")
    
    traj_paths = []
    for dir in replicate_dirs:
        traj_found = False
        for file in os.listdir(dir):
            if file.endswith(traj_suffix):
                traj_paths.append(os.path.join(dir, file))
                traj_found = True
                break
        if not traj_found:
            print(f"Warning: No trajectory file found in {dir}")
            
    print(f"\nFound {len(traj_paths)} trajectory files:")
    for traj in traj_paths:
        print(f"  - {traj}")
    
    if len(traj_paths) == 0:
        raise ValueError(f"No trajectory files found in replicate directories")

    
    # Load trajectories and concatenate
    print("\nLoading trajectories with MDAnalysis...")
    u = mda.Universe(top_path, traj_paths)
    print(f"Total number of frames: {len(u.trajectory)}")
    print(f"Number of atoms: {u.trajectory.n_atoms}")
    
    # Write concatenated trajectory
    print(f"\nWriting concatenated trajectory with stride {stride}...")
    frames_to_write = len(range(0, len(u.trajectory), stride))
    with mda.Writer(output_path, n_atoms=u.trajectory.n_atoms) as W:
        for i, ts in enumerate(u.trajectory[::stride]):
            # if (i + 1) % 100 == 0:  # Progress update every 100 frames
            print(f"Progress: {i+1}/{frames_to_write} frames written")
            W.write(u.atoms)
    
    print(f"\nConcatenation complete! Output saved to: {output_path}")
    print(f"Final trajectory contains {frames_to_write} frames")

if __name__ == "__main__":
    print("Starting MD trajectory concatenation script...")



    dir_names = ["BPTI_60", "BRD4", "HOIP", "LXRa", "MBP"]

    protein_names = ["BPTI", "BRD4", "HOIP", "LXR", "MBP"]

    RW_10_clusters = 10

    # had to rename LXRa_test to LXR_test

    # for dir_name, protein in zip(dir_names[1:], protein_names[1:]):


    #     MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_test"
    #     concatenate_trajs(MD_dir)

    #     for n in range(RW_10_clusters):
    #         MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_10_c{n}"
    #         print("Starting MD trajectory concatenation script...")
    #         concatenate_trajs(MD_dir)



    indexes = [1]

    for index in indexes:  
        dir_name = dir_names[index]
        protein = protein_names[index]

        MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_test"
        concatenate_trajs(MD_dir)

        for n in range(RW_10_clusters):
            MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_10_c{n}"
            print("Starting MD trajectory concatenation script...")
            concatenate_trajs(MD_dir)

        

        
    print("Script execution completed successfully!")

