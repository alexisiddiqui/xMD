import os
import MDAnalysis as mda

def validate_and_map_topologies(u_sim, u_ref):
    """Maps atoms between simulation and reference topologies, excluding hydrogens"""
    mapping = []
    
    for sim_res, ref_res in zip(u_sim.residues, u_ref.residues):
        # Verify residue sequence
        if sim_res.resname != ref_res.resname or sim_res.resid != ref_res.resid:
            raise ValueError(f"Residue mismatch: {sim_res.resname}{sim_res.resid} vs {ref_res.resname}{ref_res.resid}")
        
        # Map non-hydrogen atoms
        sim_atoms = {atom.name: atom for atom in sim_res.atoms if not atom.name.startswith('H')}
        ref_atoms = {atom.name: atom for atom in ref_res.atoms if not atom.name.startswith('H')}
        
        for ref_name, ref_atom in ref_atoms.items():
            if ref_name in sim_atoms:
                mapping.append((sim_atoms[ref_name].ix, ref_atom.ix))
    
    return [x[0] for x in mapping]

def strip_and_save_trajectory(traj_dir, output_topology, stride=10):
    """Strips hydrogens and saves trajectory using reference topology atom ordering"""
    # Find topology and trajectory files
    top_file = None
    traj_file = None
    
    for file in os.listdir(traj_dir):
        if file.endswith('-nojump.pdb'):
            top_file = os.path.join(traj_dir, file)
        elif file.endswith('-nojump.xtc'):
            traj_file = os.path.join(traj_dir, file)
    
    if not (top_file and traj_file):
        raise ValueError(f"Missing topology or trajectory file in {traj_dir}")
    
    # Load trajectories and map atoms
    u_sim = mda.Universe(top_file, traj_file)
    u_ref = mda.Universe(output_topology)
    atom_indices = validate_and_map_topologies(u_sim, u_ref)
    
    # Prepare output paths
    output_dir = os.path.join(os.path.dirname(traj_dir), "stripped")
    os.makedirs(output_dir, exist_ok=True)
    output_traj = os.path.join(output_dir, f"{os.path.basename(traj_dir)}_stripped.xtc")
    
    # Write stripped trajectory
    print(f"Writing stripped trajectory to {output_traj}")
    with mda.Writer(output_traj, n_atoms=len(atom_indices)) as W:
        for ts in u_sim.trajectory[::stride]:
            W.write(u_sim.atoms[atom_indices])
    
    return output_traj

def combine_trajectories(traj_paths, output_topology, output_name="combined_trajectory.xtc"):
    """Combines multiple trajectories using a single reference topology"""
    # Create output directory
    output_dir = os.path.dirname(traj_paths[0])
    output_path = os.path.join(output_dir, output_name)
    
    # Get number of atoms from reference topology
    u_ref = mda.Universe(output_topology)
    n_atoms = len([atom for atom in u_ref.atoms if not atom.name.startswith('H')])
    
    # Combine trajectories
    print(f"\nCombining {len(traj_paths)} trajectories...")
    with mda.Writer(output_path, n_atoms=n_atoms) as W:
        for traj in traj_paths:
            u = mda.Universe(output_topology, traj)
            for ts in u.trajectory:
                W.write(u.atoms)
    
    return output_path

if __name__ == "__main__":
    # Configuration
    dir_names = ["BPTI_60", "BRD4", "HOIP", "LXRa", "MBP"][:1]
    protein_names = ["BPTI", "BRD4", "HOIP", "LXR", "MBP"][:1]
    af2_pdbs = {
        "BPTI": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/BPTI_60/P00974_60_1_af_sample_127_10000_protonated_max_plddt_1050.pdb",
        "BRD4": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/BRD4/BRD4_APO_484_1_af_sample_127_10000_protonated_max_plddt_2399.pdb",
        "HOIP": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/HOIP/HOIP_apo697_1_af_sample_127_10000_protonated_max_plddt_1969.pdb",
        "LXR": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/LXRa/LXRa200_1_af_sample_127_10000_protonated_max_plddt_476.pdb",
        "MBP": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/MBP/MBP_wt_1_af_sample_127_10000_protonated_max_plddt_1791.pdb"
    }

    # Process each protein
    for dir_name, protein in zip(dir_names, protein_names):
        print(f"\nProcessing {protein}...")
        
        # Strip trajectories from each cluster
        stripped_trajs = []
        base_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}"
        
        for n in range(10):  # Process 10 clusters
            cluster_dir = os.path.join(base_dir, f"{protein}_10_c{n}/R_0")  # Assuming R_0 is the replicate directory
            print(f"\nStripping trajectory from cluster {n}...")
            stripped_traj = strip_and_save_trajectory(cluster_dir, af2_pdbs[protein])
            stripped_trajs.append(stripped_traj)
        
        # Combine stripped trajectories
        output_name = f"{protein}_all_clusters_combined.xtc"
        combined_traj = combine_trajectories(stripped_trajs, af2_pdbs[protein], output_name)
        print(f"\nCombined trajectory saved to: {combined_traj}")