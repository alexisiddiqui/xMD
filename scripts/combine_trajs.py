import os
import MDAnalysis as mda
def get_atom_name_variants():
    """Define comprehensive atom name mappings for different file formats and residues"""
    return {
        # Backbone atoms
        'N': ['N', 'NH1', 'NH2'],
        'CA': ['CA'],
        'C': ['C', 'CH3'],
        'O': ['O', 'O1', 'OC1', 'OT1', 'OXT', 'OW', 'OH'],
        
        # Terminal oxygens
        'OC1': ['O', 'O1', 'OC1', 'OT1', 'OXT'],
        'OC2': ['OC2', 'O2', 'OT2'],
        
        # ILE-specific
        'CD1': ['CD1', 'CD'],
        'CD': ['CD1', 'CD'],
        'CG1': ['CG1', 'CG'],
        'CG2': ['CG2'],
        
        # LEU-specific
        'CD1': ['CD1', 'CD'],
        'CD2': ['CD2', 'CD'],
        
        # ARG-specific
        'NH1': ['NH1', 'NH'],
        'NH2': ['NH2', 'NH'],
        
        # TYR-specific
        'CD1': ['CD1', 'CD'],
        'CD2': ['CD2', 'CD'],
        'OH': ['OH', 'OY'],
        
        # PHE-specific
        'CD1': ['CD1', 'CD'],
        'CD2': ['CD2', 'CD'],
        
        # Common hydrogens
        'H': ['H', 'HN', 'H1', 'HT1'],
        'HA': ['HA', 'HA1', 'HA2', 'HA3'],
        'HB': ['HB', 'HB1', 'HB2', 'HB3']
    }

def validate_and_map_topologies(u_sim, u_ref):
    """Validates topology matching with expanded residue-specific atom naming"""
    def get_residue_specific_variants(residue_name, atom_name):
        """Get atom name variants specific to residue type"""
        variants = get_atom_name_variants()
        
        # Handle special cases for terminal oxygens
        if atom_name in ['O', 'OC1', 'OC2']:
            return variants.get(atom_name, [atom_name])
            
        # Handle CD1/CD2 differently for different residue types
        if atom_name in ['CD1', 'CD2', 'CD']:
            if residue_name in ['ILE']:
                return ['CD1', 'CD'] if atom_name in ['CD1', 'CD'] else ['CD2']
            elif residue_name in ['LEU']:
                return ['CD1', 'CD'] if atom_name == 'CD1' else ['CD2', 'CD']
            elif residue_name in ['PHE', 'TYR']:
                return ['CD1', 'CD'] if atom_name == 'CD1' else ['CD2', 'CD']
        
        return variants.get(atom_name, [atom_name])

    # Compare residue sequences
    sim_residues = [(r.resname, r.resid) for r in u_sim.residues]
    ref_residues = [(r.resname, r.resid) for r in u_ref.residues]
    
    if sim_residues != ref_residues:
        mismatches = [f"Position {i}: Sim={sim} Ref={ref}" 
                     for i, (sim, ref) in enumerate(zip(sim_residues, ref_residues)) 
                     if sim != ref]
        raise ValueError(f"Residue sequences mismatch:\n" + "\n".join(mismatches[:5]))
    
    mapping = []
    missing_atoms = []
    
    for sim_res, ref_res in zip(u_sim.residues, u_ref.residues):
        sim_atoms = {atom.name: atom for atom in sim_res.atoms}
        ref_atoms = {atom.name: atom for atom in ref_res.atoms}
        
        is_first = sim_res == u_sim.residues[0]
        is_last = sim_res == u_sim.residues[-1]
        
        for ref_name, ref_atom in ref_atoms.items():
            mapped = False
            variants = get_residue_specific_variants(ref_res.resname, ref_name)
            
            # Try mapping using variants
            for variant in variants:
                if variant in sim_atoms:
                    mapping.append((sim_atoms[variant].ix, ref_atom.ix))
                    mapped = True
                    break
            
            # Record missing atoms
            if not mapped and not (is_last and ref_name in ['O', 'OC1', 'OC2', 'OXT']):
                missing_atoms.append((ref_res.resname, ref_res.resid, ref_name))
    
    # Sort mapping by reference index
    mapping.sort(key=lambda x: x[1])
    
    # Report missing atoms
    if missing_atoms:
        by_residue = {}
        for res_name, res_id, atom_name in missing_atoms:
            key = (res_name, res_id)
            by_residue.setdefault(key, []).append(atom_name)
        
        print("\nMissing atoms by residue:")
        for (res_name, res_id), atoms in sorted(by_residue.items()):
            print(f"  {res_name} {res_id}: {', '.join(sorted(atoms))}")
    
    return [x[0] for x in mapping]

def concatenate_stripped_trajs(MD_dir, output_topology, stride=10, output_dir=None, traj_index=6, traj_suffix='-nojump.xtc', replicate_dir_pref="R_", top_suffix="-nojump.pdb"):
    """
    Concatenates multiple molecular dynamics trajectories.
    """
    print(f"\nInitiating trajectory concatenation from directory: {MD_dir}")
    print(f"Using stride: {stride}")
    
    if output_dir is None:
        output_dir = os.path.dirname(MD_dir)
    
    output_dir = os.path.join(output_dir, "stripped")
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.basename(MD_dir)
    output_traj_name = base_name + "_concatenated_stripped.xtc"
    output_top_name = base_name + "_concatenated_stripped.pdb"
    output_traj_path = os.path.join(output_dir, output_traj_name)
    output_top_path = os.path.join(output_dir, output_top_name)
    
    replicate_dirs = []
    for dir in os.listdir(MD_dir):
        if os.path.isdir(os.path.join(MD_dir, dir)) and replicate_dir_pref in dir:
            replicate_dirs.append(os.path.join(MD_dir, dir))
    
    if len(replicate_dirs) == 0:
        raise ValueError(f"No replicate directories found in {MD_dir}")
    replicate_dirs.sort()
    
    top_dir = replicate_dirs[0]
    top_file = None

    for file in os.listdir(top_dir):
        if top_suffix in file:
            top_file = file
            top_path = os.path.join(top_dir, top_file)
            break
    
    if top_file is None:
        raise ValueError(f"No topology file found with suffix '{top_suffix}' in {top_dir}")
    
    traj_paths = []
    for dir in replicate_dirs:
        for file in os.listdir(dir):
            if file.endswith(traj_suffix):
                traj_paths.append(os.path.join(dir, file))
                break
    
    if len(traj_paths) == 0:
        raise ValueError(f"No trajectory files found in replicate directories")

    print("\nLoading trajectories...")
    u_sim = mda.Universe(top_path, traj_paths)
    u_ref = mda.Universe(output_topology)
    
    print("\nValidating residue sequences and creating atom mapping...")
    try:
        atom_indices = validate_and_map_topologies(u_sim, u_ref)
    except ValueError as e:
        print("Error: Topology validation failed!")
        print(str(e))
        raise
    
    # Create a new AtomGroup with the mapped atoms
    reordered_atoms = u_sim.atoms[atom_indices]
    
    # Create a new Universe for the reordered atoms
    temp_pdb = "temp_reordered.pdb"
    reordered_atoms.write(temp_pdb)
    u_reordered = mda.Universe(temp_pdb)
    
    # Set segment ID
    u_reordered.segments[0].segid = u_ref.segments[0].segid
    u_reordered.atoms.write(output_top_path)
    
    # Clean up temporary file
    os.remove(temp_pdb)

    print(f"Topology file saved with {len(atom_indices)} atoms")
    
    frames_to_write = len(range(0, len(u_sim.trajectory), stride))
    print(f"\nWriting trajectory with {frames_to_write} frames...")
    
    with mda.Writer(output_traj_path, n_atoms=len(atom_indices)) as W:
        for i, ts in enumerate(u_sim.trajectory[::stride]):
            print(f"\rProgress: {i+1}/{frames_to_write} frames written", end="")
            W.write(reordered_atoms)
    
    print("\nTrajectory writing complete!")
    return output_traj_path, output_top_path

def combine_trajs(traj_paths, top_path, stride=None, output_dir=None, file_suffix="combined"):
    """
    Concatenates multiple trajectories that have already been stripped and reordered into a single trajectory file.
    Uses the corresponding topology files for each trajectory.
    
    Parameters:
    -----------
    traj_paths : list
        List of paths to the trajectory files to combine
    top_path : list
        List of path to the reference topology files
    stride : int, optional
        Step size for reading frames (default: None, uses all frames)
    output_dir : str, optional
        Directory for saving output (default: directory of first trajectory)
    file_suffix : str, optional
        Suffix to add to the output filename (default: "combined")
    
    Returns:
    --------
    tuple
        Paths to the output trajectory and topology files
    """
    from tqdm.auto import tqdm
    print(f"\nInitiating trajectory combination for {len(traj_paths)} trajectories")
    
    # Set up output directory and filename
    if output_dir is None:
        output_dir = os.path.dirname(traj_paths[0])
    
    os.makedirs(output_dir, exist_ok=True)
    
    common_prefix = os.path.commonprefix([os.path.basename(p) for p in traj_paths])
    if not common_prefix:
        common_prefix = "trajectories"
    
    output_traj_name = f"{common_prefix}_{file_suffix}.xtc"
    output_top_name = f"{common_prefix}_{file_suffix}.pdb"
    output_traj_path = os.path.join(output_dir, output_traj_name)
    output_top_path = os.path.join(output_dir, output_top_name)

    try:
        # Load reference topology and first trajectory to get atom count
        u = mda.Universe(top_path, traj_paths)
        n_atoms = len(u.atoms)
        
        # Copy reference topology to output location
        u.atoms.write(output_top_path)
        
        # Calculate total number of frames
        total_frames = sum(len(mda.Universe(top_path, traj).trajectory[::stride] if stride else mda.Universe(top_path, traj).trajectory) 
                         for traj in traj_paths)
        
        # Write combined trajectory
        with mda.Writer(output_traj_path, n_atoms=n_atoms) as W:
            with tqdm(total=total_frames, desc="Combining trajectories") as pbar:
                for traj in traj_paths:
                    temp_u = mda.Universe(top_path, traj)
                    trajectory_slice = temp_u.trajectory[::stride] if stride else temp_u.trajectory
                    
                    for ts in trajectory_slice:
                        W.write(temp_u.atoms)
                        pbar.update(1)
        
        print(f"\nTrajectory combination complete!")
        print(f"Output trajectory: {output_traj_path}")
        print(f"Output topology: {output_top_path}")
        
        return output_traj_path, output_top_path
        
    except Exception as e:
        print("\nError during trajectory combination:")
        print(str(e))
        raise


def compare_atom_composition(top_path1, top_path2, description1="First topology", description2="Second topology"):
    """
    Performs a detailed comparison of atoms between two topology files.
    Ignores segment IDs when comparing atoms.
    
    Parameters:
    -----------
    top_path1 : str
        Path to the first topology file
    top_path2 : str
        Path to the second topology file
    description1 : str, optional
        Description of the first topology file
    description2 : str, optional
        Description of the second topology file
    
    Returns:
    --------
    tuple
        Sets of atoms unique to each topology
    """
    # Load universes
    top1_universe = mda.Universe(top_path1)
    top2_universe = mda.Universe(top_path2)
    
    # Create detailed atom identifiers excluding segid
    def get_atom_details(atom):
        return (atom.resname, atom.resid, atom.name)
    
    # Get sets of atoms from both files
    top1_atoms = {get_atom_details(atom): i for i, atom in enumerate(top1_universe.atoms)}
    top2_atoms = {get_atom_details(atom): i for i, atom in enumerate(top2_universe.atoms)}
    
    # Find differences
    top1_only = set(top1_atoms.keys()) - set(top2_atoms.keys())
    top2_only = set(top2_atoms.keys()) - set(top1_atoms.keys())
    
    # Print detailed comparison
    print("\nDetailed atom comparison:")
    print(f"Total atoms in {description1}: {len(top1_atoms)}")
    print(f"Total atoms in {description2}: {len(top2_atoms)}")
    print(f"Common atoms: {len(set(top1_atoms.keys()) & set(top2_atoms.keys()))}")
    
    if top1_only:
        print(f"\nAtoms present in {description1} but missing in {description2}:")
        for atom in sorted(top1_only):
            print(f"  - Resname: {atom[0]}, Resid: {atom[1]}, Atom: {atom[2]}")
    
    if top2_only:
        print(f"\nAtoms present in {description2} but missing in {description1}:")
        for atom in sorted(top2_only):
            print(f"  - Resname: {atom[0]}, Resid: {atom[1]}, Atom: {atom[2]}")
    
    # Group missing atoms by residue
    if top1_only:
        print(f"\nMissing atoms in {description2} grouped by residue:")
        by_residue = {}
        for atom in top1_only:
            key = (atom[0], atom[1])  # resname, resid
            if key not in by_residue:
                by_residue[key] = []
            by_residue[key].append(atom[2])  # atom name
        
        for res, atoms in sorted(by_residue.items()):
            print(f"  Residue {res[0]} {res[1]}: Missing atoms {', '.join(sorted(atoms))}")
    
    return top1_only, top2_only


if __name__ == "__main__":
    print("Starting MD trajectory concatenation script...")

    dir_names = ["BPTI_60", "BRD4", "HOIP", "LXRa", "MBP"]
    protein_names = ["BPTI", "BRD4", "HOIP", "LXR", "MBP"]

    af2_pdbs = {"BPTI": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/BPTI_60/P00974_60_1_af_sample_127_10000_protonated_max_plddt_1050.pdb",
                "BRD4": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/BRD4/BRD4_APO_484_1_af_sample_127_10000_protonated_max_plddt_2399.pdb",
                "HOIP": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/HOIP/HOIP_apo697_1_af_sample_127_10000_protonated_max_plddt_1969.pdb",
                "LXR": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/LXRa/LXRa200_1_af_sample_127_10000_protonated_max_plddt_476.pdb",
                "MBP": "/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/MBP/MBP_wt_1_af_sample_127_10000_protonated_max_plddt_1791.pdb"}

    RW_10_clusters = 10
    indexes = [4]

    # for index in indexes:  
    #     dir_name = dir_names[index]
    #     protein = protein_names[index]
    
    #     # Process test directory
    #     MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_test"
    #     test_traj_path, test_top_path = concatenate_stripped_trajs(MD_dir, af2_pdbs[protein])
        
    #     # Initialize lists to store paths
    #     traj_paths = []
    #     top_paths = []
        
    #     # Process cluster directories
    #     for n in range(RW_10_clusters):
    #         MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_10_c{n}"
    #         print(f"\nProcessing cluster {n}...")
    #         traj_path, top_path = concatenate_stripped_trajs(MD_dir, af2_pdbs[protein])
    #         traj_paths.append(traj_path)
    #         top_paths.append(top_path)

    #     # Compare first cluster's topology with reference
    #     print("\nComparing first cluster topology with reference structure...")
    #     top_only, traj_only = compare_atom_composition(
    #         af2_pdbs[protein], 
    #         top_paths[0],
    #         description1="Reference structure",
    #         description2="Stripped topology"
    #     )

    #     # Combine trajectories using corresponding topology files
    #     print("\nStarting MD trajectory combination...")
    #     combined_traj_path, combined_top_path = combine_trajs(
    #         traj_paths=traj_paths,
    #         top_path=af2_pdbs[protein],

    #     )
        
    # create full length trajectories for PCA
    output_dir = '/homes/hussain/hussain-simulation_hdx/projects/xMD/data/full_length_regular_MD'

    indexes = [1]

    for index in indexes:  
        dir_name = dir_names[index]
        protein = protein_names[index]
    
        # Process test directory
        MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_test"
        test_traj_path, test_top_path = concatenate_stripped_trajs(MD_dir, af2_pdbs[protein], stride=1, output_dir=output_dir)
        
        # Initialize lists to store paths
        traj_paths = []
        top_paths = []
        
        # Process cluster directories
        for n in range(RW_10_clusters):
            MD_dir = f"/homes/hussain/hussain-simulation_hdx/projects/xMD/data/MD/{dir_name}/{protein}_10_c{n}"
            print(f"\nProcessing cluster {n}...")
            traj_path, top_path = concatenate_stripped_trajs(MD_dir, af2_pdbs[protein], stride=1, output_dir=output_dir)
            traj_paths.append(traj_path)
            top_paths.append(top_path)

        # Compare first cluster's topology with reference
        print("\nComparing first cluster topology with reference structure...")
        top_only, traj_only = compare_atom_composition(
            af2_pdbs[protein], 
            top_paths[0],
            description1="Reference structure",
            description2="Stripped topology"
        )

        # Combine trajectories using corresponding topology files
        print("\nStarting MD trajectory combination...")
        combined_traj_path, combined_top_path = combine_trajs(
            traj_paths=traj_paths,
            top_path=af2_pdbs[protein],
            output_dir=output_dir
            
        )
        
