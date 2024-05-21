# %%
# xMD testing
# import pandas as pd
import os
import sys
from xMD.xMD import xMD
from xMD.MD_Settings import GROMACS_Settings
import subprocess

settings = GROMACS_Settings()
settings.gpu = True
settings.topology = "clean_top"


# amber_14_url = "https://fch.upol.cz/ff_ol/amber14sb_OL21.ff.tar.gz"

# subprocess.run(["wget", amber_14_url])
# # subprocess.run(["tar", "-xvf", "amber14sb_OL15.ff_corrected-Na-cation-params.tar.gz", "-C", os.getcwd()])
# subprocess.run(["tar", "-xvf", "amber14sb_OL21.ff.tar.gz", "-C", os.getcwd()])
# subprocess.run(["rm  amber14sb*.tar.gz"], shell=True)


# %%
# # AMPC - https://journals.asm.org/doi/10.1128/AAC.02073-20
# ampc = "6T3D"
# # KPC-2 - BLDB: http://dx.doi.org/10.1021/ACS.JMEDCHEM.7B00158
# kpc2 = "5UL8"
# # OXA-10 - BLDB: https://www.pnas.org/doi/full/10.1073/pnas.241442898
# oxa10 = "1K55"

# structures = pd.DataFrame({"PDBID": [ampc, kpc2, oxa10], 
#                            "Name": ["AmpC", "KPC2", "OXA10"]})


amber14sb_ff_path = os.path.join(os.getcwd())

# Set the GMXLIB environment variable
os.environ["GMXLIB"] = amber14sb_ff_path



# %%
# settings = GROMACS_Settings()
settings.suffix = "APO_md"
settings.search = "APO"
# settings.config = os.path.join(settings.config, "APO_MD60") 
print(settings.config)
settings.topology = os.path.join(settings.topology,"LXRa200_1_af_sample_127_10000_protonated")
print(settings.topology)
# make sure to turn on MPI for HPC 
settings.gmx_mpi_on = True





for i in range(1,5+1):
# specify ARGS: -P, -R, -N
    try:
# specify ARGS: -P, -R, -N
        md = xMD(settings, 'LXRa_test', "LXRa", i)

        # restrants = "/home/alexi/Documents/xMD/clean_top/LXRa200_1_af_sample_127_10000_protonated/LXRa200_1_af_sample_127_10000_protonated_solv_ions.gro"
        # md.check_args()
        restrants = os.path.join(settings.topology,"LXRa200_1_af_sample_127_10000_protonated_solv_ions.gro")

        md.create_directory_structure(overwrite=True)
        md.run_experiment(search="LXRa", config_files=['2_equil.mdp', '3_equil.mdp', '4_equil.mdp', '5_equil.mdp', '6_equil.mdp', '7_relax.mdp', '8_prod.mdp'], restraints=restrants)

        save_path = md.save_experiment()

    except:
        print(f"Error in replicate {i}")
        












# # %%
# def parse_gromacs_energy_log(filepath):
#     """
#     Parse a GROMACS log file to extract energy information at each step into a DataFrame.
#     """
#    # read the file
#     with open(filepath, 'r') as f:
#         log = f.readlines()
    
#     # find the start of the energy table

#         for i, line in enumerate(log):
#             #split by any whitespace
#             split = line.split()
#             # find:            Step           Time
#             if len(split)==2 and split[0] == "Step" and split[1] == "Time":
#                 start = i
#                 break
#             #find end

#         for i, line in enumerate(log):
#             #split by any whitespace
#             split = line.split()
#             if len(split)>1 and split[0] == "Statistics":
#                 end = i
#                 break

#         # extract the table
#         table = log[start+1:end]
#     # extract the headers
#     df = pd.DataFrame()

#     return df


# log_file_path = "data/MD/1K55/APO_MD60_genvel/R_1/APO_md_1K55_0.log"
# # Parse the log file again with the updated function
# energy_df_updated = parse_gromacs_energy_log(log_file_path)


# # %%



# # %%



