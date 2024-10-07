
BPTI_sciript="Run_BPTI_AF10K_RW_10.py"
HOIP_script="Run_HOIP_AF10K_RW_10.py"
LXRa_script="Run_LXRa_AF10K_RW_10.py"
MBP_script="Run_MBP_AF10K_RW_10.py"
BRD4_script="Run_BRD4_AF10K_RW_10.py"
BRD4a_script="Run_BRD4a_AF10K_RW_10.py"
BRD4b_script="Run_BRD4b_AF10K_RW_10.py"


# scripts=($BPTI_sciript $MBP_script $LXRa_script)

# # Loop through each script
# for script in "${scripts[@]}"; do
#   # Loop from 1 to 5
#   sbatch /data/chem-cat/lina4225/_home/_data/xMD/config/xMD_af10K_arg_allrep.sh $script 
  
# done


# scripts=($HOIP_script)


# # Loop through each script
# for script in "${scripts[@]}"; do
#   # Loop from 0 to 9
#   for i in {0..9}; do
#     sbatch /data/chem-cat/lina4225/_home/_data/xMD/config/xMD_af10K_arg.sh $script "$i"
#   done
# done



scripts=($BRD4_script $BRD4a_script $BRD4b_script)


# Loop through each script
for script in "${scripts[@]}"; do
  # Loop from 0 to 9
  for i in {0..5}; do
    sbatch /data/chem-cat/lina4225/_home/_data/xMD/config/xMD_af10K_arg.sh $script "$i"
  done
done