
BPTI_sciript="Run_BPTI_AF10K_HighConf.py"
HOIP_script="Run_HOIP_AF10K_HighConf.py"
LXRa_script="Run_LXRa_AF10K_HighConf.py"
MBP_script="Run_MBP_AF10K_HighConf.py"
BRD4_script="Run_BRD4_AF10K_HighConf.py"



scripts=($LXRa_script)

# Loop through each script
for script in "${scripts[@]}"; do
  # Loop from 1 to 5
  for i in {1..1}; do
    sbatch /data/chem-cat/lina4225/_home/_data/xMD/config/xMD_af10K_arg.sh $script "$i"
  done
done