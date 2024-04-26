
BPTI_sciript="Run_BPTI_AF10K_HighConf.py"
HOIP_script="Run_HOIP_AF10K_HighConf.py"
LXRa_script="Run_LXRa_AF10K_HighConf.py"
MBP_script="Run_MBP_AF10K_HighConf.py"



for script in $BPTI_sciript $HOIP_script $LXRa_script $MBP_script
do
    SBATCH /home/alexi/Documents/xMD/config/xMD_af10K_arg.sh $script
done