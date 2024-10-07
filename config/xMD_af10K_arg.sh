#!/bin/bash
#SBATCH --job-name=xMD-10K

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=48:00:00
#SBATCH --partition=medium
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=1
#SBATCH --clusters=htc

#SBATCH --output=/home/lina4225/_log_files/slurm_%j.out  # Writes standard output to this file. %j is jobnumber                             
#SBATCH --error=/home/lina4225/_log_files/slurm_%j.out   # Writes error messages to this file. %j is jobnumber


#SBATCH --mail-user=lina4225@ox.ac.uk                
#SBATCH --mail-type=begin               # Instead only email when job begins...
#S BATCH --mail-type=end                 # ... and ends
#SBATCH --mail-type=fail                # ... and fails
echo $HOME
nvidia-smi
echo $1

echo $2


HOME_DIR=${TMPDIR}

SOURCE_DIR=${DATA}/xMD
LOCAL_DIR=${TMPDIR}/xMD

# RM_LOCALCOLAB_DIR=$LOCALCOLAB_DIR-bak
# module load GROMACS/2021.5-foss-2021b-CUDA-11.4.1
module load GROMACS/2021.3-foss-2021a-CUDA-11.3.1

module load Anaconda3/2022.10

# conda activate RIN_test

source activate RIN_test

module list


# # Example project structure is:
# # my-project/code
# # my-project/data
# # my-project/results
# rm -rf ${LOCAL_DIR}
# # rm -rf ${INSTALL_DIR}

# # rm -rf ${INSTALL_DIR}

# # Create directory if needed
# mkdir -p ${LOCAL_DIR}

# # Rsync code and data from vols
# while [ 1 ]  
# do
#     rsync --partial --progress --archive ${SOURCE_DIR}/ ${LOCAL_DIR}
#     if [ "$?" = "0" ] ; then
#         echo "rsync completed normally"
#         break
#     else
#         echo "rsync failed, Retrying in 180s..."
#         sleep 180
#     fi
# done



if [ $# -eq 0 ]; then
    echo "Usage: $0 <xMD_python_script>"
    exit 1
fi

cd ${SOURCE_DIR}

# Pass both argument to the Python script
python "$1" "$2"



# # Copy results back to vols
# while [ 1 ]  
# do
#     rsync --partial --progress --archive ${LOCAL_DIR}/ ${SOURCE_DIR}
#     if [ "$?" = "0" ] ; then
#         echo "rsync completed normally"
#         break
#     else
#         echo "rsync failed, Retrying in 180s..."
#         sleep 180
#     fi
# done

# # Clean up anything that's no longer needed
# rm -rf ${LOCAL_DIR}
# pip cache purge
# conda clean --all -y
# conda deactivate

# # rm -rf ${INSTALL_DIR}