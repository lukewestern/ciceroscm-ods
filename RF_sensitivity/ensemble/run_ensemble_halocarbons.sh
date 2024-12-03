#!/bin/sh
#SBATCH --job-name=ensemble_mpgases
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:40:0
#SBATCH --mem=5G
#SBATCH --account=chem007981
#SBATCH --array=1-100

source ~/.bashrc
module --force purge
eval "$(conda shell.bash hook)"
conda activate sklearn_env

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo $current_date_time

python ensemble_halocarbons.py $SLURM_ARRAY_TASK_ID

current_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo $current_date_time
