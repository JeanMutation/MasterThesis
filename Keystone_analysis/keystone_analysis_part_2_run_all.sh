#PBS -l nodes=1:ppn=10
#PBS -l walltime=47:00:00
#PBS -l mem=120gb
#PBS -S /bin/bash
#PBS -N xsubmit

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

# Load conda environment
source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate python
echo "Conda environment activated"

mkdir -p Networks_random
mkdir -p Results
mkdir -p Logs

echo "Running Python script"
python keystone_analysis_part_2_random_cluster.py > Logs/keystone_analysis_part_2_random_cluster.log 2> Logs/keystone_analysis_part_2_random_cluster.err &

echo 'Running network analysis'

python keystone_analysis_part_2_full_nw_analysis.py > Logs/keystone_analysis_full_nw_analysis.log 2> Logs/keystone_analysis_full_nw_analysis.err &
python keystone_analysis_part_2_comparison_nw_analysis.py > Logs/keystone_analysis_comparison_nw_analysis.log 2> Logs/keystone_analysis_comparison_nw_analysis.err &
python keystone_analysis_part_2_dropout_nw_analysis.py > Logs/keystone_analysis_dropout_nw_analysis.log 2> Logs/keystone_analysis_dropout_nw_analysis.err &

wait

echo 'All Python scripts have completed'

conda deactivate

cd $PBS_O_WORKDIR
source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate R-env

echo 'Running R script'

Rscript keystone_analysis_part_3.R > Logs/keystone_analysis.R.log 2> Logs/keystone_analysis.R.err

echo 'All scripts have completed'

echo "End time: $(date)"
echo "Job script completed"