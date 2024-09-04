#PBS -l nodes=1:ppn=10
#PBS -l walltime=01:00:00
#PBS -l mem=40gb
#PBS -S /bin/bash
#PBS -N xminimap_racon

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

mkdir -p read_length_dist
mkdir -p Logs

echo "platform,length" > read_length_dist/length.csv

source miniconda3/etc/profile.d/conda.sh
conda activate bioawk
echo "Running bioawk"

bioawk -c fastx '{print "PacBio_HiFi," length($seq)}' raw_seqs/raw_reads_mapped_and_filtered_albugo.fasta >> read_length_dist/length.csv

conda deactivate

source miniconda3/etc/profile.d/conda.sh
conda activate renv-stats
echo "Running R"

Rscript vis_read_lenght_dist.R > Logs/vis_read_lenght_dist.R.log 2> Logs/vis_read_lenght_dist.R.err

conda deactivate

source miniconda3/etc/profile.d/conda.sh
conda activate assembly-stats
echo "Running assembly-stats"

assembly-stats raw_seqs/raw_reads_mapped_and_filtered.fasta >> read_length_dist/N50_stat

conda deactivate
