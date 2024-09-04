#PBS -l nodes=1:ppn=20
#PBS -l walltime=24:00:00
#PBS -l mem=120gb
#PBS -S /bin/bash
#PBS -N xhifiasm

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

mkdir -p hifiasm_assembly

# Load conda environment
source miniconda3/etc/profile.d/conda.sh
conda activate canu
echo "Conda environment activated"

canu -o hifiasm_assembly/Nc14_asm -t 20 raw_seqs/raw_reads_mapped_and_filtered.fasta
conda deactivate

source miniconda3/etc/profile.d/conda.sh
conda activate assembly-stats
echo "Running assembly-stats"

assembly-stats hifiasm_assembly/Nc14_asm.contigs.fasta >> hifiasm_assembly/N50_stat

conda deactivate 
