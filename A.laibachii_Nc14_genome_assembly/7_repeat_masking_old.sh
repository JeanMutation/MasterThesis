#PBS -l nodes=1:ppn=10
#PBS -l walltime=10:00:00
#PBS -l mem=40gb
#PBS -S /bin/bash
#PBS -N xrepeat_model

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

out_dir='repeatmodeler_2'

mkdir -p $out_dir

sample="ref/Nc14_A.laibachii.dna.toplevel.fa"

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate repeatmodeler
echo "Running repeatmodeler..."

RepeatMasker -species "Albugo laibachii" -s -parallel 10 -xsmall -alignments $sample -dir $out_dir

BuildDatabase -name $out_dir/Hifi $sample

RepeatModeler -database $out_dir/Hifi -threads 10 -LTRStruct -ninja_dir /home/tu/tu_tu/tu_zxoyf37/bin/NINJA-0.95-cluster_only/NINJA

RepeatMasker -lib $out_dir/Hifi-families.fa -s -parallel 10 -xsmall -alignments $out_dir/Nc14_A.laibachii.dna.toplevel.fa.masked -dir $out_dir

conda deactivate
echo "repeatmodeler deactivated"
