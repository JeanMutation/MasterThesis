headers_Arabidopsis.txt#PBS -l nodes=1:ppn=10
#PBS -l walltime=01:00:00
#PBS -l mem=40gb
#PBS -S /bin/bash
#PBS -N xminimap_before_assembly

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

mkdir -p mapping_before_assembly

# Convert the fastq file into fasta file
source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate seqtk
seqtk seq -a raw_seqs/CA445-001P0001.ccs.fastq.gz > raw_seqs/raw_reads_complete.fasta
conda deactivate

# Run minimap2
source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate minimap2

minimap2 -t 10 ref/Nc14_A.laibachii.dna.toplevel.fa raw_seqs/raw_reads_complete.fasta > mapping_before_assembly/mapped_raw_reads_to_Nc14.paf
minimap2 -t 10 ref/Nc2_A.candida_genome.fna raw_seqs/raw_reads_complete.fasta > mapping_before_assembly/mapped_raw_reads_to_Nc2.paf

minimap2 -t 10 ref/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa raw_seqs/raw_reads_complete.fasta > mapping_before_assembly/mapped_raw_reads_to_Arabidopsis.paf

conda deactivate

awk '{print $1}' mapping_before_assembly/mapped_raw_reads_to_Nc14.paf > mapping_before_assembly/headers_Nc14.txt
awk '{print $1}' mapping_before_assembly/mapped_raw_reads_to_Nc2.paf > mapping_before_assembly/headers_Nc2.txt

awk '{print $1}' mapping_before_assembly/mapped_raw_reads_to_Arabidopsis.paf > mapping_before_assembly/headers_Arabidopsis.txt

cat mapping_before_assembly/headers_Nc14.txt mapping_before_assembly/headers_Nc2.txt | sort | uniq > mapping_before_assembly/unique_headers_Albugo.txt

sort mapping_before_assembly/headers_Arabidopsis.txt -o mapping_before_assembly/headers_Arabidopsis_sorted.txt

# Find the intersection between unique_headers_Albugo.txt and headers_Arabidopsis.txt
comm -12 mapping_before_assembly/unique_headers_Albugo.txt mapping_before_assembly/headers_Arabidopsis_sorted.txt > mapping_before_assembly/intersection_headers.txt

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate python_bio

python filter_fasta.py

conda deactivate
