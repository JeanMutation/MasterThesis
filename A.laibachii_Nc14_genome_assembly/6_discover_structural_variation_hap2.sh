#PBS -l nodes=1:ppn=10
#PBS -l walltime=01:00:00
#PBS -l mem=40gb
#PBS -S /bin/bash
#PBS -N xstructural_variation

cd $PBS_O_WORKDIR

# Debugging information
echo "Starting job script"
echo "Current working directory: $PBS_O_WORKDIR"
echo "Node: $(hostname)"
echo "Job ID: $PBS_JOBID"
echo "User: $USER"
echo "Start time: $(date)"

out_dir="structural_variation_hap2"
mkdir -p $out_dir

# Make the scaffolds
source miniconda3/etc/profile.d/conda.sh
conda activate svim
echo "Running svim..."

#Conduct read-based SV calling using SVIM
svim reads --cores 10 --aligner minimap2 $out_dir raw_seqs/raw_reads_mapped_and_filtered_albugo.fasta hifiasm_assembly/Nc14_asm.bp.hap2.p_ctg.fa

conda deactivate
echo "SVIM deactivated"

source miniconda3/etc/profile.d/conda.sh
conda activate minimap2
echo "Running minimap2..."

minimap2 -a -x asm5 --cs -r2k -t 10 ref/Nc14_A.laibachii.dna.toplevel.fa hifiasm_assembly/Nc14_asm.bp.hap2.p_ctg.fa | samtools sort -m4G -@ 10 -O BAM -o $out_dir/minimap_assembly_to_ref.bam
samtools index $out_dir/minimap_assembly_to_ref.bam

conda deactivate
echo "minimap2 deactivated"

source miniconda3/etc/profile.d/conda.sh
conda activate svim
echo "Running svim..."

# Call SVs between the reference genome and yours
svim-asm haploid $out_dir $out_dir/minimap_assembly_to_ref.bam ref/Nc14_A.laibachii.dna.toplevel.fa

conda deactivate
echo "SVIM deactivated"
