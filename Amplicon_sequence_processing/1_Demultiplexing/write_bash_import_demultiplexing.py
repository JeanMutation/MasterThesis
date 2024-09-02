base_content = """
#PBS -l nodes=1:ppn=8
#PBS -l walltime=24:00:00
#PBS -l mem=112gb
#PBS -S /bin/bash
#PBS -N RunRUN_NR_demux

cd $PBS_O_WORKDIR

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-env

mv data_structure/RunRUN_NR/rawdata/Undetermined_S0_L001_I1_001.fastq.gz data_structure/RunRUN_NR/rawdata/barcodes.fastq.gz 
mv data_structure/RunRUN_NR/rawdata/Undetermined_S0_L001_R1_001.fastq.gz data_structure/RunRUN_NR/rawdata/forward.fastq.gz 
mv data_structure/RunRUN_NR/rawdata/Undetermined_S0_L001_R2_001.fastq.gz data_structure/RunRUN_NR/rawdata/reverse.fastq.gz 

qiime tools import \
  --type EMPPairedEndSequences \
  --input-path data_structure/RunRUN_NR/rawdata \
  --output-path sequences_in_work/RunRUN_NR/RUN_NR_EMP_paired_seqs.qza

qiime demux emp-paired \
  --m-barcodes-file data_structure/RunRUN_NR/RUN_NR_Mapfile_Bac.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs sequences_in_work/RunRUN_NR/RUN_NR_EMP_paired_seqs.qza \
  --o-per-sample-sequences sequences_in_work/RunRUN_NR/RUN_NR_bac_demux_sequences.qza \
  --o-error-correction-details sequences_in_work/RunRUN_NR/RUN_NR_bac_demux_error_correction_details.qza

qiime demux emp-paired \
  --m-barcodes-file data_structure/RunRUN_NR/RUN_NR_Mapfile_Fun.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs sequences_in_work/RunRUN_NR/RUN_NR_EMP_paired_seqs.qza \
  --o-per-sample-sequences sequences_in_work/RunRUN_NR/RUN_NR_fun_demux_sequences.qza \
  --o-error-correction-details sequences_in_work/RunRUN_NR/RUN_NR_fun_demux_error_correction_details.qza

qiime demux emp-paired \
  --m-barcodes-file data_structure/RunRUN_NR/RUN_NR_Mapfile_Oom.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs sequences_in_work/RunRUN_NR/RUN_NR_EMP_paired_seqs.qza \
  --o-per-sample-sequences sequences_in_work/RunRUN_NR/RUN_NR_oom_demux_sequences.qza \
  --o-error-correction-details sequences_in_work/RunRUN_NR/RUN_NR_oom_demux_error_correction_details.qza

qiime demux emp-paired \
  --m-barcodes-file data_structure/RunRUN_NR/RUN_NR_Mapfile_Euk.tsv \
  --m-barcodes-column BarcodeSequence \
  --p-rev-comp-mapping-barcodes \
  --i-seqs sequences_in_work/RunRUN_NR/RUN_NR_EMP_paired_seqs.qza \
  --o-per-sample-sequences sequences_in_work/RunRUN_NR/RUN_NR_euk_demux_sequences.qza \
  --o-error-correction-details sequences_in_work/RunRUN_NR/RUN_NR_euk_demux_error_correction_details.qza

qiime demux summarize \
  --i-data sequences_in_work/RunRUN_NR/RUN_NR_bac_demux_sequences.qza \
  --o-visualization download_demux_vis/RUN_NR_bac_demux_sequences_summary.qzv

qiime demux summarize \
  --i-data sequences_in_work/RunRUN_NR/RUN_NR_fun_demux_sequences.qza \
  --o-visualization download_demux_vis/RUN_NR_fun_demux_sequences_summary.qzv

qiime demux summarize \
  --i-data sequences_in_work/RunRUN_NR/RUN_NR_oom_demux_sequences.qza \
  --o-visualization download_demux_vis/RUN_NR_oom_demux_sequences_summary.qzv

qiime demux summarize \
  --i-data sequences_in_work/RunRUN_NR/RUN_NR_euk_demux_sequences.qza \
  --o-visualization download_demux_vis/RUN_NR_euk_demux_sequences_summary.qzv
"""

import string

modifications = [
    {"name": f"Import_and_demux_Run{letter}", "modification": letter, "placeholder": "RUN_NR"}
    for letter in string.ascii_uppercase[:15]
]

# Generate and write the scripts
for mod in modifications:
    script_content = base_content.replace(mod["placeholder"], mod["modification"])
    script_filename = mod["name"] + ".sh"
    with open(script_filename, "w") as f:
        f.write(script_content)

print("Scripts for demultiplexing generated successfully!")

generated_files = [mod["name"] + ".sh" for mod in modifications]

bash_script_content = "#!/bin/bash\n\n"

for file_name in generated_files:
    bash_script_content += f"dos2unix {file_name}\n"

# Write the bash script to a file
with open("run_dos2unix.sh", "w") as f:
    f.write(bash_script_content)

print("Scripts for format change generated successfully!")



# Write bash script to combine all files

