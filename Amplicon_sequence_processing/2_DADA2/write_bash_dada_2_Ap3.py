base_content = """
#PBS -l nodes=1:ppn=8
#PBS -l walltime=47:00:00
#PBS -l mem=112gb
#PBS -S /bin/bash
#PBS -N RunRUN_NR_KINGDOM_dada2

cd $PBS_O_WORKDIR

source /home/tu/tu_tu/tu_zxoyf37/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-env

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs sequences_in_work/RunRUN_NR/RUN_NR_KINGDOM_demux_sequences.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 20 \
  --p-trunc-len-f FOR \
  --p-trunc-len-r REV \
  --p-n-threads 8 \
  --o-table sequences_in_work/RunRUN_NR_ap3/RUN_NR_KINGDOM_denoised_table.qza \
  --o-representative-sequences sequences_in_work/RunRUN_NR_ap3/RUN_NR_KINGDOM_denoised_representativeSeqs.qza \
  --o-denoising-stats sequences_in_work/RunRUN_NR_ap3/RUN_NR_KINGDOM_denoising_stats.qza
  
qiime feature-table summarize \
  --i-table sequences_in_work/RunRUN_NR_ap3/RUN_NR_KINGDOM_denoised_table.qza \
  --o-visualization download_dada2_vis_ap3/RUN_NR_KINGDOM_denoised_table.qzv \
  --m-sample-metadata-file data_structure/RunRUN_NR/RUN_NR_Mapfile_Bac.tsv

qiime feature-table tabulate-seqs \
  --i-data sequences_in_work/RunRUN_NR_ap3/RUN_NR_KINGDOM_denoised_representativeSeqs.qza \
  --o-visualization download_dada2_vis_ap3/RUN_NR_KINGDOM_denoised_representativeSeqs.qzv

qiime metadata tabulate \
  --m-input-file sequences_in_work/RunRUN_NR_ap3/RUN_NR_KINGDOM_denoising_stats.qza \
  --o-visualization download_dada2_vis_ap3/RUN_NR_KINGDOM_denoising_stats.qzv
"""

import string
import pandas as pd

cutoffs = pd.read_csv("DADA2/Cutoffs_for_dada2.CSV", delimiter=";")

kingdom_short = ('bac', 'fun', 'oom', 'euk')
letters = string.ascii_uppercase[:15]

modifications = []
for index, row in cutoffs.iterrows():
  for short_handle_kingdom in kingdom_short:
      for letter in letters:
          columns_to_keep = ['Kingdom', 'Direction', letter]
          cutoff_pre = cutoffs[columns_to_keep]
          cutoff_pre = cutoff_pre.loc[cutoff_pre['Kingdom'] == short_handle_kingdom]
          cutoff_forward = cutoff_pre.loc[cutoff_pre['Direction'] == 'Forward']
          cutoff_reverse = cutoff_pre.loc[cutoff_pre['Direction'] == 'Reverse']
          cutoff_forward = cutoff_forward[letter].iloc[0]
          cutoff_reverse = cutoff_reverse[letter].iloc[0]
          modification = {
              "name": f"dada_2_{short_handle_kingdom}_{letter}",
              "modification": short_handle_kingdom,
              "placeholder": "KINGDOM",
              "modification2": letter,
              "placeholder2": "RUN_NR",
              "for_value": cutoff_forward,
              "rev_value": cutoff_reverse
          }
          modifications.append(modification)

# Generate and write the scripts
for mod in modifications:
    script_content = base_content.replace(mod["placeholder"], mod["modification"])
    script_content = script_content.replace(mod["placeholder2"], mod["modification2"])
    script_content = script_content.replace("FOR", str(mod["for_value"]))
    script_content = script_content.replace("REV", str(mod["rev_value"]))
    script_filename = mod["name"] + ".sh"
    with open('DADA2/Approach3/' + script_filename, "w") as f:
        f.write(script_content)

print("Scripts for dada2 generated successfully!")

generated_files = [mod["name"] + ".sh" for mod in modifications]

bash_script_content = "#!/bin/bash\n\n"

for file_name in generated_files:
    bash_script_content += f"dos2unix {file_name}\n"

# Write the bash script to a file
with open("DADA2/Approach3/run_dos2unix.sh", "w") as f:
    f.write(bash_script_content)

print("Script for format change generated successfully!")


# Generate file for job submission
bash_script_2_content = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N submit_otu_clustering_jobs

cd /beegfs/work/tu_zxoyf37/qiime2_run

"""

for file_name in generated_files:
    bash_script_2_content += f"qsub -q short {file_name}\n"

with open("DADA2/Approach3/submit_all_jobs.sh", "w") as f:
    f.write(bash_script_2_content)

print("Script for submission generated successfully!")