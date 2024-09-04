import pandas as pd
import string

mapfile = pd.read_csv('Approach4/Mapfiles/Mapfile.txt',sep='\t')
orientation_map = pd.read_csv('Approach4/Mapfiles/Mapfile_complete_and_corrected.csv', sep=';')

primers = (mapfile['#SampleID'].apply(lambda x: x.split('.')[0]).unique())
letters = string.ascii_uppercase[:15]

base_content_further_process = """
#PBS -l nodes=1:ppn=10
#PBS -l walltime=47:00:00
#PBS -l mem=480gb
#PBS -S /bin/bash
#PBS -N Further_proc_PRIMER

cd $PBS_O_WORKDIR

source miniconda3/etc/profile.d/conda.sh
conda activate qiime2-env

qiime vsearch cluster-features-de-novo \
  --i-table sequences_in_work_ap4/all_PRIMER_denoised_table_ap6.qza \
  --i-sequences sequences_in_work_ap4/all_PRIMER_rep_seqs_ap6.qza \
  --p-perc-identity 0.97 \
  --p-threads 10 \
  --o-clustered-table sequences_in_work_ap4/all_PRIMER_otu_table_ap6.qza \
  --o-clustered-sequences sequences_in_work_ap4/all_PRIMER_otu_seq_ap6.qza

qiime feature-table summarize \
  --i-table sequences_in_work_ap4/all_PRIMER_otu_table_ap6.qza \
  --o-visualization download_otu_frequency_per_run_ap6/all_PRIMER_denoised_table_ap6.qzv \
  --m-sample-metadata-file metadata/Mapfile_KINGDOM.tsv

qiime feature-table tabulate-seqs \
  --i-data sequences_in_work_ap4/all_PRIMER_otu_seq_ap6.qza \
  --o-visualization download_otu_frequency_per_run_ap6/all_PRIMER_otu_seq_ap6.qzv

qiime tools export \
  --input-path sequences_in_work_ap4/all_PRIMER_otu_table_ap6.qza \
  --output-path download_export_otu_frequency_ap6/exported_PRIMER

qiime feature-classifier classify-sklearn \
  --i-classifier classifier/CLASSIFIER \
  --i-reads sequences_in_work_ap4/all_PRIMER_otu_seq_ap6.qza \
  --p-n-jobs 8 \
  --p-reads-per-batch 10000 \
  --p-pre-dispatch 1*n_jobs \
  --o-classification sequences_in_work_ap4/all_PRIMER_otu_taxonomy_ap6.qza

qiime metadata tabulate \
  --m-input-file sequences_in_work_ap4/all_PRIMER_otu_taxonomy_ap6.qza  \
  --o-visualization download_otu_taxonomy_ap6/all_PRIMER_otu_taxonomy_ap6.qzv
"""

base_content_dos2unix_script = '#!/bin/bash\n\n'

base_content_submission_script = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N submit_all_further

cd /beegfs/work/tu_zxoyf37/qiime2_run/


"""

dict_kingdom_classifier = {}

for primer in primers:
    if primer.startswith('B'):
        kingdom = 'bac'
        classifier = 'silva-classifier.qza'
    elif primer.startswith('P'):
        kingdom = 'euk'
        classifier = 'pr2_classifier.qza'
    elif primer.startswith('F'):
        kingdom = 'fun'
        classifier = 'unite_classifier.qza'
    elif primer.startswith('O'):
        kingdom = 'oom'
        classifier = 'unite_classifier.qza'
    
    dict_kingdom_classifier[primer] = {'primer': kingdom, 'kingdom': kingdom, 'classifier': classifier}

for primer, data in dict_kingdom_classifier.items():
    replaced_content = base_content_further_process.replace('CLASSIFIER', data['classifier']).replace('PRIMER', data['primer']).replace('KINGDOM', str(data['kingdom']))
    script_filename = f'further_process_{primer}.sh'
    with open(f'Approach6/Further_process/{script_filename}', "w") as f:
        f.write(replaced_content)
    base_content_dos2unix_script += f"dos2unix {script_filename}\n"
    base_content_submission_script += f'qsub -q smp {script_filename}\n'

with open("Approach6/Further_process/run_dos2unix.sh", "w") as f:
    f.write(base_content_dos2unix_script)

with open("Approach6/Further_process/submit_all_further.sh", "w") as f:
    f.write(base_content_submission_script)


print('All further process files are created!')

base_content = """
#!/bin/bash\n\n


"""

for primer, data in dict_kingdom_classifier.items():
    base_content += f'biom convert -i download_export_otu_frequency_ap5/exported_{data["primer"]}/feature-table.biom \
                    -o download_export_otu_frequency_ap5/otu_table_{data["primer"]}_ap5.txt --to-tsv \n'

with open("Approach6/Further_process/Export_to_txt.sh", "w") as f:
    f.write(base_content)

print('Script created successfully!')
