import pandas as pd
import os
import shutil



table = pd.read_csv('input_flashweave_complete/dataframe.tsv', sep='\t')
table.set_index(table.columns[0], inplace=True)

dir_for_run = 'keystone_analysis_run_34_data'
dir_dataframes = 'Dataframes'
dir_scripts = 'Scripts'

if not os.path.exists(dir_for_run):
    os.makedirs(dir_for_run)

subdirs = [dir_dataframes, dir_scripts]

for subdir in subdirs:
    full_path = os.path.join(dir_for_run, subdir)
    if not os.path.exists(full_path):
        os.makedirs(full_path)

print(f"Created directories: {dir_for_run}, {os.path.join(dir_for_run, dir_dataframes)}, and {os.path.join(dir_for_run, dir_scripts)}")

# Creating a dataframe with each one OTU missing and a dictionary with relevant information on the dataframes
file_info_dict = {}

for column in table.columns:
    dropped_df = table.drop(column, axis=1)
    otu_name = column
    filename = f'{column}.tsv'
    file_info_dict[filename] = otu_name
    dropped_df.to_csv(f'{dir_for_run}/{dir_dataframes}/{filename}', sep='\t', index=True)

table.to_csv(f'{dir_for_run}/{dir_dataframes}/all_otus.tsv', sep='\t', index=True)
file_info_dict['all_otus.tsv'] = 'all_otus'

## Create all files necessary to run the networks

fastspar_script = '''
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:19:00
#PBS -l mem=40gb
#PBS -S /bin/bash
#PBS -N xOTUNAME

cd $PBS_O_WORKDIR

# Activate Julia
source juliaup/etc/profile.d/julia.sh

# Run Julia commands
julia -e 'using FlashWeave;
            data_path="FOLDER_DF/FILENAME";
            netw_results = learn_network(data_path, sensitive=true, heterogeneous=false);
            save_network("Networks/wl_nw_OTUNAME.gml", netw_results);'

'''

base_content_submission_script = """
#PBS -l nodes=1:ppn=1
#PBS -l walltime=00:10:00
#PBS -l mem=32gb
#PBS -S /bin/bash
#PBS -N xsubmit_all_further

cd $PBS_O_WORKDIR

mkdir -p Networks

"""

for filename, otu_name in file_info_dict.items():
    replaced_content = fastspar_script.replace('OTUNAME', otu_name).replace('FILENAME', filename).replace('FOLDER_DF', dir_dataframes)
    script_filename = f'flashweave_{otu_name}.sh'
    with open(f'{dir_for_run}/{dir_scripts}/{script_filename}', "w") as f:
        f.write(replaced_content)
    base_content_submission_script += f'qsub -q tiny {dir_scripts}/{script_filename}\n'

with open(f"{dir_for_run}/{dir_scripts}/submit_all_flashweave.sh", "w") as f:
    f.write(base_content_submission_script)



current_directory = os.getcwd()
destination_directory = dir_for_run

# Liste der Dateien, die kopiert werden sollen
files_to_copy = ['keystone_analysis_part_2_run_all.sh', 
                 'keystone_analysis_part_2_random_cluster.py', 
                 'keystone_analysis_full_nw_analysis.py',
                 'keystone_analysis_dropout_nw_analysis.py',
                 'keystone_analysis_comparison_nw_analysis.py',
                 'keystone_analysis_part_3.R',
                 'taxonomy.tsv']

# Sicherstellen, dass das Zielverzeichnis vorhanden ist
os.makedirs(destination_directory, exist_ok=True)

# Kopiere jede Datei aus der Liste in das Zielverzeichnis
for file_name in files_to_copy:
    source_file_path = os.path.join(current_directory, file_name)
    destination_file_path = os.path.join(destination_directory, file_name)
    shutil.copy(source_file_path, destination_file_path)
