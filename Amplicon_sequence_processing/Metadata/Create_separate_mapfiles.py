import pandas as pd

# Read in the dataframe
mapfile = pd.read_csv('Mapfile_Euro.txt',sep='\t')

run = input('From which run is the mapfile? :')

if mapfile['#SampleID'].str.startswith(run).all():
    print('All Sample IDs are good')
else:
    # Count the number of rows where the condition is not met
    count_not_start_with_run = (~mapfile['#SampleID'].str.startswith(run)).sum()
    print(f'{count_not_start_with_run} from {len(mapfile)} rows do not start with the letter {run}.')
    change = input('Do you want to change the Sampe IDs? (y/n)')
    if change == 'y':
        mapfile['#SampleID'] = run + '.' + mapfile['#SampleID']

mapfile_bac = mapfile[mapfile['#SampleID'].str.contains('Bac')]
mapfile_fun = mapfile[mapfile['#SampleID'].str.contains('F')]
mapfile_oom = mapfile[mapfile['#SampleID'].str.contains('O')]
mapfile_euk = mapfile[mapfile['#SampleID'].str.contains('PrV')]

filename_bac = f"{run}_Mapfile_Bac.tsv"
filename_fun = f"{run}_Mapfile_Fun.tsv"
filename_oom = f"{run}_Mapfile_Oom.tsv"
filename_euk = f"{run}_Mapfile_Euk.tsv"
filename_all = f"{run}_Mapfile_Euro_all.tsv"

# Save each dataframe to a tab-separated .tsv file
mapfile_bac.to_csv(filename_bac, sep='\t', index=False)
mapfile_fun.to_csv(filename_fun, sep='\t', index=False)
mapfile_oom.to_csv(filename_oom, sep='\t', index=False)
mapfile_euk.to_csv(filename_euk, sep='\t', index=False)

mapfile.to_csv(filename_all, sep='\t', index=False)

print("Dataframes saved successfully.")