import pandas as pd

kingdom_short = input('Name the short handle of your Kingdom you need the metadata for: ')

tables = []

for letter in range(ord('A'), ord('O')+1):
    filename = f"Metadatas/All_Metadatas/{chr(letter)}_Mapfile_{kingdom_short}.tsv"
    try:
        table = pd.read_csv(filename, sep='\t')
        tables.append(table)
    except FileNotFoundError:
        print(f"File {filename} not found. Skipping...")

combined_table = pd.concat(tables, ignore_index=True)

combined_table.to_csv(f"Metadatas/Mapfile_{kingdom_short}.tsv", sep='\t', index=False)

print("Tables combined successfully!")