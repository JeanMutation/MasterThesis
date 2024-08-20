from Bio import SeqIO
import matplotlib.pyplot as plt

# Input files
unique_headers_file = "mapping_before_assembly/unique_headers_Albugo.txt"
headers_Arabidopsis = "mapping_before_assembly/headers_Arabidopsis.txt"
input_fasta_file = "raw_seqs/raw_reads_complete.fasta"
output_fasta_file_included = "raw_seqs/raw_reads_mapped_and_filtered_albugo.fasta"
output_fasta_file_arabidopsis = "raw_seqs/raw_reads_mapped_and_filtered_arabidopsis.fasta"
output_fasta_file_excluded = "raw_seqs/raw_reads_not_mapped.fasta"

# Read unique headers into a set
with open(unique_headers_file, 'r') as f:
    unique_headers = set(line.strip() for line in f)

with open(headers_Arabidopsis, 'r') as f:
    headers_Arabidopsis = set(line.strip() for line in f)

# Initialize counters
mapped_count = 0
arabidopsis_count = 0
unmapped_count = 0

# Filter the FASTA file and count the reads
with open(input_fasta_file, 'r') as input_handle, \
     open(output_fasta_file_included, 'w') as output_handle_included, \
     open(output_fasta_file_included, 'w') as output_handle_arabidopsis, \
     open(output_fasta_file_excluded, 'w') as output_handle_excluded:

    for record in SeqIO.parse(input_handle, "fasta"):
        if record.id in unique_headers:
            SeqIO.write(record, output_handle_included, "fasta")
            mapped_count += 1
        elif record.id in headers_Arabidopsis:
            SeqIO.write(record, output_handle_arabidopsis, "fasta")
            arabidopsis_count += 1
        else:
            SeqIO.write(record, output_handle_excluded, "fasta")
            unmapped_count += 1

# Print counts for debugging
print(f"Mapped reads: {mapped_count}")
print(f"Mapped reads: {arabidopsis_count}")
print(f"Unmapped reads: {unmapped_count}")

# Plotting the bar graph
labels = ['Mapped Reads Albugo','Mapped Reads Arabidopsis', 'Unmapped Reads']
counts = [mapped_count, arabidopsis_count, unmapped_count]

plt.figure(figsize=(8, 6))
plt.bar(labels, counts, color=['blue', 'green', 'red'])
plt.xlabel('Read Type')
plt.ylabel('Count')
plt.title('Mapped reads overview')
plt.savefig('read_counts_bargraph.png')
plt.show()
