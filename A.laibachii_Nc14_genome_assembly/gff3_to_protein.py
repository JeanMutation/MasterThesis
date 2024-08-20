from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF

def read_gff(gff_file):
    """Read the GFF3 file and extract gene coordinates."""
    genes = {}
    with open(gff_file) as handle:
        for rec in GFF.parse(handle):
            for feature in rec.features:
                if feature.type == "gene":
                    gene_id = feature.id
                    # Extract relevant part of seqid for matching
                    seqid = rec.id.split()[0]
                    genes[gene_id] = {
                        "seqid": seqid,
                        "start": feature.location.start.position,
                        "end": feature.location.end.position,
                        "strand": feature.strand
                    }
    return genes

def extract_sequences(fasta_file, genes):
    """Extract nucleotide sequences from the masked FASTA file using gene coordinates."""
    seq_records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    gene_sequences = {}
    
    for gene_id, info in genes.items():
        seq_record = None
        for key in seq_records.keys():
            if info["seqid"] in key:
                seq_record = seq_records[key]
                break

        if not seq_record:
            #print(f"Warning: Sequence ID {info['seqid']} not found in FASTA file.")
            continue
        
        gene_seq = seq_record.seq[info["start"]:info["end"]]
        if info["strand"] == -1:
            gene_seq = gene_seq.reverse_complement()
        gene_sequences[gene_id] = gene_seq
    
    return gene_sequences

def translate_sequences(gene_sequences):
    """Translate nucleotide sequences into protein sequences."""
    protein_sequences = {}
    
    for gene_id, nucleotide_seq in gene_sequences.items():
        protein_seq = nucleotide_seq.translate(to_stop=True)
        protein_sequences[gene_id] = protein_seq
    
    return protein_sequences

def write_protein_fasta(protein_sequences, output_file):
    """Write protein sequences to a FASTA file."""
    seq_records = []
    
    for gene_id, protein_seq in protein_sequences.items():
        seq_record = SeqRecord(protein_seq, id=gene_id, description="")
        seq_records.append(seq_record)
    
    with open(output_file, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")

def main(gff_file, fasta_file, output_file):
    genes = read_gff(gff_file)
    gene_sequences = extract_sequences(fasta_file, genes)
    protein_sequences = translate_sequences(gene_sequences)
    write_protein_fasta(protein_sequences, output_file)

def modify_fasta_headers(input_fasta, output_fasta):
    """Modify FASTA headers by replacing spaces with underscores."""
    with open(input_fasta) as in_handle, open(output_fasta, "w") as out_handle:
        for record in SeqIO.parse(in_handle, "fasta"):
            # Replace spaces with underscores in the full header
            record.description = record.description.replace(' ', '_')
            record.id = record.description  # Set the ID to the modified description
            SeqIO.write(record, out_handle, "fasta")

if __name__ == "__main__":
    gff_file = "Files/old_assembly_annotated.gff3"
    fasta_file = "Files/Nc14_A.laibachii.dna.toplevel.fa.masked.masked"
    fasta_file_transformed = "Files/Nc14_fasta_with_fitting_headers.fa.masked.masked"
    output_file = "Files/protein_old_assembly.fasta"
    
    fasta_file = modify_fasta_headers(fasta_file, fasta_file_transformed)


    main(gff_file, fasta_file_transformed, output_file)