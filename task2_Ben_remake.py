import gzip
from Bio import SeqIO
from tqdm import tqdm

# Define file paths
genome_fasta = "GCF_000001405.26_GRCh38_genomic.fna.gz"
annotation_gff = "GCF_000001405.26_GRCh38_genomic.gff.gz"
target_gene = "HPSE" # From task 1

# Function to load the reference genome
def load_reference_genome(genome_fasta):
    genome = {}
    with gzip.open(genome_fasta, "rt") as handle:
        for record in tqdm(SeqIO.parse(handle, "fasta"), desc="Loading genome"):
            genome[record.id] = record.seq
    return genome

# Function to parse the GFF file, find the first exon of the gene, and return reference ID
def find_first_exon_location(gff_file, target_gene):
    exon_location = None
    found_gene = False
    reference_id = None  # Store reference/transcript ID
    with gzip.open(gff_file, "rt") as handle:
        for line in tqdm(handle, desc="Scanning GFF for first exon of target gene"):
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]
            
            # Check for the gene
            if feature_type == "gene" and (f"Name={target_gene}" in attributes or f"gene_name={target_gene}" in attributes):
                found_gene = True 
                chrom = parts[0]
                strand = parts[6]
                
                # Extract the reference or transcript ID from attributes
                if "ID=" in attributes:
                    reference_id = attributes.split("ID=")[1].split(";")[0]
                elif "transcript_id=" in attributes:
                    reference_id = attributes.split("transcript_id=")[1].split(";")[0]
                else:
                    print("Reference ID not found in attributes.")
            
            # Look for the first exon of the gene
            if found_gene and feature_type == "exon":
                exon_start = int(parts[3])
                exon_end = int(parts[4])
                exon_location = (chrom, exon_start, exon_end, strand)
                print(f"First exon found from {exon_start} to {exon_end} on strand {strand}")
                break  # Stop after finding the first exon
    return exon_location, reference_id

# Function to extract the gene sequence from genome
def extract_gene_sequence(genome, chrom, start, end, strand):
    seq = genome[chrom][start:end]
    if strand == "-":
        seq = seq.reverse_complement()
    return seq

# Load genome and find first exon for HPSE gene
genome = load_reference_genome(genome_fasta)
exon_location, reference_id = find_first_exon_location(annotation_gff, target_gene)

if exon_location and reference_id:
    chrom, start, end, strand = exon_location
    exon_seq = extract_gene_sequence(genome, chrom, start, end, strand)

    # Print the reference sequence ID and first exon sequence
    print(f"Reference sequence: {reference_id}")
    print(f"First exon sequence for {target_gene}:\n{exon_seq}\n")

    # Save the exon sequence in FASTA format, including the reference sequence ID in the header
    output_file = f"LPL_first_exon_{target_gene}_{chrom}.fasta"
    with open(output_file, "w") as f:
        f.write(f">{reference_id}_first_exon_{chrom}:{start}-{end}({strand})\n{exon_seq}\n") 
    print(f"First exon sequence saved to {output_file}")
else:
    print(f"Gene {target_gene} not found or reference ID missing.")

