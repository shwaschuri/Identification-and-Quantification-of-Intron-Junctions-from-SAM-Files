#Used ChatGPT for comments and to debug functions


import sys  # For system-specific parameters and functions
import argparse  # For command-line argument parsing
import logging  # For logging events during script execution
import gffutils  # For working with GFF3 files
import os  # For interacting with the operating system
from Bio.Seq import Seq  # For sequence manipulation
import matplotlib.pyplot as plt  # For creating visualizations

# Set up argument parser to handle input parameters from the command line
parser = argparse.ArgumentParser(description='Filter VCF file')
parser.add_argument('--vcfFile', required=True, help='Input VCF file')  # Input VCF file
parser.add_argument('--gffFile', required=True, help='GFF file')  # GFF annotation file
parser.add_argument('--fastafile', required=True, help='FASTA file')  # Reference genome FASTA file
parser.add_argument('--output', required=True, help='Prefix for output files')  # Output file prefix
args = parser.parse_args()

# Assign input parameters to variables
vcf_file = args.vcfFile
gff_file = args.gffFile
fasta_file = args.fastafile
output = args.output

# Configure logging to write logs to a file
logging.basicConfig(level=logging.INFO, filename='VCFlog.log', filemode='w',
                    format='%(asctime)s-%(levelname)s-%(message)s')

# Log the current working directory and input FASTA file
logging.info(f"Current working directory: {os.getcwd()}")
logging.info(f"fasta file: {fasta_file}")

# Function to count variants with low quality (QUAL <= 20)
def count_low_quality_variants(vcf_file):
    low_qual_count = 0  # Initialize count of low-quality variants

    try:
        with open(vcf_file, "r") as file:
            for line in file:
                if line.startswith("#"):
                    continue  # Skip header lines
                fields = line.strip().split("\t")  # Split the line into fields
                qual = fields[5]  # QUAL field is the 6th column in VCF

                if qual != "." and float(qual) <= 20:  # Check if QUAL is valid and <= 20
                    low_qual_count += 1  # Increment the count
        logging.info(f"Processed file: {vcf_file}")
        logging.info(f"Number of variants with QUAL <= 20: {low_qual_count}")
    except Exception as e:
        logging.error(f"Error processing file {vcf_file}: {e}")

# Function to categorize variants based on their location and effect
def categorize_variants(vcf_file, gff_file):
    try:
        db_path = "annotations1.db"  # Path for the GFF database
        if not os.path.exists(db_path):
            # Create a new GFF database if it doesn't exist
            gffutils.create_db(gff_file, dbfn=db_path, force=True, keep_order=True)
        db = gffutils.FeatureDB(db_path)  # Connect to the database

        counts = {'non_coding': 0, 'synonymous': 0, 'non_synonymous': 0}  # Initialize counts

        with open(vcf_file, "r") as vcf:
            for line in vcf:
                if line.startswith("#"):
                    continue  # Skip header lines
                fields = line.strip().split("\t")
                chrom = fields[0]  # Chromosome
                pos = int(fields[1])  # Position
                ref = fields[3]  # Reference base
                alt = fields[4]  # Alternate base
                qual = fields[5]  # Quality score

                if qual != "." and float(qual) > 20:  # Process variants with QUAL > 20
                    region_features = db.region(seqid=chrom, start=pos, end=pos, completely_within=False)
                    is_coding = False  # Track whether the variant is in a coding region

                    for feature in region_features:
                        if feature.featuretype == "CDS":  # Check if the region is a CDS
                            is_coding = True
                            cds_seq = feature.sequence(fasta=fasta_file)  # Get coding sequence

                            # Calculate position within the CDS
                            coding_pos = pos - feature.start if feature.strand == '+' else feature.end - pos

                            if coding_pos < 0 or coding_pos >= len(cds_seq):
                                continue  # Skip if position is outside the CDS

                            # Translate reference and alternate codons
                            codon_start = (coding_pos // 3) * 3
                            ref_codon = cds_seq[codon_start:codon_start + 3]
                            alt_codon = ref_codon[:coding_pos % 3] + alt + ref_codon[(coding_pos % 3) + 1:]

                            ref_aa = str(Seq(ref_codon).translate())
                            alt_aa = str(Seq(alt_codon).translate())

                            if ref_aa == alt_aa:
                                counts['synonymous'] += 1  # Synonymous mutation
                            else:
                                counts['non_synonymous'] += 1  # Non-synonymous mutation
                            break

                    if not is_coding:
                        counts['non_coding'] += 1  # Non-coding region

        logging.info("Variant categorization completed.")
        logging.info(f"Variant counts: {counts}")
        return counts
    except Exception as e:
        logging.error(f"Error processing files: {e}")
        sys.exit(1)

# Function to create a bar plot of variant categories
def plot_variant_categories(variant_counts, output_path="variant_categories_barplot.png"):
    try:
        categories = ['Non-coding', 'Synonymous', 'Non-synonymous']
        counts = [variant_counts.get('non_coding', 0),
                  variant_counts.get('synonymous', 0),
                  variant_counts.get('non_synonymous', 0)]

        total = sum(counts)
        proportions = [count / total for count in counts]  # Calculate proportions

        plt.figure(figsize=(8, 6))
        plt.bar(categories, proportions, color=['blue', 'green', 'red'], alpha=0.7)
        plt.xlabel("Variant Categories", fontsize=12)
        plt.ylabel("Proportion", fontsize=12)
        plt.title("Proportion of Variants by Category", fontsize=14)
        plt.ylim(0, 1)
        for i, value in enumerate(proportions):
            plt.text(i, value + 0.02, f"{value:.2%}", ha='center', fontsize=10)

        plt.tight_layout()
        plt.savefig(output_path)  # Save the plot
        plt.close()
        logging.info(f"Bar plot saved to {output_path}")
    except Exception as e:
        logging.error(f"Error generating bar plot: {e}")

# Function to generate a detailed variant table
def generate_variant_table(vcf_file, gff_file, fasta_file, output_file):
    try:
        db_path = "annotations1.db"
        if not os.path.exists(db_path):
            gffutils.create_db(gff_file, dbfn=db_path, force=True, keep_order=True)
        db = gffutils.FeatureDB(db_path)

        with open(vcf_file, "r") as vcf, open(output_file, "w") as output:
            output.write("Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein Location\tRef AA\n")

            for line in vcf:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                qual = fields[5]

                if qual != "." and float(qual) > 20:
                    region_features = db.region(seqid=chrom, start=pos, end=pos, completely_within=False)
                    is_coding = False

                    for feature in region_features:
                        if feature.featuretype == "CDS":
                            is_coding = True
                            transcript_id = feature.attributes.get("ID", ["NA"])[0]
                            cds_seq = feature.sequence(fasta=fasta_file)
                            coding_pos = pos - feature.start if feature.strand == '+' else feature.end - pos

                            if coding_pos < 0 or coding_pos >= len(cds_seq):
                                continue

                            codon_start = (coding_pos // 3) * 3
                            ref_codon = cds_seq[codon_start:codon_start + 3]
                            alt_codon = ref_codon[:coding_pos % 3] + alt + ref_codon[(coding_pos % 3) + 1:]

                            if len(ref_codon) == 3 and len(alt_codon) == 3:
                                ref_aa = str(Seq(ref_codon).translate())
                                alt_aa = str(Seq(alt_codon).translate())
                                protein_location = str((coding_pos // 3) + 1)
                                variant_type = "Synonymous" if ref_aa == alt_aa else "Non-synonymous"

                                output.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{variant_type}\t{transcript_id}\t{protein_location}\t{ref_aa}\n")

        logging.info(f"Variant table generated successfully in {output_file}")
    except Exception as e:
        logging.error(f"Error generating variant table: {e}")
        sys.exit(1)

# Main script execution
logging.info("Script started.")
logging.info(f"VCF file: {vcf_file}")
logging.info(f"GFF file: {gff_file}")

# Process VCF and generate outputs
count_low_quality_variants(vcf_file)
variant_counts = categorize_variants(vcf_file, gff_file)
plot_variant_categories(variant_counts, output_path=f"{output}_variant_categories.png")
generate_variant_table(vcf_file, gff_file, fasta_file, f"{output}_variant_table.txt")
