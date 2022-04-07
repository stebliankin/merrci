from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--file", required=True, type=str, help="Input large fasta file from Kraken.")
parser.add_argument("--out", required=True, type=str, help="Output to individual fasta for each genome.")
args=parser.parse_args()

INPUT_FASTA = args.file
OUT_FASTA = args.out
#
# INPUT_FASTA = "/disk/castalia-scratch/references/kraken-split/fasta-clean/library.fna"
# OUT_DIR = "/disk/castalia-scratch/references/kraken-individual"

# INPUT_FASTA = "../data/library.fna"
# OUT_DIR = "../data/individual_fasta"

#>kraken:taxid|648|NZ_CP025705.1 Aeromonas caviae strain R25-6 chromosome, complete genome
# OUT:
#>1307_NZ_CP017667_1_Streptococcus_suis_strain_1081_chromosome_complete_genome

with open(OUT_FASTA, "w") as f:
    for record in SeqIO.parse(INPUT_FASTA, "fasta"):
        sequence_attributes = record.description
        tax_id = sequence_attributes.split("|")[1]
        sq_name = sequence_attributes.split("|")[2].replace(" ","_").replace(",","").replace(".","_")
        updated_record = record
        updated_record.id = tax_id + "_" + sq_name
        updated_record.description=""
        SeqIO.write(updated_record, f, "fasta")
