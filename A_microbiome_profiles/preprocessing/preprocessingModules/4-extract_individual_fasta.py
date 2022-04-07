from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("--INPUT_FASTA", required=True, type=str, help="Input large fasta file from Kraken.")
parser.add_argument("--OUT_DIR", required=True, type=str, help="Output to individual fasta for each genome.")
args=parser.parse_args()

INPUT_FASTA = args.INPUT_FASTA
OUT_DIR = args.OUT_DIR
#
# INPUT_FASTA = "/disk/castalia-scratch/references/kraken-split/fasta-clean/library.fna"
# OUT_DIR = "/disk/castalia-scratch/references/kraken-individual"

# INPUT_FASTA = "../data/library.fna"
# OUT_DIR = "../data/individual_fasta"

#>kraken:taxid|648|NZ_CP025705.1 Aeromonas caviae strain R25-6 chromosome, complete genome

if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)

for record in SeqIO.parse(INPUT_FASTA, "fasta"):
    sequence_attributes = record.description.split("_")
    print(sequence_attributes)
    GENOME_NAME = sequence_attributes[1] + "_" + sequence_attributes[2] + "_" + sequence_attributes[4]

    OUT_FASTA=OUT_DIR + "/" + GENOME_NAME + ".fna"

    with open(OUT_FASTA, "w") as f:
        SeqIO.write(record, f, "fasta")
