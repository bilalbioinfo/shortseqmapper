#!/usr/bin/python3
## author: Bilal Sharif

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq
from sys import argv
import random
import gzip

input_filename = argv[1]

divergence_rate = float(argv[2])
output_filename = input_filename.split('.fa.gz')[0] + f"_div.fa.gz"

seq_translate = {"A":"T","T":"A","C":"G","G":"C"}

with gzip.open(input_filename, "rt") as input_file, gzip.open(output_filename, "wt") as output_file:
    for record in SeqIO.parse(input_file, "fasta"):
        mutated_seq = MutableSeq(record.seq)
        for i in range(len(mutated_seq)):
            if random.random() < divergence_rate:
                mutated_seq[i] = seq_translate[mutated_seq[i]]
        mutated_record = SeqRecord(mutated_seq , id=record.id, description=record.description)
        SeqIO.write(mutated_record, output_file, "fasta")

print(f"The mutated sequence with divergence {divergence_rate} is written to {output_filename}")
