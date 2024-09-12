#!/usr/bin/env python3

import argparse
import pysam
import re
from tqdm import tqdm

"""
Authors:
    Bilal Sharif: bilal.bioinfo@gmail.com
    Benjamin Guinet: benjamin.guinet95@gmail.com
Usage:
    python3 adna_sslib_damage_removal.py -b input.bam -o output.bam
"""

### Parse the arguments
parser = argparse.ArgumentParser(description="Process bam file to remove C->T substitutions from single stranded ancient DNA libraries.")
parser.add_argument("-b", "--bamfile", required=True, help="Path to the input BAM file.")
parser.add_argument("-o", "--output", required=True, help="Path to the output BAM file.")

args = parser.parse_args()

bamfile = args.bamfile
output_bam_file = args.output


### Define the functions
def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reverse_complement_seq = ''.join(complement_dict.get(base, base) for base in reversed(dna_sequence))
    return reverse_complement_seq


def reconstruct_sequences(cigar_string, md_string, sequence):
    # Regular expressions to parse CIGAR and MD strings
    cigar_pattern = re.compile(r'(\d+)([MIDSHNX=])')
    md_pattern = re.compile(r'(\d+)|(\^[A-Za-z]+)|([A-Za-z])')
    # Parse the CIGAR string into a list of operations
    cigar_operations = []
    for match in cigar_pattern.finditer(cigar_string):
        length, operation = int(match.group(1)), match.group(2)
        if operation != 'H':  # Skip hard clipping (H)
            cigar_operations.append([length, operation])
    # Initialize variables for reconstruction
    seq_reference = ''
    seq_query = ''
    md_pos = 0
    seq_pos = 0
    # Parse the MD string and reconstruct the sequences
    md_tokens = []
    for match in md_pattern.finditer(md_string):
        md_tokens.append(match.group())
    md_index = 0
    for length, operation in cigar_operations:
        if operation == 'M':  # Match or mismatch
            while length > 0:
                if md_index < len(md_tokens):
                    token = md_tokens[md_index]
                    if token.isdigit():  # Match segment
                        match_length = min(length, int(token))
                        seq_reference += sequence[seq_pos:seq_pos + match_length]
                        seq_query += sequence[seq_pos:seq_pos + match_length]
                        seq_pos += match_length
                        length -= match_length
                        if match_length == int(token):
                            md_index += 1
                        else:
                            md_tokens[md_index] = str(int(token) - match_length)
                    else:  # Mismatch
                        seq_reference += token
                        seq_query += sequence[seq_pos]
                        seq_pos += 1
                        length -= 1
                        md_index += 1
        elif operation == 'I':  # Insertion
            seq_query += sequence[seq_pos:seq_pos + length]
            seq_reference += '-' * length
            seq_pos += length
        elif operation == 'D':  # Deletion
            deletion_length = length
            seq_reference += md_tokens[md_index][1:]  # Skip the '^'
            seq_query += '-' * deletion_length
            md_index += 1
    return [seq_reference, seq_query]


with pysam.AlignmentFile(bamfile, "rb") as bam, pysam.AlignmentFile(output_bam_file, "wb", header=bam.header) as out_bam:
    ## Check if the bam file is sorted
    if bam.header["HD"]["SO"] == "unsorted":
        print("Bamfile is unsorted. Only sorted bamfiles are supported. Exiting...")
        exit()
    # Initialize tqdm progress bar with the total number of reads
    total_reads = bam.mapped
    for read in tqdm(bam.fetch(), total=total_reads, desc="Processing reads"):
            ## skip the read if mapping quality is zero or read is unmapped
            if read.mapping_quality == 0 or read.is_unmapped:
                continue
            ## setting up the variables
            reverse = read.is_reverse
            read_seq = read.query_sequence
            ref_pos = read.reference_start
            MD_tag = read.get_tag("MD")
            cigar = read.cigarstring
            ## reconstruct the reference sequence from the read sequence and MD tag
            try:
                ref_reconstructed, read_seq_updated = reconstruct_sequences(cigar, MD_tag, read_seq)
            except:
                print("Error in reconstructing the reference sequence from the read sequence and MD tag. Exiting...")
                print("if you get this error and can't fix it, could you send the following information to bilal.bioinfo@gmail.com")
                print(MD_tag, cigar, ref_pos, ref_reconstructed, read_seq, read_seq_updated, reverse)
                exit()

            ## remove C->T substitution
            if not reverse:
                for base in range(len(ref_reconstructed)):
                    if ref_reconstructed[base] == "C" and read_seq_updated[base] == "T":
                        read_seq_updated = read_seq_updated[:base] + "N" + read_seq_updated[base+1:]
            else:
                for base in range(len(ref_reconstructed)):
                    if ref_reconstructed[base] == "G" and read_seq_updated[base] == "A":
                        read_seq_updated = read_seq_updated[:base] + "N" + read_seq_updated[base+1:]

            if len(ref_reconstructed) != len(read_seq_updated):
                print("Length of reference resconstructed and read sequence do not match. Exiting...")
                print("if you get this error and can't fix it, could you send the following information to bilal.bioinfo@gmail.com")
                print(MD_tag, cigar, ref_pos, ref_reconstructed, read_seq, read_seq_updated, reverse)
                exit()

            ## update the read sequence in the read object and remove - from the read sequence
            read.query_sequence = read_seq_updated.replace("-", "")

            ## write the read to the output bam file
            out_bam.write(read)

