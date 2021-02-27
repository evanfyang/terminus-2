#!/usr/bin/python

import sys
import subprocess

from tqdm import tqdm
from arguments import parse_arguments
from extraction import extract_telomeric_reads

def main():
    interleaved_fastq_path, r1_fastq_path, r2_fastq_path = parse_arguments()

    if interleaved_fastq_path is None:
        # find all telomeric reads in forward read FASTQ file and write to file
        print("Extracting telomeric reads from R1 FASTQ dataset...")
        r1_tel_path, r1_fastq_len = extract_telomeric_reads(r1_fastq_path, "forward")
        print()

        # find all telomeric reads in reverse read FASTQ file and write to file
        print("Extracting telomeric reads from R2 FASTQ dataset...")
        r2_tel_path, r2_fastq_len = extract_telomeric_reads(r2_fastq_path, "reverse")
        print()
    else:
        # find all telomeric reads in forward read FASTQ file and write to file
        print("Extracting telomeric reads from interleaved FASTQ dataset...")
        r1_tel_path, r1_fastq_len = extract_telomeric_reads(interleaved_fastq_path, "interleaved")
        print()    

if __name__ == "__main__":
    main()