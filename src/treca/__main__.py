#!/usr/bin/python

import sys
import subprocess

from treca.arguments import parse_arguments
from treca.extraction import SequenceData

def main():
    interleaved_fastq_path, r1_fastq_path, r2_fastq_path, out_directory = parse_arguments()

    if interleaved_fastq_path is None:
        r1_sequence_data = SequenceData(r1_fastq_path, "forward", out_directory)
        r2_sequence_data = SequenceData(r2_fastq_path, "reverse", out_directory)
        
        r1_sequence_data.extract_telomeric_reads()
        r2_sequence_data.extract_telomeric_reads()
    else:
        interleaved_sequence_data = SequenceData(interleaved_fastq_path, "interleaved", out_directory)
        interleaved_sequence_data.extract_telomeric_reads()

if __name__ == "__main__":
    main()