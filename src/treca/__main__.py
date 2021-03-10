#!/usr/bin/python

import sys
import subprocess

from treca.arguments import parse_arguments
from treca.data import SequenceData

def main():
    interleaved_fastq_path, r1_fastq_path, r2_fastq_path, out_directory = parse_arguments()

    if interleaved_fastq_path is None:
        seperated_sequence_data = SequenceData(read_type="seperated", r1_read_filepath=r1_fastq_path, r2_read_filepath=r2_fastq_path, out_directory=out_directory)
        seperated_sequence_data.extract_telomeric_reads()
        seperated_sequence_data.cluster_telomeric_reads()
        seperated_sequence_data.extract_subtelomeric_reads()
        seperated_sequence_data.cluster_subtelomeric_reads()
    else:
        interleaved_sequence_data = SequenceData(read_type="interleaved", interleaved_read_filepath=interleaved_fastq_path, out_directory=out_directory)
        interleaved_sequence_data.extract_telomeric_reads()
        interleaved_sequence_data.cluster_telomeric_reads()
        interleaved_sequence_data.extract_subtelomeric_reads()
        interleaved_sequence_data.cluster_subtelomeric_reads()
if __name__ == "__main__":
    main()