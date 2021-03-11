#!/usr/bin/python

import sys
import argparse
from treca.validation import validate_file

def parse_arguments():
    # define parser for arguments
    parser = argparse.ArgumentParser(prog="treca", description="TRECA: Telomeric Read Extraction, Clustering, and Assembly")
    parser.add_argument("-i", metavar="interleaved.fastq", help="specify path to FASTQ file containing interleaved paired-end reads")
    parser.add_argument("-r1", metavar="forward_reads.fastq", help="specify path to FASTQ file containing R1 pair-end reads")
    parser.add_argument("-r2", metavar="reverse_reads.fastq", help="specify path to FASTQ file containing R2 pair-end reads")
    parser.add_argument("-o", metavar="out_directory", default=None, help="specify path to directory to store output files")
    parser.add_argument("-n", metavar="n_ratio", default=0.5, help="throw out reads with a ratio of N's greater than specified")
    args = vars(parser.parse_args())

    out_directory = args["o"]
    if out_directory is not None and out_directory[-1] == "/":
        out_directory = out_directory[0:-1]

    # print help page if improper input is revieced from command line
    if not ((args["i"] and args["r1"] is None and args["r2"] is None) or (args["i"] is None and args["r1"] and args["r2"])):
        print("error: on input: arguments must specify either '-i' or 'r1' and 'r2'", file=sys.stderr)
        parser.print_help()
        exit(1)
    
    if args["i"] and args["r1"] is None and args["r2"] is None:
        # validate interleved FASTQ file
        interleaved_fastq_path = args["r1"]
        interleaved_extension = interleaved_fastq_path.split("/")[-1].split(".")[-1]
        validate_file(interleaved_fastq_path, interleaved_extension, "interleaved")

        return interleaved_fastq_path, None, None, out_directory
    
    if args["i"] is None and args["r1"] and args["r2"]:
        # validate forward read FASTQ file
        r1_fastq_path = args["r1"]
        r1_extension = r1_fastq_path.split("/")[-1].split(".")[-1]
        validate_file(r1_fastq_path, r1_extension, "forward")

        # validate reverse read FASTQ file
        r2_fastq_path = args["r2"]
        r2_extension = r2_fastq_path.split("/")[-1].split(".")[-1]
        validate_file(r2_fastq_path, r2_extension, "reverse")

        return None, r1_fastq_path, r2_fastq_path, out_directory

if __name__ == "__main__":
    parse_arguments()