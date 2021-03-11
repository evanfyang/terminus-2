#!/usr/bin/python

import os
import re

# define constants
HEADER = 1
SEQUENCE = 2

FORWARD = 1
REVERSE = 2

FA_LINES_PER_READ = 2
FQ_LINES_PER_READ = 4
class SequenceData:
    def __init__(self, read_type, interleaved_read_filepath=None, r1_read_filepath=None, r2_read_filepath=None, out_directory=None):
        self._read_type = read_type

        self._interleaved_reads_filepath = interleaved_read_filepath 
        self._r1_reads_filepath = r1_read_filepath 
        self._r2_reads_filepath = r2_read_filepath

        self._out_directory = out_directory

        self._r1_telomeric_reads_filepath = None
        self._r2_telomeric_reads_filepath = None

        self._r1_clustered_telomeric_reads_info_filepath = None
        self._r2_clustered_telomeric_reads_info_filepath = None

        self._r1_clustered_telomeric_reads_filepath = None
        self._r2_clustered_telomeric_reads_filepath = None

        self._r1_subtelomeric_reads_filepath = None
        self._r2_subtelomeric_reads_filepath = None

        self._r1_clustered_subtelomeric_reads_info_filepath = None
        self._r1_clustered_subtelomeric_reads_info_filepath = None

        self._r1_clustered_subtelomeric_reads_filepath = None
        self._r1_clustered_subtelomeric_reads_filepath = None
    
    def get_read_type(self):
        return self._read_type

    def get_interleaved_reads_filepath(self):
        return self._interleaved_reads_filepath

    def get_r1_reads_filepath(self):
        return self._r1_reads_filepath

    def get_r2_reads_filepath(self):
        return self._r2_reads_filepath
    
    def get_out_directory(self):
        return self._out_directory
    
    def get_r1_telomeric_reads_filepath(self):
        return self._r1_telomeric_reads_filepath
    
    def get_r2_telomeric_reads_filepath(self):
        return self._r2_telomeric_reads_filepath
    
    def get_r1_clustered_telomeric_reads_info_filepath(self):
        return self._r1_clustered_telomeric_reads_info_filepath
    
    def get_r2_clustered_telomeric_reads_info_filepath(self):
        return self._r2_clustered_telomeric_reads_info_filepath

    def get_r1_clustered_telomeric_reads_filepath(self):
        return self._r1_clustered_telomeric_reads_filepath
    
    def get_r2_clustered_telomeric_reads_filepath(self):
        return self._r2_clustered_telomeric_reads_filepath

    def get_r1_subtelomeric_reads_filepath(self):
        return self._r1_subtelomeric_reads_filepath
    
    def get_r2_subtelomeric_reads_filepath(self):
        return self._r2_subtelomeric_reads_filepath
    
    def get_r1_clustered_subtelomeric_reads_info_filepath(self):
        return self._r1_clustered_subtelomeric_reads_info_filepath
    
    def get_r2_clustered_subtelomeric_reads_info_filepath(self):
        return self._r2_clustered_subtelomeric_reads_info_filepath

    def get_r1_clustered_subtelomeric_reads_filepath(self):
        return self._r1_clustered_subtelomeric_reads_filepath
    
    def get_r2_clustered_subtelomeric_reads_filepath(self):
        return self._r2_clustered_subtelomeric_reads_filepath
    
    def extract_telomeric_reads(self):
        if self._read_type == "interleaved":
            print("Extracting telomeric reads from an interleaved FASTQ dataset...")
            if self._out_directory == None: 
                path_prefix = "/".join(self._interleaved_reads_filepath.split("/")[0:-1]) + "/telomeric_reads/"
            else:
                path_prefix = self._out_directory + "/telomeric_reads/"
            
            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)

            self._r1_telomeric_reads_filepath = path_prefix + self._interleaved_reads_filepath.split("/")[-1].split(".")[0] + "_r1_telomeric_reads.fasta"
            r1_telomeric_reads_output = open(self._r1_telomeric_reads_filepath, "w+")

            self._r2_telomeric_reads_filepath = path_prefix + self._interleaved_reads_filepath.split("/")[-1].split(".")[0] + "_r2_telomeric_reads.fasta"
            r2_telomeric_reads_output = open(self._r2_telomeric_reads_filepath, "w+")

            self._process_reads(self._interleaved_reads_filepath, r1_telomeric_reads_output, r2_telomeric_reads_output)
            
            print("Done.\n")
            print("R1 telomeric reads saved to " + self._r1_telomeric_reads_filepath)
            print("R2 telomeric reads saved to " + self._r2_telomeric_reads_filepath)
            print()
        if self._read_type == "seperated":
            print("Extracting telomeric reads from R1 FASTQ dataset...")
            
            if self._out_directory == None: 
                path_prefix = "/".join(self._r1_reads_filepath.split("/")[0:-1]) + "/telomeric_reads/"
            else:
                path_prefix = self._out_directory + "/telomeric_reads/"

            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            self._r1_telomeric_reads_filepath = path_prefix + self._r1_reads_filepath.split("/")[-1].split(".")[0] + "_r1_telomeric_reads.fasta"
            r1_telomeric_reads_output = open(self._r1_telomeric_reads_filepath, "w+")

            self._process_reads(input_filepath=self._r1_reads_filepath, r1_output=r1_telomeric_reads_output, r2_output=None)
            print("Done.\n")
            
            print("Extracting telomeric reads from R2 FASTQ dataset...")

            if self._out_directory == None: 
                path_prefix = "/".join(self._r2_reads_filepath.split("/")[0:-1]) + "/telomeric_reads/"
            else:
                path_prefix = self._out_directory + "/telomeric_reads/"

            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            self._r2_telomeric_reads_filepath = path_prefix + self._r2_reads_filepath.split("/")[-1].split(".")[0] + "_r2_telomeric_reads.fasta"
            r2_telomeric_reads_output = open(self._r2_telomeric_reads_filepath, "w+")
            
            self._process_reads(input_filepath=self._r2_reads_filepath, r1_output=None, r2_output=r2_telomeric_reads_output)
            print("Done.\n")

            print("R1 telomeric reads saved to '" + self._r1_telomeric_reads_filepath + "'")
            print("R2 telomeric reads saved to '" + self._r2_telomeric_reads_filepath + "'")
            print()
    
    def cluster_telomeric_reads(self):   
        # install and make wcd for clustering         # install and make wcd for clustering
        print("Installing and preparing wcdest for clustering...\n")
        os.system("git submodule update --init --recursive")
        os.system("../lib/wcdest/code/configure >/dev/null 2>&1")
        os.system("make clean -C ../lib/wcdest/code/ >/dev/null 2>&1")
        os.system("make -C ../lib/wcdest/code/ >/dev/null 2>&1")
        os.system("make install -C ../lib/wcdest/code/ >/dev/null 2>&1")
        print("\nDone.\n")

        if self._out_directory == None: 
            path_prefix = "/".join(self._r1_telomeric_reads_filepath.split("/")[0:-2]) + "/clusters/cluster_info/telomeric_reads/"
        else:
            path_prefix = self._out_directory + "/clusters/cluster_info/telomeric_reads/"
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        self._r1_clustered_telomeric_reads_info_filepath = path_prefix + self._r1_telomeric_reads_filepath.split("/")[-1].split(".")[0] + ".ans"
        self._r2_clustered_telomeric_reads_info_filepath = path_prefix + self._r2_telomeric_reads_filepath.split("/")[-1].split(".")[0] + ".ans"

        # cluster forward telomeric reads and write results to file 
        print("Clustering telomeric reads in '" + self._r1_telomeric_reads_filepath + "'...")
        os.system("./wcdest/code/src/wcd -o " + self._r1_clustered_telomeric_reads_info_filepath + " -c " + self._r1_telomeric_reads_filepath + " >/dev/null 2>&1")
        self._r1_clustered_telomeric_reads_filepath = self._get_clustered_tels(self._r1_telomeric_reads_filepath, self._r1_clustered_telomeric_reads_info_filepath, "forward", "tel")
        print("Clustering results saved to '" + self._r1_clustered_telomeric_reads_filepath + "'")
        print("Done.\n")
        
        # cluster reverse telomeric reads and write results to file 
        print("Clustering telomeric reads in '" + self._r2_telomeric_reads_filepath + "'...")
        os.system("./wcdest/code/src/wcd -o " + self._r2_clustered_telomeric_reads_info_filepath + " -c " + self._r2_telomeric_reads_filepath + " >/dev/null 2>&1")
        self._r2_clustered_telomeric_reads_filepath = self._get_clustered_tels(self._r2_telomeric_reads_filepath, self._r2_clustered_telomeric_reads_info_filepath, "forward", "tel")
        print("Clustering results saved to '" + self._r2_clustered_telomeric_reads_filepath + "'")
        print("Done.\n")
    
    def extract_subtelomeric_reads(self):
        # given forward telomeric clusters, retrieve the subtelomeric sequence reads
        #  from the reverse read file  
        print("Extracting R1 subtelomeric reads from R2 read file...")
        for tel_file in os.listdir(self._r1_clustered_telomeric_reads_filepath):
            if tel_file.endswith(".fasta"): 
                self._r1_subtelomeric_reads_filepath = self._extract_paired_ends(self._r2_telomeric_reads_filepath, self._r1_clustered_telomeric_reads_filepath + tel_file, "reverse")
        print("Done.\n")

        # given reverse telomeric clusters, retrieve the subtelomeric sequence reads
        #  from the forward read file  
        print("Extracting R2 subtelomeric reads from R1 read file...")
        for tel_file in os.listdir(self._r2_clustered_telomeric_reads_filepath):
            if tel_file.endswith(".fasta"): 
                self._r2_subtelomeric_reads_filepath = self._extract_paired_ends(self._r1_telomeric_reads_filepath, self._r2_clustered_telomeric_reads_filepath + tel_file, "reverse")
        print("Done.\n")

        print("R1 subtelomeric reads saved to '" + self._r1_subtelomeric_reads_filepath + "'")
        print("R2 subtelomeric reads saved to '" + self._r2_subtelomeric_reads_filepath + "'")
        print("\n")
    
    def cluster_subtelomeric_reads(self):
        print("Clustering subtelomeric reads in " + self._r1_subtelomeric_reads_filepath + "...")
        # prepare directories for writing results of clustering 
        if self._out_directory == None: 
            path_prefix = "/".join(self._r1_telomeric_reads_filepath.split("/")[0:-2]) + "/clusters/cluster_info/subtelomeric_reads/r1/"
        else:
            path_prefix = self._out_directory + "/clusters/cluster_info/subtelomeric_reads/r1/"
    
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        # cluster every forward subtelomeric read file obtained from a corresponding 
        # reverse telomeric read cluster
        for subtel_file in os.listdir(self._r1_subtelomeric_reads_filepath):
            if subtel_file.endswith(".fasta"): 
                self._r1_clustered_subtelomeric_reads_info_filepath = path_prefix + subtel_file.split(".")[0] + ".ans"
                os.system("./wcdest/code/src/wcd -o " + self._r1_clustered_subtelomeric_reads_info_filepath + " -c " + self._r1_subtelomeric_reads_filepath + subtel_file + " >/dev/null 2>&1")
                self._r1_clustered_subtelomeric_reads_filepath = self._get_clustered_tels(self._r1_telomeric_reads_filepath, self._r1_clustered_subtelomeric_reads_info_filepath, "forward", "subtel")
        print("Clustering results saved to '" + self._r1_clustered_subtelomeric_reads_filepath + "'")
        print("Done.\n")

        # cluster every reverse subtelomeric read file obtained from a corresponding 
        # forward telomeric read cluster
        print("Clustering subtelomeric reads in " + self._r2_subtelomeric_reads_filepath + "...")
        # prepare directories for writing results of clustering 
        path_prefix = "/".join(self._r2_telomeric_reads_filepath.split("/")[0:-2]) + "/clusters/cluster_info/subtelomeric_reads/r2/"
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        # cluster every forward subtelomeric read file obtained from a corresponding 
        # reverse telomeric read cluster
        for subtel_file in os.listdir(self._r2_subtelomeric_reads_filepath):
            if subtel_file.endswith(".fasta"): 
                self._r2_clustered_subtelomeric_reads_info_filepath = path_prefix + subtel_file.split(".")[0] + ".ans"
                os.system("./wcdest/code/src/wcd -o " + self._r2_clustered_subtelomeric_reads_info_filepath + " -c " + self._r2_subtelomeric_reads_filepath + subtel_file + " >/dev/null 2>&1")
                self._r2_clustered_subtelomeric_reads_filepath = self._get_clustered_tels(self._r2_telomeric_reads_filepath, self._r2_clustered_subtelomeric_reads_info_filepath, "forward", "subtel")
        print("Clustering results saved to '" + self._r2_clustered_subtelomeric_reads_filepath + "'")        
        print("Done.\n")
    
    def _process_reads(self, input_filepath, r1_output, r2_output):
        # define regular expressions for matching telomeric sequences
        start_tel_regex = "^(CCCTAACCCTAACCCTAA)+|^A(CCCTAACCCTAACCCTAA)+|^AA(CCCTAACCCTAACCCTAA)+|^TAA(CCCTAACCCTAACCCTAA)+|^CTAA(CCCTAACCCTAACCCTAA)+|^CCTAA(CCCTAACCCTAACCCTAA)+"
        end_tel_regex = "(TTAGGGTTAGGGTTAGGG)+$|G(TTAGGGTTAGGGTTAGGG)+$|GG(TTAGGGTTAGGGTTAGGG)+$|GGG(TTAGGGTTAGGGTTAGGG)+$|AGGG(TTAGGGTTAGGGTTAGGG)+$|TAGGG(TTAGGGTTAGGGTTAGGG)+$"

        line_counter = 0
        header_line = ""
        with open(input_filepath) as sequence_reads:
            for line in sequence_reads:
                line_counter += 1
                if line_counter % FQ_LINES_PER_READ == HEADER:
                    header_line = line
                    # replace for output to FASTA format 
                    header_line = header_line.replace("@", ">")
                if line_counter % FQ_LINES_PER_READ == SEQUENCE:
                    start_tel_exists = re.search(start_tel_regex, line)
                    end_tel_exists = re.search(end_tel_regex, line)
                    if bool(start_tel_exists) ^ bool(end_tel_exists):
                        if int(header_line.split(" ")[1].split(":")[0]) == FORWARD:
                            r1_output.write(header_line)
                            r1_output.write(line)
                        if int(header_line.split(" ")[1].split(":")[0]) == REVERSE:
                            r2_output.write(header_line)
                            r2_output.write(line)
    
    # retrieve telomeric reads based on clustering by wcdest
    def _get_clustered_tels(self, tel_path, tel_clusters_info_path, read, type):
        tel_clusters = list()
        rows = list()
        row_num = 0 # keep track of row number of each cluster
        with open(tel_clusters_info_path) as tel_clusters_input:
            # for each line in wcdest cluster file, convert to list of ints and 
            # remove newline and ending "." Adjust indexing.
            for line in tel_clusters_input:
                line = [int(i) for i in line.strip(".\n").split(" ")]
                line = [i * FA_LINES_PER_READ for i in line]
                # filter out singleton clusters
                if len(line) > 1:
                    tel_clusters.append(line)
                    rows.append(row_num)
                row_num += 1
        row_idx = 0
        # get all R1 or R2 telomeric reads 
        tel_input = open(tel_path, "r")
        tel_lines = tel_input.readlines()
        for cluster in tel_clusters:
            # specify where to save FASTA files based on read direction and the 
            # type reads that are being clustered 
            if read == "forward" and type == "tel":
                path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r1_telomeric_reads/"
                os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
                out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
                clustered_tel_output = open(out_file_path, "w+")
            if read == "reverse" and type == "tel":
                path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r2_telomeric_reads/"
                os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
                out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
                clustered_tel_output = open(out_file_path, "w+")
            if read == "forward" and type == "subtel":
                path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r1_subtelomeric_reads/"
                os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
                out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
                clustered_tel_output = open(out_file_path, "w+")
            if read == "reverse" and type == "subtel":
                path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r2_subtelomeric_reads/"
                os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
                out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
                clustered_tel_output = open(out_file_path, "w+")
            # given the list of clusters with singletons filtered out, append
            # cluster reads to file 
            for i in range(len(cluster)):
                index = cluster[i]
                header_line = tel_lines[index]
                sequence_line = tel_lines[index + 1]
                # edit header to match the read direction of fastq file 
                if read == "reverse":
                    header_line = header_line.replace(" 1:N", " 2:N")
                if read == "forward":
                    header_line = header_line.replace(" 2:N", " 1:N")
                clustered_tel_output.write(header_line)
                clustered_tel_output.write(sequence_line)
            row_idx += 1
        return path_prefix

    # retrieve paired ends of a set of reads for a file with telomeric reads from
    # a corresponding FASTQ file containing the paired end
    def _extract_paired_ends(self, fastq_path, tel_path, read):
        tel_headers = list()
        line_counter = 1
        header_line = ""
        # store all headers from telomeric read file
        with open(tel_path) as ref_tel_input:
            for line in ref_tel_input:
                if line_counter % FA_LINES_PER_READ == HEADER:
                    # edit header to match the read direction of fastq file 
                    if read == "reverse":
                        line = line.replace(" 1:N", " 2:N")
                    if read == "forward":
                        line = line.replace(" 2:N", " 1:N")
                    # edit header to match FASTQ format
                    tel_headers.append(line.replace(">", "@"))
                line_counter += 1
        # specify where to save FASTA files based on the read type
        if read == "forward":
            path_prefix = "/".join(tel_path.split("/")[0:-3]) + "/paired_ends/r1_subtelomeric_reads/"
            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            out_file_path = path_prefix + "/" + tel_path.split("/")[-1].split(".")[0] + ".fasta"
            out_file_path = out_file_path.replace("telomeric", "subtelomeric")
            fastq_output = open(out_file_path, "w+")
        if read == "reverse":
            path_prefix = "/".join(tel_path.split("/")[0:-3]) + "/paired_ends/r2_subtelomeric_reads/"
            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            out_file_path = path_prefix + "/" + tel_path.split("/")[-1].split(".")[0] + ".fasta"
            out_file_path = out_file_path.replace("telomeric", "subtelomeric")
            fastq_output = open(out_file_path, "w+")
        line_counter = 1
        header_line = ""
        # write all paired ends to file 
        with open(fastq_path) as fastq_input:
            for line in fastq_input:
                if line_counter % FQ_LINES_PER_READ == HEADER:
                    header_line = line
                if line_counter % FQ_LINES_PER_READ == SEQUENCE:
                    # write sequence from FASTQ to file if its header is specified 
                    # in the list tel_headers constructed earlier 
                    if header_line in tel_headers:
                        # edit header to match FASTA format
                        header_line = header_line.replace("@", ">")
                        fastq_output.write(header_line)
                        fastq_output.write(line)
                line_counter += 1
        return path_prefix