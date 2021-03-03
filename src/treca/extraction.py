import os
import re
import subprocess
from tqdm import tqdm

# define constants
HEADER = 1
SEQUENCE = 2
NUM_LINES = 92475176
FA_LINES_PER_READ = 2
FQ_LINES_PER_READ = 4

# define lookup dictionary for read direction
# define regular expressions for matching telomeric sequences
read_type = {"forward":1, "reverse":2}
start_tel_regex = "^(CCCTAACCCTAACCCTAA)+|^A(CCCTAACCCTAACCCTAA)+|^AA(CCCTAACCCTAACCCTAA)+|^TAA(CCCTAACCCTAACCCTAA)+|^CTAA(CCCTAACCCTAACCCTAA)+|^CCTAA(CCCTAACCCTAACCCTAA)+"
end_tel_regex = "(TTAGGGTTAGGGTTAGGG)+$|G(TTAGGGTTAGGGTTAGGG)+$|GG(TTAGGGTTAGGGTTAGGG)+$|GGG(TTAGGGTTAGGGTTAGGG)+$|AGGG(TTAGGGTTAGGGTTAGGG)+$|TAGGG(TTAGGGTTAGGGTTAGGG)+$"

# find the number of lines in a file via system call to wc
def get_num_lines(file_path):
    out = subprocess.Popen(['wc', '-l', file_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = out.communicate()
    return int(stdout.split()[0])

# add a file to the end of another file and delete the appended file
def append_file(filepath1, filepath2):
    output_file = open(filepath1, "a")
    with open(filepath2) as input_file:
        for line in input_file:
                output_file.write(line)
    os.remove(filepath2)

# get the telomeric reads from FASTQ file
def extract_telomeric_reads(fastq_path, read):
    # prepare output for forward telomeric reads
    if read == "forward":
        path_prefix = "/".join(fastq_path.split("/")[0:-1]) + "/tel_reads/"
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        out_file_path = path_prefix + fastq_path.split("/")[-1].split(".")[0] + "_telomeric_reads.fasta"
        fastq_output = open(out_file_path, "w+")
    # prepare output for reverse telomeric reads
    if read == "reverse":
        path_prefix = "/".join(fastq_path.split("/")[0:-1]) + "/tel_reads/"
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        out_file_path = path_prefix + fastq_path.split("/")[-1].split(".")[0] + "_telomeric_reads.fasta"
        fastq_output = open(out_file_path, "w+")
    line_counter = 1
    header_line = ""
    # get number of lines for tqdm progress bar
    fastq_file_length = get_num_lines(fastq_path)  
    # for each 4 lines in the FASTQ file, attempt to match the sequence line to
    # the regular expressions defined. If a match occurs, write both the 
    # sequence and header to file
    with open(fastq_path) as fastq_input:
        for line in tqdm(fastq_input, total=fastq_file_length, desc="Progress"):
            if line_counter % FQ_LINES_PER_READ == HEADER:
                header_line = line
                # replace for output to FASTA format 
                header_line = header_line.replace("@", ">")
            if line_counter % FQ_LINES_PER_READ == SEQUENCE:
                start_tel_exists = re.search(start_tel_regex, line)
                end_tel_exists = re.search(end_tel_regex, line)
                if bool(start_tel_exists) ^ bool(end_tel_exists):
                    fastq_output.write(header_line)
                    fastq_output.write(line)
            line_counter += 1
    # return path to file containing sequence reads
    return out_file_path, fastq_file_length

# retrieve telomeric reads based on clustering by wcdest
def get_clustered_tels(tel_path, tel_clusters_info_path, read, type):
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
            path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r1_tel/"
            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
            clustered_tel_output = open(out_file_path, "w+")
        if read == "reverse" and type == "tel":
            path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r2_tel/"
            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
            clustered_tel_output = open(out_file_path, "w+")
        if read == "forward" and type == "subtel":
            path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r1_subtel/"
            os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
            out_file_path = path_prefix + tel_path.split("/")[-1].split(".")[0] + "_cluster"+ str(rows[row_idx]) +".fasta"
            clustered_tel_output = open(out_file_path, "w+")
        if read == "reverse" and type == "subtel":
            path_prefix = "/".join(tel_path.split("/")[0:-2]) + "/clusters/r2_subtel/"
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
def extract_paired_ends(fastq_path, fastq_len, tel_path, read):
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
        path_prefix = "/".join(tel_path.split("/")[0:-3]) + "/paired_ends/r1_subtel/"
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        out_file_path = path_prefix + "/" + tel_path.split("/")[-1].split(".")[0] + ".fasta"
        out_file_path = out_file_path.replace("telomeric", "subtelomeric")
        fastq_output = open(out_file_path, "w+")
    if read == "reverse":
        path_prefix = "/".join(tel_path.split("/")[0:-3]) + "/paired_ends/r2_subtel/"
        os.makedirs(os.path.dirname(path_prefix), exist_ok=True)
        out_file_path = path_prefix + "/" + tel_path.split("/")[-1].split(".")[0] + ".fasta"
        out_file_path = out_file_path.replace("telomeric", "subtelomeric")
        fastq_output = open(out_file_path, "w+")
    line_counter = 1
    header_line = ""
    # write all paired ends to file 
    with open(fastq_path) as fastq_input:
        for line in tqdm(fastq_input, total=fastq_len, desc="Progress"):
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