import os

def validate_file(fastq_path, extension, read):
    read_type = {"interleaved":0, "forward":1, "reverse":2}
    # checks if the file exists
    if not os.path.exists(fastq_path):
        print("\nThe specified file does not exist:", fastq_path)
        print("\nExiting...\n")
        exit(0)
    # checks if the file is in the correct format
    if not ((extension == "fastq") or (extension == "fq")):
        print("\nThe specified file is in the incorrect format:", fastq_path.split("/")[-1])
        print("\nSupported file formats for sequence reads: [*.fastq, *.fq]")
        print("\nExiting...\n")
        exit(0)
    # checks if the file are specified in the correct order
    f = open(fastq_path, "r+")
    line = f.readline()
    if not read_type[read] == int(line.split(" ")[1].split(":")[0]):
        if read == "forward":
            print("\nExpected a R1 (forward) read file but got a R2 (reverse) read file:", fastq_path.split("/")[-1])
        if read == "reverse":
            print("\nExpected a R2 (reverse) read file but got a R1 (forward) read file:", fastq_path.split("/")[-1])
        print("\nExiting...\n")
        exit(0)