import csv
from housekeeping import flatten, print_elements
import os
import re
import sys

def clean_fasta(fasta):
    '''
    Filter a fasta file to only retain sequences that don't include N bases.
    '''
    names, seq = read_fasta(fasta)
    seq = [i.upper() for i in seq]
    names = [names[i] for i in range(len(names)) if "N" not in seq[i]]
    seq = [i for i in seq if "N" not in i]
    write_to_fasta(names, seq, fasta) 

def read_families(file):
    '''
    Read a families file (one family of paralogous genes per line, the member genes separated by commas) into a list,
    with each sublist containing the identifiers of the genes belonging to one family.
    '''
    families = []
    with open(file) as families_file:
        for line in families_file:
            current_family = line.rstrip("\n")
            current_family = current_family.split(",")
            current_family = [i for i in current_family if i != ""]
            families.append(current_family)
    return(families)

def read_fasta(input_file):
    '''
    Given a fasta file return a first lists containing the sequence identifiers and a second list containing teh sequences (in the same order).
    '''
    file_to_read = open(input_file)
    input_lines = file_to_read.readlines()
    file_to_read.close()
    input_lines = [i.rstrip("\n") for i in input_lines]
    names = [i.lstrip(">") for i in input_lines if i[0] == ">"]
    sequences = [i for i in input_lines if i[0] != ">"]
    if len(sequences) != len(names):
        print("Problem extracting data from fasta file!")
        print(len(sequences))
        print(len(names))
        raise Exception
    if len(sequences) == 0:
        print("No sequences were extracted!")
        raise Exception
    return(names, sequences)

def read_genbank(input_file):
    '''
    Parse a GenBank file into a dictionary.
    '''
    output = {}
    with open(input_file) as file:
        data = "".join(file)
    version = re.search("(VERSION[ ]*)([\w\d\._]*)(  GI)", data).group(2)
    output["identifier"] = version
    definition = re.search("(DEFINITION[ ]*)([\w\d\ \(\),\.\-\n]*)(\nACCESSION)", data).group(2)
    definition = re.sub("\n", " ", definition)
    definition = re.sub(" {2,}", " ", definition)
    output["definition"] = definition
    translation = re.search("(\/translation\=\")([\w\n ]*)(\")", data).group(2)
    translation = re.sub("\n", "", translation)
    translation = re.sub(" ", "", translation)
    output["translation"] = translation
    length = int(re.search("(LOCUS[ ]*[\w\d\._]*)( *)(\d*)", data).group(3))
    output["mRNA length"] = length
    exon_lines = re.findall("[ ]*exon[ ]*\d*\.\.\d*", data)
    exon_lines = [re.findall("\d*", i) for i in exon_lines]
    exon_lines = [[int(j) for j in i if j != ""] for i in exon_lines]
    CDS_coords = re.search("([ ]*CDS[ ]*)(\d*)(\.\.)(\d*)", data)
    CDS_coords = [int(CDS_coords.group(2)), int(CDS_coords.group(4))]
    CDS_pos = []
    for line in exon_lines:
        if (line[1] < CDS_coords[0]) or (line[0] > CDS_coords[1]):
            pass
        elif line[0] < CDS_coords[0]:
            CDS_pos.append([CDS_coords[0], line[1]])
        elif line[1] > CDS_coords[1]:
            CDS_pos.append([line[0], CDS_coords[1]])
        else:
            CDS_pos.append(line)
    sequence = re.search("(ORIGIN[ \n1]*)([\w \n]*)(\/\/)", data).group(2)
    sequence = re.findall("[a-z]*", sequence)
    sequence = "".join(sequence).upper()
    output["mRNA sequence"] = sequence
    if len(sequence) != length:
        print("mRNA sequence not extracted correctly!")
        print("Observed length: {0}; expected length: {1}".format(len(sequence), length))
        raise Exception
    CDS_sequence = [sequence[(i[0] - 1):i[1]] for i in CDS_pos]
    CDS_sequence = "".join(CDS_sequence)
    output["CDS sequence"] = CDS_sequence
    if len(CDS_sequence) != (len(translation) * 3) + 3:
        print("CDS sequence not extracted correctly!")
        raise Exception
    return(output)

def read_many_fields(input_file, delimiter):
    '''
    Read a csv/tsv/... into a lists of lists with each sublist correpsonding to one line.
    '''
    file_to_read = open(input_file)
    field_reader = csv.reader(file_to_read, delimiter = delimiter)
    lines = []
    for i in field_reader:
        lines.append(i)
    file_to_read.close()
    return(lines)

def read_motifs(input_file):
    '''
    Given a motifs file (a fasta with RBPs as names and pipe-separated motifs as sequences),
    read it into a dictionary with the RBPs as keys and lists of motifs as values.
    '''
    names, motifs = read_fasta(input_file)
    motif_dict = {}
    for i in names:
        motif_dict[i] = motifs[names.index(i)].split("|")
    return(motif_dict)

def read_names(input_file):
    '''
    Read in a list of strings.
    '''
    file_to_read = open(input_file)
    names = file_to_read.readlines()
    names = [i.rstrip("\n") for i in names]
    file_to_read.close()
    return(names)

def read_numbers(input_file):
    '''
    Read in a list of floats.
    '''
    file_to_read = open(input_file)
    numbers = file_to_read.readlines()
    numbers = [float(i.rstrip("\n")) for i in numbers]
    file_to_read.close()
    return(numbers)

def write_all(input_string, file_name):
    '''
    Write a string to file.
    '''
    output_file = open(file_name, "w")
    output_file.write(input_string)
    output_file.close()

def write_names(names, file_name):
    '''
    Write a list of strings to file.
    '''
    with open(file_name, "w") as file:
        #this is so you wouldn't get an empty line at the end
        for name in names[:-1]:
            file.write(name)
            file.write("\n")
        file.write(names[-1])

def write_to_csv(input_list, file_name, delimiter, flat = False, verbose = False):
    '''
    Write a list of lists to a csv/tsv/...
    '''
    output_file = open(file_name, "w")
    output_writer = csv.writer(output_file, delimiter = delimiter, lineterminator = os.linesep)
    counter = 0
    if flat:
        for i in input_list:
            output_writer.writerow([i])
            counter = counter + 1
    else:
        for i in input_list:
            output_writer.writerow(i)
            counter = counter + 1
    output_file.close()
    if verbose:
        print("Wrote {0} records to file.".format(counter))

def write_to_fasta(names, seq, fasta_name):
    '''
    Write a set of sequence identifiers and sequences to fasta file.
    '''
    with open(fasta_name, "w") as file:
        for i in range(len(names)):
            file.write(">{0}\n".format(names[i]))
            file.write("{0}\n".format(seq[i]))
