'''
Author: Rosina Savisaar.
Write the median motif lengths of a series of motif sets to file.
'''

from housekeeping import parse_arguments
import numpy as np
import read_and_write as rw

def main():

    description = "Write the median motif lengths of a series of motif sets to file."
    args = parse_arguments(description, ["input_file", "output_file"])
    [input_file, output_file] = [args.input_file, args.output_file]
    
    #parse motifs from FASTA
    names, motifs = rw.read_fasta(input_file)
    motifs = [i.split("|") for i in motifs]
    motif_lengths = [[len(j) for j in i] for i in motifs]
    #write down and print out motif lengths
    with open(output_file, "w") as file:
        for pos, lengths_list in enumerate(motif_lengths):
            file.write("{0}\t{1}\n".format(names[pos], np.median(lengths_list)))
            print(np.median(lengths_list))

if __name__ == "__main__":
    main()
