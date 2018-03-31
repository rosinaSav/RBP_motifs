'''
Author: Rosina Savisaar.
Replace true RBP target motifs with a simulant set
that preserves certain properties of the original set (depending on options).
'''

from housekeeping import flatten
import nucleotide_comp as nc
import read_and_write as rw
import sys

def main():
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    mode = sys.argv[3]
    motifs = rw.read_motifs(input_file)
    comp_dict = {"A": 0.25, "T": 0.25, "C": 0.25, "G": 0.25}
    new_motifs = {}
    for i in motifs:
        if mode == "structure":
            #systematically replace each instance of a particular base with another base
            #(say, all the As with Cs, all the Gs with Ts...)
            new_motifs[i] = nc.reassign_bases(motifs[i].copy())
        elif mode == "random":
            #generate random motifs, with all bases equiprobable
            #if all of the motifs always had the same lengths, you could just generate them in one go.
            #but here you have to go one by one, which is why n is 1.
            new_motifs[i] = [nc.kmers_from_nc(len(j), 1, comp_dict) for j in motifs[i]]
            #kmers from nc retruns a list of motifs. In this case, there is only one motif.
            new_motifs[i] = flatten(new_motifs[i])
        elif mode == "genome_comp":
            #generate random motifs, sampling from the hg38 mononucleotide composition
            new_motifs[i] = [nc.kmers_from_nc(len(j), 1, genome_comp = True) for j in motifs[i]]
            new_motifs[i] = flatten(new_motifs[i])
            print(new_motifs[i])
        else:
            print("Invalid mode!")
            print(mode)
            sys.exit()
        print(motifs[i])
        print(new_motifs[i])
        print("\n")
    names = [i for i in sorted(list(new_motifs.keys()))]
    seq = [new_motifs[i] for i in sorted(list(new_motifs.keys()))]
    seq = ["|".join(i) for i in seq]
    rw.write_to_fasta(names, seq, output_file)
                

if __name__ == "__main__":
    main()
