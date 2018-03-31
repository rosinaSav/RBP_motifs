'''
Author: Rosina Savisaar.
Calculate the combined density of a set of motif sets.
'''

from bedtools_games import Feature_Set
from housekeeping import flatten, list_to_dict, make_dir, parse_arguments
import nucleotide_comp as nc
import numpy as np
import read_and_write as rw

def main():

    description = "Calculate the combined density of a set of motif sets."
    args = parse_arguments(description, ["motifs_file_name", "summary_file_name", "dataset_name", "output_folder_name", "output_file_name", "n_sim", "features_file_name", "genome", "families_file_name", "fasta_name", "ND_column", "seed", "output_suffix", "negative_ND", "new_filters", "upper_quarter", "lower_quarter", "full_set", "newer_filters", "two_seqs"], ints = [5, 10, 11], flags = [13, 14, 15, 16, 17, 18, 19])
    [motifs_file_name, summary_file_name, dataset_name, output_folder_name, output_file_name, n_sim, features_file_name, genome, families_file_name, fasta_name, ND_column, seed, output_suffix, negative_ND, new_filters, upper_quarter, lower_quarter, full_set, newer_filters, two_seqs] = [args.motifs_file_name, args.summary_file_name, args.dataset_name, args.output_folder_name, args.output_file_name, args.n_sim, args.features_file_name, args.genome, args.families_file_name, args.fasta_name, args.ND_column, args.seed, args.output_suffix, args.negative_ND, args.new_filters, args.upper_quarter, args.lower_quarter, args.full_set, args.newer_filters, args.two_seqs]

    #make a dictionary with RBPs as keys and ND/p values as values.
    if summary_file_name != "None":
        summary_data = rw.read_many_fields(summary_file_name, "\t")
        #because some of the files are tab-separated, while others are comma-separated and have a header row
        if len(summary_data[0]) == 1:
            summary_data = rw.read_many_fields(summary_file_name, ",")
            summary_data = summary_data[1:]

        summary_dict = list_to_dict(summary_data, 0, ND_column, floatify = True)

    #make a dictionary with RBPs as keys and lists of associated motifs as values        
    motifs = rw.read_motifs(motifs_file_name)

    #if you only want to be using a subset of the motifs
    if not full_set:
        #motifs with negative ND
        if negative_ND:
            motifs = [motifs[RBP] for RBP in motifs if summary_dict[RBP] < 0]
        #the most significantly enriched motifs
        elif upper_quarter:
            motifs = [motifs[RBP] for RBP in motifs if summary_dict[RBP] < 0.1]
        #the most significantly depleted motifs
        elif lower_quarter:
            motifs = [motifs[RBP] for RBP in motifs if summary_dict[RBP] > 0.9]
        #motifs with positive ND
        else:
            motifs = [motifs[RBP] for RBP in motifs if summary_dict[RBP] >= 0]

    #shove all the remaining motifs into a great big flattened and uniquified bag
    motifs = list(set(flatten(list(motifs.values()))))

    print(len(motifs))
    make_dir(output_folder_name)

    #if you want to average over families
    if features_file_name != "None":
        fs = Feature_Set(features_file_name, genome)
        fs.set_dataset(dataset_name)
        families = rw.read_families(families_file_name)
        fs.add_families(families)
    else:
        fs = None

    #generate 100 1000 bp long random sequences based on the hg38 mononucleotide composition and use that as your sequence fasta
    if fasta_name == "random":
        names = [i for i in range(100)]
        seqs = nc.kmers_from_nc(1000, 100, genome_comp = True)
        fasta_name = "RBP/random_sequences_from_genome_comp.fasta"
        rw.write_to_fasta(names, seqs, fasta_name)

    with open(output_file_name, "w") as output_file:
        #generate n_sim sets of simulant motifs (constraining the space of simulants based on different sets of filters)
        if new_filters:
            simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True, seed = seed)
        elif newer_filters:
            simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True, seed = seed, concat = False, no_duplicates = True)
        else:
            current_simulants = nc.make_simulants(motifs, n_sim, seed = seed)
        #calculate the density parameters of the motifs in the sequence fasta
        output_dict = nc.get_sequence_set_density(fasta_name, None, motifs, simulants, n_sim,
                                                   "{0}/overall_density_{1}.csv".format(output_folder_name, output_suffix),
                                                   "{0}/overall_sim_density_{1}.csv".format(output_folder_name, output_suffix),
                                                   "{0}/overall_positions.csv_{1}".format(output_folder_name, output_suffix),
                                                   "{0}/overall_sim_positions_{1}".format(output_folder_name, output_suffix),
                                                   concat = False, positions = False, feature_set = fs, verbose = True, two_seqs = two_seqs)
        record = [str(output_dict["median density"]), str(np.mean(output_dict["simulated densities"])), str(output_dict["median ND"]), str(output_dict["effective p"]), str(output_dict["Z"]), str(output_dict["depletion p"]), str(len(motifs)), str(output_dict["simulant sd"])]
        #write to output file
        output_file.write("\t".join(record))
        print(record)
    

if __name__ == "__main__":
    main()
    
