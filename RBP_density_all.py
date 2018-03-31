'''
Author: Rosina Savisaar.
Determine the density of RBP target motif occurrences in a set of CDSs.
Also generate an empirical distribution of densities for simulant motifs.
'''

import argparse
from bedtools_games import Feature_Set
from housekeeping import make_dir
import nucleotide_comp as nc
import numpy as np
import os
import read_and_write as rw
import sys

def main():

    parser = argparse.ArgumentParser(description="Calculate the density of a series of RBP motifs in any type of sequence.")
    parser.add_argument("RBP_file_name", type = str, help = "name of file with RBP motifs")
    parser.add_argument("output_folder_name", type = str, help = "name of folder that will contain analysis results")
    parser.add_argument("output_file_name", type = str, help = "name of file that will contain analysis results")
    parser.add_argument("input_file_name", type = str, help = "name of fasta file with the sequences")
    parser.add_argument("n_sim", type = int, help = "number of simulants")
    parser.add_argument("features_file_name", type = str, help = "name of GTF file")
    parser.add_argument("genome", type = str, help = "genome name")
    parser.add_argument("dataset_name", type = str, help = "dataset name")
    parser.add_argument("families_file_name", type = str, help = "families file name")
    parser.add_argument("--simulants_within", dest = "simulants_within", action = "store_true", help = "Should simulants be generated only from dinucleotides within each particular motif?")
    parser.add_argument("--sequence_control", dest = "sequence_control", action = "store_true", help = "Should shuffled sequences be used as control?")
    parser.add_argument("--remove_stops", dest = "remove_stops", action = "store_true", help = "Should simulant motifs not incldue motifs that contain stop codon sequences? (boolean)")
    parser.add_argument("--markov", dest = "markov", action = "store_true", help = "Should simulants be generated using a Markov model?")
    parser.add_argument("--new_filters", dest = "new_filters", action = "store_true", help = "Should simulants be generated using the old method but capping mononucleotide runs and removing existing motifs?")
    parser.add_argument("--no_concat", dest = "no_concat", action = "store_true", help = "Should a density be calculated for each gene?")
    parser.add_argument("--newer_filters", dest = "newer_filters", action = "store_true", help = "Like new_filters, but also not allowing duplicates in the simulants and without concatenation.")
    parser.add_argument("--two_seqs", dest = "two_seqs", action = "store_true", help = "Set to true if the sequence fasta has two sequences separated by a pipe in each line.")

    args = parser.parse_args()
    [RBP_file_name, output_folder_name, output_file_name, input_file_name, n_sim, features_file_name, genome, dataset_name, families_file_name, simulants_within, sequence_control, remove_stops, markov, new_filters, no_concat, newer_filters, two_seqs] = [args.RBP_file_name,
                                                                                                                                                                    args.output_folder_name,
                                                                                                                                                                    args.output_file_name,
                                                                                                                                                                    args.input_file_name,
                                                                                                                                                                    args.n_sim,
                                                                                                                                                                    args.features_file_name,
                                                                                                                                                                    args.genome,
                                                                                                                                                                    args.dataset_name,
                                                                                                                                                                    args.families_file_name,
                                                                                                                                                                                      args.simulants_within,
                                                                                                                                                                                    args.sequence_control,
                                                                                                                                                                                                     args.remove_stops,
                                                                                                                                                                                                             args.markov,
                                                                                                                                                                                                                          args.new_filters,
                                                                                                                                                                                                                                     args.no_concat,
                                                                                                                                                                                                                                                    args.newer_filters,
                                                                                                                                                                                                                                                              args.two_seqs]   
    make_dir(output_folder_name)

    #if you want to average over families
    if features_file_name != "None":
        fs = Feature_Set(features_file_name, genome)
        fs.set_dataset(dataset_name)
        families = rw.read_families(families_file_name)
        fs.add_families(families)
    else:
        fs = None

    #if concat, sum motif hit base counts across sequences and divide by the total sequence length,
    #otherwise produce a density estimate separately for each sequence and use the median as the final statistic
    if no_concat:
        concat = False
    else:
        concat = True

    #read in RBP motifs
    RBPs, motifs = rw.read_fasta(RBP_file_name)
    with open(output_file_name, "w") as output_file:
        for pos, RBP in enumerate(RBPs):
            curr_motifs = motifs[pos].split("|")
            #if, as control, you want to shuffle the codons within sequences
            if sequence_control:
                current_simulants = 3
                output_suffix = "_sequence_control"
            #if, as control, you want to calculate the density of simulant motifs
            else:
                #generate simulant motifs, applying different sets of filters onto the simulant motifs
                output_suffix = ""
                if simulants_within:
                    current_simulants = nc.make_simulants_within(curr_motifs, n_sim)
                elif markov:
                    current_simulants = nc.make_simulants_markov(curr_motifs, n_sim, remove_stops = remove_stops, remove_existing = True)
                elif new_filters:
                    current_simulants = nc.make_simulants(curr_motifs, n_sim, remove_stops = remove_stops, remove_existing = True, cap_runs = True)
                elif newer_filters:
                    current_simulants = nc.make_simulants(curr_motifs, n_sim, remove_stops = remove_stops, remove_existing = True, cap_runs = True, no_duplicates = True, concat = False)                   
                else:
                    current_simulants = nc.make_simulants(curr_motifs, n_sim, remove_stops = remove_stops)
            #get raw density, normalized density, p, Z... for current RBP
            current_dict = nc.get_sequence_set_density(input_file_name, None, curr_motifs, current_simulants, n_sim,
                                                       "{0}/{1}_{2}_density.csv".format(output_folder_name, RBP, output_suffix),
                                                       "{0}/{1}_{2}_sim_density.csv".format(output_folder_name, RBP, output_suffix),
                                                       "{0}/{1}_{2}_positions.csv".format(output_folder_name, RBP, output_suffix),
                                                       "{0}/{1}_{2}_sim_positions".format(output_folder_name, RBP, output_suffix),
                                                       concat = concat, positions = False, feature_set = fs, two_seqs = two_seqs)
            if concat:
                current_record = [RBP, str(current_dict["density"]), str(np.mean(current_dict["simulated densities"])), str(current_dict["ND"]), str(current_dict["effective p"]), str(current_dict["Z"]), str(current_dict["depletion p"]), str(len(curr_motifs)), str(current_dict["simulant sd"])]
            else:
                current_record = [RBP, str(current_dict["median density"]), str(np.mean(current_dict["simulated densities"])), str(current_dict["median ND"]), str(current_dict["effective p"]), str(current_dict["Z"]), str(current_dict["depletion p"]), str(len(curr_motifs)), str(current_dict["simulant sd"])]
            output_file.write("\t".join(current_record))
            output_file.write("\n")
            print(current_record)


if __name__ == "__main__":
    main()
