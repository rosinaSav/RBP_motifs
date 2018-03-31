'''
Author: Rosina Savisaar.
Calculate the combined density of a set of motif sets.
'''

from bedtools_games import Feature_Set
import conservation
from housekeeping import flatten, list_to_dict, make_dir, parse_arguments
import nucleotide_comp as nc
import numpy as np
import os
import random
import read_and_write as rw

def main():

    description = "Calculate the combined density of a set of motif sets."
    args = parse_arguments(description, ["motifs_file_name", "summary_file_name", "dataset_name", "correspondances_file_name", "alignment_folder_name", "output_folder_name", "output_file_name", "n_sim", "features_file_name", "genome", "families_file_name", "fasta_name", "ND_column", "output_suffix", "validity_folder_name", "negative_ND", "new_filters", "upper_quarter", "lower_quarter", "full_set", "gene_families", "newer_filters", "baseml"], ints = [7, 12], flags = [15, 16, 17, 18, 19, 20, 21, 22])
    [motifs_file_name, summary_file_name, dataset_name,  correspondances_file_name, alignment_folder_name, output_folder_name, output_file_name, n_sim, features_file_name, genome, families_file_name, fasta_name, ND_column, output_suffix, validity_folder_name, negative_ND, new_filters, upper_quarter, lower_quarter, full_set, gene_families, newer_filters, baseml] = [args.motifs_file_name, args.summary_file_name, args.dataset_name,  args.correspondances_file_name, args.alignment_folder_name, args.output_folder_name, args.output_file_name, args.n_sim, args.features_file_name, args.genome, args.families_file_name, args.fasta_name, args.ND_column, args.output_suffix, args.validity_folder_name, args.negative_ND, args.new_filters, args.upper_quarter, args.lower_quarter, args.full_set, args.gene_families, args.newer_filters, args.baseml]

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
        #which RBPs fulfill the necessary information content criteria?
        validity = rw.read_many_fields("{0}/sufficient_information_fraction05.csv".format(validity_folder_name), "\t")
        validity = list_to_dict(validity, 0, 1)
        #motifs with negative ND
        if negative_ND:
            motifs = [motifs[RBP] for RBP in motifs if (summary_dict[RBP] < 0) and (validity[RBP] == "True")]
        #the most significantly enriched motifs
        elif upper_quarter:
            motifs = [motifs[RBP] for RBP in motifs if (summary_dict[RBP] < 0.1) and (validity[RBP] == "True")]
        #the most significantly depleted motifs
        elif lower_quarter:
            motifs = [motifs[RBP] for RBP in motifs if (summary_dict[RBP] > 0.9) and (validity[RBP] == "True")]
        #motifs with positive ND
        else:
            motifs = [motifs[RBP] for RBP in motifs if (summary_dict[RBP] >= 0) and (validity[RBP] == "True")]

    #shove all the remaining motifs into a great big flattened and uniquified bag
    motifs = list(set(flatten(list(motifs.values()))))

    make_dir(output_folder_name)

    #prepare a Feature_Set object (a genome gtf associated to a particular genome and to a set of transcript identifiers)
    if features_file_name != "None":
        fs = Feature_Set(features_file_name, genome)
        fs.set_dataset(dataset_name)
        transcripts = fs.get_transcripts()
        CDS = fs.get_CDS()
        #paralogous families
        families = rw.read_families(families_file_name)
        #the families file might use gene identifiers, whereas the Feature_Set object uses transcript identifiers
        if gene_families:
            families = fs.convert_families_to_ENST(families, transcripts)
        fs.add_families(families)
        #pick a random member from each paralogous family
        picked_trans = fs.pick_random_members()
        names = rw.read_fasta(fasta_name)[0]
        if picked_trans[0] not in names:
            picked = [fs.convert_between_ENST_and_ENSG(i, transcripts, "ENSG") for i in picked_trans]
        else:
            picked = picked_trans
        print(len(picked))
    else:
        picked = None

    if baseml:
        method = "baseml"
    else:
        method = "gy"

    #write the input data for the conservation analysis into a file
    input_dict_file_name = "temp_data/temp_{0}.txt".format(random.random())
    conservation.input_dict_for_dS(correspondances_file_name, alignment_folder_name, fasta_name, input_dict_file_name, picked = picked)
    with open(output_file_name, "w") as file:
        file.write(",".join(["real_dS", "mean_sim_dS", "norm_dS", "p", "motif_number"]))
        file.write("\n")
        #make n_sim simulant sets for the motifs, filtering the simulants based on different sets of criteria
        if new_filters:
            simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True, seed = 1)
        elif newer_filters:
            simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True, seed = 1, no_duplicates = True, concat = False)               
        else:
            simulants = nc.make_simulants(motifs, n_sim, seed = 100)
        #file where the simulants dS values will be stored
        sim_output_file_name = "{0}/{1}_sim_ds.csv".format(output_folder_name, output_suffix)
        #calculate dS within motifs and simulants
        output_dict = conservation.dS_from_hits(motifs, alignment_folder_name, input_dict_file_name, n_sim = n_sim, simulants = simulants, sim_output_file_name = sim_output_file_name, method = method)
        print(output_dict)
        print("\n")
        #write to output file
        if output_dict != None:
            file.write(",".join([str(output_dict["dS"]), str(output_dict["mean simulated dS"]), str(output_dict["normalized dS"]), str(output_dict["effective p"]), str(len(motifs))]))
        else:
            file.write(",".join([str(None), str(None), str(None), str(None), str(None)]))
    os.remove(input_dict_file_name)

if __name__ == "__main__":
    main()
    
