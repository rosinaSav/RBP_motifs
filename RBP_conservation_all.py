import argparse
from bedtools_games import Feature_Set
import conservation
from housekeeping import flatten, list_to_dict, print_elements
import nucleotide_comp as nc
import os
import random
import read_and_write as rw

def main():
    parser = argparse.ArgumentParser(description="Calculate the conservation level of a series of RBP motifs.")
    parser.add_argument("features_file_name", type = str, help = "name of GTF file with genome features")
    parser.add_argument("dataset_name", type = str, help = "dataset name")
    parser.add_argument("genome", type = str, help = "genome assembly name")
    parser.add_argument("RBP_file_name", type = str, help = "name of file with RBP motifs")
    parser.add_argument("correspondances_file_name", type = str, help = "name of file with correspondances between genes in dataset and orthologs")
    parser.add_argument("fasta_file_name", type = str, help = "name of fasta file with the sequences")
    parser.add_argument("families_file_name", type = str, help = "name of file that contains families")
    parser.add_argument("output_file_name", type = str, help = "file for output data")
    parser.add_argument("output_folder_name", type = str, help = "folder that will contain simulated dS scores")
    parser.add_argument("alignment_folder_name", type = str, help = "name of folder that contains alignments")
    parser.add_argument("n_sim", type = int, help = "number of simulants")
    parser.add_argument("--valid_file", nargs = "?", const = "False")
    parser.add_argument("--gene_families", action = "store_true", help = "does the families file use gene identifiers?")
    parser.add_argument("--markov", dest = "markov", action = "store_true", help = "Should simulants be generated using a Markov model?")
    parser.add_argument("--new_filters", dest = "new_filters", action = "store_true", help = "Should simulants be generated using the old method but capping mononucleotide runs and removing existing motifs?")
    parser.add_argument("--newer_filters", dest = "newer_filters", action = "store_true", help = "Like new_filters but without concatenation and without allowing duplicates within simulant sets.")
    parser.add_argument("--goldman_yang", dest = "goldman_yang", action = "store_true", help = "Should Goldman & Yang's method be used for calculating dS?")
    parser.add_argument("--baseml", dest = "baseml", action = "store_true", help = "Should baseml be used instead of codeml?")
    args = parser.parse_args()
    [features_file_name, dataset_name, genome, RBP_file_name, correspondances_file_name, output_folder_name, fasta_file_name, families_file_name, output_file_name, output_folder_name, alignment_folder_name, n_sim, valid_file, gene_families, markov, new_filters, newer_filters, goldman_yang, baseml] = [args.features_file_name, args.dataset_name, args.genome, args.RBP_file_name, args.correspondances_file_name, args.output_folder_name, args.fasta_file_name, args.families_file_name, args.output_file_name, args.output_folder_name, args.alignment_folder_name, args.n_sim, args.valid_file, args.gene_families, args.markov, args.new_filters, args.newer_filters, args.goldman_yang, args.baseml]   

    #pick a random member from each paralogous family
    if features_file_name != "None":
        fs = Feature_Set(features_file_name, genome)
        fs.set_dataset(dataset_name)
        families = rw.read_families(families_file_name)
        #if the families file uses gene identifiers rather than transcript identifiers
        if gene_families:
            families = fs.convert_families_to_ENST(families, transcripts)
        fs.add_families(families)
        picked_trans = fs.pick_random_members()
        #if the fasta uses gene identifiers but the feature set uses transcript identifiers
        names = rw.read_fasta(fasta_file_name)[0]
        if picked_trans[0] not in names:
            transcripts = fs.get_transcripts()
            picked = []
            for i in picked_trans:
                picked.append(fs.convert_between_ENST_and_ENSG(i, transcripts, "ENSG"))
        else:
            picked = picked_trans
        print(len(picked))
    else:
        picked = None

    motif_dict = rw.read_motifs(RBP_file_name)

    #valid_file says which proteins pass information content criteria. Only analyze the ones that do.
    if not valid_file:
        validity = rw.read_many_fields("{0}/sufficient_information_fraction05.csv".format(output_folder_name), "\t")
        validity = list_to_dict(validity, 0, 1)
    elif valid_file == "None":
        validity = {i: "True" for i in motif_dict}
    else:
        validity = rw.read_many_fields(valid_file, "\t")        
        validity = list_to_dict(validity, 0, 1)
    protein_names = sorted([name for name in list(motif_dict.keys()) if validity[name] == "True"])

    #whether to use PAML codeml or yn00.
    if baseml:
        method = "baseml"
    elif goldman_yang:
        method = "gy"
    else:
        method = "yn"

    #write the input data for the conservation analysis to file
    input_dict_file_name = "temp_data/temp_{0}.txt".format(random.random())
    conservation.input_dict_for_dS(correspondances_file_name, alignment_folder_name, fasta_file_name, input_dict_file_name, picked = picked)
    with open(output_file_name, "w") as file:
        file.write(",".join(["protein_name", "real_dS", "mean_sim_dS", "norm_dS", "p", "motif_number"]))
        file.write("\n")
        for protein in protein_names:
            print(protein)
            motifs = motif_dict[protein]
            #use one of several different methods to generate simulant motifs
            if markov:
                simulants = nc.make_simulants_markov(motifs, n_sim)
            elif new_filters:
                simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True)
            elif newer_filters:
                simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True, no_duplicates = True, concat = False, seed = 1)               
            else:
                simulants = nc.make_simulants(motifs, n_sim)
            sim_output_file_name = "{0}/{1}_sim_ds.csv".format(output_folder_name, protein)
            #determine the conservation parameters of the current protein
            output_dict = conservation.dS_from_hits(motifs, alignment_folder_name, input_dict_file_name, n_sim = n_sim, simulants = simulants, sim_output_file_name = sim_output_file_name, method = method)
            print(output_dict)
            print("\n")
            if output_dict != None:
                file.write(",".join([protein, str(output_dict["dS"]), str(output_dict["mean simulated dS"]), str(output_dict["normalized dS"]), str(output_dict["effective p"]), str(len(motifs))]))
            else:
                file.write(",".join([protein, str(None), str(None), str(None), str(None), str(None)]))
            file.write("\n")
    os.remove(input_dict_file_name)

if __name__ == "__main__":
    main()
