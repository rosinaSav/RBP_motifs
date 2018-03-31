'''
Author: Rosina Savisaar.
Determine the synonymous rate of evolution (dS) of RBP target motif occurrences in a set of exons.
Also generate an empirical distribution of dS values for simulant motifs.
'''

import argparse
from bedtools_games import Feature_Set
import conservation
from housekeeping import flatten, list_to_dict, print_elements
import nucleotide_comp as nc
import os
import random
import read_and_write as rw
import time

def do_dS_calc(protein_names, RBP_dict, uf_fasta, df_fasta, c_fasta, n_sim, output_folder_name, correspondances_file_name, alignment_folder_name, output_file_name, regions_dict, markov, new_filters, goldman_yang):
    '''
    The core of the script that does the actual conservation calculations. 
    '''
    with open(output_file_name, "w") as file:
        #write the header of the output file
        file.write(",".join(["protein_name", "motif_number", "upstream_ds", "norm_upstream_ds", "upstream_p", "core_ds", "norm_core_ds", "core_p", "downstream_ds", "norm_downstream_ds", "downstream_p"]))
        file.write("\n")
        temp_file_name = "temp_data/temp_{0}.txt".format(random.random())
        #write the input data for the conservation analysis to file
        conservation.input_dict_for_dS(correspondances_file_name, alignment_folder_name, None, temp_file_name, map_from_regions = regions_dict)
        #loop over the RBPs
        for protein_name in protein_names:
            start_time = time.time()
            print(protein_name)
            #fetch the current motifs and generate simulants
            motifs = RBP_dict[protein_name]
            if markov:
                current_simulants = nc.make_simulants_markov(motifs, n_sim)
            elif new_filters:
                current_simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True)
            else:
                current_simulants = nc.make_simulants(motifs, n_sim)
            to_print = []
            #loop over the three different exonic subregions (5' flank, core, 3' flank)
            for pos, i in enumerate([(uf_fasta, "upstream"), (c_fasta, "core"), (df_fasta, "downstream")]):
                sim_output_file_name = "{0}/{1}_{2}_sim_ds.csv".format(output_folder_name, protein_name, i[1])
                #should you use PAML codeml or PAML yn00?
                if goldman_yang:
                    method = "gy"
                else:
                    method = "yn"
                #get conservation for that particular RBP in that particular exonic subregion
                current_dict = conservation.dS_from_hits(motifs, alignment_folder_name, temp_file_name, n_sim = n_sim, simulants = current_simulants, sim_output_file_name = sim_output_file_name, map_from_regions = regions_dict, region_name = pos, method = method)
                if current_dict != None:
                    to_print.append("{0},{1},{2}".format(current_dict["dS"], current_dict["normalized dS"], current_dict["effective p"]))
                else:
                    to_print.append(",".join([str(None), str(None), str(None)]))
            #print to output file
            to_print = ",".join(to_print)
            to_print = ",".join(["{0},{1},{2}".format(protein_name, len(motifs), to_print)])
            print(to_print)
            to_print = to_print + "\n"
            file.write(to_print)
            current_time = time.time()
            print("Time spent:")
            print(str(current_time - start_time))
            print("\n")
    os.remove(temp_file_name)

def main():
    parser = argparse.ArgumentParser(description="Calculate the conservation of a series of RBP motifs in exon cores and flanks.")
    parser.add_argument("features_file_name", type = str, help = "name of GTF file with genome features")
    parser.add_argument("dataset_name", type = str, help = "dataset name")
    parser.add_argument("genome", type = str, help = "genome assembly name")
    parser.add_argument("RBP_file_name", type = str, help = "name of file with RBP motifs")
    parser.add_argument("correspondances_file_name", type = str, help = "name of file with correspondances between genes in dataset and orthologs")
    parser.add_argument("fasta_file_prefix", type = str, help = "prefix for fasta files with the sequences")
    parser.add_argument("output_file_name", type = str, help = "file for output data")
    parser.add_argument("output_folder_name", type = str, help = "folder that will contain simulated dS scores")
    parser.add_argument("alignment_folder_name", type = str, help = "name of folder that contains alignments")
    parser.add_argument("n_sim", type = int, help = "number of simulants")
    parser.add_argument("--markov", dest = "markov", action = "store_true", help = "Should simulants be generated using a Markov model?")
    parser.add_argument("--new_filters", dest = "new_filters", action = "store_true", help = "Should simulants be generated using the sampling method but removing existing motifs and capping mononucleotide runs?")
    parser.add_argument("--goldman_yang", dest = "goldman_yang", action = "store_true", help = "Should the Goldman & Yang method (rather tahn Yang & Nielsen) be used for calculating dS?")
    parser.add_argument("--validity", dest = "validity", action = "store_true", help = "Should RBPs be filtered based on information content?")

    args = parser.parse_args()
    [features_file_name, dataset_name, genome, RBP_file_name, correspondances_file_name, output_folder_name, fasta_file_prefix, output_file_name, output_folder_name, alignment_folder_name, n_sim, markov, new_filters, goldman_yang, validity] = [args.features_file_name, args.dataset_name, args.genome, args.RBP_file_name, args.correspondances_file_name, args.output_folder_name, args.fasta_file_prefix, args.output_file_name, args.output_folder_name, args.alignment_folder_name, args.n_sim, args.markov, args.new_filters, args.goldman_yang, args.validity]   

    #make dictionary with RBPs as keys and lists of associated motifs as values
    motif_dict = rw.read_motifs(RBP_file_name)

    uf_fasta = "{0}_uf.fasta".format(fasta_file_prefix)
    df_fasta = "{0}_df.fasta".format(fasta_file_prefix)
    c_fasta = "{0}_c.fasta".format(fasta_file_prefix)

    fs = Feature_Set(features_file_name, genome)
    fs.set_dataset(dataset_name)
    transcripts = fs.get_transcripts()
    gene_name_dict = fs.get_gene_name_dict(transcripts)
    CDS = fs.get_CDS()

    #prepare the dictionary that is going to be necessray for mapping between exonic subregions and full CDSs
    regions_dict = {}
    regions_dict["gene name dict"] = gene_name_dict
    regions_dict["CDS"] = CDS
    regions_bed_file_names = [i[:-6] + ".bed" for i in [uf_fasta, c_fasta, df_fasta]]
    regions_dict["regions bed file"] = regions_bed_file_names
    regions_dict["fastas"] = [uf_fasta, c_fasta, df_fasta]

    #leave only those RBPs that pass the information content cutoff
    if validity:
        validity_5 = rw.read_many_fields("{0}/sufficient_information_fraction05_fiveprime.csv".format(output_folder_name), "\t")
        validity_core = rw.read_many_fields("{0}/sufficient_information_fraction05_core.csv".format(output_folder_name), "\t")
        validity_3 = rw.read_many_fields("{0}/sufficient_information_fraction05_threeprime.csv".format(output_folder_name), "\t")
        validity_5 = list_to_dict(validity_5, 0, 1)
        validity_core = list_to_dict(validity_core, 0, 1)
        validity_3 = list_to_dict(validity_3, 0, 1)
        protein_names = sorted([name for name in list(motif_dict.keys()) if validity_5[name] == "True" and validity_3[name] == "True" and validity_core[name] == "True"])
    else:
        protein_names = sorted(list(motif_dict.keys())) 

    #run conservation analysis
    do_dS_calc(protein_names, motif_dict, uf_fasta, df_fasta, c_fasta, n_sim, output_folder_name, correspondances_file_name, alignment_folder_name, output_file_name, regions_dict, markov, new_filters, goldman_yang)

if __name__ == "__main__":
    main()
