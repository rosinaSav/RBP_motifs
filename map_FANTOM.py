from bedtools_games import Feature_Set, intersect_bed
from housekeeping import list_to_dict, parse_arguments, print_elements
import numpy as np
import re
import read_and_write as rw
import scipy.stats
import sys

def main():

    description = "Take an output file from prepare_FANTOM.py and make a file with the expression data for each gene."
    args = parse_arguments(description, ["intronless_folder_name", "introncontaining_folder_name", "features_file_name", "genome", "intronless_families", "introncontaining_families", "promoters_file_name", "cage_file_name", "output_file_name", "TPM_threshold","filter_genes"], flags = [10], ints = [9])
    [intronless_folder_name, introncontaining_folder_name, features_file_name, genome, intronless_families, introncontaining_families, promoters_file_name, cage_file_name, output_file_name, TPM_threshold, filter_genes] = [args.intronless_folder_name, args.introncontaining_folder_name, args.features_file_name, args.genome, args.intronless_families, args.introncontaining_families, args.promoters_file_name, args.cage_file_name, args.output_file_name, args.TPM_threshold, args.filter_genes]

    suffix = ""
    if filter_genes:
        #if you only want to use the genes included in your dataset

        #make a Feature_Set object (a set of Ensembl transcript identifiers, along with associated gene annotations)
        fs = Feature_Set(features_file_name, genome)
        #probably shouldn't be hard-coded
        fs.set_dataset("multi-exon")
        #associate a clustering of genes into paralogous families
        ic_families = rw.read_families(introncontaining_families)
        fs.add_families(ic_families)

        transcripts = fs.get_transcripts()
        #a dictionary where the keys are gene identifiers and the values are associated transcript identifiers
        gene_name_dict = fs.get_gene_name_dict(transcripts)

        #write a bed file with promoter coordinates
        with open(promoters_file_name, "w") as file:
            for idn in transcripts:
                #you're also converting from base 1 to base 0
                #note that you're putting down gene identifiers rather than transcript identifiers
                if transcripts[idn][6] == "+":
                    current_line = ["chr" + transcripts[idn][0], transcripts[idn][2] - 500 - 1, transcripts[idn][2] + 500, transcripts[idn][5], ".", transcripts[idn][6]]
                else:
                    current_line = ["chr" + transcripts[idn][0], transcripts[idn][3] - 500 - 2, transcripts[idn][3] + 500 - 1, transcripts[idn][5], ".", transcripts[idn][6]]
                file.write("\t".join([str(i) for i in current_line]))
                file.write("\n")


        #get CAGE data, as formatted by prepare_FANTOM.py (and CrossMapped)
        cage_data = rw.read_many_fields(cage_file_name, "\t")
        #keep only those CAGE peaks that overlap with the promoters
        overlapping_peaks = intersect_bed(cage_data, promoters_file_name, write_both = True,
                      force_strand = True, no_name_check = False, no_dups = False, bed_input = True)

        #make a dictionary where the keys are gene identifiers and the values are the coordinates of all the peaks that overlap the associated promoter
        peaks_dict = {idn: [] for idn in gene_name_dict}
        for peak in overlapping_peaks:
            gene_ID = peak[5][3]
            peaks_dict[gene_ID].append(peak)

        #sum the TPMs associated to each gene within samples and across peaks
        sum_dict = {}
        np.set_printoptions(suppress = True)
        for i in peaks_dict:
            if len(peaks_dict[i]) > 0:
                current_mat = np.array([[float(k) for k in j[3].split("|")] for j in peaks_dict[i]])
                sums = np.mean(current_mat, axis = 0)
                sum_dict[i] = sums

        #take the various expression measures
        final_dict = {}
        for gene in sum_dict:
            expressed = len([i for i in sum_dict[gene] if i > TPM_threshold])
            fraction = expressed/len(sum_dict[gene])
            maximum = np.max(sum_dict[gene])
            median_expr = np.median(sum_dict[gene])
            median_if_expressed = np.median([i for i in sum_dict[gene] if i > TPM_threshold])
            final_dict[gene] = [fraction, maximum, median_expr, median_if_expressed]

        if filter_genes:
            gene_name_dict = fs.get_gene_name_dict(transcripts)
            final_dict_ENST = {fs.convert_between_ENST_and_ENSG(i, gene_name_dict, "ENST"): final_dict[i] for i in final_dict}                       
            final_dict_ENST = fs.average_over_families_2d(final_dict_ENST)
            #this is an embarrassing hack but hopefully nobody will ever read this
            #this is so you could average over families!
            final_dict_ENST = {fs.convert_between_ENST_and_ENSG(i, transcripts, "ENSG", families = True): final_dict_ENST[i] for i in final_dict_ENST}                       
        else:
            final_dict_ENST = final_dict
        
        with open(output_file_name, "w") as file:
            file.write("gene\tbreadth\tmax\tmedian\tmedian_expr\n")
            for i in sorted(list(final_dict_ENST.keys())):
                if final_dict_ENST[i] != None:
                    file.write("\t".join([i] + [str(j) for j in final_dict_ENST[i]]))
                    file.write("\n")

    else:
        with open(cage_file_name) as file, open(output_file_name, "w") as output_file:
            for pos, line in enumerate(file):
                line = line.rstrip("\n")
                line = line.split("\t")
                line = [float(i) for i in line[3].split("|")]
                expressed = len([i for i in line if i > TPM_threshold])
                fraction = expressed/len(line)
                maximum = np.max(line)
                median_expr = np.median(line)
                to_write = [pos, fraction, maximum, median_expr]
                to_write = [str(i) for i in to_write]
                output_file.write("\t".join(to_write))
                output_file.write("\n")                   
    

if __name__ == "__main__":
    main()
