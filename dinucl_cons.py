from bedtools_games import Feature_Set
import conservation
from housekeeping import flatten, make_dir, parse_arguments
import nucleotide_comp as nc
import read_and_write as rw

def main():
    description = "Calculate the conservation of a set of motifs separately for each dinucleotide."
    args = parse_arguments(description, ["features_file_name", "dataset_name", "genome", "RBP_file_name", "correspondances_file_name", "fasta_file_name", "families_file_name", "output_file_name", "alignment_folder_name", "flanks"], flags = [9])
    [features_file_name, dataset_name, genome, RBP_file_name, correspondances_file_name, fasta_file_name, families_file_name, output_file_name, alignment_folder_name, flanks] = [args.features_file_name, args.dataset_name, args.genome, args.RBP_file_name, args.correspondances_file_name, args.fasta_file_name, args.families_file_name, args.output_file_name, args.alignment_folder_name, args.flanks]

    #prepare an object for storing the genome annotations associated to the sequences in the sequence file
    fs = Feature_Set(features_file_name, genome)
    fs.set_dataset(dataset_name)
    #make a dictionary with RBPs as keys and lists of associated motifs as values
    motif_dict = rw.read_motifs(RBP_file_name)
    transcripts = fs.get_transcripts()
    gene_name_dict = fs.get_gene_name_dict(transcripts)

    #if working with full CDSs
    if not flanks:
        #pick a random memebr from each paralogous family
        families = rw.read_families(families_file_name)
        families = fs.convert_families_to_ENST(families, transcripts)
        fs.add_families(families)
        picked_trans = fs.pick_random_members()
        picked = []
        for i in picked_trans:
            for j in gene_name_dict:
                if gene_name_dict[j][0][4] == i:
                    picked.append(j)
        print(len(picked))
        map_from_regions = None
    #if working with exon subregions
    #my exon subregions file already has regions from only one transcript per paralogous family
    else:
        picked = None
        CDS = fs.get_CDS()
        bed_file_name = "{}.bed".format(fasta_file_name)
        fasta_file_name = "{0}.fasta".format(fasta_file_name)
        map_from_regions = conservation.map_regions_to_CDS(fasta_file_name, bed_file_name, fs, gene_name_dict, CDS)       

    #generate all possible DNA dinucleotides
    dinucl = nc.generate_all_kmers(2)

    motifs = flatten(list(motif_dict.values()))

    with open(output_file_name, "w") as file:
        file.write("dinucleotide\tmotif rate\tmotif frequency\tnonmotif rate\tnonmotif frequency\n")
        #calculate the rate of evolution wihtin vs outside of motifs separately for each dinucleotide
        freqs_dict = conservation.cons_by_dinucl(fasta_file_name, motifs, correspondances_file_name, alignment_folder_name, dinucl, picked = picked, map_from_regions = map_from_regions)
        for dint in sorted(list(freqs_dict.keys())):
            if (freqs_dict[dint]["subst. in motifs"] != None) and (freqs_dict[dint]["subst. in non-motifs"] != None):
                to_write = [dint, freqs_dict[dint]["subst. in motifs"], freqs_dict[dint]["frequency in motifs"], freqs_dict[dint]["subst. in non-motifs"], freqs_dict[dint]["frequency in non-motifs"]]
                to_write = "\t".join([str(i) for i in to_write])
                file.write(to_write)
                file.write("\n")
        #get an over-all estimate by taking a weighted avergae (weighted by dinucleotide frequency) of the frequencies of all the different dinucleotides
        output_dict = conservation.weight_cons_by_dinucl(freqs_dict, dinucl)
        print(output_dict)

if __name__ == "__main__":
    main()
