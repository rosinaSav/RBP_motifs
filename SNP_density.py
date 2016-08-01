from bedtools_games import bed_to_CDS_indices, convert_coords, Feature_Set, intersect_bed_return_bed, uniquify_lines, write_features_to_bed
import conservation
from housekeeping import flatten, list_to_dict, parse_arguments, print_elements, run_in_parallel
import my_stats as ms
import nucleotide_comp as nc
import read_and_write as rw

def write_relative_coords(CDS, input_file_name, output_file_name, seqs, names, remove_C = False):
    '''
    Write a file where for each CDS analyzed, you have a list of all the SNPs that map to it (in coordinates relative to the CDS sequence).
    '''
    relative_coords_dict = {i: [] for i in list(CDS.keys())}
    counter = 0
    with open(input_file_name) as file:
        error_counter = 0
        #loop over bed file with SNPs
        for line in file:
            if counter % 1000 == 0:
                print(counter)
            counter = counter + 1
            bed_record = line.split("\t")
            info = bed_record[6].split("$")
            #reported ancestral base
            ref_base = info[1]
            #this is so you could filter out certain SNPs based on the ancestral allele
            ancestral_allowed = nc._canon_bases_.copy()
            if remove_C:
                ancestral_allowed.remove("C")
                ancestral_allowed.remove("G")
            #filter out anything that isn't simple SNPs
            if ref_base in ancestral_allowed:
                var_base = [i for i in info[2].split(",") if i in nc._canon_bases_]
                if var_base:
                    #the SNP bed files contain information on which transcript they overlap
                    trans = bed_record[3]
                    #convert bed record to indices relative to the CDS
                    #mismatch_warning_only = True means that if the SNP coordinates don't map onto the CDS,
                    #return an error message but don't crash. This happens when, say, the SNP was on a patch in hg19 and
                    #as a result, the hg38 CDS coordinates mapped onto hg19 but then when converting the SNP back into hg38,
                    #the conversion didn't work.
                    current_index = bed_to_CDS_indices(bed_record, CDS[trans], mismatch_warning_only = True)[0]
                    if current_index != ("error"):
                        #check if the reported ancestral allele matches the base that is at that position in the CDS
                        #it's expected that sometimes it won't but it should most of the time if everything's worked out properly
                        strand = CDS[trans][0][0][6]
                        if strand == "-":
                            ref_base = nc.rev_comp(ref_base)
                        CDS_base = seqs[names.index(trans)][current_index]
                        if ref_base != CDS_base:
                            print("PROBLEM!")
                            print(ref_base)
                            print(CDS_base)
                            print(var_base)
                            print("\n")
                        relative_coords_dict[trans].append(str(current_index))
                    else:
                        error_counter = error_counter + 1
                        relative_coords_dict[trans].append("NA")
    with open(output_file_name, "w") as file:
        for trans in relative_coords_dict:
            to_write = "\t".join([trans, ",".join(relative_coords_dict[trans])])
            file.write(to_write)
            file.write("\n")
    relative_coords_dict = {}
    print("Number of conversion errors:")
    print(error_counter)

def main():
    description = "Calculate the SNP density of a set of motifs."
    args = parse_arguments(description, ["features_file", "genome", "summary_file_name", "dataset", "families_file", "sequence_file", "motifs_file", "n_sim", "output_file_name", "protein", "merge", "flanks", "remove_C"], flags = [10, 11, 12], ints = [7])
    [features_file, genome, summary_file_name, dataset, families_file, sequence_file, motifs_file, n_sim, output_file_name, protein, merge, flanks, remove_C] = [args.features_file, args.genome, args.summary_file_name, args.dataset, args.families_file, args.sequence_file, args.motifs_file, args.n_sim, args.output_file_name, args.protein, args.merge, args.flanks, args.remove_C]

    print(output_file_name)

    #make a dictionary with RBPs as keys and lists of associated motifs as values
    motifs = rw.read_motifs(motifs_file)
    #if merge, consider all motifs, otherwise only the most depleted ones
    if not merge:
        summary_data = rw.read_many_fields(summary_file_name, "\t")
        #because some of the files are tab-separated, while others are comma-separated and have a header row
        if len(summary_data[0]) == 1:
            summary_data = rw.read_many_fields(summary_file_name, ",")
            summary_data = summary_data[1:]
        #summary_dict has RBPs as keys and the p-values from the single-RBP analsysi as values
        summary_dict = list_to_dict(summary_data, 0, 4, floatify = True)
        motifs = {i: motifs[i] for i in motifs if summary_dict[i] > 0.9}

    #flatten and uniquify motifs
    motifs = list(set(flatten(list(motifs.values()))))

    fs = Feature_Set(features_file, genome)
    fs.set_dataset(dataset)
    CDS = fs.get_CDS()
    transcripts = fs.get_transcripts()
    gene_name_dict = fs.get_gene_name_dict(transcripts)
    #if doing full CDSs, pick a random member from each paralogous family
    if not flanks:
        families = rw.read_families(families_file)
        fs.add_families(families)
        picked = fs.pick_random_members()
        map_from_regions = None
    #the exonic subregions fastas already have regions from only one gene per paralogous family so you don't need to worry about that here
    else:
        picked = None
        bed_file_name = "{}.bed".format(sequence_file)
        sequence_file = "{0}.fasta".format(sequence_file)
        map_from_regions = conservation.map_regions_to_CDS(sequence_file, bed_file_name, fs, gene_name_dict, CDS)
        for idn in map_from_regions:
            map_from_regions[idn]["idn"] = fs.convert_between_ENST_and_ENSG(map_from_regions[idn]["idn"], gene_name_dict, "ENST")

    #read in the CDS sequences
    names, seqs = rw.read_fasta(sequence_file)
    if not flanks:
        #the full CDS fasta uses gene identifiers but you want to be using transcript identifiers in the analysis
        names = [fs.convert_between_ENST_and_ENSG(i, gene_name_dict, "ENST") for i in names]

    #normally, CDS coordinates are stored along with phase information on each CDS region
    #remove that extra information
    flat_CDS = [[j[0] for j in i] for i in list(CDS.values())]
    CDS_bed = "RBP/multi-exon_CDS.bed"
    CDS_bed_hg19 = "RBP/multi-exon_CDS_hg19.bed"
    write_features_to_bed(flat_CDS, CDS_bed, modify_chr_ids = True)
    
    SNP_file_name_hg19 = "RBP/multi-exon_SNPs_hg19_tabix.bed"
    SNP_file_name_hg38 = "RBP/multi-exon_SNPs_hg38_tabix.bed"
    #convert CDS coordinates to hg19 because that's what 1000Genomes uses
    convert_coords(CDS_bed, CDS_bed_hg19, "hg38", "hg19")
    #get the SNPs overlapping the CDS regions
    conservation.tabix(CDS_bed_hg19, SNP_file_name_hg19)
    #convert them back to hg38
    convert_coords(SNP_file_name_hg19, SNP_file_name_hg38, "hg19", "hg38")

    if remove_C:
        CDS_SNP_file_name = "RBP/multi-exon_SNPs_hg38_relative_no_C.bed"
    else:
        CDS_SNP_file_name = "RBP/multi-exon_SNPs_hg38_relative.bed"

    #convert the SNP coordinates to coordinates relative to the CDS that they overlap
    write_relative_coords(CDS, SNP_file_name_hg38, CDS_SNP_file_name, seqs, names, remove_C = remove_C)

    SNP_dict = rw.read_many_fields(CDS_SNP_file_name, "\t")
    SNP_dict = list_to_dict(SNP_dict, 0, 1)
    SNP_dict = {i: [int(j) for j in SNP_dict[i].split(",") if (SNP_dict[i] and ("NA" not in SNP_dict[i]))] for i in SNP_dict}

    #what fraction of fourfold degenarte sites within motifs overlap SNPs?
    real_fraction = conservation.get_SNP_density(picked, motifs, seqs, names, SNP_dict, map_from_regions)
    print(real_fraction)
    sim_fractions = []
    print("Simulants:")
    #and in simulants?
    for sim in range(n_sim):
        if sim%100 == 0:
            print(sim)
        #n_sim is passed as 1 to make_simulants because you're only making one set at a time
        current_simulants = nc.make_simulants(motifs, 1, cap_runs = True, remove_existing = True, no_duplicates = True, concat = False)[0]
        sim_fractions.append(conservation.get_SNP_density(picked, current_simulants, seqs, names, SNP_dict, map_from_regions))
        print("{0}:{1}".format(sim, sim_fractions[-1]))

    norm_fraction = ms.normalize(real_fraction, sim_fractions)
    p = ms.calc_eff_p(real_fraction, sim_fractions, greater = False)
    sim_output_file_name = "{0}_sim.csv".format(output_file_name[:-4])
    with open(output_file_name, "w") as file:
        file.write("real_fraction\tnorm_fraction\tp\n")
        file.write("\t".join([str(real_fraction), str(norm_fraction), str(p)]))
        file.write("\n")
    with open(sim_output_file_name, "w") as file:
        for fraction in sim_fractions:
            file.write("{0}\n".format(fraction))
   
if __name__ == "__main__":
    main()
