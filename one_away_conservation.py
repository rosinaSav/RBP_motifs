'''
Author: Rosina Savisaar.
Calculate the conservation of k-mers that are a single point mutation away
from being part of a set of motifs.
'''

from bedtools_games import Feature_Set
from Bio import SeqIO
import conservation
from housekeeping import flatten, list_to_dict, parse_arguments, run_in_parallel
import my_stats as ms
import nucleotide_comp as nc
import numpy as np
import read_and_write as rw

def get_mutation_to_motif(picked, correspondance_dict, alignment_folder_name, CDS, names, motifs, neighbours, neighbour_lengths):
    '''
    Determine the number of fourfold degenerate sites (summed across all sequences) that are a single base substitution away from one of the motifs.
    Determine what fraction of these actually have a base that would generate the motif in the focal species at the orthologous position in the ortholog species.
    '''
    site_number = 0
    mutation_score = 0
    for pos, gene in enumerate(picked):
        orth_idn = correspondance_dict[gene]
        if "_" not in orth_idn:
            orth_idn = orth_idn + "_0"
        phy_file_name = "{0}/{1}_{2}.phy".format(alignment_folder_name, gene, orth_idn)
        aligned_sequences = [str(i.seq) for i in SeqIO.parse(phy_file_name, "phylip-sequential")]
        current_CDS = CDS[names.index(gene)]
        mutation_score, site_number = conservation.mutation_to_motif(current_CDS, aligned_sequences, motifs, neighbours, neighbour_lengths, mutation_score, site_number) 
    return(site_number, mutation_score)

def main():

    description = "Calculate the conservation of k-mers that are a single point mutation away from being part of a set of motifs."
    args = parse_arguments(description, ["motifs_file_name", "summary_file_name", "output_folder_name", "p_column", "alignment_folder_name", "correspondances_file_name", "output_file_name", "dataset_name", "features_file_name", "n_sim", "output_suffix", "sequences_file_name", "families_file_name", "genome", "by_RBP"], ints = [3, 9], flags = [14])
    [motifs_file_name, summary_file_name, output_folder_name, p_column, alignment_folder_name, correspondances_file_name, output_file_name,  dataset_name, features_file_name, n_sim, output_suffix, sequences_file_name, families_file_name, genome, by_RBP] = [args.motifs_file_name, args.summary_file_name, args.output_folder_name, args.p_column, args.alignment_folder_name, args.correspondances_file_name, args.output_file_name, args.dataset_name, args.features_file_name, args.n_sim, args.output_suffix, args.sequences_file_name, args.families_file_name, args.genome, args.by_RBP]

    RBPs = rw.read_motifs(motifs_file_name)

    #only leave those RBPs hat pass information content criteria
    validity = rw.read_many_fields("{0}/sufficient_information_fraction05.csv".format(output_folder_name), "\t")
    validity = list_to_dict(validity, 0, 1)
    RBPs = {i: RBPs[i] for i in RBPs if validity[i] == "True"}

    #if you're not doing this by RBP, pool motifs from the most significantly depleted sets
    if not by_RBP:
        summary_data = rw.read_many_fields(summary_file_name, "\t")
        if len(summary_data[0]) == 1:
            summary_data = rw.read_many_fields(summary_file_name, ",")    
        summary_dict = list_to_dict(summary_data, 0, p_column, floatify = True)            
        RBPs = {i: RBPs[i] for i in RBPs if summary_dict[i] > 0.9}
        motifs = list(set(flatten(list(RBPs.values()))))
        RBPs = {"all": motifs}

    #randomly pick one gene from each paralogous family
    fs = Feature_Set(features_file_name, genome)
    fs.set_dataset(dataset_name)
    transcripts = fs.get_transcripts()
    families = rw.read_families(families_file_name)
    families = fs.convert_families_to_ENST(families, transcripts)
    fs.add_families(families)
    picked_from_families = fs.pick_random_members()
    gene_name_dict = fs.get_gene_name_dict(transcripts)
    picked = [fs.convert_between_ENST_and_ENSG(i, gene_name_dict, "ENSG") for i in picked_from_families]

    names, CDS = rw.read_fasta(sequences_file_name)

    #make a dictionary where the keys are genes from the focal species and the values are orthologs from another species
    correspondances = rw.read_many_fields(correspondances_file_name, ",")
    correspondance_dict = {}
    for i in correspondances:
        correspondance_dict[i[0]] = i[1]

    output_dict = {}

    #loop over the RBPs
    for protein in sorted(RBPs):

        #fetch the current motifs
        print(protein)
        motifs = RBPs[protein]
        print("There are {0} motifs.".format(len(motifs)))
        #generate all unique motifs that are a single base substitution away from one of the motifs but are not actually in the set
        neighbours = nc.get_neighbours(motifs)
        print("There are {0} neighbours.".format(len(neighbours)))            

        #make simulants for the motifs. don't allow simulants to be part of the set of neighbours.
        simulants = nc.make_simulants(motifs, n_sim, remove_existing = True, cap_runs = True, exclude = neighbours, no_duplicates = True, concat = False)

        neighbour_lengths = [len(i) for i in neighbours]        
        neighbours = nc.motif_to_regex(neighbours)

        #determine the true frequency at which fourfold degenarte sites that are a single substitution away from a motif in human actually contain the base that
        #would give rise to the motif in the orthologous species
        site_number = 0
        mutation_score = 0
        motifs = [list(i) for i in motifs]
        true_result = run_in_parallel(picked, ["foo", correspondance_dict, alignment_folder_name, CDS, names, motifs, neighbours, neighbour_lengths], get_mutation_to_motif) 
        for i in true_result:
            current = i.get()
            site_number = site_number + current[0]
            mutation_score = mutation_score + current[1]
        if site_number > 0:
            real_fraction = mutation_score/site_number
        else:
            real_fraction = None
        print("Real fraction:")
        print(real_fraction)

        neighbours = ""      
        sim_site_numbers = np.zeros((n_sim))
        sim_mutation_scores = np.zeros((n_sim))

        #obtain this estimate also for each simulant set
        #I'm doing this in this awkward manner because I don't have enough RAM to hold all the simulated neighbours in memory at once
        for sim in range(n_sim):
            if sim%10 == 0:
                print(sim)
            current_simulants = simulants[sim]
            current_neighbours = nc.get_neighbours(current_simulants)
            current_neighbour_lengths = [len(i) for i in current_neighbours]        
            current_neighbours = nc.motif_to_regex(current_neighbours)
            current_simulants = [list(i) for i in current_simulants]
            current_result = run_in_parallel(picked, ["foo", correspondance_dict, alignment_folder_name, CDS, names, current_simulants, current_neighbours, current_neighbour_lengths], get_mutation_to_motif)
            for i in current_result:
                current = i.get()
                sim_site_numbers[sim] = sim_site_numbers[sim] + current[0]
                sim_mutation_scores[sim] = sim_mutation_scores[sim] + current[1]

        #normalize the real fraction, calculate p
        sim_fractions = np.divide(sim_mutation_scores, sim_site_numbers)
        sim_fractions = [i for i in sim_fractions if i != np.inf]
        p = ms.calc_eff_p(real_fraction, sim_fractions, greater = False)
        norm_fraction = ms.normalize(real_fraction, sim_fractions) 

        output_dict[protein] = [protein, mutation_score, site_number, real_fraction, np.mean(sim_fractions), p, norm_fraction]
        print(output_dict[protein])
        
    with open(output_file_name, "w") as output_file:
        #write header to output file
        output_file.write("protein\tmutation score\tsite number\treal fraction\tmean sim fraction\tp\tnormalized fraction\n")
        #write the rest of the output data
        for protein in sorted(list(output_dict.keys())):
            to_write = output_dict[protein]
            to_write = [str(i) for i in to_write]
            output_file.write("\t".join(to_write))
            output_file.write("\n")

    
if __name__ == "__main__":
    main()   
