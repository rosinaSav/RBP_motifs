from bedtools_games import Feature_Set
from housekeeping import list_to_dict, parse_arguments, print_elements
import nucleotide_comp as nc
import read_and_write as rw

def main():
    description = "Check whether a set of motifs has enough information (= the density is high enough) for it to be included in subsequent analysis."
    args = parse_arguments(description, ["motifs_file_name", "summary_file_name", "dataset_name", "sim_file_folder", "sim_file_suffix", "fraction", "threshold", "output_file_name", "features_file_name", "genome", "families_file_name", "fasta_name", "density_column"], floats = [5], ints = [6, 12])
    [motifs_file_name, summary_file_name, dataset_name, sim_file_folder, sim_file_suffix, fraction, threshold, output_file_name, features_file_name, genome, families_file_name, fasta_name, density_column] = [args.motifs_file_name,
                                                                                                                                                                              args.summary_file_name,
                                                                                                                                                                              args.dataset_name,
                                                                                                                                                                              args.sim_file_folder,
                                                                                                                                                                              args.sim_file_suffix,
                                                                                                                                                                              args.fraction, args.threshold,
                                                                                                                                                                              args.output_file_name,
                                                                                                                                                                              args.features_file_name,
                                                                                                                                                                              args.genome,
                                                                                                                                                                              args.families_file_name,
                                                                                                                                                                                     args.fasta_name,
                                                                                                                                                                                     args.density_column]

    proteins, motifs = rw.read_fasta(motifs_file_name)

    #make a dictionary with protein names as keys and the raw densities as values
    density_data = rw.read_many_fields(summary_file_name, "\t")
    #the output file can have different delimiters
    if len(density_data[0]) == 1:
        density_data = rw.read_many_fields(summary_file_name, ",")
    #to get rid of a header file in case there is one
    try:
        density_data[0][1] = float(density_data[0][1])
    except ValueError:
        density_data = density_data[1:]
    #you can't hard-code which column in the density output file has the raw density because this will be different in cores/flanks vs full CDS
    density_dict = list_to_dict(density_data, 0, density_column)
    for i in density_dict:
        density_dict[i] = float(density_dict[i])

    #if you don't provide a fasta file, the total sequence length will be determined from the dataset information
    if fasta_name == "None":
        fs = Feature_Set(features_file_name, genome)
        fs.set_dataset(dataset_name)
        CDS = fs.get_CDS()
        lengths = fs.get_lengths(CDS, CDS = True)
        transcripts = fs.get_transcripts()
        gene_name_dict = fs.get_gene_name_dict(transcripts)

        #you can't just loop over _lengths_ as you normally would with a dictionary because you
        #would end up cycling back onto identifiers you've already converted and try to convert an ENSG to an ENSG
        for identifier in sorted(list(lengths.keys())):
            lengths[fs.convert_between_ENST_and_ENSG(identifier, gene_name_dict, "ENSG")] = lengths.pop(identifier)

        families = rw.read_families(families_file_name)
        fs.add_families(families)

        #those proteins will be considered that have motifs that overlap with more than 100 basepairs (or where at least 50% of the simulants do)
        averaged_lengths = fs.average_over_families(lengths)
        summed_lengths = sum(list(averaged_lengths.values()))
    else:
        lengths = {}
        names, seq = rw.read_fasta(fasta_name)
        for pos, name in enumerate(names):
            lengths[name] = len(seq[pos])
        summed_lengths = sum(list(lengths.values()))
        fs = None

    threshold_density = threshold/summed_lengths

    output_dict = nc.check_amount_of_information(density_dict, lengths, sim_file_folder, sim_file_suffix, fraction, threshold_density, fs = fs)
    sorted_proteins = sorted(list(output_dict.keys()))
    with open(output_file_name, "w") as file:
        for i in sorted_proteins:
            file.write("\t".join([i, str(output_dict[i])]))
            file.write("\n")
    

if __name__ == "__main__":
    main()
