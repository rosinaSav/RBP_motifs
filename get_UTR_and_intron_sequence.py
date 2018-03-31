'''
Author: Rosina Savisaar.
Write UTR and intron sequences to file (aligments to macaque).
'''

from bedtools_games import bed_from_coords, coords_from_bed, Feature_Set
import ensembl_ops as eo
from housekeeping import flatten, parse_arguments, remove_file, run_process
import random
import read_and_write as rw

def main():
    description = "Write UTR and intron sequences to file."
    args = parse_arguments(description, ["features_file", "genome", "dataset", "families_file", "output_prefix"])
    [features_file, genome, dataset, families_file, output_prefix] = [args.features_file, args.genome, args.dataset, args.families_file, args.output_prefix]

    #prepare a Feature_Set object
    fs = Feature_Set(features_file, genome)
    fs.set_dataset(dataset)

    #add paralogous families
    families = rw.read_families(families_file)
    fs.add_families(families)

    #get the coordinates of various gene features
    CDSs = fs.get_CDS()
    exons = fs.get_exons()
    UTRs = fs.get_UTRs_new(exons, CDSs)
    #get the exon number of all the genes in the dataset
    exon_numbers = fs.get_exon_numbers(exons)

    #pick random members from paralogous families
    picked = fs.pick_random_members()
    #only keep multi-exon genes
    picked = [i for i in picked if exon_numbers[i] > 1]
    UTRs = {i: UTRs[i] for i in picked}
    exons = {i: exons[i] for i in picked}

    #create a second Feature_Set object that would contain all the genes in the relevant release of Ensembl
    global_fs = Feature_Set(features_file, genome)
##    global_fs.create_dataset("hg38_all", filter_trans = False)
    global_fs.set_dataset("hg38_all")
    #get the coordinates of all the exons in the second Feature_Set object
    all_exons = global_fs.get_exons()
    all_exons = flatten(list(all_exons.values()))

    introns_bed = "temp_data/introns.bed"
    temp_introns_bed = "temp_data/temp_introns.bed"
    exons_bed = "temp_data/exons.bed"

    five_utrs = flatten([UTRs[i][5] for i in UTRs])
    three_utrs = flatten([UTRs[i][3] for i in UTRs])
    #get the coordinates of full introns/exon proximal intronic regions
    introns = list(fs.get_introns_from_coords(exons).values())
    introns_upstr = flatten(list(fs.get_introns_from_coords(exons, upstream = 100).values()))
    introns_downstr = flatten(list(fs.get_introns_from_coords(exons, downstream = 100).values()))

    #only keep one randomly chosen intron per gene
    introns_temp = []
    for trans in introns:
        chosen_index = random.randrange(0, len(trans))
        introns_temp.append(trans[chosen_index])
    introns = introns_temp

    #write the full intron coordinates, as well as the coordinates of the exons from the second feature set, to bed files
    bed_from_coords(introns, temp_introns_bed)
    bed_from_coords(all_exons, exons_bed)

    #make a bed file that has only those full introns that don't overlap with any exons
    bedtools_output = run_process(["subtractBed", "-a", temp_introns_bed, "-b", exons_bed, "-A"], file_for_output = introns_bed)

    #read in the coordinates from that bed file
    introns = coords_from_bed(introns_bed, "intron")

    #remove temporary files
    remove_file(introns_bed)
    remove_file(temp_introns_bed)
    remove_file(exons_bed)

    #things that probably shouldn't be hard-coded
    query_species = "homo_sapiens"
    other_species = "macaca_mulatta"
    version = 85

    temp_coords_file = "temp_data/temp_coords.txt"

    #for all the sequence region types, get a fasta file that has the sequence from that region in both human and in macaque (separated by a pipe)
    print("5' UTRs")
    eo.get_pairwise_alignment(five_utrs, temp_coords_file, query_species, other_species, version, "{0}_five_utrs.txt".format(output_prefix))

    print("3' UTRs")
    eo.get_pairwise_alignment(three_utrs, temp_coords_file, query_species, other_species, version, "{0}_three_utrs.txt".format(output_prefix))

    print("introns")
    eo.get_pairwise_alignment(introns, temp_coords_file, query_species, other_species, version, "{0}_introns.txt".format(output_prefix))

    print("upstream introns")
    eo.get_pairwise_alignment(introns_upstr, temp_coords_file, query_species, other_species, version, "{0}_introns_upstr.txt".format(output_prefix))

    print("downstream introns")
    eo.get_pairwise_alignment(introns_downstr, temp_coords_file, query_species, other_species, version, "{0}_introns_downstr.txt".format(output_prefix))



if __name__ == "__main__":
    main()
