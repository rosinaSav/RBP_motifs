'''
Author: Rosina Savisaar.
Prepare a clean dataset of protein-coding genes.
'''

import argparse
from bedtools_games import Feature_Set
from conservation import check_ORF_integrity, keep_conserved_pc
import csv
from housekeeping import flatten, make_dir, print_elements
import re
import read_and_write as rw

def main():  
    parser = argparse.ArgumentParser(description="Prepare a clean dataset of protein-coding genes.")
    parser.add_argument("features_file_name", type = str, help = "name of GTF file with genome features")
    parser.add_argument("ortholog_features_file_name", type = str, help = "name of GTF file with genome features for the orthologous genome")
    parser.add_argument("genome", type = str, help = "genome assembly name")
    parser.add_argument("ortholog_genome", type = str, help = "ortholog genome assembly name")
    parser.add_argument("dataset_name", type = str, help = "dataset name")
    parser.add_argument("ortholog_dataset_name", type = str, help = "ortholog dataset name")
    parser.add_argument("orthologs_file_name", type = str, help = "csv with orthologous pairs")
    parser.add_argument("dS_threshold", type = float, help = "csv with orthologus pair")
    parser.add_argument("alignment_folder", type = str, help = "folder where phy alignment files will be stored")
    parser.add_argument("raw_orth_seq_file", type = str, help = "file with the raw ortholog CDS sequences (downloaded via ensembl biomart)")

    args = parser.parse_args()
    [features_file_name, ortholog_features_file_name, genome, ortholog_genome, dataset_name, ortholog_dataset_name, orthologs_file_name, dS_threshold, alignment_folder, raw_orth_seq_file] = [args.features_file_name,
                                                                                                                                                          args.ortholog_features_file_name, args.genome, args.ortholog_genome,
                                                                                                                                                                                               args.dataset_name,
                                                                                                                                                                                               args.ortholog_dataset_name,
                                                                                                                                                                                               args.orthologs_file_name,
                                                                                                                                                                                               args.dS_threshold,
                                                                                                                                                                                               args.alignment_folder,
                                                                                                                                                                                               args.raw_orth_seq_file]
    make_dir(alignment_folder)
    trans_id_pattern = re.compile("ENS\w*T\d*")
    ids_to_keep = []
    #loop over an ensembl GTF file
    with open(features_file_name) as features_file:
        #skip the metadata
        for i in range(5):
            features_file.readline()
            
        for i in features_file:
            #only consider features that have been localized to chromosomes and that are from protein-coding genes
            if "PATCH" not in i and "gene_biotype \"protein_coding\"" in i and i[0] in "123456789XY" and i[1] in "0123456789XY\t":
                trans_id_obj = re.search(trans_id_pattern, i)
                if trans_id_obj:
                    trans_id = trans_id_obj.group(0)
                    #store the transcript ID
                    ids_to_keep.append(trans_id)

    #make a list of the unique transcript IDs you got in the previous step
    ids_to_keep = list(set(ids_to_keep))

    #create a feature set object from the transcript IDs,
    #that is to say, make a file that has all the associated gene feature annotations
    fs = Feature_Set(features_file_name, genome)
    #the dataset only needs to be created if it didn't exist previously
##    fs.create_dataset(dataset_name, input_list = ids_to_keep)
    fs.set_dataset(dataset_name)
    print("Created dataset with {0} transcripts.".format(len(fs.names)))
    #this file will have the mappings between genes from the focal species and genes from the orthologus species
    final_pairs_file_name = "general/{0}_{1}_pc_pairs.csv".format(genome, ortholog_genome)

    CDS = fs.get_CDS()
    CDS = {i: CDS[i] for i in CDS if CDS[i]}
    #write the full ORF sequences of the genes to FASTA, filtering based on reading frame integrity. Also check that
    #there are no premature termination codons.
    fs.write_full_CDS(CDS, check_ORF = True, bare_name = True, PTC_check = True)

    ids_to_keep = rw.read_fasta("{0}_{1}_full_CDS.fasta".format(fs.features_file_name[:-4], fs.dataset))[0]

    print("{0} transcripts pass the check for ORF integrity.".format(len(ids_to_keep)))

    transcripts = fs.get_transcripts()
    transcripts = {i: transcripts[i] for i in ids_to_keep}

    #for genes with several associated transcript IDs, only keep the longest.   
    gene_name_dict = fs.get_gene_name_dict(transcripts)
    ids_to_keep = []
    for gene in gene_name_dict:
        current_CDS = [CDS[j] for j in gene_name_dict[gene]]
        current_lengths = [sum([j[0][3] - j[0][2] + 1 for j in k]) for k in current_CDS]
        id_to_keep = gene_name_dict[gene][current_lengths.index(max(current_lengths))]
        ids_to_keep.append(id_to_keep)

    print("After only keeping one transcript per gene (the longest), {0} transcripts remain.".format(len(ids_to_keep)))

    #this is a file that has the orthologs of your gens from Ensmebl biomart
    orth_data = rw.read_many_fields(orthologs_file_name, ",")
    #make a dictionary for the gene-to-ortholog mapping
    pairs_dict = {}

    for line in orth_data:
        if line[1] not in pairs_dict:
            pairs_dict[line[1]] = []
        pairs_dict[line[1]].append(line[2])

    #only keep genes for which there is an ortholog in the comparator species
    #transcript identifiers
    ids_to_keep = [i for i in ids_to_keep if i in pairs_dict]

    #gene identifiers
    orth_ids_to_keep = list(pairs_dict.values())
    orth_ids_to_keep = list(set(flatten(orth_ids_to_keep)))

    #create a feature set for the other species based on the genes that are orthologous to the genes in your focal set
    orth_fs = Feature_Set(ortholog_features_file_name, ortholog_genome)
##    orth_fs.create_dataset(ortholog_dataset_name, input_list = orth_ids_to_keep, input_type = "gene")
    orth_fs.set_dataset(ortholog_dataset_name)
    orth_CDS = orth_fs.get_CDS()
    orth_CDS = {i: orth_CDS[i] for i in orth_CDS if orth_CDS[i]}
    #write the ortholog ORFs to FASTA. Filter based on reading frame integrity and PTC content.
    orth_fs.write_full_CDS(orth_CDS, check_ORF = True, bare_name = True, PTC_check = True)
    orth_full_CDS_file = "{0}_{1}_full_CDS.fasta".format(ortholog_features_file_name[:-4], ortholog_dataset_name)

    #in some cases, if the genome assembly for the ortholog is not very good, it can take forever to get the sequences using faidx.
    #In that case, you can get the sequences via biomart. Uncomment the code below!
##    rw.write_names(list(orth_CDS.keys()), "general/{0}_trans_IDs.txt".format(ortholog_dataset_name))
##    with open(raw_orth_seq_file) as file:
##        raw_orth_seq = "".join(file)
##    raw_orth_seq = re.sub("([A-Z])\n([A-Z])", "\\1\\2", raw_orth_seq)
##    raw_orth_seq = raw_orth_seq.split("\n")
##    raw_orth_seq = [i for i in raw_orth_seq if len(i) > 0]
##    raw_orth_names = [i for i in raw_orth_seq if i[0] == ">"]
##    raw_orth_seq = [i for i in raw_orth_seq if i[0] != ">"]

##    with open(orth_full_CDS_file, "w") as file:
##        for pos, seq in enumerate(raw_orth_seq):
##            ORF_check = check_ORF_integrity(seq, PTC_check = True)
##            if ORF_check[0]:
##                file.write("{0}\n".format(raw_orth_names[pos]))
##                file.write("{0}\n".format(seq))
##            else:
##                print(pos)
##                print(ORF_check[1])
##                print(raw_orth_names[pos])
##                print(seq)
##                print("\n")            

    #read in the full ORF sequences from both species
    CDS_names, CDS_seq = rw.read_fasta("{0}_{1}_full_CDS.fasta".format(fs.features_file_name[:-4], fs.dataset))
    orth_CDS_names, orth_CDS_seq = rw.read_fasta(orth_full_CDS_file)

    orth_transcripts = orth_fs.get_transcripts()

    orth_gene_name_dict = orth_fs.get_gene_name_dict(orth_transcripts)

    final_pairs = {}
    counter = 0
    #loop over the remaining genes
    for i in ids_to_keep:
        if counter%1000 == 0:
            print(counter)
        counter = counter + 1
        #get the IDs of the orthologous genes in the ortholog species
        orth_ids = pairs_dict[i]
        #get all the associated transcript identifiers
        orth_ids_trans = [[orth_gene_name_dict[j][k] for k in range(len(orth_gene_name_dict[j]))] for j in orth_ids if j in orth_gene_name_dict]
        orth_ids_trans = flatten(orth_ids_trans)
        CDS = CDS_seq[CDS_names.index(i)]
        orth_CDS = []
        ids_to_remove = []
        #get all the ortholog ORF sequences
        for j in orth_ids_trans:
            try:
                current_CDS = orth_CDS_seq[orth_CDS_names.index(j)]
                orth_CDS.append(current_CDS)
            #this is because some of the transcripts produced from the gene might be non-coding or have a wonky ORF and therefore not appear in the CDS fasta
            except ValueError:
                ids_to_remove.append(j)
        orth_ids_trans = [j for j in orth_ids_trans if j not in ids_to_remove]
        #check that the sequence from the focal species aligns to an ortholog with dN/dS below 0.5 and dS below the specified threshold
        if orth_ids_trans:
            conservation_check = keep_conserved_pc(i, orth_ids_trans, CDS, orth_CDS, dS_threshold, alignment_folder)
            if conservation_check[0]:
                #also store which ortholog transcript gave the lowest dS in the alignment
                final_pairs[i] = conservation_check[1]
            
    print("After filtering by conservation, {0} transcripts remain.".format(len(list(final_pairs.values()))))
    #write the final retained ortholog gene pairs to file
    with open(final_pairs_file_name, "w") as file:
        output_writer = csv.writer(file, delimiter = ",")
        for i in final_pairs:
            output_writer.writerow([i, final_pairs[i]])

    print("Wrote ortholog pairs to {0}.".format(final_pairs_file_name))

    #write the remaining ORF sequences to fasta
    CDS_seq = [i for pos, i in enumerate(CDS_seq) if CDS_names[pos] in final_pairs]
    CDS_names = [i for i in CDS_names if i in final_pairs]
    rw.write_to_fasta(CDS_names, CDS_seq, "general/filtered_{0}_wo_low_omega.fasta".format(dataset_name))

    #create a feature set with the remaining genes
    filtered_fs = Feature_Set(features_file_name, genome)
    filtered_fs.create_dataset("filtered_{0}".format(dataset_name), input_list = list(final_pairs.keys()))
    print("All done.")

if __name__ == "__main__":
    main()
