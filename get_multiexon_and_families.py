from bedtools_games import Feature_Set
import conservation
from housekeeping import parse_arguments
import read_and_write as rw

def main():
    description = "Pick out the multi-exon genes from a dataset and generate families."
    args = parse_arguments(description, ["features_file", "genome", "dataset", "fasta"])
    [features_file, genome, dataset, fasta] = [args.features_file, args.genome, args.dataset, args.fasta]

    fs = Feature_Set(features_file, genome)
    fs.set_dataset(dataset)
    exons = fs.get_exons()
    exon_numbers = fs.get_exon_numbers(exons)

    output_fasta_name = "{0}_multiexon.fasta".format(fasta[:-6])
    
    multi_exon = [i for i in exon_numbers if exon_numbers[i] > 1]

    fs_new = Feature_Set(features_file, genome)
    fs_new.create_dataset("{0}_multiexon".format(dataset), input_list = multi_exon)
    fs_new.set_dataset("{0}_multiexon".format(dataset))
    
    names, seqs = rw.read_fasta(fasta)
    seqs = [seqs[pos] for pos, i in enumerate(names) if i in multi_exon]
    names = [i for i in names if i in multi_exon]
    rw.write_to_fasta(names, seqs, output_fasta_name)

    transcripts = fs_new.get_transcripts()
    gene_name_dict = fs_new.get_gene_name_dict(transcripts)

    conservation.find_families(output_fasta_name, "general/{0}_multiexon".format(dataset))
        

if __name__ == "__main__":
    main()
