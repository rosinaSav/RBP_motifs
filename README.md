Source code for paper on RBP motifs. housekeeping.py, nucleotide_comp.py, plotting.py, read_and_write.py, 
cython_func.pyx, my_stats.py, bedtools_games.py, ensembl_ops.py and conservation.py are modules and not all of the code
that they contain is relevant to the paper. 
The R scripts are necessary for the module my_stats.py. 
The perl script are necessary for the module ensembl_ops.py.

check_amount_of_information_RBP.py: filtering RBP motif sets based on the information content of hits.

dinucl_cons.py: motif rate of evolution, calculated separately for each dinucleotide.

get_clean_pc_sequences.py: preparing a clean set of protein-coding sequences for a given species. This script was used for generating the mouse dataset. The human dataset had been assembled for a previous publication using different code, however, the filtering steps were the same.

get_multiexon_and_families.py: given a dataset of clean ORFs, making a subdataset that only contains ORFs from intron-containing genes. Also clustering into paralogous families.

get_UTR_and_intron_sequence.py: obtaining pairwise alignments of UTR and intronic sequences.

map_FANTOM.py: mapping FANTOM5 CAGE peaks to gene promoters and calculating a series of expression parameters for each gene.

one_away_conservation.py: selection against gaining avoided motifs.

overall_conservation.py: conservation parameters of the full set of RBP motifs.

overall_density.py: density parameters of the full set of RBP motifs.

prepare_FANTOM.py: filtering a FANTOM5 osc file and writing data to BED format for CrossMapping.

RBP_conservation.py: conservation parameters of individual RBP motif sets, determined in exonic subregions.

RBP_density_all.py: density parameters of individual RBP motif sets, determined in full CDSs.

RBP_motif_lengths.py: calculating the median motif length of each set of motifs.

RBP_motifs.py: obtaining the RBP motifs

reassign_bases_in_RBP.py: generating random motifs based on the hg38 mononucleotide composition.
