***DISCLAIMER***
The code in this repository has not been packaged to be run out of the box. For instance, there are many software dependencies that haven't been explicitly documented, and the documentation regarding how to run the scripts is also insufficiently documented for external use. The purpose of the repository is to serve as an explicit record of the analyses performed in the paper. It is thus primarily a supplement to the methods. However, if anybody did wish to reuse any of the code and was having issues or questions, please do contact me (rosinasavisaar@gmail.com) and we'll get it to work.
****************

Source code for the Savisaar and Hurst 2016 Molecular Biology and Evolution paper. housekeeping.py, nucleotide_comp.py, plotting.py, read_and_write.py, 
cython_func.pyx, my_stats.py, bedtools_games.py, ensembl_ops.py and conservation.py are custom modules that I've written for my personal use and not all of the code that they contain is relevant to the paper. 
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
