Source code for paper on RBP motifs. housekeeping.py, nucleotide_comp.py, plotting.py, read_and_write.py, 
cython_func.pyx, my_stats.py, bedtools_games.py and conservation.py are modules and not all of the code
that they contain is relevant to the paper. The R scripts are necessary for the module my_stats.py.

RBP_motifs.py: obtaining the RBP motifs

overall_density.py: density parameters of the full set of RBP motifs.

overall_conservation.py: conservation parameters of the full set of RBP motifs.

dinucl_cons.py: motif rate of evolution, calculated separately for each dinucleotide.

RBP_density_all.py: density parameters of individual RBP motif sets, determined in full CDSs.

check_amount_of_information_RBP.py: filtering RBP motif sets based on the information content of hits.

reassign_bases_in_RBP.py: generating random motifs based on the hg38 mononucleotide composition.

RBP_motif_lengths.py: calculating the median motif length of each set of motifs.

RBP_conservation.py: conservation parameters of individual RBP motif sets, determined in exonic subregions.

one_away_conservation.py: selection against gaining avoided motifs.