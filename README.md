# mimobaits
hybrid capture in mimosoid/Caesalpinioid legumes

The file mimobaits_phylotargets_Albizia964_cds.fasta contains coding sequences of the 964 low-copy nuclear genes that are described by Koenen et al. (submitted) "Hybrid capture of 964 nuclear genes resolves generic-level evolutionary relationships in the mimosoid clade (Leguminosae-Caesalpinioideae) and reveals a putative hard polytomy that spawned a large pantropical radiation". American Journal of Botany.

The file mimobaits_v1_AlbiziaInga1387targets_cds.fasta contains 1387 genes based on Albizia and Inga transcriptomes, which includes both the 964 low-copy nuclear genes mentioned above and the genes used in the study on the genus Inga by Nicholls et al. (2015) "Using targeted enrichment of nuclear genes to increase phylogenetic resolution in the neotropical rain forest genus Inga (Leguminosae: Mimosoideae)" Frontiers in Plant Science https://doi.org/10.3389/fpls.2015.00710, as well as some additional "functional" genes that have not been analyzed yet. AlbiziaInga1387_Baits_mim_Tpad2-120-40.fas is the corresponding probe sequences file for this set of targets.

The file mimobaits_v2_Albizia997targets_cds.fasta contains 997 genes included in a redesign of the baits that we are currently using in the Hughes lab at the University of Zurich. Albizia997_baits_ultrastringent_180816je.fas is the corresponding probe sequences file for this set of targets.

The python script RBH4.py (and the accompanying fasta_stuff_rbh4.py) was used to gather the genes with a reciprocal best hit algorithm for the 4 taxa to select the genes to target, including building alignments and small 4-taxon trees from these. Requirements are to have blast, mafft and raxml installed and in your path, as well as BMGE but the path to BMGE will have to be changed in the RBH4.py script. After running this script with 4 different percentage identity cut-offs in cd-hit (90, 95, 97 and 99%), the script select_rbh.sh was used to select the putative orthologs with the appropriate level of divergence among the taxa for each of the cd-hit cut-offs.

More info soon.
