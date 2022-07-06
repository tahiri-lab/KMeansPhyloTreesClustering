all: main.cpp
	g++ main.cpp fonctions.cpp hgt_int.cpp K-means.cpp Super-trees.cpp utils_tree.cpp -openmp -g -o KMPTC
	chmod -R 777 * .
	@echo '=> ============================================================================================================================='
	@echo '=> Program   : KMeansSuperTreeClustering - 2021'
	@echo '=> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Québec a Montreal)'
	@echo '=> This program clusters phylogenetic trees using the k-means partitioning algorithm.'
	@echo '=> These trees may have the same or different, but mutually overlapping, sets of leaves (the multiple supertree problem).'
	@echo '=> Phylogenetic trees must be given in the Newick format (program input).'
	@echo '=> A partitioning of the input trees in K clusters of trees is returned as output. '
	@echo '=> The optimal number of clusters can be determined either by the Calinski-Harabasz (CH) or by the Ball-Hall (BH) cluster '
	@echo '=> validity index adapted for tree clustering.'
	@echo '=> A supertree can then be inferred for each cluster of trees.'
	@echo '=> The Robinson and Foulds topological distance is used in the objective function of K-means.'
	@echo '=> The list of the program parameters is specified below.'
	@echo '=> ============================================================================================================================='
	@echo ' '
	@echo 'Please execute the following command lines (on two steps):'
	@echo '===> 1) make install'
	@echo '===> 2) For trees: ./KMPTC -tree input_file cluster_validity_index alpha Kmin Kmax'
	@echo ' '
	@echo 'input_file: the input file for the program'
	@echo 'cluster_validity_index: the cluster validity index used in K-means (1 for Calinski-Harabasz and 2 for Ball-Hall)'
	@echo 'alpha: the penalty parameter for species overlap in phylogenetic trees (must be between 0 and 1)'
	@echo 'Kmin: the minimum number of clusters in K-means.'
	@echo '            - For CH, Kmin>=2,'
	@echo '            - For BH, Kmin>=1.'
	@echo 'Kmax: the maximum number of clusters in K-means.'
	@echo '            - Kmax must be less or equal to N-1 (where N is the number of input trees).'
	@echo ''
	@echo '===> Command line execution example (input_file = data/Covid-19_trees.txt, cluster_validity_index = CH, alpha = 0.2, Kmin = 3, Kmax = 8):'
	@echo '===> ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.2 3 8'
	
install: main.cpp
	g++ main.cpp fonctions.cpp hgt_int.cpp K-means.cpp Super-trees.cpp utils_tree.cpp -openmp -g -o KMPTC
	chmod -R 777 * .
	
help:
	@echo 'Please execute the following command lines (on two steps):'
	@echo '===> 1) make install'
	@echo '===> 2) For trees: ./KMPTC -tree input_file cluster_validity_index alpha Kmin Kmax'
	@echo ' '
	@echo 'input_file: the input file for the program'
	@echo 'cluster_validity_index: the cluster validity index used in K-means (1 for Calinski-Harabasz and 2 for Ball-Hall)'
	@echo 'alpha: the penalty parameter for species overlap in phylogenetic trees (must be between 0 and 1)'
	@echo 'Kmin: the minimum number of clusters in K-means.'
	@echo '            - For CH, Kmin>=2,'
	@echo '            - For BH, Kmin>=1.'
	@echo 'Kmax: the maximum number of clusters in K-means.'
	@echo '            - Kmax must be less or equal to N-1 (where N is the number of input trees).'
	@echo ''
	@echo '===> Command line execution example (input_file = data/Covid-19_trees.txt, cluster_validity_index = CH, alpha = 0.2, Kmin = 3, Kmax = 8):'
	@echo '===> ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.2 3 8'
	
clean:
	$(RM) *.o input_.txt
	$(RM) KMPTC

