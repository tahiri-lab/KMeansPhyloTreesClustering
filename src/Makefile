all: install help
	
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
	@echo ' '
	@echo 'You also can enter step by step your argyment, for that, please execute this line :'
	@echo '===> ./KMPTC'
	
execute:
	@echo 'We gonna execute this line as example :'
	@echo '===> ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.2 3 8'
	$. ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.2 3 8

clean:
	$(RM) *.o input_.txt
	$(RM) KMPTC

