all: install help execute clean
	

install: main.cpp
	g++ -std=c++11 main.cpp fonctions.cpp hgt_int.cpp K-means.cpp Super-trees.cpp utils_tree.cpp -fopenmp -g -o KMPTC
	@if [ ! -d "../output" ]; then mkdir "../output"; fi
	@if [ ! -d "../output/stat.csv" ]; then touch "../output/stat.csv"; fi
	@if [ ! -d "../output/output.txt" ]; then touch "../output/output.txt"; fi
	chmod -R 777 * .
	@echo ' '
	@echo 'The program has been compiled correctly'
	@echo ' '
	
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
	@echo ' '
	
execute:
	@echo 'We are going to execute this line as example :'
	@echo '===> ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.2 3 8'
	@echo ' '
	$. ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.2 3 8
	@echo ' '

clean:
	$(RM) *.o input_.txt
	$(RM) KMPTC
	@echo ' '

.PHONY: pybind clean-py

# Build the pybind11 Python extension
PYEXT_NAME = KMPTC_core
PY_SRC = py_wrapper.cpp
PY_INCLUDES = $(shell python3 -m pybind11 --includes)
PY_SUFFIX = $(shell python3-config --extension-suffix)

pybind: $(PY_SRC)
	@echo 'Building Python extension $(PYEXT_NAME)'
	g++ -O3 -shared -std=c++11 -fPIC $(PY_INCLUDES) $(PY_SRC) -o $(PYEXT_NAME)$(PY_SUFFIX)

clean-py:
	$(RM) $(PYEXT_NAME)*$(PY_SUFFIX)
	@echo 'Removed Python extension files'
