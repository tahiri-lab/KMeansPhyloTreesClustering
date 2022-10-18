# KMeansPhylogeneticTreesClustering
A new fast method for clustering phylogenetic trees using K-means and inferring multiple supertrees.

# About
**Program**   : KMeansPhylogeneticTreesClustering - 2022  
**Authors**   : Benjamin Albertelli and Nadia Tahiri (University of Sherbrooke)  

This program clusters phylogenetic trees using the k-means partitioning algorithm.  
These trees may have the same (the multiple consensus tree problem) or different, but mutually overlapping, sets of leaves (the multiple supertree problem).
    
Phylogenetic trees must be given in the Newick format (program input). A partitioning of the input trees in K clusters of trees is returned as output. 
The optimal number of clusters can be determined either:

1) by the Calinski-Harabasz (CH) or 
2) by the Ball-Hall (BH) cluster validity index adapted for tree clustering.

A supertree can then be inferred for each cluster of trees.The Robinson and Foulds topological distance is used in the objective function of K-means.
The list of the program parameters is specified below.

# Requirements
This package works with : 

- boost 1.78.0_1 for boost regex library (boost has to be added to your PATH)
- git 2.35.1
- macOS Monterey version 12.5
- macOS Terminal version 2.12.7

Warning : that's not the minimal requirement, the software should work with previous versions of software.

# How to install
First copy and paste this link in your shell after you go to the folder where you want to download the software :

    $ git clone https://github.com/tahiri-lab/KMeansPhylogeneticTreesClustering.git

Then please execute these two lines.

    $ cd src
    $ make install


# Help
If you need help for execution, please execute this line in the src folder:

    $ make help

# Examples

Please execute the following command line:  

=> For trees: ./KMPTC -tree input_file cluster_validity_index α Kmin Kmax

=> input_file: the input file for the program  
=> cluster_validity_index: the cluster validity index used in K-means (1 for Calinski-Harabasz and 2 for Ball-Hall)  
=> α: is the penalty parameter for species overlap in phylogenetic trees (must be between 0 and 1)  
=> Kmin: is the minimum number of clusters in K-means.  
- For CH, Kmin>=2,  
- For BH, Kmin>=1.   
  
=> Kmax: the maximum number of clusters in K-means.  
- Kmax must be less or equal to N-1 (where N is the number of input trees).

Command line execution examples:

input_file = data/Covid-19_trees.txt, cluster_validity_index = CH, α = 0.1, Kmin = 3, Kmax = 8):

    $ ./KMPTC -tree ../data/Covid-19_trees.txt 1 0.1 3 8
    
Or by using the Makefile instruction as follows (as in the previous example):

    $ make execute
        
input_file = data/all_trees_woese.txt, cluster_validity_index = CH, α = 1, Kmin = 2, Kmax = 10):

    $ ./KMPTC -tree ../data/all_trees_woese.txt 1 1 2 10

# Input
The input data sets are located in the folder "data".  
You can also use your own data, please ensure that the file respect the needed format.

# Output
See the folder "output"  
The output is in the following files:  
1) stat.csv - for the clustering statistics.  
2) output.txt - for the cluster content.

# Clean Project

To clean the project, please execute:

    $ make clean
