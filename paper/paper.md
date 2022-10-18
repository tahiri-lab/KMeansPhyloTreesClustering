---
title: 'KMPTC: a new tool to classify phylogenetic trees with or without the same taxon sets'
tags:
  - bioinformatics
  - phylogeny
  - supertree
  - consensus tree
  - classification
  - clustering
  - Robinson and Foulds distance
authors:
  - name: Benjamin Albertelli
    affiliation: "1, 2"
  - name: Nadia Tahiri
    affiliation: 1
    orcid: 0000-0002-1818-208X
    corresponding: true
    email: Nadia.Tahiri@USherbrooke.ca
affiliations:
  - name: département d’Informatique, Université de Sherbrooke, 2500 Boulevard de l’Université, Sherbrooke, Québec J1K 2R1, Canada
    index: 1
  - name: ENSEA, 6 avenue du Ponceau 95000 Cergy, France
    index: 2
date: 18 October 2022
bibliography: paper.bib
---

# Summary 

KMeansPhylogeneticTreesClustering (KMPTC) [@tahiri2022building] is part of a research work in the field of phylogeny. Phylogeny is the study of the relationships between organisms. The goals of phylogeny are, among others, to study the evolution of living beings and to estimate the time of divergence since their last common ancestor. The KMPTC package allows users to cluster a set of phylogenetic trees into an optimal partition taken as input. The algorithm is fast and efficient, using the k-means partitioning algorithm and leaves the opportunity to the user to choose many criteria in order to have a complete control over the use of the algorithm. The algorithm has the vocation to evolve to propose new metrics for comparison between trees and an ever-increasing optimization. The repository has a properly detailed readme allowing everyone to install the package on his machine. Moreover, we gave a clear and documented example to allow the user to use it from the first use. We created a computational framework to study the theory of evolutionary history of genes to better understand the impact of biological modifications on genomes. We implemented an interactive console version (version 1.0) for each of our new algorithms. These console versions allow researchers to easily integrate them into their scientific workflows, e.g., easily run them on High Performance Computers (HPC) (e.g., Compute Canada and Compute Quebec clusters) to use the power of the servers and the parallelism of the threads.
 
Thereafter, we will develop user-friendly versions of our programs to make them essential tools for analysis of evolutionary data. All algorithms will be implemented and released as free and open-source software to make the conducted analyses reproducible and reliable. 

# Statement of need

KMeansPhylogeneticTreesClustering [@tahiri2022IFCS] is a package that aims to bring together in a single algorithm all the necessary tools to cluster phylogenetic trees. This algorithm can be used by everyone thanks to its installation and execution which detailed step by step (cf make help) but it aims above all to help the scientific community by scientific community by offering them a complete, efficient and malleable algorithm. This algorithm comes at a time when the world objective is to create the tree of life, the tree of of all bipartitioned living beings. For this, the steps of a phylogenetic analysis of molecular data consist in obtaining gene sequences (genomic extraction => PCR => sequencing), to align them, to calculate the matrix distances or similarities between species (depending on the inference method chosen) and apply a tree reconstruction method. Once the phylogenetic tree is obtained, it is necessary to be able to compare the trees since the different tree reconstruction methods give different results. Moreover, depending on the information chosen to perform the analysis, the phylogenies differ, because each gene has its own evolutionary history. It is therefore at this last step of comparison that the algorithm intervenes. There is currently no other tool as complete and available in open source for free for everyone. The algorithm therefore needs to be efficient in order to help the scientific community scientific community to accomplish the objective [@makarenkov2022]. 

# State of the field

Our research team has developed a method for inferring multiple consensus trees using k-medoids [@tahiri2018new]. Then, recently, our research team contributed to extend the classification of phylogenetic trees using convolutional neural networks [@tahiri2022invariant].

The main functionality of the algorithm allows the user to group trees by allowing him to define each of the input parameters. Let's take a concrete case:  
 
- First, the user has to download and install the package on his machine. All these steps are detailed in the README.  
  
- Then, the user has two possibilities.  
Firstly, he can run examples in order to familiarize himself with the use of the software. For that, with the help of the readme and its "example" part, it will be very easy for the user to execute an example.  
Otherwise, he can choose to simply execute the program by choosing his criteria. Note that the file containing the trees to be grouped must be given in Newick format (program input). The user can then modify the validity index of the clusters used in the K-means, the penalty parameter for the overlap of species in the phylogenetic trees, the minimum and maximum number of clusters in the K-means. 

To do this the user must run the command line ./KMPTC in the src folder. After that the algorithm will allow him to choose his own criteria.
  
- Finally, after execution, the user has access to the clustering statistics as well as to the distribution of the trees within the clusters. To do this, he must go to the files in the "output" folder. He can then find the partition found according to the criteria as well as the clustering score. 

# Acknowledgements

The authors would like to thank the Department of Computer Science, University of Sherbrooke, Quebec, Canada for providing the necessary resources to conduct this research. We also thank Compute Canada for providing access to high-performance computing facilities. This work was supported by a grant by NSERC discovery grant. 

# References
