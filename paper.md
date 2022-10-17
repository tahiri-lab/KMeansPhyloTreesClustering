---

title: ''
tags:
  - phylogeny
  - 
authors:
  - name: Benjamin Albertelli
    orcid: 
  - name: Nadia Tahiri
affiliations:
  - name:
    index:
date: 17 October 2022
bibliography: paper.bib

---

# Summary 

KMeansPhylogeneticTreesClustering (KMPTC) is part of a research work in the field of phylogeny. Phylogeny is the study of the relationships between organisms. The goals of phylogeny are, among others, to study the evolution of living beings and to estimate the time of divergence since their last common ancestor. The KMPTC package allows users to cluster a set of phylogenetic trees into an optimal partition taken as input. The algorithm is fast and efficient, using the k-means partitioning algorithm and leaves the opportunity to the user to choose many criteria in order to have a complete control on the use of over the use of the algorithm. The algorithm has the vocation to evolve to propose new metrics for comparison between trees and an ever-increasing optimization. The repository has a properly detailed readme allowing everyone to install the package on his machine. Moreover, we gave a clear and documented example to allow the user to use it from the first use. We created a computational framework to study the theory of evolutionary history of genes to better understand the impact of biological modifications on genomes. We implemented an interactive console version (version 1.0) for each of our new algorithms. These console versions allow researchers to easily integrate them into their scientific workflows, e.g., easily run them on High Performance Computers (HPC) (e.g., Compute Canada and Compute Quebec clusters) to use the power of the servers and the parallelism of the threads.
 
Thereafter, we will develop user-friendly versions of our programs to make them essential tools for analysis of evolutionary data. All algorithms will be implemented and released as free and open-source software to make the conducted analyses reproducible and reliable. 

# Statement of need

KMeansPhylogeneticTreesClustering is a package that aims to bring together in a single algorithm all the same algorithm all the necessary tools to cluster phylogenetic trees. This algorithm can be used by everyone thanks to its installation and execution which detailed step by step (cf make help) but it aims above all to help the scientific community by scientific community by offering them a complete, efficient and malleable algorithm. This algorithm comes at a time when the world objective is to create the tree of life, the tree of of all bipartitioned living beings. For this, the steps of a phylogenetic analysis of molecular data consist in obtaining gene sequences (genomic extraction => PCR => sequencing), to align them, to calculate the matrix distances or similarities between species (depending on the inference method chosen) and apply a tree reconstruction method. Once the phylogenetic tree is obtained, it is necessary to be able to compare the trees since the different tree reconstruction methods give different results. Moreover, depending on the information chosen to perform the analysis, the phylogenies differ, because each gene has its own evolutionary history. It is therefore at this last step of comparison that the algorithm intervenes. There is currently no other tool as complete and available in open source for free for everyone. The algorithm therefore needs to be efficient in order to help the scientific community scientific community to accomplish the objective. 

# State of the field



# Acknowledgements

The authors would like to thank the Department of Computer Science, University of Sherbrooke, Quebec, Canada for providing the necessary resources to conduct this research. We also thank Compute Canada for providing access to high-performance computing facilities. This work was supported by a grant by NSERC discovery grant. 

# References
