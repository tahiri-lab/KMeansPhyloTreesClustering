#include <pybind11/pybind11.h>
#include <cstdlib>
#include <string>
#include <iostream>

namespace py = pybind11;

/**
 * Lance le clustering K-means sur un fichier d'arbres phylogénétiques
 * en appelant le binaire compilé KMPTC via la ligne de commande.
 *
 * Construit et exécute la commande :
 *   ./KMPTC -tree <path> <index> <alpha> <kmin> <kmax>
 *
 * Paramètres :
 *   path  - chemin vers le fichier d'arbres au format Newick
 *   index - indice de validité de cluster (1 = CH, 2 = BH)
 *   alpha - pénalisation du chevauchement des espèces (entre 0 et 1)
 *   kmin  - nombre minimal de clusters
 *   kmax  - nombre maximal de clusters
 */
void executer_clustering(std::string path, int index, double alpha, int kmin, int kmax) {
    std::string commande = "./KMPTC -tree " + path + " " + std::to_string(index) +
                           " " + std::to_string(alpha) + " " + std::to_string(kmin) +
                           " " + std::to_string(kmax);

    std::cout << "[C++ Core] Exécution de la commande : " << commande << std::endl;
    int ret = std::system(commande.c_str());
    (void)ret;
}

PYBIND11_MODULE(KMPTC_core, m) {
    m.doc() = "Interface Python pour le moteur C++ KMPTC — clustering d'arbres phylogénétiques par K-means";
    m.def("executer_clustering", &executer_clustering,
          "Lance le clustering K-means via le binaire KMPTC compilé",
          py::arg("path"), py::arg("index"), py::arg("alpha"),
          py::arg("kmin"), py::arg("kmax"));
}
