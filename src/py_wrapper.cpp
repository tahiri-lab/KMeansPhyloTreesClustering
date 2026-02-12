#include <pybind11/pybind11.h>
#include <cstdlib>
#include <string>
#include <iostream>

namespace py = pybind11;

void executer_clustering(std::string path, int index, double alpha, int kmin, int kmax) {
    // Reconstruction de la commande exacte que tu as testée avec succès
    std::string commande = "./KMPTC -tree " + path + " " + std::to_string(index) + 
                           " " + std::to_string(alpha) + " " + std::to_string(kmin) + 
                           " " + std::to_string(kmax);
    
    std::cout << "[C++ Core] Exécution de la commande : " << commande << std::endl;
    std::system(commande.c_str());
}

PYBIND11_MODULE(KMPTC_core, m) {
    m.doc() = "Interface Python pour le moteur C++ KMPTC";
    m.def("executer_clustering", &executer_clustering, "Lance le calcul via le binaire compilé");
}
