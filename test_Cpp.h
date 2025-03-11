#include <iostream>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense>

Matrix3f m;
m << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;
std::cout << m;

std::cout << "Matrice carré avec que des zeros:\n";
Array33f a1 = Array33f::Zero();
std::cout << a1 << "\n\n";
 
 
std::cout << "Vecteur colonne\n";
ArrayXf a2 = ArrayXf::Zero(3);
std::cout << a2 << "\n\n";
 
 
std::cout << "Matrice 3x4 avec que des zeros:\n";
ArrayXXf a3 = ArrayXXf::Zero(3, 4);
std::cout << a3 << "\n";
