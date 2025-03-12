#include <iostream>
#include <vector>
#include <stdexcept>
#include <Eigen/Dense> 

// Exercice 1
// Question 1
Eigen::MatrixXd Vandermonde(const std::vector<double>& x) {
    int N = x.size();
    Eigen::MatrixXd Vander(N, N);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            Vander(i, j) = pow(x[i], j);
        }
    }
    return Vander;
}

// Question 1 - Alternative sans boucle
Eigen::MatrixXd matrice_vandermonde_mieux(const std::vector<double>& x) {
    Eigen::VectorXd X = Eigen::Map<const Eigen::VectorXd>(x.data(), x.size());
    int N = X.size();
    Eigen::MatrixXd V(N, N);
    for (int j = 0; j < N; ++j) {
        V.col(j) = X.array().pow(j);
    }
    return V;
}

// Question 2
Eigen::VectorXd solvand(const std::vector<double>& x, const std::vector<double>& y) {
    int N = x.size();
    int NN = y.size();
    if (N != NN) {
        throw std::invalid_argument("Probleme de taille entre les vecteurs");
    }
    Eigen::MatrixXd Vander = Vandermonde(x);
    Eigen::VectorXd y_eigen = Eigen::Map<const Eigen::VectorXd>(y.data(), y.size());
    Eigen::VectorXd a = Vander.colPivHouseholderQr().solve(y_eigen);
    return a;
}

// Question 3
int main() {
    std::vector<double> x = {0, 1, 2, 3};
    std::vector<double> y = {1, 2, 9, 28};

    Eigen::VectorXd rep = solvand(x, y);
    std::cout << "Les coefficients de la solution sont : " << rep.transpose() << std::endl;

    return 0;
}

// Question 4
int main() {
    // Question 4
    int size = 100;
    std::vector<double> x_random, y_random;

    generateRandomVectors(x_random, y_random, size);

    Eigen::VectorXd rep_random = solvand(x_random, y_random);
    std::cout << "Les coefficients de la solution pour les vecteurs aléatoires sont : " << rep_random.transpose() << std::endl;

    return 0;
