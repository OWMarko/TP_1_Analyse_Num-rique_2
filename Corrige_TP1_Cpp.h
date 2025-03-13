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

    RandomVectors(x_random, y_random, size);

    Eigen::VectorXd rep_random = solvand(x_random, y_random);
    std::cout << "Les coefficients sont : " << rep_random.transpose() << std::endl;

    return 0;


//Exercice 2 
//Question 1

double evaluerPolynomial(const std::vector<double>& coefficients, double z) {
    double result = 0.0;
    for (int i = coefficients.size() - 1; i >= 0; --i) {
        result = result * z + coefficients[i];
    }
    return result;
}

//Question 2
double evaluerPolynomialSansFor(const std::vector<double>& coefficients, double z, int n) {
    if (n == 0) {
        return coefficients[0];
    }
    return coefficients[n] + z * evaluerPolynomialSansFor(coefficients, z, n - 1);

//Question 3
double horner(const vector<double>& coeffs, double x) {
    double result = coeffs[0];
    for (size_t i = 1; i < coeffs.size(); ++i) {
        result = result * x + coeffs[i];
    }
    return result;
}

//Exercice 3
vector<vector<double>> diffdiv(const vector<double>& x, const vector<double>& y) {
    int n = x.size();
    vector<vector<double>> div_diff(n, vector<double>(n));
    for (int i = 0; i < n; ++i) {
        div_diff[i][0] = y[i];
    }
    for (int j = 1; j < n; ++j) {
        for (int i = 0; i < n - j; ++i) {
            div_diff[i][j] = (div_diff[i + 1][j - 1] - div_diff[i][j - 1]) / (x[i + j] - x[i]);
        }
    }
    return div_diff;


