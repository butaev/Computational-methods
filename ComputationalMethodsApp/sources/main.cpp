#include "Functions.h"
#include "Vector.h"
#include "Matrix.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

int main ()
{
    std::ifstream f;
    f.open("../Data/input.txt");
    if (!f.is_open())
    {
        std::cout << "File is not exist" << std::endl;
        return -1;
    }

    Matrix m;
    f>>m;

    std::cout << "m =" << std::endl;
    std::cout << m << std::endl;

    Vector v;
    f>>v;
    f.close();

    std::cout << "v =" << std::endl;
    std::cout << v << std::endl;

    double detm = 0.0;

    auto result = Functions::GaussS(m, v, detm);

    std::cout << "Result of Gauss method:" << std::endl;
    std::cout << result << std::endl;

    std::cout << "det(m) = " << detm << std::endl;

    Matrix orthogonalMatrix, upperTriangular;
    double eps = 1e-6;
    Functions::QRDecomposition(m, orthogonalMatrix, upperTriangular);

    std::cout << "Orthogonal matrix" << std::endl;
    std::cout << orthogonalMatrix  << std::endl;

    std::cout << "Upper triangular" << std::endl;
    std::cout << upperTriangular  << std::endl;

    std::cout << "Result of BisectionMethod:" << std::endl;
    std::cout << Functions::BisectionMethod(m, eps) << std::endl;

    system("pause");
    return 0;
}