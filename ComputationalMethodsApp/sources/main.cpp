#include "Functions.h"
#include "Vector.h"
#include "Matrix.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

int main ()
{
    std::ifstream F;
    F.open("../Data/input.txt");
    if (!F.is_open())
    {
        std::cout << "File is not exist" << std::endl;
        return -1;
    }

    Matrix A;
    F>>A;

    std::cout << "A =" << std::endl;
    std::cout << A << std::endl;

    Vector b;
    F>>b;
    F.close();

    std::cout << "b =" << std::endl;
    std::cout << b << std::endl;

    double detA = 0.0;

    auto x = Functions::GaussS(A, b, detA);

    std::cout << "x =" << std::endl;
    std::cout << x << std::endl;

    std::cout << "detA = " << detA << std::endl;

    Matrix Q, R;
    double eps = 1e-6;
    Functions::QRDecomposition(A, Q, R);

    std::cout << "Q =" << std::endl;
    std::cout << Q << std::endl;

    std::cout << "R =" << std::endl;
    std::cout << R << std::endl;

    std::cout << "Result of BisectionMethod:" << std::endl;
    std::cout << Functions::BisectionMethod(A, eps) << std::endl;

    system("pause");
    return 0;
}