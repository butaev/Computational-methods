#pragma once

#include "Vector.h"
#include "Matrix.h"
#include <fstream>

class Functions
{
public:
    static void QRDecomposition(Matrix a, Matrix& q, Matrix& r);

    // methods for finding eigenvalues of matrix
    static Vector QRAlgorithm(Matrix A, double eps);
    static Vector BisectionMethod(Matrix A, double eps);

    // methods for solution system of linear equations
    static Vector RotationMethod(Matrix A, Vector b);
    static Vector GMRES(Matrix A, Vector b, double eps, int& numIter);// generalized minimal residual method
    static Vector SimpleIteration(Matrix A, Vector b, double eps, int& numIter);
    
    static Vector TridiagMatrix(Matrix A, Vector f);
    static Vector Gauss(Matrix A, Vector b, double& detA);
    static Vector GaussS(Matrix A, Vector b, double& detA);

    static int NumLeftValues(double a, const Matrix& A);// number of eigenvalues A less than a
    static Vector SearchValues(double& a, double& b, double eps, const Matrix& A);//finding eigenvalues between a and b
};