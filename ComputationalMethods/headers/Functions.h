#pragma once

#include "Vector.h"
#include "Matrix.h"

namespace Functions
{
    void QRDecomposition(Matrix m, Matrix& q, Matrix& r);

    // methods for finding eigenvalues of matrix
    Vector QRAlgorithm(Matrix m, double eps);
    Vector BisectionMethod(Matrix m, double eps);

    // methods for solution system of linear equations
    Vector GMRES(Matrix m, Vector v, double eps, int& numIter);// generalized minimal residual method
    Vector SimpleIteration(Matrix m, Vector v, double eps, int& numIter);
    
    Vector TridiagMatrix(Matrix m, Vector v);
    Vector Gauss(Matrix m, Vector v, double& detm);
    Vector GaussS(Matrix m, Vector v, double& detm);
    Vector RotationMethod(Matrix m, Vector v);

    int NumLeftValues(double boundary, const Matrix& m);// number of eigenvalues A less than a
    Vector SearchValues(double& leftBound, double& rightBound, double eps, const Matrix& m);//finding eigenvalues between a and b, like binary search
};