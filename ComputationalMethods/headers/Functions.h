#include "Vector.h"
#include "Matrix.h"
#include <fstream>

#pragma once

class Functions
{
private:
	static int GetLeftValue(double& a, const Matrix& A);
	static Vector Serch(double& a, double& b, double eps, const Matrix& A);
public:
	static Vector Gauss(Matrix A, Vector b, double& detA);
	static Vector GaussS(Matrix A, Vector b, double& detA);
	static Vector RotationMethod(Matrix A, Vector b);
	static Vector SweepMethod(Matrix A, Vector f);
	static Vector GetValues(Matrix A, double eps);
	static Vector MinResidual(Matrix A, Vector b, double eps, int& k);
	static Vector SimpleIteration(Matrix A, Vector b, double eps, int& k);
	static Vector BisectionMethod(Matrix A, double eps);
};

