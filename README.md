README
========

**Computational methods for square matrix.**
---------------------------------------------

***Method of decomposition matrix:***

QRDecomposition is decomposition of a matrix A into a product A = QR of an orthogonal matrix Q and an upper triangular matrix R.

***Methods for finding eigenvalues of matrix to within epsilon:***

QRAlgorithm use qr decomposition for solve problem.

BisectionMethod algorithm which uses the bisection method to find eigenvalues of matrix A. In this algorithm is used two main functions:
	NumLeftValues(a, A) return number of eigenvalues A less than a,
	SearchValues(a, b, eps,  A) Return eigenvalues A between a and b.


***Methods for solution system of linear equations of form Ax = b:***

GMRES is generalized minimal residual method.

SimpleIteration is iterative method.

TridiagMatrix is tridiagonal matrix algorithm.

Gauss - Gauss method or Gaussian elimination.

GaussS - Gauss method with selection of the main element by column.

RotationMethod with help rotation matrices leads the matrix A to the upper triangular form and solve equation.

For matrices, some additional functions are implemented, such as:
matrix norms (spectral and 1 - norm), matrix inversion, calculation of maximum and minimum eigenvalues, calculation of condition number and boundary of the eigenvalues of the matrix.