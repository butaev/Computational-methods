README
========

**Computational methods for square matrix.**
---------------------------------------------

QRDecomposition(A, Q, R) is decomposition of a matrix A into a product A = QR of an orthogonal matrix Q and an upper triangular matrix R.

***Methods for finding eigenvalues of matrix to within epsilon:***

QRAlgorithm(A, eps) use qr decomposition.
BisectionMethod(A, eps) algorithm which uses the bisection method to find eigenvalues of matrix A. In this algorithm is used two main functions:
	NumLeftValues(a, A) return number of eigenvalues A less than a,
	SearchValues(a, b, eps,  A) Return eigenvalues A between a and b.


***Methods for solution system of linear equations of form Ax = b***

GMRES is generalized minimal residual method.
SimpleIteration is iterative method.
TridiagMatrix is tridiagonal matrix algorithm.
Gauss - Gauss method or Gaussian elimination.
GaussS - Gauss method with selection of the main element by column.
RotationMethod with help rotation matrices leads the matrix A to the upper triangular form and solve equation.