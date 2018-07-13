#include "Functions.h"
#include "Vector.h"
#include "Matrix.h"
#include <fstream>
#include <iostream>
#include <cmath>

void Functions::QRDecomposition(Matrix A, Matrix& q, Matrix& r)
{
    Matrix B(A.Size());
    Matrix Q(A.Size());
    Matrix a = A;

    Q.Identity();

    for (int i = 0; i < A.Size(); ++i)
    {
        for (int j = i + 1; j < A.Size(); ++j)
        {
            B.Identity();
            B[i][i] = a[i][i] / sqrt(a[i][i] * a[i][i] + a[j][i] * a[j][i]);//  cos
            B[i][j] = a[j][i] / sqrt(a[i][i] * a[i][i] + a[j][i] * a[j][i]);// sin
            B[j][i] = -B[i][j];
            B[j][j] = B[i][i];
            a = B * a;
        }
    }
    r = a;

    q = (a * A.Inverse()).Inverse();
    //std::cout<<r<<std::endl;
}

int Functions::NumLeftValues(double l, const Matrix& A)
{
	Vector d(A.Size());
	d[0] = A[0][0] - l;
	for (int i = 1; i < A.Size(); ++i)
		d[i] = (A[i][i] - l) - (A[i][i - 1] * A[i][i - 1]) / d[i - 1];
	int k = 0;
	for (int i = 0; i < A.Size(); ++i)
		if (d[i] < 0)
			++k;
	return k;
}

Vector Functions::SearchValues(double& a, double& b, double eps, const Matrix& A)
{
	Vector left, right;
	double l = (a + b) / 2.0;
	int la = NumLeftValues(a, A), lb = NumLeftValues(b, A), lc = NumLeftValues(l, A);
	Vector v(lb - la);
	
	if((b - a) < eps && (lb - la) == 1)
	{
		v[0] = l;
		return v;
	}
	
	if(lc - la == 0)
		right = SearchValues(l , b , eps, A);
	else if(lb - lc == 0)
		left = SearchValues(a, l, eps, A);
	else
	{
		left = SearchValues(a, l, eps, A);
		right = SearchValues(l , b , eps, A);
	}
	
	for(int i = 0; i < left.Size(); ++i)
		v[i] = left[i];
	for(int i = 0; i < right.Size(); ++i)
		v[i + left.Size()] = right[i];

	//std::cout<<v;

	return v;
}

Vector Functions::Gauss(Matrix A, Vector b, double& detA)
{
	double a, c;
	detA = 1.0;
	//std::cout<<A;
	int size = A.Size();
	for(int i = 0; i < A.Size(); ++i)
	{
		detA *= A[i][i];
		a = A[i][i];
		b[i] /= A[i][i];
		for(int k = i; k < size; ++k)
			A[i][k] = A[i][k] / a;

		for(int j = i + 1; j < size; ++j)
		{	
			c = A[j][i];
			b[j] -= b[i]*c;
			for(int k = i; k < size; ++k)
				A[j][k] -= A[i][k]*c;
		}
	}
	Vector x(size);
	x[size - 1] = b[size - 1];
	
	for(int i = size - 2; i >= 0; --i)
	{
		c = 0.0;
		for(int j = size - 1; j > i; --j)
			c -=  A[i][j]*x[j];
		x[i] = b[i] + c;
	}
	return x;
}

Vector Functions::GaussS(Matrix A, Vector b, double& detA)
{
	double a, c;
	detA = 1.0;

	//std::cout<<A;
	int size = A.Size();
	Vector x(size);
	for(int i = 0; i < A.Size(); ++i)
	{
		c = A[i][i];
		int j;
		for(int l = 0; l < size; ++l)
		{
			if(A[l][i] > c)
			{
				c = A[l][i];
				j = l;
			}
			else
				j = i;
		}
		if(i != j)
		{
			A.SwapLine(i,j);
			b.Swap(i,j);
			//std::cout<<A<<std::endl;
			detA *= -1;
		}
		if(A[i][i] != 0 && detA != 0.0)
			detA *= A[i][i];
		else
		{
			detA = 0.0;
			return x;
		}
		a = A[i][i];
		b[i] /= A[i][i];
		for(int k = i; k < size; ++k)
			A[i][k] = A[i][k] / a;

		for(int j = i + 1; j < size; ++j)
		{	
			c = A[j][i];
			b[j] -= b[i]*c;
			for(int k = i; k < size; ++k)
				A[j][k] -= A[i][k]*c;
		}
	}
	
	x[size - 1] = b[size - 1];
	
	for(int i = size - 2; i >= 0; --i)
	{
		c = 0.0;
		for(int j = size - 1; j > i; --j)
			c -=  A[i][j]*x[j];
		x[i] = b[i] + c;
	}
	return x;
}

Vector Functions::RotationMethod(Matrix A, Vector b)
{
	//Matrix a = A;
	Matrix B(A.Size());
	for(int i = 0; i < A.Size(); ++i)
	{
		int k = i;
		if (A[i][i] == 0)
		{
			for(int l = i + 1; l < A.Size(); ++l)
			{
				if(A[l][i] > 0)
					k = l;
				else
					if(l == A.Size() - 1)
					{
						Vector x(A.Size());
						return x;
					}
			}
		}
		if(i != k)
		{
			A.SwapLine(i, k);
			b.Swap(i, k);
		}

		for(int j = i + 1; j < A.Size(); ++j)
		{
			B.Identity();
			B[i][i] = A[i][i] / sqrt(A[i][i] * A[i][i] + A[j][i] * A[j][i]);//  cos
			B[i][j] = A[j][i] / sqrt(A[i][i] * A[i][i] + A[j][i] * A[j][i]);// sin
			B[j][i] = - B[i][j];
			B[j][j] = B[i][i];
			b = B * b;
			A = B * A;
		}
	}
	//std::cout<<A * A.Inverse()<<std::endl;
	Vector x(A.Size());
	x[A.Size() - 1] = b[A.Size() - 1] / A[A.Size() - 1][A.Size() - 1];
	double c;
	for(int i = A.Size() - 2; i >= 0; --i)
	{
		c = 0.0;
		for(int j = A.Size() - 1; j > i; --j)
			c -=  A[i][j]*x[j];
		x[i] = (b[i] + c) / A[i][i];
	}
	return x;
}

Vector Functions::TridiagMatrix(Matrix A, Vector f)
{
	int N = A.Size();
	Vector a(N);
	Vector b(N);
	Vector c(N);
	Vector alpha(N);
	Vector beta(N);
	Vector x(N);
	for(int i = 0; i < N; ++i)
	{
		c[i] = A[i][i];
		//std::cout<<c[i]<<'\n';
		if(i < A.Size() - 1)
		{
			b[i] = -A[i][i + 1];
			//std::cout<<b[i]<<'\n';
			a[i] = -A[i + 1][i];
			//std::cout<<a[i]<<'\n';
		}
	}
    alpha[0] = b[0] / c[0];
    beta[0] = f[0] / c[0];
	for(int i = 1; i < N; ++i)
	{
        alpha[i] = b[i] / (c[i] - alpha[i - 1] * a[i]);
        beta[i] = (f[i] + beta[i - 1] * a[i]) / (c[i] - alpha[i - 1] * a[i]);
	}

    x[N - 1] = beta[N - 1];

	for(int i = N - 2; i >= 0; --i)
	{
		x[i] = alpha[i] * x[i + 1] + beta[i];
	}
	return x;
}

Vector Functions::QRAlgorithm(Matrix a, double eps)
{
	int n = a.Size();
    Matrix q, r; 
    Matrix A = a;
	Vector v(n);
	while (true)
	{
		int s = 0;
		QRDecomposition(A, q, r);
		A = r * q;
		for(int j = 0; j < n - 1; ++j)
		{
			double Norm = 0.0;
			for(int i = j + 1; i < n; ++i)
				Norm += A[i][j] * A[i][j];
			if(sqrt(Norm) < eps)
			{
				++s;
				for(int k = 0; k < n; ++k)
					v[k] = A[k][k];
				if (s == n - 1)
					return v;
			}
		}
	}
}

Vector Functions::GMRES(Matrix A, Vector b, double eps, int& numIter)
{
	int N = A.Size();
	Vector x(N);

	for(int i = 0; i < N; ++i)
		x[i] = 1;

	Vector r = A * x - b;
    numIter = 0;
	while (r.Norm() > eps)
	{
		++numIter;
		Vector q = A * r;
		x = x  - ((q * r) / (q * q)) * r;
		r = A * x - b;
	}
	return x;
}

Vector Functions::SimpleIteration(Matrix A, Vector b, double eps, int& numIter)
 {
	int N = A.Size();
	Vector x(N);

	for(int i = 0; i < N; ++i)
		x[i] = 1;

	double t = 2 / (A.MinValue(eps) + A. MaxValue(eps));

	Vector r = A * x - b;
    numIter = 0;
	while (r.Norm() > eps)
	{
		++numIter;
		
		x = x  - t * r;
		r = A * x - b;
	}
	return x;}

Vector Functions::BisectionMethod(Matrix A, double eps)
{
	/*
    Matrix E(A.Size()), X(A.Size()), H(A.Size());
	Vector v(A.Size()), e(A.Size()), x(A.Size());
	double t;
	 
	E.Identity();

	for(int i = 0; i < A.Size() - 1; ++i)
	{
		Vector b(A.Size());
		e[i + 1] = 1;
		if(i > 0)
			e[i] = 0;
		for(int l = 0; l < A.Size(); ++l)
		{
			 if(l >= i + 1)
				 b[l] = A[l][i];
		}

		t = (b[i + 1] > 0 ? 1 : -1) * b.Norm();
		v = b + t * e;

		for(int j = 0; j < A.Size(); ++j)
			for(int k = 0; k < A.Size(); ++k)
				X[j][k] = v[j] * v[k];

		H = E - (2 / (v * v)) * X;
		A = H * A * H;
	}
    */
	//std::cout<<"A\n"<<A<<"\n";

    Vector bounds(2);
	bounds = A.ValuesBounds();
	//std::cout<<"bounds[0] = "<< bounds[0] <<std::endl;
    //std::cout << "bounds[1] = " << bounds[1] << std::endl;

	return Functions::SearchValues(bounds[0], bounds[1], eps, A);
}