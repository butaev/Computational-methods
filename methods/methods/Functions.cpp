#include "Functions.h"
#include "Vector.h"
#include "Matrix.h"
#include <fstream>
#include <iostream>
#include <cmath>

int Functions::GetLeftValue(double& l, const Matrix& A)
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

Vector Functions::Serch(double& a, double& b, double eps, const Matrix& A)
{
	Vector left, right;
	double l = (a + b) / 2.0;
	int la = GetLeftValue(a, A), lb = GetLeftValue(b, A), lc = GetLeftValue(l, A);
	Vector v(lb - la);
	
	if((b - a) < eps && (lb - la) == 1)
	{
		v[0] = l;
		return v;
	}
	
	// la <=> left <=> lc <=> right <=> lb
	if(lc - la == 0)
		right = Serch(l , b , eps, A);
	else if(lb - lc == 0)
		left = Serch(a, l, eps, A);
	else
	{
		left = Serch(a, l, eps, A);
		right = Serch(l , b , eps, A);
	}
	
	for(int i = 0; i < left.Length(); ++i)
		v[i] = left[i];
	for(int i = 0; i < right.Length(); ++i)
		v[i + left.Length()] = right[i];

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

Vector Functions::SweepMethod(Matrix A, Vector f)
{
	int N = A.Size();
	Vector a(N - 1);
	Vector b(N - 1);
	Vector c(N);
	Vector d(N - 1); // �����
	Vector e(N - 1); // ����
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
	d[0] = b[0] / c[0];
	e[0] = f[0] / c[0];
	for(int i = 1; i < N - 1; ++i)
	{
		d[i] = b[i - 1] / (c[i - 1] - d[i - 1] * a[i - 1]);
		//std::cout<<d[i]<<'\n';
		e[i] = (f[i] + e[i - 1] * a[i - 1]) / (c[i - 1] - d[i - 1] * a[i - 1]);
		//std::cout<<e[i]<<'\n';
	}

	x[N - 1] = (f[N - 1] + e[N - 2] * a[N - 2]) / (c[N - 2] - d[N - 2] * a[N - 2]);

	for(int i = N - 2; i >= 0; --i)
	{
		x[i] = d[i] * x[i + 1] + e[i];
	}
	return x;
}

Vector Functions::GetValues(Matrix a, double eps)
{
	int n = a.Size();
	Matrix q, r, A = a;
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

Vector Functions::MinResidual(Matrix A, Vector b, double eps, int& k)
{
	int N = A.Size();
	Vector x(N);

	for(int i = 0; i < N; ++i)
		x[i] = 1;

	Vector r = A * x - b;
	k = 0;
	while (r.Norm() > eps)
	{
		++k;
		Vector q = A * r;
		x = x  - ((q * r) / (q * q)) * r;
		r = A * x - b;
	}
	return x;
}

Vector Functions::SimpleIteration(Matrix A, Vector b, double eps, int& k)
 {
	int N = A.Size();
	Vector x(N);

	for(int i = 0; i < N; ++i)
		x[i] = 1;

	double t = 2 / (A.MinValue(eps) + A. MaxValue(eps));

	Vector r = A * x - b;
	k = 0;
	while (r.Norm() > eps)
	{
		++k;
		
		x = x  - t * r;
		r = A * x - b;
	}
	return x;}

Vector Functions::BisectionMethod(Matrix A, double eps)
{
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
	std::cout<<"A\n"<<A<<"\n";
	double a, b;

	a = A.Segment()[0];
	//std::cout<<"a\n"<<a<<std::endl;
	b = A.Segment()[1];

	return Functions::Serch(a, b, eps, A);
}