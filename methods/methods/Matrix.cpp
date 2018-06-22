#include "Matrix.h"
#include "Functions.h"
#include <iostream>
#include <cmath>

Matrix::Matrix()
{
	matrix = NULL;
	size = 0;
}

int Matrix::Size() const
{
	return size;
}

Matrix::Matrix(int s)
{
	size = s;
	matrix = new Vector[s];
	for(int i = 0; i < s; ++i)
		matrix[i] = Vector(s);
}

Matrix::Matrix(const Matrix &m)
{
	size = m.size;
	matrix = new Vector[size];
	for(int i = 0; i < size; ++i)
		matrix[i] = Vector(size);
	for (int i = 0; i < size; ++i)
		for (int j = 0; j < size; ++j)
			matrix[i][j] = m[i][j];
}

Vector Matrix::operator[](const int i) const
{
	if (i < 0 || i >= size)
		std::cout<<"ERROR"<<std::endl;
	return matrix[i];
}

Vector& Matrix::operator[](const int i)
{
	if (i < 0 || i >= size)
		std::cout<<"ERROR"<<std::endl;
	return matrix[i];
}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
	for(int i = 0; i < m.Size(); ++i)
		for(int j = 0; j < m.Size(); ++j)
		{
			os << m[i][j] << ' ';
			if (j == m.Size() -1)
				os << '\n';
		}
	
	return os;
}

std::istream& operator>>(std::istream& is, Matrix& m)
{
	int size;
	is>>size;
	m = Matrix(size);
	for(int i = 0; i < size; ++i)
		for(int j = 0; j < size; ++j)
			is >> m[i][j];

	return is;
}

Matrix& Matrix::operator=(const Matrix& m)
{
	if (this != &m)
	{
		delete[] matrix;

		size = m.Size();

		matrix = new Vector[size];
		for (int i = 0; i < size; ++i)
			matrix[i] = m[i];
	}
	
	return *this;
}

Matrix operator*(const Matrix&a, const Matrix& b)
{
	Matrix m(a.Size());
	for (int i = 0; i < b.Size(); ++i)
		for (int k = 0; k < b.Size(); ++k)
			for (int j = 0; j < b.Size(); ++j)
				m[i][k] += a[i][j]*b[j][k];
	return m;
}

Vector operator*(const Matrix& a, const Vector& b)
{
	Vector v(a.Size());
	for(int i = 0; i < a.Size(); ++i)
		for(int j = 0; j < a.Size(); ++j)
			v[i] += a[i][j]*b[j];
	return v;
}

Matrix operator*(const double x, Matrix& a)
{
	for(int i = 0; i < a.Size(); ++i)
		a[i] = x * a[i];
	return a;
}

Matrix operator+(const Matrix& a, const Matrix& b)
{
	Matrix c(a.Size());
	for(int i = 0; i < a.Size(); ++i)
		for(int j = 0; j < a.Size(); ++j)
			c[i][j] = a[i][j] + b[i][j];
	return c;
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
	Matrix c(a.Size());
	for(int i = 0; i < a.Size(); ++i)
		for(int j = 0; j < a.Size(); ++j)
			c[i][j] = a[i][j] - b[i][j];
	return c;
}

void QRDecomposition(const Matrix& A, Matrix& q, Matrix& r)
{
	Matrix B(A.Size());
	Matrix Q(A.Size());
	Matrix R(A.Size());
	Matrix a = A;
	Q.Identity();
	R.Identity();
	for(int i = 0; i < A.Size(); ++i)
	{
		for(int j = i + 1; j < A.Size(); ++j)
		{
			B.Identity();
			B[i][i] = a[i][i] / sqrt(a[i][i] * a[i][i] + a[j][i] * a[j][i]);//  cos
			B[i][j] = a[j][i] / sqrt(a[i][i] * a[i][i] + a[j][i] * a[j][i]);// sin
			B[j][i] = - B[i][j];
			B[j][j] = B[i][i];
			a = B * a;
		}
	}
	r = a;
	q = (a * A.Inverse()).Inverse();
	//std::cout<<r<<std::endl;
}

void Matrix::SwapLine(int i, int j)
{
	Vector v(size);
	for(int k = 0; k < size; k++)
		v[k] = matrix[i][k];
	for(int k = 0; k < size; k++)
		matrix[i][k] = matrix[j][k];
	for(int k = 0; k < size; k++)
		matrix[j][k] = v[k];
}

void Matrix::Identity()
{
	for(int i = 0; i < size; ++i)
		for(int j = 0; j < size; ++j)
			if(i == j)
				matrix[i][j] = 1;
			else
				matrix[i][j] = 0;
}

Matrix Matrix::Inverse() const
{
	Vector x(size);
	Matrix B(size);
	double detA;
	for(int i = 0; i < size; ++i)
	{
		x[i] = 1;
		if(i > 0)
			x[i - 1] = 0;
		for(int j = 0; j  < size; ++j)
			B[j][i] = Functions::GaussS(*this, x,  detA)[j];		
	}
	return B;
}

double Matrix::MaxValue(double eps)
{
	Vector x(size), y(size);
	for(int i = 0; i < size; ++i)
		x[i] = 1;
	y = (1 / x.Norm()) * x;
	//int i = 0;
	Matrix A = *this;
	double l1 = 0.0, l2 = 0.0;
	do
	{
		l1 = l2;
		x = A * x;
		l2 = ((A * y) * y) / (y * y);
		y = (1 / x.Norm()) * x;
	}
	while (abs(l2 - l1) > eps);
	return l2;
}

double Matrix::MinValue(double eps)
{
	Matrix A = *this;
	Matrix E(A.Size());
	E.Identity();
	Matrix B = A.MaxValue(eps) * E - A;
	//std::cout<<B<<std::endl;
	return abs(B.MaxValue(eps) - A.MaxValue(eps));
}

double Matrix::Norm()
{
	double n = 0.0;
	for(int i = 0; i < size; ++i)
		for(int j = 0; j < size; ++j)
			n += matrix[i][j] * matrix[i][j];
	return sqrt(n);
}

double Matrix::Norm1()
{
	double n = 0.0, m;
	for(int i = 0; i < size; ++i)
	{
		m = n;
		n = 0.0;
		for(int j = 0; j < size; ++j)
		{
			n += abs(matrix[i][j]);
		}
		n = (n > m ? n : m);
	}
	return n;
}

Vector Matrix::Segment() const
{
	Matrix A = *this;
	Vector c(size), f(size), v(2);
	c[0] = A[0][0];

	for(int i = 0; i < A.Size(); ++i)
	{
		c[i] = A[i][i];
		if(i < A.Size() - 1)
			f[i + 1] = -A[i][i +1];
	}
	//std::cout<<"f\n"<<f<<std::endl;
	//system("pause");
	//std::cout<<"c\n"<<c<<std::endl;
	//system("pause");
	double min1, min2, max1, max2;
	min1 = (c[0] - abs(f[1]) < c[size -1] - abs(f[size - 1]) ? c[0] - abs(f[1]) : c[size -1] - abs(f[size - 1]));
	max1 = (c[0] + abs(f[1]) > c[size -1] + abs(f[size - 1]) ? c[0] + abs(f[1]) : c[size -1] + abs(f[size - 1]));
	min2 = c[1] - abs(f[1]) - abs(f[2]);
	max2 = c[1] + abs(f[1]) + abs(f[2]);
	for(int i = 2; i < size - 1; ++i)
	{
		min2 = (c[i] - abs(f[i]) - abs(f[i + 1]) < min2 ? c[i] - abs(f[i]) - abs(f[i + 1]) : min2);
		max2 = (c[i] + abs(f[i]) + abs(f[i + 1]) > max2 ? c[i] + abs(f[i]) + abs(f[i + 1]) : max2);
	}
	v[0] = (min1 < min2 ? min1 : min2);
	v[1] = (max1 > max2 ? max1 : max2);
	return v;
}

double Matrix::ConditionNumber()
{
	return Matrix::Norm() * Matrix::Inverse().Matrix::Norm();
}

Matrix::~Matrix()
{
	delete[] matrix;
}
