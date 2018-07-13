#include "Vector.h"
#include <iostream>
#include <cmath>

Vector::Vector() : vector(std::vector<double>(0)) {}

Vector::Vector(const size_t size) : vector(std::vector<double>(size, 0)) {}

Vector::Vector(Vector && v)
{
    std::swap(vector, v.vector);
}

double Vector::operator[](const int i) const
{
	return vector[i];
}

double& Vector::operator[](const int i)
{
	return vector[i];
}
 
size_t Vector::Size() const
{
	return vector.size();
}
 
Vector operator+(const Vector a, const Vector b)
{
	Vector c(a.Size());
	for (int i = 0; i < a.Size(); ++i)
	{
		c[i] = a[i] + b[i];
	}
	return c;
}
 
Vector::Vector(const Vector &a)
{
	vector = std::vector<double>(a.vector);
}

Vector& Vector::operator=(const Vector& a)
{
    vector = a.vector;
	return *this;
}

Vector & Vector::operator=(Vector&& v)
{
    std::swap(vector, v.vector);
    return *this;
}

void Vector::Swap(const int i, const int j)
{
    std::swap(vector[i], vector[j]);
}

std::ostream& operator<<(std::ostream& os, const Vector& a)
{
	for(int i = 0; i < a.Size(); ++i)
	{
		os<<a[i]<<'\n';
	}
	return os;
}
 
std::istream& operator>>(std::istream& is, Vector& a)
{
	size_t size;
	is>>size;
	a = Vector(size);
	for(int i = 0; i < a.Size(); ++i)
		is>>a[i];
	return is;
}

double operator*(const Vector& a, const Vector& b)
{
	double s = 0.0;
	for(int i = 0; i < a.Size(); ++i)
		s += a[i] * b[i];
	return s;
}

Vector operator-(const Vector& a, const Vector& b)
{
	Vector c(a);
	for(int i = 0; i < a.Size(); ++i)
		c[i] -= b[i];
	return c;
}

double Vector::Norm() const
{
	double n = 0.0;
	for(int i = 0; i < vector.size(); ++i)
		n += vector[i] * vector[i];
	return sqrt(n);
}

Vector operator*(const double a, Vector& b)
{
	for(int i = 0; i < b.Size(); ++i)
		b[i] = a * b[i];
	return b;
}