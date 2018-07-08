#include "Vector.h"
#include <iostream>
#include <cmath>

Vector::Vector()
{
	length = 0;
	vector = NULL;
}
 
Vector::Vector(const int l)
{
	length = l;
	vector = new double[l];
	for (int i = 0; i < l; ++i)
		vector[i] = 0.0;
}
 
double Vector::operator[](const int i) const
{
    if (i < 0 || i >= length)
    {
        std::cout << "Vector index out of range" << std::endl;
    }
	return vector[i];
}

double& Vector::operator[](const int i)
{
    if (i < 0 || i >= length)
    {
        std::cout << "Vector index out of range" << std::endl;
    }
	return vector[i];
}
 
size_t Vector::Length() const
{
	return length;
}
 
Vector operator+(const Vector a, const Vector b)
{
	Vector c(a.Length());
	for (int i = 0; i < a.Length(); ++i)
	{
		c[i] = a[i] + b[i];
	}
	return c;
}
 
Vector::Vector(const Vector &a)
{
	length = a.Length();
	vector = new double [a.Length()];
	for(int i = 0; i < length; ++i)
		vector[i] = a[i];
}

Vector& Vector::operator=(const Vector& a)
{
	if (this != &a)
	{
		delete[] vector;
		length = a.Length();
		vector = new double[length];
		for (int i = 0; i < length; ++i)
			vector[i] = a[i];
	}
	
	return *this;
}

void Vector::Swap(const int i, const int j)
{
	double k;
	k = vector[i];
	vector[i] = vector[j];
	vector[j] = k;
}

std::ostream& operator<<(std::ostream& os, const Vector& a)
{
	for(int i = 0; i < a.Length(); ++i)
	{
		os<<a[i]<<'\n';
	}
	return os;
}
 
std::istream& operator>>(std::istream& is, Vector& a)
{
	int l;
	is>>l;
	a = Vector(l);
	for(int i = 0; i < a.Length(); ++i)
		is>>a[i];
	return is;
}

double operator*(const Vector& a, const Vector& b)
{
	double s = 0.0;
	for(int i = 0; i < a.Length(); ++i)
		s += a[i] * b[i];
	return s;
}

Vector operator-(const Vector& a, const Vector& b)
{
	Vector c = a;
	for(int i = 0; i < a.Length(); ++i)
		c[i] -= b[i];
	return c;
}

double Vector::Norm() const
{
	double n = 0.0;
	for(int i = 0; i < length; ++i)
		n += vector[i] * vector[i];
	return sqrt(n);
}

Vector operator*(const double a, Vector& b)
{
	for(int i = 0; i < b.Length(); ++i)
		b[i] = a * b[i];
	return b;
}

Vector::~Vector()
{
	if(vector != NULL)
		delete[] vector;
}