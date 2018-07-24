#include "Vector.h"
#include <iostream>
#include <cassert>

Vector::Vector() : Vector(0) {}

Vector::Vector(const size_t size) : vector(std::vector<double>(size, 0)) {}

Vector::Vector(std::initializer_list<double> list) : vector(std::vector<double>(list)) {}

double Vector::operator[](const size_t i) const
{
    assert(i >= 0 && i < vector.size());
    return vector[i];
}

double& Vector::operator[](const size_t i)
{
    assert(i >= 0 && i < vector.size());
    return vector[i];
}
 
size_t Vector::Size() const
{
    return vector.size();
}
 
Vector operator+(const Vector v1, const Vector v2)
{
    assert(v1.Size() == v2.Size());
    Vector resultVector(v1);
    for (size_t i = 0; i < v1.Size(); ++i)
    {
        resultVector[i] += v2[i];
    }
return resultVector;
}

void Vector::Swap(const size_t i, const size_t j)
{
    std::swap(vector[i], vector[j]);
}

std::ostream& operator<<(std::ostream& os, const Vector& v)
{
    for(size_t i = 0; i < v.Size(); ++i)
    {
        os << v[i] << std::endl;
    }
    return os;
}
 
std::istream& operator>>(std::istream& is, Vector& v)
{
    size_t size;
    is>>size;
    v = Vector(size);
    for(size_t i = 0; i < size; ++i)
        is>>v[i];
    return is;
}

double operator*(const Vector& v1, const Vector& v2)
{
    assert(v1.Size() == v2.Size());
    double result = 0.0;
    for(size_t i = 0; i < v1.Size(); ++i)
        result += v1[i] * v2[i];
    return result;
}

Vector operator-(const Vector& v1, const Vector& v2)
{
    assert(v1.Size() == v2.Size());
    Vector result(v1);
    for(size_t i = 0; i < v1.Size(); ++i)
        result[i] -= v2[i];
    return result;
}

double Vector::Norm() const
{
    double normSquare = 0.0;
    for(size_t i = 0; i < vector.size(); ++i)
        normSquare += vector[i] * vector[i];
    return sqrt(normSquare);
}

Vector operator*(const double scalar, Vector& v)
{
    for(size_t i = 0; i < v.Size(); ++i)
        v[i] = scalar * v[i];
    return v;
}