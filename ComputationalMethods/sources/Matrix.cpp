#include "Matrix.h"
#include "Functions.h"
#include <iostream>
#include <cmath>
#include <cassert>

int Matrix::Size() const
{
    return size;
}

Matrix::Matrix() : size(0), matrix(std::vector<Vector>(0)) {}

Matrix::Matrix(size_t s)
{
    size = s;
    matrix = std::vector<Vector>(s);
    for(size_t i = 0; i < s; ++i)
        matrix[i] = Vector(s);
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

Vector Matrix::operator[](const size_t i) const
{
    assert(i >= 0 &&  i < size);
    return matrix[i];
}

Vector& Matrix::operator[](const size_t i)
{
    assert(i >= 0 && i < size);
    return matrix[i];
}

Matrix operator*(const Matrix& a, const Matrix& b)
{
    assert(a.Size() == b.Size());
    Matrix m(a.Size());
    for (int i = 0; i < b.Size(); ++i)
        for (int k = 0; k < b.Size(); ++k)
            for (int j = 0; j < b.Size(); ++j)
                m[i][k] += a[i][j]*b[j][k];
    return m;
}

Vector operator*(const Matrix& a, const Vector& b)
{
    assert(a.Size() == b.Size());
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
    assert(a.Size() == b.Size());
    Matrix c(a.Size());
    for(int i = 0; i < a.Size(); ++i)
        for(int j = 0; j < a.Size(); ++j)
            c[i][j] = a[i][j] + b[i][j];
    return c;
}

Matrix operator-(const Matrix& a, const Matrix& b)
{
    assert(a.Size() == b.Size());
    Matrix c(a.Size());
    for(int i = 0; i < a.Size(); ++i)
        for(int j = 0; j < a.Size(); ++j)
            c[i][j] = a[i][j] - b[i][j];
    return c;
}


void Matrix::SwapLine(int i, int j)
{
    std::swap(matrix[i], matrix[j]);
}

void Matrix::Identity()
{
    for(size_t i = 0; i < size; ++i)
        for(size_t j = 0; j < size; ++j)
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
    for(size_t i = 0; i < size; ++i)
    {
        x[i] = 1;
        if(i > 0)
            x[i - 1] = 0;
        for(size_t j = 0; j  < size; ++j)
            B[j][i] = Functions::GaussS(*this, x,  detA)[j];
    }
    return B;
}

double Matrix::MaxValue(double eps)
{
    Vector x(size), y(size);
    for(size_t i = 0; i < size; ++i)
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
    //std::cout<< "bMaxValue = " << B.MaxValue(eps) <<std::endl;
    //std::cout << "aMaxValue = " <<A.MaxValue(eps) << std::endl;

    return abs(B.MaxValue(eps) - A.MaxValue(eps));
}

double Matrix::Norm()
{
    double n = 0.0;
    for(size_t i = 0; i < size; ++i)
        for(size_t j = 0; j < size; ++j)
            n += matrix[i][j] * matrix[i][j];
    return sqrt(n);
}

double Matrix::Norm1()
{
    double n = 0.0;
    double m;
    for(size_t i = 0; i < size; ++i)
    {
        m = n;
        n = 0.0;
        for(size_t j = 0; j < size; ++j)
        {
            n += abs(matrix[i][j]);
        }
        n = (n > m ? n : m);
    }
    return n;
}

Vector Matrix::ValuesBounds() const
{
    Matrix A = *this;
    Vector c(size), f(size + 1), v(2);
    c[0] = A[0][0];

    for(int i = 0; i < A.Size(); ++i)
    {
        c[i] = A[i][i];
        if(i < A.Size() - 1)
            f[i + 1] = -A[i][i +1];
    }
    double min1, min2, max1, max2;
    min1 = (c[0] - abs(f[1]) < c[size -1] - abs(f[size - 1])
        ? c[0] - abs(f[1]) : c[size -1] - abs(f[size - 1]));
    max1 = (c[0] + abs(f[1]) > c[size -1] + abs(f[size - 1])
        ? c[0] + abs(f[1]) : c[size -1] + abs(f[size - 1]));
    min2 = c[1] - abs(f[1]) - abs(f[2]);
    max2 = c[1] + abs(f[1]) + abs(f[2]);
    for(size_t i = 2; i < size - 1; ++i)
    {
        min2 = (c[i] - abs(f[i]) - abs(f[i + 1]) < min2
            ? c[i] - abs(f[i]) - abs(f[i + 1]) : min2);
        max2 = (c[i] + abs(f[i]) + abs(f[i + 1]) > max2
            ? c[i] + abs(f[i]) + abs(f[i + 1]) : max2);
    }
    v[0] = (min1 < min2 ? min1 : min2);
    v[1] = (max1 > max2 ? max1 : max2);
    return v;
}

double Matrix::ConditionNumber()
{
    return Matrix::Norm() * Matrix::Inverse().Matrix::Norm();
}