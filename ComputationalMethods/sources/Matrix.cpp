#include "Matrix.h"
#include "Functions.h"
#include <iostream>
#include <cassert>

size_t Matrix::Size() const
{
    return size;
}

Matrix::Matrix() : Matrix(0) {}

Matrix::Matrix(size_t s)
{
    size = s;
    matrix = std::vector<Vector>(s);
    for(size_t i = 0; i < s; ++i)
        matrix[i] = Vector(s);
}

Matrix::Matrix(std::initializer_list<Vector> list) : size(list.size()), matrix(std::vector<Vector>(list)) {}

std::ostream& operator<<(std::ostream& os, const Matrix& m)
{
    for(size_t i = 0; i < m.Size(); ++i)
        for(size_t j = 0; j < m.Size(); ++j)
        {
            os << m[i][j] << ' ';
            if (j == m.Size() -1)
                os << std::endl;
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

Matrix operator*(const Matrix& m1, const Matrix& m2)
{
    assert(m1.Size() == m2.Size());
    Matrix m(m1.Size());
    for (size_t i = 0; i < m1.Size(); ++i)
        for (size_t k = 0; k < m1.Size(); ++k)
            for (size_t j = 0; j < m1.Size(); ++j)
                m[i][k] += m1[i][j]*m2[j][k];
    return m;
}

Vector operator*(const Matrix& m, const Vector& v)
{
    assert(m.Size() == v.Size());
    Vector vector(m.Size());
    for(size_t i = 0; i < m.Size(); ++i)
        for(size_t j = 0; j < m.Size(); ++j)
            vector[i] += m[i][j]*v[j];
    return vector;
}

Matrix operator*(const double scalar, Matrix& m)
{
    for(size_t i = 0; i < m.Size(); ++i)
        m[i] = scalar * m[i];
    return m;
}

Matrix operator+(const Matrix& m1, const Matrix& m2)
{
    assert(m1.Size() == m2.Size());
    Matrix resultMatrix(m1);

    for(size_t i = 0; i < m1.Size(); ++i)
        for(size_t j = 0; j < m1.Size(); ++j)
            resultMatrix[i][j] += m2[i][j];

    return resultMatrix;
}

Matrix operator-(const Matrix& m1, const Matrix& m2)
{
    assert(m1.Size() == m2.Size());
    Matrix resultMatrix(m1);
    for(size_t i = 0; i < m1.Size(); ++i)
        for(size_t j = 0; j < m1.Size(); ++j)
            resultMatrix[i][j] -= m2[i][j];
    return resultMatrix;
}


void Matrix::SwapLine(size_t i, size_t j)
{
    std::swap(matrix[i], matrix[j]);
}

void Matrix::Identity()
{
    Matrix identityMatrix(size);

    for(size_t i = 0; i < size; ++i)
        identityMatrix[i][i] = 1;
    *this = identityMatrix;
}

Matrix Matrix::Inverse() const
{
    Vector vector(size);
    Matrix resultMatrix(size);
    double det;
    for(size_t i = 0; i < size; ++i)
    {
        vector[i] = 1;
        if(i > 0)
            vector[i - 1] = 0;
        for(size_t j = 0; j  < size; ++j)
            resultMatrix[j][i] = Functions::GaussS(*this, vector,  det)[j];
    }
    return resultMatrix;
}

double Matrix::MaxValue(double eps)
{
    Vector vector1(size), vector2(size);
    for(size_t i = 0; i < size; ++i)
        vector1[i] = 1;
    vector2 = (1 / vector1.Norm()) * vector1;

    Matrix matrix1 = *this;
    double coef1 = 0.0, coef2 = 0.0;
    do
    {
        coef1 = coef2;
        vector1 = matrix1 * vector1;
        coef2 = ((matrix1 * vector2) * vector2) / (vector2 * vector2);
        vector2 = (1 / vector1.Norm()) * vector1;
    }
    while (abs(coef2 - coef1) > eps);

    return coef2;
}

double Matrix::MinValue(double eps)
{
    Matrix identityMatrix(size);
    identityMatrix.Identity();
    Matrix matrix2 = this->MaxValue(eps) * identityMatrix - *this;

    return abs(matrix2.MaxValue(eps) - this->MaxValue(eps));
}

double Matrix::Norm()
{
    double normSquare = 0.0;
    for(size_t i = 0; i < size; ++i)
        for(size_t j = 0; j < size; ++j)
            normSquare += matrix[i][j] * matrix[i][j];
    return sqrt(normSquare);
}

double Matrix::Norm1()
{
    double currNorm = 0.0;
    double prevNorm = 0.0;
    for(size_t i = 0; i < size; ++i)
    {
        prevNorm = currNorm;
        currNorm = 0.0;
        for(size_t j = 0; j < size; ++j)
        {
            currNorm += abs(matrix[i][j]);
        }
        currNorm = (currNorm > prevNorm ? currNorm : prevNorm);
    }
    return currNorm;
}

Vector Matrix::ValuesBounds() const
{
    Vector resultVector(2), vector1(size), vector2(size + 1);
    vector1[0] = matrix[0][0];

    for(size_t i = 0; i < size; ++i)
    {
        vector1[i] = matrix[i][i];
        if(i < size - 1)
            vector2[i + 1] = -matrix[i][i +1];
    }
    double min1, min2, max1, max2;
    min1 = (vector1[0] - abs(vector2[1]) < vector1[size -1] - abs(vector2[size - 1])
        ? vector1[0] - abs(vector2[1]) : vector1[size -1] - abs(vector2[size - 1]));
    max1 = (vector1[0] + abs(vector2[1]) > vector1[size -1] + abs(vector2[size - 1])
        ? vector1[0] + abs(vector2[1]) : vector1[size -1] + abs(vector2[size - 1]));
    min2 = vector1[1] - abs(vector2[1]) - abs(vector2[2]);
    max2 = vector1[1] + abs(vector2[1]) + abs(vector2[2]);
    for(size_t i = 2; i < size - 1; ++i)
    {
        min2 = (vector1[i] - abs(vector2[i]) - abs(vector2[i + 1]) < min2
            ? vector1[i] - abs(vector2[i]) - abs(vector2[i + 1]) : min2);
        max2 = (vector1[i] + abs(vector2[i]) + abs(vector2[i + 1]) > max2
            ? vector1[i] + abs(vector2[i]) + abs(vector2[i + 1]) : max2);
    }
    resultVector[0] = (min1 < min2 ? min1 : min2);
    resultVector[1] = (max1 > max2 ? max1 : max2);
    return resultVector;
}

double Matrix::ConditionNumber()
{
    return Matrix::Norm() * Matrix::Inverse().Matrix::Norm();
}