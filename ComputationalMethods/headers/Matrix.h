#pragma once

#include "Vector.h"
#include <vector>

class Matrix
{
public:
    Matrix();
    Matrix(size_t s);
    Matrix(std::initializer_list<Vector> list);

    size_t Size() const;
    void SwapLine(size_t i, size_t j);
    void Identity();
    Matrix Inverse() const;
    double MaxValue(double eps);
    double MinValue(double eps);
    double Norm();
    double Norm1();
    double ConditionNumber();
    Vector ValuesBounds() const;

    Vector operator[](const size_t i) const;
    Vector& operator[](const size_t i);

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    friend std::istream& operator>>(std::istream& is, Matrix& m);
    friend Matrix operator*(const Matrix& m1, const Matrix& m2);
    friend Vector operator*(const Matrix& m, const Vector& v);
    friend Matrix operator*(const double scalar, Matrix& m);
    friend Matrix operator+(const Matrix& m1, const Matrix& m2);
    friend Matrix operator-(const Matrix& m1, const Matrix& m2);

private:
    std::vector<Vector> matrix;
    size_t size;
};