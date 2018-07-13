#pragma once

#include "Vector.h"
#include <fstream>
#include <vector>

class Matrix
{
public:
    Matrix();
    Matrix(size_t s);
    Matrix(const Matrix &m);
    Matrix(Matrix && m);
	
    int Size() const;
    void SwapLine(int i, int j);
    void Identity();
    Matrix Inverse() const;
    double MaxValue(double eps);
    double MinValue(double eps);
    double Norm();
    double Norm1();
    double ConditionNumber();
    Vector ValuesBounds() const;

    Vector operator[](const int i) const;
    Vector& operator[](const int i);
    Matrix& operator=(const Matrix& m);
    Matrix& operator=(Matrix && m);

    friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
    friend std::istream& operator>>(std::istream& is, Matrix& m);
	friend Matrix operator*(const Matrix& a, const Matrix& b);
	friend Vector operator*(const Matrix& a, const Vector& b);
	friend Matrix operator*(const double x, Matrix& a);
	friend Matrix operator+(const Matrix& a, const Matrix& b);
	friend Matrix operator-(const Matrix& a, const Matrix& b);

private:
    std::vector<Vector> matrix = std::vector<Vector>(0);
    size_t size = 0;
};