#pragma once

#include "Vector.h"
#include <fstream>

class Matrix
{
private:
	Vector *matrix;
	int size;
public:
	Matrix();
	Matrix(int s);
	Matrix(const Matrix &m);
	Vector operator[](const int i) const;
	Vector& operator[](const int i);
	friend std::ostream& operator<<(std::ostream& os, const Matrix& m);
	friend std::istream& operator>>(std::istream& is, Matrix& m);
	int Size() const;
	Matrix& operator=(const Matrix& m);
	friend Matrix operator*(const Matrix& a, const Matrix& b);
	friend Vector operator*(const Matrix& a, const Vector& b);
	friend Matrix operator*(const double x, Matrix& a);
	friend Matrix operator+(const Matrix& a, const Matrix& b);
	friend Matrix operator-(const Matrix& a, const Matrix& b);
	friend void QRDecomposition(const Matrix& a, Matrix& q, Matrix& r);
	void SwapLine(int i, int j);
	void Identity();
	Matrix Inverse() const;
	double MaxValue(double eps);
	double MinValue(double eps);
	double Norm();
	double Norm1();
	double ConditionNumber();
	Vector Segment() const;
	~Matrix();
};