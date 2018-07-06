#pragma once

#include <fstream>

class Vector
{
public:
    Vector();
    ~Vector();
    Vector(const int l);
    Vector(const Vector &a);

    void Swap(const int i, const int j);
    size_t Length() const;
    double Norm() const;

    double operator[](const int i) const;
    double& operator[](const int i);
    Vector& operator=(const Vector& a);

    friend std::ostream& operator<<(std::ostream& os, const Vector& a);
    friend std::istream& operator>>(std::istream& is, Vector& a);
    friend double operator*(const Vector& a, const Vector& b);
    friend Vector operator*(const double a, Vector& b);
    friend Vector operator+(const Vector a, const Vector b);
    friend Vector operator-(const Vector& a, const Vector& b);

private:
    size_t length;
    double *vector;
};
