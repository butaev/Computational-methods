#pragma once

#include <vector>
#include <fstream>

class Vector
{
public:
    Vector();
    Vector(const size_t size);
    Vector(const Vector &a);
    Vector(Vector && v);

    void Swap(const int i, const int j);
    size_t Size() const;
    double Norm() const;

    double operator[](const int i) const;
    double& operator[](const int i);
    Vector& operator=(const Vector& a);
    Vector& operator=(Vector&& v);

    friend std::ostream& operator<<(std::ostream& os, const Vector& a);
    friend std::istream& operator>>(std::istream& is, Vector& a);
    friend double operator*(const Vector& a, const Vector& b);
    friend Vector operator*(const double a, Vector& b);
    friend Vector operator+(const Vector a, const Vector b);
    friend Vector operator-(const Vector& a, const Vector& b);

private:
    std::vector<double> vector = std::vector<double>(0);
};
