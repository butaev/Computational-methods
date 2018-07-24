#pragma once

#include <vector>

class Vector
{
public:
    Vector();
    Vector(const size_t size);
    Vector(std::initializer_list<double> list);

    void Swap(const size_t i, const size_t j);
    size_t Size() const;
    double Norm() const;

    double operator[](const size_t i) const;
    double& operator[](const size_t i);

    friend std::ostream& operator<<(std::ostream& os, const Vector& v);
    friend std::istream& operator>>(std::istream& is, Vector& v);
    friend double operator*(const Vector& v1, const Vector& v2);
    friend Vector operator*(const double scalar, Vector& v);
    friend Vector operator+(const Vector v1, const Vector v2);
    friend Vector operator-(const Vector& v1, const Vector& v2);

private:
    std::vector<double> vector;
};
