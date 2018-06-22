#pragma once

#include <fstream>

class Vector
{
private:
        int length;
        double *vector;
public:
        Vector();
        Vector(const int l);
        Vector(const Vector &a);
        double operator[](const int i) const;
        double& operator[](const int i);
        friend Vector operator+(const Vector a, const Vector b);
        Vector& operator=(const Vector& a);
		void Swap(const int i, const int j);
		//////////////
		friend std::ostream& operator<<(std::ostream& os, const Vector& a);
		friend std::istream& operator>>(std::istream& is, Vector& a);
		friend double operator*(const Vector& a, const Vector& b);
		friend Vector operator*(const double a, Vector& b);
		friend Vector operator-(const Vector& a, const Vector& b);
		//////////////

        int Length() const;
		double Norm() const;


        ~Vector();
};