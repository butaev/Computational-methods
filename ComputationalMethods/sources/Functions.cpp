#include "Functions.h"
#include "Vector.h"
#include "Matrix.h"
#include <iostream>
#include <cmath>

void Functions::QRDecomposition(Matrix m, Matrix& q, Matrix& r)
{
    Matrix rotationMatrix(m.Size());
    Matrix upperTriangular = m;

    for (size_t i = 0; i < m.Size(); ++i)
    {
        for (size_t j = i + 1; j < m.Size(); ++j)
        {
            rotationMatrix.Identity();
            rotationMatrix[i][i] = upperTriangular[i][i] / sqrt(upperTriangular[i][i] * upperTriangular[i][i]
                + upperTriangular[j][i] * upperTriangular[j][i]);//  cos
            rotationMatrix[i][j] = upperTriangular[j][i] / sqrt(upperTriangular[i][i] * upperTriangular[i][i]
                + upperTriangular[j][i] * upperTriangular[j][i]);// sin
            rotationMatrix[j][i] = -rotationMatrix[i][j];
            rotationMatrix[j][j] = rotationMatrix[i][i];
            upperTriangular = rotationMatrix * upperTriangular;
        }
    }
    r = upperTriangular;

    q = (upperTriangular * m.Inverse()).Inverse();
}

int Functions::NumLeftValues(double boundary, const Matrix& m)
{
    Vector vector(m.Size());
    vector[0] = m[0][0] - boundary;
    for (size_t i = 1; i < m.Size(); ++i)
        vector[i] = (m[i][i] - boundary) - (m[i][i - 1] * m[i][i - 1]) / vector[i - 1];
    int numValues = 0;
    for (size_t i = 0; i < m.Size(); ++i)
        if (vector[i] < 0)
            ++numValues;
    return numValues;
}

Vector Functions::SearchValues(double& leftBound, double& rightBound, double eps, const Matrix& m)
{
    Vector leftValues, rightValues;
    double middleBound = (leftBound + rightBound) / 2.0;
    int numLeftValLeft = NumLeftValues(leftBound, m), numLeftValRight = NumLeftValues(rightBound, m),
        numLeftValMiddle = NumLeftValues(middleBound, m);
    Vector values(numLeftValRight - numLeftValLeft);
    
    if((rightBound - leftBound) < eps && (numLeftValRight - numLeftValLeft) == 1)
    {
        values[0] = middleBound;
        return values;
    }
    
    if(numLeftValMiddle - numLeftValLeft == 0)
        rightValues = SearchValues(middleBound, rightBound, eps, m);
    else if(numLeftValRight - numLeftValMiddle == 0)
        leftValues = SearchValues(leftBound, middleBound, eps, m);
    else
    {
        leftValues = SearchValues(leftBound, middleBound, eps, m);
        rightValues = SearchValues(middleBound, rightBound, eps, m);
    }
    
    for(size_t i = 0; i < leftValues.Size(); ++i)
        values[i] = leftValues[i];
    for(size_t i = 0; i < rightValues.Size(); ++i)
        values[i + leftValues.Size()] = rightValues[i];

    return values;
}

Vector Functions::Gauss(Matrix m, Vector v, double& detm)
{
    double coef1 = 0.0, coef2 = 0.0;
    detm = 1.0;
    size_t size = m.Size();

    for(size_t i = 0; i < m.Size(); ++i)
    {
        detm *= m[i][i];
        coef1 = m[i][i];
        v[i] /= m[i][i];
        for(size_t k = i; k < size; ++k)
            m[i][k] = m[i][k] / coef1;

        for(size_t j = i + 1; j < size; ++j)
        {    
            coef2 = m[j][i];
            v[j] -= v[i] * coef2;
            for(size_t k = i; k < size; ++k)
                m[j][k] -= m[i][k] * coef2;
        }
    }

    Vector resultVector(size);
    resultVector[size - 1] = v[size - 1];
    
    for(int i = size - 2; i >= 0; --i)
    {
        coef1 = 0.0;
        for(int j = size - 1; j > i; --j)
            coef1 -=  m[i][j] * resultVector[j];
        resultVector[i] = v[i] + coef1;
    }

    return resultVector;
}

Vector Functions::GaussS(Matrix m, Vector v, double& detm)
{
    double coef1 = 0.0, coef2 = 0.0;
    detm = 1.0;
    int size = m.Size();

    Vector resultVector(size);
    for(size_t i = 0; i < m.Size(); ++i)
    {
        coef2 = m[i][i];
        int j;
        for(int l = 0; l < size; ++l)
        {
            if(m[l][i] > coef2)
            {
                coef2 = m[l][i];
                j = l;
            }
            else
                j = i;
        }
        if(i != j)
        {
            m.SwapLine(i,j);
            v.Swap(i,j);
            detm *= -1;
        }
        if(m[i][i] != 0 && detm != 0.0)
            detm *= m[i][i];
        else
        {
            detm = 0.0;
            return resultVector;
        }

        coef1 = m[i][i];
        v[i] /= m[i][i];
        for(int k = i; k < size; ++k)
            m[i][k] = m[i][k] / coef1;

        for(int j = i + 1; j < size; ++j)
        {    
            coef2 = m[j][i];
            v[j] -= v[i] * coef2;
            for(int k = i; k < size; ++k)
                m[j][k] -= m[i][k] * coef2;
        }
    }
    
    resultVector[size - 1] = v[size - 1];
    
    for(int i = size - 2; i >= 0; --i)
    {
        coef1 = 0.0;
        for(int j = size - 1; j > i; --j)
            coef1 -=  m[i][j] * resultVector[j];
        resultVector[i] = v[i] + coef1;
    }

    return resultVector;
}

Vector Functions::RotationMethod(Matrix m, Vector v)
{
    Matrix rotationMatrix(m.Size());

    for(size_t i = 0; i < m.Size(); ++i)
    {
        int k = i;
        if (m[i][i] == 0)
        {
            for(size_t l = i + 1; l < m.Size(); ++l)
            {
                if(m[l][i] > 0)
                    k = l;
                else
                    if(l == m.Size() - 1)
                    {
                        Vector resultVector(m.Size());
                        return resultVector;
                    }
            }
        }

        if(i != k)
        {
            m.SwapLine(i, k);
            v.Swap(i, k);
        }

        for(size_t j = i + 1; j < m.Size(); ++j)
        {
            rotationMatrix.Identity();
            rotationMatrix[i][i] = m[i][i] / sqrt(m[i][i] * m[i][i] + m[j][i] * m[j][i]);//  cos
            rotationMatrix[i][j] = m[j][i] / sqrt(m[i][i] * m[i][i] + m[j][i] * m[j][i]);// sin
            rotationMatrix[j][i] = -rotationMatrix[i][j];
            rotationMatrix[j][j] = rotationMatrix[i][i];
            v = rotationMatrix * v;
            m = rotationMatrix * m;
        }
    }

    Vector resultVector(m.Size());
    resultVector[m.Size() - 1] = v[m.Size() - 1] / m[m.Size() - 1][m.Size() - 1];

    double coef = 0.0;
    for(int i = m.Size() - 2; i >= 0; --i)
    {
        coef = 0.0;
        for(int j = m.Size() - 1; j > i; --j)
            coef -=  m[i][j] * resultVector[j];
        resultVector[i] = (v[i] + coef) / m[i][i];
    }

    return resultVector;
}

Vector Functions::TridiagMatrix(Matrix m, Vector f)
{
    size_t size = m.Size();
    Vector leftDiagonal(size), rightDiagonal(size), mainDiagonal(size);
    Vector alpha(size), beta(size);
    Vector resultVector(size);

    for(size_t i = 0; i < size; ++i)
    {
        mainDiagonal[i] = m[i][i];
        if(i < size - 1)
        {
            rightDiagonal[i] = -m[i][i + 1];
            leftDiagonal[i] = -m[i + 1][i];
        }
    }

    alpha[0] = rightDiagonal[0] / mainDiagonal[0];
    beta[0] = f[0] / mainDiagonal[0];
    for(size_t i = 1; i < size; ++i)
    {
        alpha[i] = rightDiagonal[i] / (mainDiagonal[i] - alpha[i - 1] * leftDiagonal[i]);
        beta[i] = (f[i] + beta[i - 1] * leftDiagonal[i]) / (mainDiagonal[i] - alpha[i - 1] * leftDiagonal[i]);
    }

    resultVector[size - 1] = beta[size - 1];

    for(int i = size - 2; i >= 0; --i)
    {
        resultVector[i] = alpha[i] * resultVector[i + 1] + beta[i];
    }

    return resultVector;
}

Vector Functions::QRAlgorithm(Matrix m, double eps)
{
    size_t size = m.Size();
    Matrix orthogonalMatrix, upperTriangular; //upperTriangular
    Matrix transformationMatrix = m;
    Vector resultVector(size);

    while (true)
    {
        size_t stopInd = 0;
        QRDecomposition(m, orthogonalMatrix, upperTriangular);
        transformationMatrix = upperTriangular * orthogonalMatrix;
        
        for(size_t j = 0; j < size - 1; ++j)
        {
            double norm = 0.0;

            for(size_t i = j + 1; i < size; ++i)
                norm += transformationMatrix[i][j] * transformationMatrix[i][j];

            if(sqrt(norm) < eps)
            {
                ++stopInd;
                for(size_t k = 0; k < size; ++k)
                    resultVector[k] = transformationMatrix[k][k];
            
                if (stopInd == size - 1)
                    return resultVector;
            }
        }
    }
}

Vector Functions::GMRES(Matrix m, Vector v, double eps, int& numIter)
{
    size_t size = m.Size();
    Vector resultVector(size);

    for(size_t i = 0; i < size; ++i)
        resultVector[i] = 1;

    Vector errorVector = m * resultVector - v;
    numIter = 0;

    while (errorVector.Norm() > eps)
    {
        ++numIter;
        Vector vector = m * errorVector;
        resultVector = resultVector  - ((vector * errorVector) / (vector * vector)) * errorVector;
        errorVector = m * resultVector - v;
    }

    return resultVector;
}

Vector Functions::SimpleIteration(Matrix m, Vector v, double eps, int& numIter)
 {
    size_t size = m.Size();
    Vector resultVector(size);

    for(size_t i = 0; i < size; ++i)
        resultVector[i] = 1;

    double coef = 2 / (m.MinValue(eps) + m. MaxValue(eps));

    Vector errorVector = m * resultVector - v;
    numIter = 0;
    while (errorVector.Norm() > eps)
    {
        ++numIter;
        
        resultVector = resultVector  - coef * errorVector;
        errorVector = m * resultVector - v;
    }

    return resultVector;
}

Vector Functions::BisectionMethod(Matrix m, double eps)
{
    size_t size = m.Size();
    Matrix identityMatrix(size), matrix1(size), matrix2(size);
    Vector vector1(size), vector2(size);
    double coef = 0.0;
     
    identityMatrix.Identity();

    // matrix transformation to Hessenberg form
    for(size_t i = 0; i < size - 1; ++i)
    {
        Vector vector3(size);
        vector2[i + 1] = 1;
        if(i > 0)
            vector2[i] = 0;

        for(size_t j = i + 1; j < size; ++j)
        {
                 vector3[j] = m[j][i];
        }

        coef = (vector3[i + 1] > 0 ? 1 : -1) * vector3.Norm();
        vector1 = vector3 + coef * vector2;

        for(size_t j = 0; j < size; ++j)
            for(size_t k = 0; k < size; ++k)
                matrix1[j][k] = vector1[j] * vector1[k];

        if (!(vector1 * vector1))
            continue;
        matrix2 = identityMatrix - (2 / (vector1 * vector1)) * matrix1;
        m = matrix2 * m * matrix2;
    }

    Vector bounds(2);
    bounds = m.ValuesBounds();

    return Functions::SearchValues(bounds[0], bounds[1], eps, m);
}