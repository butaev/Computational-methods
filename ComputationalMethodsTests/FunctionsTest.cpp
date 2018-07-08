#include "pch.h"
#include "Functions.h"

class FunctionTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        x = Vector(3);
        x[0] = 1;
        x[1] = -1;

        A = Matrix(2);
        A[0][0] = 1, A[0][1] = 2, A[1][0] = 0, A[1][1] = 2;
    }

    // virtual void TearDown() {}
    Matrix A;
    Vector x;
};

TEST_F(FunctionTest, QRDecomposition)
{
    Matrix q, r;
    Functions::QRDecomposition(A, q, r);

    EXPECT_LE(std::abs(q[0][0] - 1), 1e-6);
    EXPECT_LE(std::abs(q[0][1] - 0), 1e-6);
    EXPECT_LE(std::abs(q[1][0] - 0), 1e-6);
    EXPECT_LE(std::abs(q[1][1] - 1), 1e-6);

    EXPECT_LE(std::abs(r[0][0] - 1), 1e-6);
    EXPECT_LE(std::abs(r[0][1] - 2), 1e-6);
    EXPECT_LE(std::abs(r[1][0] - 0), 1e-6);
    EXPECT_LE(std::abs(r[1][1] - 2), 1e-6);
}

TEST_F(FunctionTest, NumLeftValues) {
    int numberVal = Functions::NumLeftValues(3, A);
    EXPECT_LE(std::abs(numberVal - 2), 1e-6);
}

TEST_F(FunctionTest, QRAlgorithm) {
    Vector eigenvalues(A.Size());
    eigenvalues = Functions::QRAlgorithm(A, 1e-6);
    EXPECT_LE(std::abs(eigenvalues[0] - 1), 1e-6);
    EXPECT_LE(std::abs(eigenvalues[1] - 2), 1e-6);
}

TEST_F(FunctionTest, BisectionMethod) {
    Vector eigenvalues(A.Size());
    eigenvalues = Functions::BisectionMethod(A, 1e-6);
    EXPECT_LE(std::abs(eigenvalues[0] - 1), 1e-6);
    EXPECT_LE(std::abs(eigenvalues[1] - 2), 1e-6);
}

TEST_F(FunctionTest, RotationMethod)
{
    Vector b = A * x;
    Vector xx = Functions::RotationMethod(A, b);
    EXPECT_LE(std::abs(xx[0] - x[0]), 1e-6);
    EXPECT_LE(std::abs(xx[1] - x[1]), 1e-6);
}

TEST_F(FunctionTest, GMRES)
{
    Vector b = A * x;
    int numIter = 0;
    Vector xx = Functions::GMRES(A, b, 1e-6, numIter);
    std::cout << "Number of iteration: " << numIter << std::endl;
    EXPECT_LE(std::abs(xx[0] - x[0]), 1e-6);
    EXPECT_LE(std::abs(xx[1] - x[1]), 1e-6);
}

TEST_F(FunctionTest, SimpleIteration)
{
    Vector b = A * x;
    int numIter = 0;
    Vector xx = Functions::SimpleIteration(A, b, 1e-6, numIter);
    std::cout << "Number of iteration: " << numIter << std::endl;
    EXPECT_LE(std::abs(xx[0] - x[0]), 1e-6);
    EXPECT_LE(std::abs(xx[1] - x[1]), 1e-6);
}

TEST_F(FunctionTest, TridiagMatrix)
{
    Vector b = A * x;
    Vector xx = Functions::TridiagMatrix(A, b);
    EXPECT_LE(std::abs(xx[0] - x[0]), 1e-6);
    EXPECT_LE(std::abs(xx[1] - x[1]), 1e-6);
}

TEST_F(FunctionTest, Gauss)
{
    Vector b = A * x;
    double detA = 0;
    Vector xx = Functions::Gauss(A, b, detA);
    EXPECT_LE(std::abs(xx[0] - x[0]), 1e-6);
    EXPECT_LE(std::abs(xx[1] - x[1]), 1e-6);
    EXPECT_LE(std::abs(detA - 2), 1e-6);
}

TEST_F(FunctionTest, GaussS)
{
    Vector b = A * x;
    double detA = 0;
    Vector xx = Functions::GaussS(A, b, detA);
    EXPECT_LE(std::abs(xx[0] - x[0]), 1e-6);
    EXPECT_LE(std::abs(xx[1] - x[1]), 1e-6);
    EXPECT_LE(std::abs(detA - 2), 1e-6);
}

TEST_F(FunctionTest, SearchValues)
{
    double a = 0.5;
    double b = 1.5;
    Vector values = Functions::SearchValues(a, b, 1e-6, A);
    EXPECT_LE(std::abs((int)values.Length() - 1), 1e-6);
    EXPECT_LE(std::abs(values[values.Length() - 1] - 1), 1e-6);

}