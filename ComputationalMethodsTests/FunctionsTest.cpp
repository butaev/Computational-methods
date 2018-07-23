#include "pch.h"
#include "Functions.h"

const double eps = 1e-6;

class FunctionTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        v = Vector(2);
        v[0] = 1, v[1] = -1;

        m = Matrix(2);
        m[0][0] = 1, m[0][1] = 2, m[1][0] = 0, m[1][1] = 2;
    }

    Matrix m;
    Vector v;
};

TEST_F(FunctionTest, QRDecomposition)
{
    Matrix q, r;
    Functions::QRDecomposition(m, q, r);

    EXPECT_LE(std::abs(q[0][0] - 1), eps);
    EXPECT_LE(std::abs(q[0][1] - 0), eps);
    EXPECT_LE(std::abs(q[1][0] - 0), eps);
    EXPECT_LE(std::abs(q[1][1] - 1), eps);

    EXPECT_LE(std::abs(r[0][0] - 1), eps);
    EXPECT_LE(std::abs(r[0][1] - 2), eps);
    EXPECT_LE(std::abs(r[1][0] - 0), eps);
    EXPECT_LE(std::abs(r[1][1] - 2), eps);
}

TEST_F(FunctionTest, NumLeftValues) {
    EXPECT_LE(std::abs(Functions::NumLeftValues(3, m) - 2), eps);
}

TEST_F(FunctionTest, QRAlgorithm) {
    Vector eigenvalues = Functions::QRAlgorithm(m, eps);
    EXPECT_LE(std::abs(eigenvalues[0] - 1), eps);
    EXPECT_LE(std::abs(eigenvalues[1] - 2), eps);
}

TEST_F(FunctionTest, BisectionMethod) {
    Vector eigenvalues = Functions::BisectionMethod(m, eps);
    EXPECT_LE(std::abs(eigenvalues[0] - 1), eps);
    EXPECT_LE(std::abs(eigenvalues[1] - 2), eps);
}

TEST_F(FunctionTest, RotationMethod)
{
    Vector b = m * v;
    Vector x = Functions::RotationMethod(m, b);
    EXPECT_LE(std::abs(v[0] - x[0]), eps);
    EXPECT_LE(std::abs(v[1] - x[1]), eps);
}

TEST_F(FunctionTest, GMRES)
{
    Vector b = m * v;
    int numIter = 0;
    Vector x = Functions::GMRES(m, b, eps, numIter);
    EXPECT_LE(std::abs(x[0] - v[0]), eps);
    EXPECT_LE(std::abs(x[1] - v[1]), eps);
}

TEST_F(FunctionTest, SimpleIteration)
{
    Vector b = m * v;
    int numIter = 0;
    Vector x = Functions::SimpleIteration(m, b, eps, numIter);
    EXPECT_LE(std::abs(v[0] - x[0]), eps);
    EXPECT_LE(std::abs(v[1] - x[1]), eps);
}

TEST_F(FunctionTest, TridiagMatrix)
{
    Vector b = m * v;
    Vector x = Functions::TridiagMatrix(m, b);
    EXPECT_LE(std::abs(v[0] - x[0]), eps);
    EXPECT_LE(std::abs(v[1] - x[1]), eps);
}

TEST_F(FunctionTest, Gauss)
{
    Vector b = m * v;
    double detm = 0;
    Vector x = Functions::Gauss(m, b, detm);
    EXPECT_LE(std::abs(v[0] - x[0]), eps);
    EXPECT_LE(std::abs(v[1] - x[1]), eps);
    EXPECT_LE(std::abs(detm - 2), eps);
}

TEST_F(FunctionTest, GaussS)
{
    Vector b = m * v;
    double detm = 0;
    Vector x = Functions::GaussS(m, b, detm);
    EXPECT_LE(std::abs(v[0] - x[0]), eps);
    EXPECT_LE(std::abs(v[1] - x[1]), eps);
    EXPECT_LE(std::abs(detm - 2), eps);
}

TEST_F(FunctionTest, SearchValues)
{
    double leftBound = 0.5;
    double rightBound = 1.5;
    Vector values = Functions::SearchValues(leftBound, rightBound, eps, m);
    EXPECT_EQ(values.Size(), 1);
    EXPECT_LE(std::abs(values[values.Size() - 1] - 1), eps);
}