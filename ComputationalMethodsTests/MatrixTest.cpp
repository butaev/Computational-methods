#include "pch.h"
#include "Matrix.h"

class MatrixTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        a = Matrix(2);
        a[0][0] = 1;
        a[1][1] = -1;

        b = Matrix(2);
        b[0][0] = 1;
        b[0][1] = 2;
        b[1][0] = 3;
        b[1][1] = 4;
    }

    // virtual void TearDown() {}
    Matrix a;
    Matrix b;
};

TEST_F(MatrixTest, Addition) {
    auto summ = a + b;
    EXPECT_LE(std::abs(summ[0][0] - 2), 1e-6);
    EXPECT_LE(std::abs(summ[0][1] - 2), 1e-6);
    EXPECT_LE(std::abs(summ[1][0] - 3), 1e-6);
    EXPECT_LE(std::abs(summ[1][1] - 3), 1e-6);
}

TEST_F(MatrixTest, Subtraction) {
    auto summ = a - b;
    EXPECT_LE(std::abs(summ[0][0] - 0), 1e-6);
    EXPECT_LE(std::abs(summ[0][1] + 2), 1e-6);
    EXPECT_LE(std::abs(summ[1][0] + 3), 1e-6);
    EXPECT_LE(std::abs(summ[1][1] + 5), 1e-6);
}

TEST_F(MatrixTest, DotProduct) {
    Matrix dotProduct = a * b;
    EXPECT_LE(std::abs(dotProduct[0][0] - 1), 1e-6);
    EXPECT_LE(std::abs(dotProduct[0][1] - 2), 1e-6);
    EXPECT_LE(std::abs(dotProduct[1][0] + 3), 1e-6);
    EXPECT_LE(std::abs(dotProduct[1][1] + 4), 1e-6);
}

TEST_F(MatrixTest, MatrixVectorProduct)
{
    Vector v(2);
    v[0] = 1;
    v[1] = -1;
    Vector dotProduct = a * v;
    EXPECT_LE(std::abs(dotProduct[0] - 1), 1e-6);
    EXPECT_LE(std::abs(dotProduct[1] - 1), 1e-6);
}

TEST_F(MatrixTest, dotScalarProduct)
{
    Matrix scalarProduct = 2 * b;
    EXPECT_LE(std::abs(scalarProduct[0][0] - 2), 1e-6);
    EXPECT_LE(std::abs(scalarProduct[0][1] - 4), 1e-6);
    EXPECT_LE(std::abs(scalarProduct[1][0] - 6), 1e-6);
    EXPECT_LE(std::abs(scalarProduct[1][1] - 8), 1e-6);
}

TEST_F(MatrixTest, SwapLine)
{
    b.SwapLine(0, 1);
    EXPECT_LE(std::abs(b[0][0] - 3), 1e-6);
    EXPECT_LE(std::abs(b[0][1] - 4), 1e-6);
    EXPECT_LE(std::abs(b[1][0] - 1), 1e-6);
    EXPECT_LE(std::abs(b[1][1] - 2), 1e-6);
}

TEST_F(MatrixTest, Identity)
{
    b.Identity();
    EXPECT_LE(std::abs(b[0][0] - 1), 1e-6);
    EXPECT_LE(std::abs(b[0][1] - 0), 1e-6);
    EXPECT_LE(std::abs(b[1][0] - 0), 1e-6);
    EXPECT_LE(std::abs(b[1][1] - 1), 1e-6);
}

TEST_F(MatrixTest, Inverse)
{
    a.Inverse();
    EXPECT_LE(std::abs(a[0][0] - 1), 1e-6);
    EXPECT_LE(std::abs(a[0][1] - 0), 1e-6);
    EXPECT_LE(std::abs(a[1][0] - 0), 1e-6);
    EXPECT_LE(std::abs(a[1][1] + 1), 1e-6);
}

TEST_F(MatrixTest, MaxValueSymmetricPositiveMatrix)
{
    Matrix m(2);
    m[0][0] = 3;
    m[0][1] = 0;
    m[1][0] = 0;
    m[1][1] = 5;
    double maxValue = m.MaxValue(1e-6);
    EXPECT_LE(std::abs(maxValue - 5), 1e-6);
}

TEST_F(MatrixTest, MinValueSymmetricPositiveMatrix)
{
    Matrix m(2);
    m[0][0] = 3;
    m[0][1] = 0;
    m[1][0] = 0;
    m[1][1] = 5;
    double minValue = m.MinValue(1e-6);
    EXPECT_LE(std::abs(minValue - 3), 1e-6);
}


TEST_F(MatrixTest, SpectralNorm)
{
    double norm = b.Norm();
    double trueNorm = std::sqrt(b[0][0] * b[0][0] + b[0][1] * b[0][1] + b[1][0] * b[1][0] + b[1][1] * b[1][1]);
    EXPECT_LE(std::abs(norm - trueNorm), 1e-6);
}

TEST_F(MatrixTest, 1Norm)
{
    double norm1 = b.Norm1();
    double trueNorm1 = std::max(std::abs(b[0][0]) + std::abs(b[0][1]), std::abs(b[1][0]) + std::abs(b[1][1]));
    EXPECT_LE(std::abs(norm1 - trueNorm1), 1e-6);
}

TEST_F(MatrixTest, ConditionNumber)
{
    double conditionNumber = b.ConditionNumber();
    double trueMatrixNorm = std::sqrt(b[0][0] * b[0][0] + b[0][1] * b[0][1] + b[1][0] * b[1][0] + b[1][1] * b[1][1]);
    double trueInverseMatrixNorm = std::sqrt(4 + 0.25 + 1 + 2.25);
    double trueConditionNumber = trueMatrixNorm * trueInverseMatrixNorm;
    EXPECT_LE(std::abs(conditionNumber - trueConditionNumber), 1e-6);
}

TEST_F(MatrixTest, ValuesBounds)
{
    Matrix M(2);
    M[0][0] = 2, M[0][1] = -1, M[1][0] = -1, M[1][1] = 2;
    //std::cout << "M = " << M << std::endl;
    Vector valuesBounds(2);
    valuesBounds = M.ValuesBounds();

    //std::cout << "valuesBounds = " << valuesBounds << std::endl;

    EXPECT_LE(std::abs(valuesBounds[0] - 1), 1e-6);
    EXPECT_LE(std::abs(valuesBounds[1] - 3), 1e-6);
}