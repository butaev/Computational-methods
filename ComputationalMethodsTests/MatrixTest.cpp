#include "pch.h"
#include "Matrix.h"

const double eps = 1e-6;

class MatrixTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        a = { {1, 0}, {0, -1} };
        b = { {1, 2}, {3, 4} };
    }
    
    Matrix a;
    Matrix b;
};

TEST_F(MatrixTest, Addition) {
    auto summ = a + b;
    EXPECT_LE(std::abs(summ[0][0] - 2), eps);
    EXPECT_LE(std::abs(summ[0][1] - 2), eps);
    EXPECT_LE(std::abs(summ[1][0] - 3), eps);
    EXPECT_LE(std::abs(summ[1][1] - 3), eps);
}

TEST_F(MatrixTest, Subtraction) {
    auto summ = a - b;
    EXPECT_LE(std::abs(summ[0][0] - 0), eps);
    EXPECT_LE(std::abs(summ[0][1] + 2), eps);
    EXPECT_LE(std::abs(summ[1][0] + 3), eps);
    EXPECT_LE(std::abs(summ[1][1] + 5), eps);
}

TEST_F(MatrixTest, DotProduct) {
    Matrix dotProduct = a * b;
    EXPECT_LE(std::abs(dotProduct[0][0] - 1), eps);
    EXPECT_LE(std::abs(dotProduct[0][1] - 2), eps);
    EXPECT_LE(std::abs(dotProduct[1][0] + 3), eps);
    EXPECT_LE(std::abs(dotProduct[1][1] + 4), eps);
}

TEST_F(MatrixTest, MatrixVectorProduct)
{
    Vector v = { 1, -1 };
    Vector dotProduct = a * v;
    EXPECT_LE(std::abs(dotProduct[0] - 1), eps);
    EXPECT_LE(std::abs(dotProduct[1] - 1), eps);
}

TEST_F(MatrixTest, dotScalarProduct)
{
    Matrix scalarProduct = 2 * b;
    EXPECT_LE(std::abs(scalarProduct[0][0] - 2), eps);
    EXPECT_LE(std::abs(scalarProduct[0][1] - 4), eps);
    EXPECT_LE(std::abs(scalarProduct[1][0] - 6), eps);
    EXPECT_LE(std::abs(scalarProduct[1][1] - 8), eps);
}

TEST_F(MatrixTest, SwapLine)
{
    b.SwapLine(0, 1);
    EXPECT_LE(std::abs(b[0][0] - 3), eps);
    EXPECT_LE(std::abs(b[0][1] - 4), eps);
    EXPECT_LE(std::abs(b[1][0] - 1), eps);
    EXPECT_LE(std::abs(b[1][1] - 2), eps);
}

TEST_F(MatrixTest, Identity)
{
    b.Identity();
    EXPECT_LE(std::abs(b[0][0] - 1), eps);
    EXPECT_LE(std::abs(b[0][1] - 0), eps);
    EXPECT_LE(std::abs(b[1][0] - 0), eps);
    EXPECT_LE(std::abs(b[1][1] - 1), eps);
}

TEST_F(MatrixTest, Inverse)
{
    a.Inverse();
    EXPECT_LE(std::abs(a[0][0] - 1), eps);
    EXPECT_LE(std::abs(a[0][1] - 0), eps);
    EXPECT_LE(std::abs(a[1][0] - 0), eps);
    EXPECT_LE(std::abs(a[1][1] + 1), eps);
}

TEST_F(MatrixTest, MaxValueSymmetricPositiveMatrix)
{
    Matrix m = { {3, 0}, {0, 5} };
    EXPECT_LE(std::abs(m.MaxValue(eps) - 5), eps);
}

TEST_F(MatrixTest, MinValueSymmetricPositiveMatrix)
{
    Matrix m = { { 3, 0 },{ 0, 5 } };
    EXPECT_LE(std::abs(m.MinValue(eps) - 3), eps);
}


TEST_F(MatrixTest, SpectralNorm)
{
    double norm = std::sqrt(b[0][0] * b[0][0] + b[0][1] * b[0][1] + b[1][0] * b[1][0] + b[1][1] * b[1][1]);
    EXPECT_LE(std::abs(norm - b.Norm()), eps);
}

TEST_F(MatrixTest, 1Norm)
{
    double norm1 = std::max(std::abs(b[0][0]) + std::abs(b[0][1]), std::abs(b[1][0]) + std::abs(b[1][1]));
    EXPECT_LE(std::abs(norm1 - b.Norm1()), eps);
}

TEST_F(MatrixTest, ConditionNumber)
{
    double matrixNorm = std::sqrt(b[0][0] * b[0][0] + b[0][1] * b[0][1] + b[1][0] * b[1][0] + b[1][1] * b[1][1]);
    double inverseMatrixNorm = std::sqrt(4 + 0.25 + 1 + 2.25);
    double conditionNumber = matrixNorm * inverseMatrixNorm;
    EXPECT_LE(std::abs(conditionNumber - b.ConditionNumber()), eps);
}

TEST_F(MatrixTest, ValuesBounds)
{
    Matrix m = { {2, -1}, {-1, 2} };
    Vector valuesBounds = m.ValuesBounds();

    EXPECT_LE(std::abs(valuesBounds[0] - 1), eps);
    EXPECT_LE(std::abs(valuesBounds[1] - 3), eps);
}