#include "pch.h"
#include "Vector.h"

const double eps = 1e-6;

class VectorTest : public ::testing::Test {
protected:
    virtual void SetUp() {
        a = Vector(2);
        a[0] = 1;
        a[1] = -1;

        b = Vector(2);
        b[0] = 1;
        b[1] = 2;
    }

    Vector a;
    Vector b;
};

TEST_F(VectorTest, Addition) {
    auto summ = a + b;
    EXPECT_LE(std::abs(summ[0] - 2), eps);
    EXPECT_LE(std::abs(summ[1] - 1), eps);
}

TEST_F(VectorTest, Subtraction) {
    auto summ = a - b;
    EXPECT_LE(std::abs(summ[0] - 0), eps);
    EXPECT_LE(std::abs(summ[1] + 3), eps);
}

TEST_F(VectorTest, DotProduct) {
    double dotProduct = a * b;
    EXPECT_LE(std::abs(dotProduct + 1), eps);
}

TEST_F(VectorTest, ScalarProduct)
{
    Vector scalarProduct = 2 * a;
    EXPECT_LE(std::abs(scalarProduct[0] - 2), eps);
    EXPECT_LE(std::abs(scalarProduct[1] + 2), eps);
}

TEST_F(VectorTest, SwapElement)
{
    a.Swap(0, 1);
    EXPECT_LE(std::abs(a[0] + 1), eps);
    EXPECT_LE(std::abs(a[1] - 1), eps);
}

TEST_F(VectorTest, Norm)
{
    EXPECT_LE(std::abs(a.Norm() - std::sqrt(a[0] * a[0] + a[1] * a[1])), eps);
}

