#include "pch.h"
#include "Vector.h"

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

	// virtual void TearDown() {}
	Vector a;
	Vector b;
};

TEST_F(VectorTest, Addition) {
	auto summ = a + b;
	EXPECT_LE(std::abs(summ[0] - 2), 1e-6);
	EXPECT_LE(std::abs(summ[1] - 1), 1e-6);
}

TEST_F(VectorTest, DotProduct) {
	double dotProduct = a * b;
	EXPECT_LE(std::abs(dotProduct + 1), 1e-6);
}