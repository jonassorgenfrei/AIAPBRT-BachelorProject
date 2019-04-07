#include "tests/gtest/gtest.h"
#include "pbrt.h"
#include "geometry.h"

using namespace pbrt;

TEST(Bounds2, IteratorBasic) {
	Bounds2i b{ {0, 1}, {2, 3} };
	Point2i e[] = { {0, 1}, {1, 1}, {0, 2}, {1, 2} };
	int offset = 0;
	for (auto p : b) {
		EXPECT_LT(offset, sizeof(e) / sizeof(e[0]));
		EXPECT_EQ(e[offset], p) << "offset = " << offset;
		++offset;
	}
}