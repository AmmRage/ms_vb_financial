#include "gtest/gtest.h"
#include "src/financial.h"

using namespace  MS_VB_Financial;

TEST(Financial, IRR)
{
    std::vector<double> vect1{ -1, 1 };
    EXPECT_EQ(0, Financial::IRR(vect1, 0));
    std::vector<double> vect2{ -70000.0, 22000.0, 25000.0, 28000.0, 31000.0 };
    EXPECT_EQ(0.177435884421108, Financial::IRR(vect2, 0.1));
    std::vector<double> vect3{ -10000.0, 6000.0, -2000.0, 7000.0, 1000.0 };
    EXPECT_EQ(0.08667204742917164, Financial::IRR(vect3, 0.1));
    std::vector<double> vect4{ -30000.0, -10000.0, 25000.0, 12000.0, 15000.0 };
    EXPECT_EQ(0.10928101434575987, Financial::IRR(vect4, 0.1));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}