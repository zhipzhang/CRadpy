#include "Utils.h"
#include "Eigen/Core"
#include "gtest/gtest.h"
#include "constants.h"




using CRad::Utils;

testing::AssertionResult CheckVectorPostive(const Eigen::VectorXd& vec)
{
    for(const auto& i: vec)
    {
        if( i < 0)
        {
            return testing::AssertionFailure();
        }
    }
    return testing::AssertionSuccess();
}

class UtilsTest: public testing::Test{
    protected:
        UtilsTest()
        {
            x = Utils::GetLogspaceVec(1.0, 1000, 30);
            y = Utils::GetLogspaceVec(1.0, 1000, 30);
            x2 = pow(x.array(),2).matrix();
        }
        Eigen::VectorXd x;
        Eigen::VectorXd y;
        Eigen::VectorXd x2;
};


TEST_F(UtilsTest, MethodGetLogspaceVec)
{
    Eigen::VectorXd res(30);
    res <<  1.        ,    1.268961  ,    1.61026203,    2.04335972,
          2.5929438 ,    3.29034456,    4.17531894,    5.29831691,
          6.72335754,    8.53167852,   10.82636734,   13.73823796,
         17.43328822,   22.12216291,   28.07216204,   35.6224789 ,
         45.20353656,   57.3615251 ,   72.78953844,   92.36708572,
        117.21022975,  148.73521073,  188.73918221,  239.502662  ,
        303.91953823,  385.66204212,  489.39009185,  621.01694189,
        788.04628157, 1000.;
    std::cout << "x is " <<x << std::endl;
    for (int i = 0; i < x.size(); i++)
    {
        EXPECT_NEAR(res[i], x[i], 1e-5) << "Difference at index" << i;
    }
}
TEST_F(UtilsTest, MethodBlackBodyPositive)
{
    double tempreture = 2.7;
    auto max_energy  = tempreture * constants::kb * 1e6;
    auto min_energy  = tempreture * constants::kb * 1e-12;
    auto energy = Utils::GetLogspaceVec(min_energy, max_energy, 200);
    auto density = Utils::BlackBody(energy, tempreture);

     EXPECT_TRUE(CheckVectorPostive(density)) ;
}

TEST_F(UtilsTest, MethodInterpolate)
{
    double x0 = x[10];
    double y2_interpolate = Utils::interpolate(x, x2, x0);

    EXPECT_NEAR(y2_interpolate , x2[10],  1e-5);
}

TEST_F(UtilsTest, MethodIntegrate)
{
    auto linear_func = [](double x)->double {return x;};
    double result = Utils::integrate_1d(linear_func, 0, 1);
    EXPECT_NEAR(result, 0.5, 1e-5);
}

