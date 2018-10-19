#include "NaturalSpline.h"
#include <gtest/gtest.h>
#include <memory>
#include <cstdio>

TEST(zeros_test,NaturalSpline)
{
	auto x = std::make_unique<double[]>(10);
	auto y = std::make_unique<double[]>(10);

	int i=0;
	for (double xi : Lerp<double>(-5.,+5.,10))
		x[i++] = xi;

	y[0] = 0.0;
	y[1] = 1.0;
	y[2] = 1.5;
	y[3] = -2.0;
	y[4] = 0.5;
	y[5] = -1.0;
	y[6] = 1.;
	y[7] = 0.5;
	y[8] = 0.1;
	y[9] = 0.0;

	auto spline = NaturalSpline<double>::Interpolate(x.get(),y.get(),10,3);
	auto Dspline = spline.Diff(1);

	//for (double xi : Lerp<double>(-5.,+5.,100))
	//	printf("%g %g %g\n", xi, spline(xi), Dspline(xi));
	
	EXPECT_FALSE(spline.IsMonotonic());
	EXPECT_TRUE(spline.IsMonotonic(-5,-3.15));
	EXPECT_TRUE(spline.IsMonotonic(-3.1,-1.65));
	EXPECT_FALSE(spline.IsMonotonic(-5,-1.65));
	EXPECT_TRUE(spline.IsMonotonic(-1.55,-0.5));
	EXPECT_FALSE(spline.IsMonotonic(-1.55,0.5));

	double zeros[3];

	EXPECT_EQ(Dspline.RealRoots(0,zeros),0);
	EXPECT_EQ(Dspline.RealRoots(1,zeros),1);
	EXPECT_NEAR(zeros[0],-3.12389,1e-4);
	EXPECT_EQ(Dspline.RealRoots(2,zeros),0);
	EXPECT_EQ(Dspline.RealRoots(3,zeros),1);
	EXPECT_NEAR(zeros[0],-1.64096,1e-4);
	EXPECT_EQ(Dspline.RealRoots(4,zeros),1);
	EXPECT_NEAR(zeros[0],-0.468429,1e-4);
	EXPECT_EQ(Dspline.RealRoots(5,zeros),1);
	EXPECT_NEAR(zeros[0],0.57857,1e-4);
	EXPECT_EQ(Dspline.RealRoots(6,zeros),1);
	EXPECT_NEAR(zeros[0],1.94933,1e-4);
	EXPECT_EQ(Dspline.RealRoots(7,zeros),0);
	EXPECT_EQ(Dspline.RealRoots(8,zeros),0);

#if 0
	for (int k=0; k<spline.GetNumSegments(); ++k)
	{
		double zeros[3];
		int nz = Dspline.RealRoots(k,zeros);
		printf("segment %i: %i real roots:",k,nz);
		for (int j=0; j<nz; ++j)
			printf(" %g (->%g)", zeros[j], Dspline(zeros[j]));
		printf("\n");
	}
#endif
}

