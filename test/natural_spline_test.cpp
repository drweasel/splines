#include "NaturalSpline.h"
#include <gtest/gtest.h>
#include <random>
#include <fstream>

TEST(natural_spline_test,NaturalSpline)
{
	std::random_device rd;
	std::mt19937_64 PRNG(rd());
	std::uniform_real_distribution<double> randu(0.,1.);
	std::normal_distribution<double> randn(0.,1.);

	int N = 100;
	std::unique_ptr< double[] > x = std::make_unique<double[]>(N);
	std::unique_ptr< double[] > y = std::make_unique<double[]>(N);

	std::ofstream samples("samples.csv");
	int k=0;
	double shift = randn(PRNG);
	double shift2 = randn(PRNG);
	for (double t : Lerp<double>(0.,2*3.14159265,N))
	{
		x[k] = t;
		y[k] = std::cos(2.*x[k]+shift) + std::sin(x[k]-shift2) + 0.1*randn(PRNG);

		samples << x[k] << ' ' << y[k] << std::endl;
		k++;
	}
	samples.close();

	auto spline = NaturalSpline<double>::CrossValidatedFit(
			x.get(),y.get(),N,1,3,false,20);
	fprintf(stderr,"obtained a spline with %i segments\n", spline.GetNumSegments());
	spline.Dump();

	for (int s=1; s<spline.GetNumSegments(); ++s)
	{
		printf("test segment %i for C0,C1 ...\n", s);
		EXPECT_TRUE(spline.IsContinuous(s,0,1e-6));
	}

	EXPECT_TRUE(bool(spline));

	double l = spline.Length(spline.GetAbscissa(0),spline.GetAbscissa(spline.GetNumSegments()),1e-6);
	fprintf(stderr,"length = %.16f\n", l);

	std::ofstream values("values.csv");
	for (double x : Lerp<double>(0.,2*3.14159265,200))
	{
		values << x << ' ' << spline(x) << std::endl;
	}
	values.close();

	fprintf(stderr,"S=dlmread('samples.csv'); P=dlmread('values.csv'); "
		"clf; plot(S(:,1),S(:,2),'r.'); hold on; plot(P(:,1),P(:,2),'b-');\n");

	//std::function< double(double) > f = std::bind(&XXX<double>::Eval,XXX<double>(3),std::placeholders::_1);
	//auto g = std::bind(&NaturalSpline<double>::GetAbscissa,spline,std::placeholders::_1);
	//fprintf(stderr,"f=%f g=%f\n", f(10.),g(1));

	double appInt, appErr;
	int nEvals;

	std::tie(appInt,appErr,nEvals) = ::Integrate(
		[&](double x)->double{ return spline(x); },
		spline.GetAbscissa(0),
		spline.GetAbscissa(spline.GetNumSegments()),
		1e-9
		);

	printf("integral: %.16f (exakt: %.16f) error: %g, #evals: %i\n", appInt,
		spline.Integrate(
			spline.GetAbscissa(0),
			spline.GetAbscissa(spline.GetNumSegments())
		),
		appErr, nEvals
		);
}

