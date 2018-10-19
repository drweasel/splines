#include "QR.h"
#include "Polynomials.h"
#include <gtest/gtest.h>
#include <cstdio>
#include <complex>

TEST(polynomials_test,Polynomials)
{
	double H[] = {
1, 9, 3, 6, 5,
7, 2, 2, 5, 3,
0, 4, 3, 1, 2,
0, 0, 8, 7,-7,
0, 0, 0, 9, 2
};

	std::complex<double> W[5];
	HessenbergEig(H,5,W);
	for (int i=0; i<5; ++i)
	{
		printf("%g %+g i\n", W[i].real(), W[i].imag());
	}

	printf("roots:\n");
	double a[] = { 5, 8, 42, 3, 1, 2, 3.5, 99 };
	constexpr int N = (int)(sizeof(a)/sizeof(a[0]));
	std::complex<double> w[N-1];
	double work[(N-1)*(N-1)];

	PolynomialRoots(a,N,w,work);
	for (int i=0; i<N-1; ++i)
	{
		printf("x = %g %+g i\n", w[i].real(), w[i].imag());

		std::complex<double> x = w[i], X(1.,0.), y(0.,0.);
		for (int j=0; j<N; ++j)
		{
			y += a[j]*X;
			X *= x;
		}
		printf("\t->  y(x) = %g %+g i\n", y.real(), y.imag());
		EXPECT_LT(std::abs(y),1e-6);
	}
}

