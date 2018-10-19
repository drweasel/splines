#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "Curve.h"
#include "LeastSquares.h"
#include "Lerp.h"
#include "LinearPartition.h"
#include <memory>
#include <random>

/**
 * Natural Spline.
 *
 * @author Michael Weitzel <mich@elweitzel.de>
 */
template< typename T > class NaturalSpline
	: public Curve<T>
{
private:
	static void computeCoeffs1(
		int n,
		T (*C)[4],
		const T* x,
		const T* y
		);
	static void computeCoeffs2(
		int n,
		T (*C)[4],
		const T* x,
		const T* y
		);
	static void computeCoeffs3(
		int n,
		T (*C)[4],
		const T* x,
		const T* y
		);

	static NaturalSpline fit1(
		const T* x,
		const T* y,
		unsigned int N,
		const T* xs,
		unsigned int S,
		bool periodic
		);
	static NaturalSpline fit2(
		const T* x,
		const T* y,
		unsigned int N,
		const T* xs,
		unsigned int S,
		unsigned int C,
		bool periodic
		);
	static NaturalSpline fit3(
		const T* x,
		const T* y,
		unsigned int N,
		const T* xs,
		unsigned int S,
		unsigned int C,
		bool periodic
		);

	// distance of r to line segment (p,q)
	static T dist_to_line(
		const T& px, const T& py,
		const T& qx, const T& qy,
		const T& rx, const T& ry
		)
	{
		T dx = qx-px, dy = qy-py;
		T len = std::sqrt(dx*dx + dy*dy);
		return std::abs(dx*(py-ry) - (px-rx)*dy) / len;
	}

	// Ramer-Douglas-Peucker algorithm (TODO: testing!)
	static unsigned int approximate(
		const T* x,
		const T* y,
		bool* keep,
		unsigned int N,
		T eps
		)
	{
		if (N < 3)
			return 0;
		T px = x[0], py = y[0];
		T qx = x[N-1], qy = y[N-1];

		unsigned int max_k = 1;
		T max_d(0);
		for (unsigned int k=1; k<N-1; ++k)
		{
			T d = dist_to_line(px,py,qx,qy,x[k],y[k]);
			if (d > max_d) { max_k = k; max_d = d; }
			keep[k] = false;
		}
		if (max_d > eps)
		{
			keep[max_k] = true;
			unsigned int c_l = approximate(x,y,keep,max_k+1,eps);
			unsigned int c_r = approximate(x+max_k,y+max_k,keep+max_k,N-max_k,eps);
			return c_r+c_r;
		}
		return 0;
	}

	// Frontend for the RDP algorithm (TODO: testing!)
	static std::tuple< std::unique_ptr< T[] >,unsigned int > approximate(
		const T* x,
		const T* y,
		unsigned int N,
		T eps
		)
	{
		bool* keep = new bool[N];
		keep[0] = keep[N-1] = true;
		unsigned int nk;
		nk = approximate(x,y,keep,N,eps) + 2;
		auto xs = std::make_unique<T[]>(nk);
		for (unsigned int k=0,j=0; j<nk; ++k)
			if (keep[k])
				xs[j++] = x[k];
		return {std::move(xs),nk};
	}

	NaturalSpline(int n, const T* X, const T C[][4])
		: Curve<T>(n,X,C) { }

public:
	NaturalSpline()
		: Curve<T>() { }

	NaturalSpline(
		const NaturalSpline& cpy
		); // untested

	NaturalSpline(
		NaturalSpline&& cpy
		);

	~NaturalSpline();

	NaturalSpline& operator=(const NaturalSpline& cpy);
	NaturalSpline& operator=(NaturalSpline&& cpy);

	NaturalSpline Diff(
		unsigned int order = 1
		) const;

	bool IsMonotonic() const;
	bool IsMonotonic(
		const T& lo,
		const T& hi
		) const;

	/**
	 * Interpolation.
	 *
	 * @param	x	data (x-values)
	 * @param	y	data (y-values)
	 * @param	N	number of (x,y) value pairs
	 * @param	D	polynomial degree
	 */
	static NaturalSpline Interpolate(
		const T* x,
		const T* y,
		int N,
		int D
		);

	/**
	 * Regression Spline.
	 *
	 * @param	x	data (x-values)
	 * @param	y	data (y-values)
	 * @param	N	number of (x,y) value pairs
	 * @param	xs	segment start points (S)
	 * @param	S	number of segments (w/ start points xs)
	 * @param	C	continuity of segment transitions (C0,C1,C2,C3)
	 * @param	D	polynomial degree of segments
	 * @param	periodic	flag; periodic spline?
	 */
	static NaturalSpline Fit(
		const T* x,
		const T* y,
		unsigned int N,
		const T* xs,
		unsigned int S,
		unsigned int C,
		unsigned int D,
		bool periodic
		);

	/**
	 * Cross-Validated Regression Spline.
	 *
	 * This method uses cross validation to determine the optimal
	 * number of spline segments.
	 *
	 * @param	x	data (x-values)
	 * @param	y	data (y-values)
	 * @param	N	number of (x,y) value pairs
	 * @param	C	continuity of segment transitions (C0,C1,C2,C3)
	 * @param	D	polynomial degree of segments
	 * @param	periodic	flag; periodic spline?
	 * @param	Smax	upper limit for number of segments (-1 = auto)
	 */
	static NaturalSpline CrossValidatedFit(
		const T* x,
		const T* y,
		unsigned int N,
		unsigned int C,
		unsigned int D,
		bool periodic,
		unsigned int Smax = -1
		);

	void Dump() const;

}; // class NaturalSpline

// Template implementations

template<typename T> void NaturalSpline<T>::computeCoeffs1(
	int n,
	T (*C)[4],
	const T* x,
	const T* y
	)
{
	// s_j(x) = a_j + b_j(x-x_j)
	//
	// a_j = y_j
	// b_j = (y_{j+1}-y_j)/(x_{j+1}-x_j)
	//
	for (int j=0; j<n; ++j)
	{
		C[j][0] = y[j]; // a_j
		C[j][1] = (y[j+1]-y[j])/(x[j+1]-x[j]); // b_j
		C[j][2] = 0; // c_j
		C[j][3] = 0; // d_j
	}
}

template<typename T> void NaturalSpline<T>::computeCoeffs2(
	int n,
	T (*C)[4],
	const T* x,
	const T* y
	)
{
	// s_j(x) = a_j + b_j(x-x_j) + c_j(x-x_j)^2
	//
	// a_j = y_j;
	// b_j = z_j;
	// c_j = (z_{j+1}-z_j)/(2(x_{j+1}-x_j))
	// d_j = 0
	//
	// z_{j+1} = -z_j+2(y_{j+1}-y_j)/(x_{j+1}-x_j)
	//     z_0 = 0 -- or something else

	// z_jp1 is the slope of the spline at x[0]; different options:
	// z_jp1 = 0 => "natural spline"
	// z_jp1 = (y[1]-y[0])/(x[1]-x[0]) => slope of linear spline (degree 1)
	// z_jp1 = ds_0(x[0])/dx => slope of the first parabola at x[0]; preferred
	T x0sq = x[0]*x[0], x1sq = x[1]*x[1], x2sq = x[2]*x[2];
	T det = (x[0]*(x2sq-x1sq)-x[1]*x2sq+x1sq*x[2]+x0sq*(x[1]-x[2]));
	T a1 = -(x0sq*(y[2]-y[1])-x1sq*y[2]+x2sq*y[1]+(x1sq-x2sq)*y[0])/det;
	T a2 = (x[0]*(y[2]-y[1])-x[1]*y[2]+x[2]*y[1]+(x[1]-x[2])*y[0])/det;

	T z_jp1 = 2*a2*x[0] + a1; // = d s_0(x[0]) / d x

	for (int j=0; j<n; ++j)
	{
		T z_j = z_jp1;
		z_jp1 = -z_j + 2*(y[j+1]-y[j])/(x[j+1]-x[j]);

		C[j][0] = y[j];
		C[j][1] = z_j;
		C[j][2] = (z_jp1-z_j)/(2*(x[j+1]-x[j]));
		C[j][3] = 0; // cubic coeffcient is zero for degree 2
	}
}

template<typename T> void NaturalSpline<T>::computeCoeffs3(
	int n,
	T (*C)[4],
	const T* x,
	const T* y
	)
{
	int i, j;
	T* g, *h, *xc;
	g = new T[n-1];
	h = new T[n];
	xc = new T[n-1];

	for (j=0; j<n; ++j)
	{
		C[j][0] = y[j];
		h[j] = x[j+1] - x[j];
	}

	for (j=2; j<=n; ++j)
		g[j-2] = 3*((y[j]-y[j-1])/h[j-1] - (y[j-1]-y[j-2])/h[j-2]);

	// "Thomas" algorithm for the tridiagonal system of equations
	// a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i
	// a := lower secondary diagonal
	// b := main diagonal
	// c := upper secondary diagonal
	//
	// | b1 c1  0 ...      |
	// | a2 b2 c2  0 ...   |
	// |  0 a3 b3 c3  0 ...|
	// |        ...        |
	// |          ...c(n-1)|
	// |  0 ...   am     bm|
	//
	// a(i) = h[i], b(i) := 2*(h[i]+h[i+1]), c(i) := h[i+1], m := n-1
	//
	// Forward elimination
	xc[0] = 2*(h[0]+h[1]);
	for (i=1; i<n-1; ++i)
	{
		T m = h[i] / xc[i-1];
		xc[i] = 2*(h[i]+h[i+1]) - m * h[i];
		g[i] -= m * g[i-1];
	}
	// Back substitution
	xc[n-2] = g[n-2] / xc[n-2];
	for (i=n-3; i>=0; --i)
		xc[i] = (g[i]-h[i+1]*xc[i+1]) / xc[i];

	// copying of coeffcients
	C[0][2] = 0;
	for (i=1; i<n; ++i)
		C[i][2] = xc[i-1];

	for (j=1; j<n; ++j)
	{
		C[j-1][1] = (y[j]-y[j-1])/h[j-1] - h[j-1]*(2*C[j-1][2]+C[j][2])/3;
		C[j-1][3] = (C[j][2]-C[j-1][2])/(3*h[j-1]);
	}
	C[n-1][1] = (y[n]-y[n-1])/h[n-1] - h[n-1]*(2*C[n-1][2])/3;
	C[n-1][3] = -C[n-1][2]/(3*h[n-1]);

	delete[] g;
	delete[] h;
	delete[] xc;
}

template<typename T> NaturalSpline<T> NaturalSpline<T>::fit1(
	const T* x,
	const T* y,
	unsigned int N,
	const T* xs,
	unsigned int S,
	bool periodic
	)
{
	// number of unknowns
	const int W = 2*S;
	// periodic?
	const int p = periodic ? 1 : 0;

	if (!(N >= W - (S-1+p)))
	{
		// to few samples for given number of segments
		return NaturalSpline();
	}

	// overdetermined system A.c=b
	T* A = new T[N*W]();
	T* b = new T[W>N ? W : N](); // rhs, solution
	// equality constraints (continuity) B.c=d
	T* B = new T[(S-1+p)*W]();
	T* d = new T[S-1+p]();

	for (unsigned int s=0,i=0; s<S; ++s)
	{
		T xl = xs[s];
		T xh = xs[s+1];
		unsigned int j = 2*s;
		for (; i<N && x[i]<xh; ++i)
		{
			T xm = x[i] - xl;
			A[i*W+j+0] = xm;
			A[i*W+j+1] = T(1);
			b[i] = y[i];
		}

		// equality constraints
		if (s+1 < S)
		{
			// segment transition s -> s+1
			unsigned int k = s;
			xl = xs[s+1] - xs[s];
			xh = T(0); // =xs[s+1] - xs[s+1];
			// C0 continuity
			B[k*W+j+0] = xl;
			B[k*W+j+1] = T(1);
			B[k*W+j+2] = -(xh);
			B[k*W+j+3] = -T(1);
			d[k] = T(0);
		} // if (s+1 < S)
	}

	if (periodic)
	{
		unsigned int k = S-1;
		unsigned int j = W-2;
		T xl = x[N-1] - xs[S-1];
		T xh = T(0);
		// C0 continuity
		B[k*W+j+0] = xl;
		B[k*W+j+1] = T(1);
		B[k*W+0] = -(xh);
		B[k*W+1] = -T(1);
		d[k] = T(0);
	}

	// solve overdetermined system A.x=b subject to equality constraints B.x=d
	if (!LSESolve(A,b,B,d,N,W,S-1+p))
	{
		delete[] A; delete[] b;
		delete[] B; delete[] d;
		return NaturalSpline();
	}

	T (*c)[4] = new T[S][4]();
	T* X = new T[S+1]();
	for (unsigned int i=0; i<S; ++i)
	{
		X[i] = xs[i];
		for (unsigned int j=0; j<2; ++j)
			c[i][j] = b[2*i+(1-j)];
		c[i][2] = T(0);
		c[i][3] = T(0);
	}
	X[S] = x[N-1];

	delete[] A; delete[] b;
	delete[] B; delete[] d;
	return NaturalSpline(S,X,c);
}

template<typename T> NaturalSpline<T> NaturalSpline<T>::fit2(
	const T* x,
	const T* y,
	unsigned int N,
	const T* xs,
	unsigned int S,
	unsigned int C,
	bool periodic
	)
{
	if (C>2)
	{
		// invalid continuity
		return NaturalSpline();
	}

	// number of unknowns
	const int W = 3*S;
	// periodic?
	const int p = periodic ? 1 : 0;

	if (!(N >= W - (C+1)*(S-1+p)))
	{
		// too few samples for given number of segments / continuity
		return NaturalSpline();
	}

	// overdetermined system A.c=b
	T* A = new T[N*W]();
	T* b = new T[W>N ? W : N](); // rhs, Lösung
	// equality constraints (continuity) B.c=d
	T* B = new T[((C+1)*(S-1+p))*W]();
	T* d = new T[(C+1)*(S-1+p)]();

	for (unsigned int s=0,i=0; s<S; ++s)
	{
		T xl = xs[s];
		T xh = xs[s+1];
		unsigned int j = 3*s;
		for (; i<N && x[i]<xh; ++i)
		{
			T xm = x[i] - xl;
			A[i*W+j+0] = xm*xm;
			A[i*W+j+1] = xm;
			A[i*W+j+2] = T(1);
			b[i] = y[i];
		}

		// equality constraints
		if (s+1 < S)
		{
			// segment transition s -> s+1
			unsigned int k = (C+1)*s;
			xl = xs[s+1] - xs[s];
			xh = T(0); // =xs[s+1] - xs[s+1];
			switch (C)
			{
			case 2:
				// C2 continuity
				B[k*W+j+0] = T(2);
				B[k*W+j+1] = B[k*W+j+2] = T(0);
				B[k*W+j+3] = T(-2);
				B[k*W+j+4] = B[k*W+j+5] = T(0);
				d[k] = T(0);
				k++;
			case 1:
				// C1 continuity
				B[k*W+j+0] = T(2)*(xl);
				B[k*W+j+1] = T(1);
				B[k*W+j+2] = T(0);
				B[k*W+j+3] = T(-2)*(xh);
				B[k*W+j+4] = T(-1);
				B[k*W+j+5] = T(0);
				d[k] = T(0);
				k++;
			case 0:
				// C0 continuity
				B[k*W+j+0] = xl*xl;
				B[k*W+j+1] = xl;
				B[k*W+j+2] = T(1);
				B[k*W+j+3] = -(xh*xh);
				B[k*W+j+4] = -(xh);
				B[k*W+j+5] = -T(1);
				d[k] = T(0);
			}
		} // if (s+1 < S)
	}

	if (periodic)
	{
		unsigned int k = (C+1)*(S-1);
		unsigned int j = W-3;
		T xl = x[N-1] - xs[S-1];
		T xh = T(0);
		switch (C)
		{
		case 2:
			// C2 continuity
			B[k*W+j+0] = T(2);
			B[k*W+j+1] = B[k*W+j+2] = T(0);
			B[k*W+0] = T(-2);
			B[k*W+1] = B[k*W+2] = T(0);
			d[k] = T(0);
			k++;
		case 1:
			// C1 continuity
			B[k*W+j+0] = T(2)*(xl);
			B[k*W+j+1] = T(1);
			B[k*W+j+2] = T(0);
			B[k*W+0] = T(-2)*(xh);
			B[k*W+1] = T(-1);
			B[k*W+2] = T(0);
			d[k] = T(0);
			k++;
		case 0:
			// C0 continuity
			B[k*W+j+0] = xl*xl;
			B[k*W+j+1] = xl;
			B[k*W+j+2] = T(1);
			B[k*W+0] = -(xh*xh);
			B[k*W+1] = -(xh);
			B[k*W+2] = -T(1);
			d[k] = T(0);
		}
	}

	// solve overdetermined system A.x=b subject to equality constraints B.x=d
	if (!LSESolve(A,b,B,d,N,W,(C+1)*(S-1+p)))
	{
		delete[] A; delete[] b;
		delete[] B; delete[] d;
		return NaturalSpline();
	}

	T (*c)[4] = new T[S][4]();
	T* X = new T[S+1]();
	for (unsigned int i=0; i<S; ++i)
	{
		X[i] = xs[i];
		for (unsigned int j=0; j<3; ++j)
			c[i][j] = b[3*i+(2-j)];
		c[i][3] = T(0);
	}
	X[S] = x[N-1];

	delete[] A; delete[] b;
	delete[] B; delete[] d;
	return NaturalSpline(S,X,c);
}

template<typename T> NaturalSpline<T> NaturalSpline<T>::fit3(
	const T* x,
	const T* y,
	unsigned int N,
	const T* xs,
	unsigned int S,
	unsigned int C,
	bool periodic
	)
{
	if (C>3)
	{
		// invalid continuity
		return NaturalSpline();
	}

	// number of unknowns
	const int W = 4*S;
	// periodic?
	const int p = periodic ? 1 : 0;

	if (!(N >= W - (C+1)*(S-1+p)))
	{
		// too few samples for given number of segments / continuity
		return NaturalSpline();
	}

	// overdetermined system A.c=b
	T* A = new T[N*W]();
	T* b = new T[W>N ? W : N](); // rhs, solution
	// equality constraints (continuity) B.c=d
	T* B = new T[((C+1)*(S-1+p))*W]();
	T* d = new T[(C+1)*(S-1+p)]();

	for (unsigned int s=0,i=0; s<S; ++s)
	{
		T xl = xs[s];
		T xh = xs[s+1];
		unsigned int j = 4*s;
		for (; i<N && x[i]<xh; ++i)
		{
			T xm = x[i] - xl;
			A[i*W+j+0] = xm*xm*xm;
			A[i*W+j+1] = xm*xm;
			A[i*W+j+2] = xm;
			A[i*W+j+3] = T(1);
			b[i] = y[i];
		}

		// equality constraints
		if (s+1 < S)
		{
			// segment transition s -> s+1
			unsigned int k = (C+1)*s;
			xl = xs[s+1] - xs[s];
			xh = T(0); // =xs[s+1] - xs[s+1];
			switch (C)
			{
			case 3:
				// C3 continuity
				B[k*W+j+0] = T(6);
				B[k*W+j+1] = B[k*W+j+2] = B[k*W+j+3] = T(0);
				B[k*W+j+4] = T(-6);
				B[k*W+j+5] = B[k*W+j+6] = B[k*W+j+7] = T(0);
				d[k] = T(0);
				k++;
			case 2:
				// C2 continuity
				B[k*W+j+0] = T(6)*(xl);
				B[k*W+j+1] = T(2);
				B[k*W+j+2] = B[k*W+j+3] = T(0);
				B[k*W+j+4] = T(-6)*(xh);
				B[k*W+j+5] = T(-2);
				B[k*W+j+6] = B[k*W+j+7] = T(0);
				d[k] = T(0);
				k++;
			case 1:
				// C1 continuity
				B[k*W+j+0] = T(3)*(xl*xl);
				B[k*W+j+1] = T(2)*(xl);
				B[k*W+j+2] = T(1);
				B[k*W+j+3] = T(0);
				B[k*W+j+4] = T(-3)*(xh*xh);
				B[k*W+j+5] = T(-2)*(xh);
				B[k*W+j+6] = T(-1);
				B[k*W+j+7] = T(0);
				d[k] = T(0);
				k++;
			case 0:
				// C0 continuity
				B[k*W+j+0] = xl*xl*xl;
				B[k*W+j+1] = xl*xl;
				B[k*W+j+2] = xl;
				B[k*W+j+3] = T(1);
				B[k*W+j+4] = -(xh*xh*xh);
				B[k*W+j+5] = -(xh*xh);
				B[k*W+j+6] = -(xh);
				B[k*W+j+7] = -T(1);
				d[k] = T(0);
			}
		} // if (s+1 < S)
	}

	if (periodic)
	{
		unsigned int k = (C+1)*(S-1);
		unsigned int j = W-4;
		T xl = x[N-1] - xs[S-1];
		T xh = T(0);
		switch (C)
		{
		case 3:
			// C3 continuity
			B[k*W+j+0] = T(6);
			B[k*W+j+1] = B[k*W+j+2] = B[k*W+j+3] = T(0);
			B[k*W+0] = T(-6);
			B[k*W+1] = B[k*W+2] = B[k*W+3] = T(0);
			d[k] = T(0);
			k++;
		case 2:
			// C2 continuity
			B[k*W+j+0] = T(6)*(xl);
			B[k*W+j+1] = T(2);
			B[k*W+j+2] = B[k*W+j+3] = T(0);
			B[k*W+0] = T(-6)*(xh);
			B[k*W+1] = T(-2);
			B[k*W+2] = B[k*W+3] = T(0);
			d[k] = T(0);
			k++;
		case 1:
			// C1 continuity
			B[k*W+j+0] = T(3)*(xl*xl);
			B[k*W+j+1] = T(2)*(xl);
			B[k*W+j+2] = T(1);
			B[k*W+j+3] = T(0);
			B[k*W+0] = T(-3)*(xh*xh);
			B[k*W+1] = T(-2)*(xh);
			B[k*W+2] = T(-1);
			B[k*W+3] = T(0);
			d[k] = T(0);
			k++;
		case 0:
			// C0 continuity
			B[k*W+j+0] = xl*xl*xl;
			B[k*W+j+1] = xl*xl;
			B[k*W+j+2] = xl;
			B[k*W+j+3] = T(1);
			B[k*W+0] = -(xh*xh*xh);
			B[k*W+1] = -(xh*xh);
			B[k*W+2] = -(xh);
			B[k*W+3] = -T(1);
			d[k] = T(0);
		}
	}

	// solve overdetermined system A.x=b subject to equality constraints B.x=d
	if (!LSESolve(A,b,B,d,N,W,(C+1)*(S-1+p)))
	{
		delete[] A; delete[] b;
		delete[] B; delete[] d;
		return NaturalSpline();
	}

	T (*c)[4] = new T[S][4]();
	T* X = new T[S+1]();
	for (unsigned int i=0; i<S; ++i)
	{
		X[i] = xs[i];
		for (unsigned int j=0; j<4; ++j)
			c[i][j] = b[4*i+(3-j)];
	}
	X[S] = x[N-1];

	delete[] A; delete[] b;
	delete[] B; delete[] d;
	return NaturalSpline(S,X,c);
}

template<typename T> NaturalSpline<T>::NaturalSpline(
	const NaturalSpline& cpy
	)
{
	Curve<T>::n_ = cpy.n_;
	Curve<T>::X_ = new T[Curve<T>::n_+1]();
	Curve<T>::C_ = new T[Curve<T>::n_][4]();
	Curve<T>::lastSeg_ = cpy.lastSeg_;
	for (int s=0; s<Curve<T>::n_; ++s)
	{
		const_cast<T&>(Curve<T>::X_[s]) = cpy.X_[s];
		for (unsigned int j=0; j<4; ++j)
			const_cast<T&>(Curve<T>::C_[s][j]) = cpy.C_[s][j];
	}
	const_cast<T&>(Curve<T>::X_[Curve<T>::n_]) = cpy.X_[Curve<T>::n_];
}

template<typename T> NaturalSpline<T>::NaturalSpline(
	NaturalSpline&& cpy
	)
{
	Curve<T>::n_ = cpy.n_;
	Curve<T>::X_ = cpy.X_;
	Curve<T>::C_ = cpy.C_;
	Curve<T>::lastSeg_ = cpy.lastSeg_;
	cpy.n_ = 0;
	cpy.X_ = nullptr;
	cpy.C_ = nullptr;
	cpy.lastSeg_ = 0;
}

template<typename T> NaturalSpline<T>& NaturalSpline<T>::operator=(
	const NaturalSpline& cpy
	)
{
	return operator=(std::move(NaturalSpline(cpy)));
}

template<typename T> NaturalSpline<T>& NaturalSpline<T>::operator=(
	NaturalSpline&& cpy
	)
{
	delete[] Curve<T>::X_;
	delete[] Curve<T>::C_;
	Curve<T>::n_ = cpy.n_;
	Curve<T>::X_ = cpy.X_;
	Curve<T>::C_ = cpy.C_;
	Curve<T>::lastSeg_ = cpy.lastSeg_;
	cpy.n_ = 0;
	cpy.X_ = nullptr;
	cpy.C_ = nullptr;
	cpy.lastSeg_ = 0;
	return *this;
}

template<typename T> NaturalSpline<T>::~NaturalSpline()
{
	delete[] Curve<T>::X_;
	delete[] Curve<T>::C_;
}

template<typename T> NaturalSpline<T> NaturalSpline<T>::Diff(
	unsigned int order
	) const
{
	T* x = new T[Curve<T>::n_+1];
	T (*C)[4] = new T[Curve<T>::n_][4]();

	for (int k=0; k<Curve<T>::n_; ++k)
	{
		x[k] = Curve<T>::X_[k];
		switch (order)
		{
		case 0:
			C[k][0] = Curve<T>::C_[k][0];
			C[k][1] = Curve<T>::C_[k][1];
			C[k][2] = Curve<T>::C_[k][2];
			C[k][3] = Curve<T>::C_[k][3];
			break;
		case 1:
			C[k][0] = Curve<T>::C_[k][1];
			C[k][1] = T(2) * Curve<T>::C_[k][2];
			C[k][2] = T(3) * Curve<T>::C_[k][3];
			C[k][3] = T(0);
			break;
		case 2:
			C[k][0] = T(2) * Curve<T>::C_[k][2];
			C[k][1] = T(6) * Curve<T>::C_[k][3];
			C[k][2] = C[k][3] = T(0);
			break;
		case 3:
			C[k][0] = T(6) * Curve<T>::C_[k][3];
			C[k][1] = C[k][2] = C[k][3] = T(0);
			break;
		default:
			C[k][0] = C[k][1] = C[k][2] = C[k][3] = T(0);
			break;
		}
	}
	x[Curve<T>::n_] = Curve<T>::X_[Curve<T>::n_];
	return NaturalSpline(Curve<T>::n_,x,C);
}

template<typename T> bool NaturalSpline<T>::IsMonotonic() const
{
	T zeros[3];
	NaturalSpline<T> Dspl = Diff(1);
	for (int k=0; k<Curve<T>::n_; ++k)
		if (Dspl.RealRoots(k,zeros) != 0)
			return false;
	return true;
}

template<typename T> bool NaturalSpline<T>::IsMonotonic(
	const T& lo,
	const T& hi
	) const
{
	T zeros[3];
	NaturalSpline<T> Dspl = Diff(1);
	int lo_s = Curve<T>::findSegment(lo);
	int hi_s = Curve<T>::findSegment(hi);

	int nz = Dspl.RealRoots(lo_s,zeros);
	while (nz-- > 0)
		if (zeros[nz]>=lo && zeros[nz]<=hi)
			return false;
	if (lo_s == hi_s)
		return true;
	nz = Dspl.RealRoots(hi_s,zeros);
	while (nz-- > 0)
		if (zeros[nz]>=lo && zeros[nz]<=hi)
			return false;
	for (int s=lo_s+1; s<hi_s; ++s)
		if (Dspl.RealRoots(s,zeros) != 0)
			return false;
	return true;
}

template<typename T>
NaturalSpline<T> NaturalSpline<T>::Interpolate(
	const T* x,
	const T* y,
	int N,
	int D
	)
{
	int n = N-1;
	T* X = new T[N];
	for (int i=0; i<N; ++i)
		X[i] = x[i];

	if (D > 3)
		D = 3;
	if (D < 1)
		D = 1;
	if (N == 3 && D>2)
		D = 2;
	else if (N == 2 && D>1)
		D = 1;

	T (*C)[4] = new T[n][4]();
	switch (D)
	{
	case 1:
		computeCoeffs1(n,C,x,y);
		break;
	case 2:
		computeCoeffs2(n,C,x,y);
		break;
	default:
		computeCoeffs3(n,C,x,y);
	}
	return NaturalSpline(n,X,C);
}

template<typename T> NaturalSpline<T> NaturalSpline<T>::Fit(
	const T* x,
	const T* y,
	unsigned int N,
	const T* xs,
	unsigned int S,
	unsigned int C,
	unsigned int D,
	bool periodic
	)
{
	switch (D)
	{
	case 1:
		return fit1(x,y,N,xs,S,periodic);
	case 2:
		return fit2(x,y,N,xs,S,C,periodic);
	case 3:
		return fit3(x,y,N,xs,S,C,periodic);
	}
	return NaturalSpline();
}

template<typename T> NaturalSpline<T> NaturalSpline<T>::CrossValidatedFit(
	const T* x,
	const T* y,
	unsigned int N,
	unsigned int C,
	unsigned int D,
	bool periodic,
	unsigned int Smax // = -1
	)
{
	std::default_random_engine PRNG(47110815);
	std::uniform_int_distribution<int> randu(0,N-1);
	auto random = std::bind(randu,PRNG);
	auto P = std::make_unique<int[]>(N);

	auto kFoldCrossValidation = [&](int k, int S)
	{
		// random permutation
		std::iota(P.get(),P.get()+N,0);
		std::random_shuffle(P.get(),P.get()+N,random);

		// generate k partitions, sort permutation indices within partitions
		LinearPartition part(N,k);
		for (auto interval : part)
			std::sort(P.get()+interval.low,P.get()+interval.high+1);

		// S segments (S+1 segment boundaries)
		auto xs = Lerp<T>(x[0],x[N-1],S+1).Render();

		T CVscore(0);
		for (int validation_part=0; validation_part<k; ++validation_part)
		{
			// assemble data for fitting
			int training_set_size = N-int(part[validation_part].size());
			auto Pt = std::make_unique<int[]>(training_set_size);
			auto xt = std::make_unique<T[]>(training_set_size);
			auto yt = std::make_unique<T[]>(training_set_size);
			int i=0;
			for (int p=0; p<k; ++p)
			{
				if (p == validation_part)
					continue; // omit partition that is used for validation
				for (int j : part[p])
					Pt[i++] = P[j];
			}
			std::sort(Pt.get(),Pt.get()+training_set_size);
			for (int j=0; j<training_set_size; ++j)
			{
				xt[j] = x[Pt[j]];
				yt[j] = y[Pt[j]];
			}

			auto spline = NaturalSpline<T>::Fit(xt.get(),yt.get(),
					training_set_size,xs.get(),S,C,D,periodic);
			if (!spline)
			{
				// fitting failed
				return std::numeric_limits<T>::infinity();
			}

			T MSE(0);
			for (int j : part[validation_part])
			{
				T e = spline(x[P[j]]) - y[P[j]];
				MSE += e*e;
			}
			MSE /= T(part[validation_part].size());
			CVscore += MSE;
		}
		CVscore /= T(k);
		return CVscore;
	}; // kFoldCrossValidation(k,S)

	if (Smax == -1)
	{
		// limit the number of segments to 75% of maximum
		Smax = 3*(N/(D+1))/4;
	}
	// 3-fold cross validation
	const int k = 3;
	T minScore(-1);
	int Sopt(0);
	for (int S=1; S<=Smax; ++S)
	{
		T score = kFoldCrossValidation(k,S);
		if (minScore<0 || score<minScore)
		{
			minScore = score;
			Sopt = S;
		}
	}
	auto xs = Lerp<T>(x[0],x[N-1],Sopt+1).Render();
	return Fit(x,y,N,xs.get(),Sopt,C,D,periodic);
}

template<typename T> void NaturalSpline<T>::Dump() const
{
	for (unsigned int s=0; s<Curve<T>::GetNumSegments(); ++s)
	{
		printf("s_%u(x) = ", s);
		for (unsigned int d=0,D=Curve<T>::GetDegree(s); d<=D; ++d)
		{
			printf("%+g*(x-%g)^%u", Curve<T>::GetCoeff(s,d),
				Curve<T>::GetAbscissa(s), d);
		}
		printf("\n");
	}
}

// vim: fenc=utf-8 noet:
