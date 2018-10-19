#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "Polynomials.h"
#include "Integration.h"
#include <limits>
#include <tuple>

/**
 * Segmented Polynomial Curve (up to polynomials of degree 3).
 * The segments are not necessarily continuous.
 *
 * N.B.: this class uses pre-allocated fields, i.e. memory is
 * neither allocated nor released.
 *
 * @author	Michael Weitzel <mich@elweitzel.de>
 */
template<typename T> class Curve
{
protected:
	/** Number of segments */
	unsigned int n_;
	/** Abscissas of the interpolated points */
	const T* X_;
	/** Coefficients of the polynomials s_j(x) for segment j=1,...,n
		with s_j(x) = ... */
	const T (*C_)[4];
	/** Last segment evaluated */
	mutable unsigned int lastSeg_;

protected:
	inline unsigned int findSegment(const T& x) const;

public:
	/**
	 * C'tor.
	 * Beware: a two-dim array is passed here. The type of the passed
	 * array decays into T(*)[4]. Since C++11, these arrays can be
	 * allocated like this "T (*name)[4] = new T[42][4];" (or, with
	 * default init.: "T (*name)[4] = new T[42][4]();"), and deallocated
	 * simply using "delete[] name".
	 *
	 * @param	n	Number of segments
	 * @param	X	n+1 segment boundaries
	 * @param	C	n x 4 coefficients, C[s][0]+C[s][1]*x¹+C[s][2]*x²+C[s][3]*x³
	 */
	inline Curve(unsigned int n, const T* X, const T C[][4])
		: n_(n)
		, X_(X)
		, C_(C)
		, lastSeg_(0) { }

public:
	inline Curve() : n_(0), X_(nullptr), C_(nullptr), lastSeg_(0) { }

	inline T Eval(const T& x) const;
	inline T Eval(const T& x, unsigned int diff) const;
	inline T operator()(const T& x) const { return Eval(x); }
	inline T operator()(const T& x, unsigned int diff) const { return Eval(x,diff); }

	inline unsigned int GetNumSegments() const { return n_; }
	inline unsigned int GetNumSegments(const T& lo, const T& hi) const
	{ return findSegment(hi)-findSegment(lo)+1; }
	inline const T& GetAbscissa(unsigned int i) const { return X_[i]; }
	inline const T& GetCoeff(unsigned int i, unsigned int j) const { return C_[i][j]; }

	inline operator bool() const { return n_>0; }

	/** Analytical definite integration */
	T Integrate(T lo, T hi) const;
	/** Numerical computation of curve length with given accuracy */
	T Length(const T& lo, const T& hi, const T& accuracy) const;
	/** Real roots of segment s */
	unsigned int RealRoots(unsigned int s, T zeros[3]) const;
	/** Degree of segment s */
	unsigned int GetDegree(unsigned int s) const;
	/** Test for continuity of segment s-1 and segment s */
	bool IsContinuous(unsigned int s, unsigned int diff, const T tol) const;

};

// Template implementations

template<typename T> inline unsigned int Curve<T>::findSegment(const T& x) const
{
	if (x >= X_[lastSeg_] && x<X_[lastSeg_+1])
		return lastSeg_;
	unsigned int mid;
	if (x >= X_[n_])
		mid = n_-1;
	else if (x <= X_[0])
		mid = 0;
	else
	{
		unsigned int lo = 0;
		unsigned int hi = n_;
		for (;;)
		{
			mid = lo+(hi-lo)/2;
			if (x>=X_[mid] && x<X_[mid+1]) break;
			if (x<X_[mid]) hi=mid; else lo=mid;
		}
	}
	lastSeg_ = mid;
	return mid;
}

template< typename T > inline T Curve<T>::Eval(
	const T& x
	) const
{
	unsigned int mid = findSegment(x);
	T xd = x - X_[mid];
	return C_[mid][0] + xd*(C_[mid][1] + xd*(C_[mid][2] + xd*C_[mid][3]));
}

template< typename T > T Curve<T>::Eval(
	const T& x,
	unsigned int diff
	) const
{
	if (diff > 4)
		return T(0);

	unsigned int mid = findSegment(x);
	T xd = x - X_[mid];
	switch (diff)
	{
	case 0: return C_[mid][0] + xd*(C_[mid][1] + xd*(C_[mid][2] + xd*C_[mid][3]));
	case 1: return C_[mid][1] + xd*(2*C_[mid][2] + 3*C_[mid][3]*xd);
	case 2: return 2*C_[mid][2] + 6*C_[mid][3]*xd;
	case 3: return 6*C_[mid][3];
	default: return 0;
	}
}

template<typename T> T Curve<T>::Integrate(
	T lo,
	T hi
	) const
{
	T s(1);
	unsigned int j,l,h,m_lo,m_hi;

	if (lo>hi) { T t(lo); lo = hi; hi = t; s=T(-1); }

	// search index of the control point <= lo
	if (lo < X_[0])
		m_lo = 0;
	else if (lo >= X_[n_])
		m_lo = n_-1;
	else
	{
		for (l=0,h=n_;;)
		{
			m_lo = l+(h-l)/2;
			if (lo>=X_[m_lo] && lo<X_[m_lo+1]) break;
			if (lo<X_[m_lo]) h=m_lo; else l=m_lo;
		}
	}
	// search index of the control point <= hi
	if (hi < X_[0])
		m_hi = 0;
	else if (hi >= X_[n_])
		m_hi = n_-1;
	else
	{
		for (l=m_lo,h=n_;;)
		{
			m_hi = l+(h-l)/2;
			if (hi>=X_[m_hi] && hi<X_[m_hi+1]) break;
			if (hi<X_[m_hi]) h=m_hi; else l=m_hi;
		}
	}

	T t(0);
	T dx;
	// everything in one segment:
	if (m_lo == m_hi)
	{
		dx = hi-lo;
		t += C_[m_lo][0]*dx;
		t += C_[m_lo][1]*(dx*dx)/2.;
		t += C_[m_lo][2]*((dx*dx)*dx)/3.;
		t += C_[m_lo][3]*(((dx*dx)*dx)*dx)/4.;
		return s*t;
	}
	// upper remainder [X_[m_hi],hi]
	dx = hi - X_[m_hi];
	t += C_[m_hi][0]*dx;
	t += C_[m_hi][1]*(dx*dx)/2.;
	t += C_[m_hi][2]*((dx*dx)*dx)/3.;
	t += C_[m_hi][3]*(((dx*dx)*dx)*dx)/4.;
	// full segments [X_[m_lo+1],X_[m_hi]]
	for (j=m_lo+1; j<m_hi; ++j)
	{
		dx = X_[j+1]-X_[j];
		t += C_[j][0]*dx;
		t += C_[j][1]*(dx*dx)/2.;
		t += C_[j][2]*((dx*dx)*dx)/3.;
		t += C_[j][3]*(((dx*dx)*dx)*dx)/4.;
	}
	// lower remainder [lo,X_[m_lo+1]]
	dx = X_[m_lo+1] - lo;
	t += C_[m_lo][0]*dx;
	t += C_[m_lo][1]*(dx*dx)/2.;
	t += C_[m_lo][2]*((dx*dx)*dx)/3.;
	t += C_[m_lo][3]*(((dx*dx)*dx)*dx)/4.;
	return s*t;
}

template<typename T> T Curve<T>::Length(
	const T& lo,
	const T& hi,
	const T& accuracy
	) const
{
	return std::get<0>(::Integrate(
		[this](const T& x)->T
		{
			T dp = Eval(x,1);
			T ds = std::sqrt(1 + dp*dp);
			return ds;
		}, lo, hi, accuracy));
}

template<typename T> unsigned int Curve<T>::RealRoots(
	unsigned int s,
	T zeros[3]
	) const
{
	T work[3*3];
	std::complex<T> w[3];
	unsigned int d = GetDegree(s);
	if (!PolynomialRoots(C_[s],d+1,w,work))
		return 0; // Fehler!

	unsigned int j=0;
	T lo = X_[s];
	T hi = X_[s+1] + (s==n_-1 ? 100*std::numeric_limits<T>::epsilon() : T(0));
	for (int k=0; k<d; ++k)
	{
		if (std::abs(w[k].imag()) > 100*std::numeric_limits<T>::epsilon())
			continue;
		T z = w[k].real() + lo;
		if (z >= lo && z < hi)
			zeros[j++] = z;
	}
	return j;
}

template<typename T> unsigned int Curve<T>::GetDegree(unsigned int s) const
{
	unsigned int deg = 3;
	while (deg>0 && GetCoeff(s,deg) == T(0))
		deg--;
	return deg;
}

template<typename T> bool Curve<T>::IsContinuous(
	unsigned int s,
	unsigned int diff,
	const T tol
	) const
{
	// first/last X_ are never continuous
	if (s<=0 || s>=n_ || diff >= GetDegree(s))
		return false;

	T xd = X_[s] - X_[s-1];
	switch (diff)
	{
	case 0:
		return std::abs(C_[s][0] -
			(C_[s-1][0] + xd*(C_[s-1][1] + xd*(C_[s-1][2] + xd*C_[s-1][3])))
			) <= tol;
	case 1:
		return std::abs(C_[s][1] -
			(C_[s-1][1] + xd*(2*C_[s-1][2] + 3*C_[s-1][3]*xd))
			) <= tol;
	case 2:
		return std::abs(2*C_[s][2] -
			(2*C_[s-1][2] + 6*C_[s-1][3]*xd)
			) <= tol;
	case 3:
		return std::abs(6*C_[s][3] -
			(6*C_[s-1][3])
			) <= tol;
	}
	return false;
}

// vim: fenc=utf-8 noet:
