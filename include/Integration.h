#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include <cmath>
#include <functional>
#include <tuple>
#include <limits>

enum class QuadratureMethod
{
	AdaptiveSimpson,
	AdaptiveBoole,
	AdaptiveGaussLegendre
};

template<typename T,typename FuncT> class Integrator
{
	static std::tuple<T, T, unsigned int> adaptiveSimpson(FuncT f,
		const T & sab,
		const T & fa,
		const T & fc,
		const T & fb,
		const T & a,
		const T & b,
		const T & tol
		)
	{
		T c = (a + b) / 2;
		T fd = f((a + c) / 2);
		T fe = f((c + b) / 2);
		T sac = (c - a) * (fa + 4 * fd + fc) / 6;
		T scb = (b - c) * (fc + 4 * fe + fb) / 6;
		T estErr = std::abs(sab - sac - scb);

		if (estErr < 10 * tol)
			return {sac + scb, estErr / 10, 2};
		else
		{
			T appInt1, appInt2, appErr1, appErr2;
			int nEval1, nEval2;
			std::tie(appInt1, appErr1, nEval1) =
				adaptiveSimpson(f, sac, fa, fd, fc, a, c, tol / 2);
			std::tie(appInt2, appErr2, nEval2) =
				adaptiveSimpson(f, scb, fc, fe, fb, c, b, tol / 2);
			return {appInt1 + appInt2, appErr1 + appErr2, 2 + nEval1 + nEval2};
		}
	}

	static std::tuple<T, T, unsigned int> adaptiveBoole(FuncT f,
		const T & sab,
		const T & fa,
		const T & fc,
		const T & fd,
		const T & fe,
		const T & fb,
		const T & a,
		const T & b,
		const T & tol
		)
	{
		T h = (b - a) / 8;
		T f1 = f(a + h);
		T f3 = f(a + 3 * h);
		T f5 = f(a + 5 * h);
		T f7 = f(a + 7 * h);
		T sad = 2 * h * (7 * (fa + fd) + 32 * (f1 + f3) + 12 * fc) / 45;
		T sdb = 2 * h * (7 * (fd + fb) + 32 * (f5 + f7) + 12 * fe) / 45;

		T estErr = std::abs(sab - sad - sdb);
		if (estErr < 42 * tol) return {sad + sdb, estErr / 42, 4};

		T appInt1, appInt2, appErr1, appErr2;
		int nEval1, nEval2;
		std::tie(appInt1, appErr1, nEval1) =
			adaptiveBoole(f, sad, fa, f1, fc, f3, fd, a, a + 4 * h, tol / 2);
		std::tie(appInt2, appErr2, nEval2) =
			adaptiveBoole(f, sdb, fd, f5, fe, f7, fb, a + 4 * h, b, tol / 2);

		return {appInt1 + appInt2, appErr1 + appErr2, 4 + nEval1 + nEval2};
	}

	static std::tuple<T, T, unsigned int> adaptiveGaussLegendre(
		FuncT f,
		const T & sab,
		const T & a,
		const T & b,
		const T & tol
		)
	{
		constexpr T w = std::sqrt(T(3) / T(5));
		T h2 = (b - a) / 4;
		T fml = f(a + h2 - w * h2);
		T fcl = f(a + h2);
		T fpl = f(a + h2 + w * h2);
		T fmr = f(b - h2 - w * h2);
		T fcr = f(b - h2);
		T fpr = f(b - h2 + w * h2);
		T sac = h2 * (5 * fml + 8 * fcl + 5 * fpl) / 9;
		T scb = h2 * (5 * fmr + 8 * fcr + 5 * fpr) / 9;

		T estErr = std::abs(sab - sac - scb);
		if (estErr < 42 * tol) return {sac + scb, estErr / 42, 6};

		T appInt1, appInt2, appErr1, appErr2;
		int nEval1, nEval2;
		std::tie(appInt1, appErr1, nEval1) =
			adaptiveGaussLegendre(f, sac, a, a + 2 * h2, tol / 2);
		std::tie(appInt2, appErr2, nEval2) =
			adaptiveGaussLegendre(f, scb, a + 2 * h2, b, tol / 2);
		return {appInt1 + appInt2, appErr1 + appErr2, 6 + nEval1 + nEval2};
	}

public:
	static std::tuple<T, T, unsigned int> Run(FuncT f,
		const T & a,
		const T & b,
		const T & tol,
		QuadratureMethod formula
		)
	{
		switch (formula)
		{
		case QuadratureMethod::AdaptiveSimpson:
			return AdaptiveSimpson(f,a,b,tol);
		case QuadratureMethod::AdaptiveBoole:
			return AdaptiveBoole(f,a,b,tol);
		case QuadratureMethod::AdaptiveGaussLegendre:
			return AdaptiveGaussLegendre(f,a,b,tol);
		}
		return {T(0),std::numeric_limits<T>::infinity(),0};
	}

	/**
	 * Approximation of the definite integral or an unary function
	 * with specified error tolerance using adaptive quadrature
	 * based on Simpson's rule.
	 */
	static std::tuple<T, T, unsigned int> AdaptiveSimpson(
		FuncT f,
		const T & a,
		const T & b,
		const T & tol
		)
	{
		T fa = f(a);
		T fc = f((a + b) / 2);
		T fb = f(b);

		T sab = (b - a) * (fa + 4 * fc + fb) / 6;
		T appInt, appErr;
		int nEval;
		std::tie(appInt, appErr, nEval) =
			adaptiveSimpson(f, sab, fa, fc, fb, a, b, tol);
		return {appInt, appErr, 3 + nEval};
	}

	/**
	 * Approximation of the definite integral or an unary function
	 * with specified error tolerance using adaptive quadrature
	 * based on Boole's rule (closed Newton-Cotes quadrature rule
	 * with n=4).
	 */
	static std::tuple<T, T, unsigned int> AdaptiveBoole(
		FuncT f,
		const T & a,
		const T & b,
		const T & tol
		)
	{
		T h = (b - a) / 4;
		T fa = f(a);
		T fb = f(b);
		T fc = f(a + h);
		T fd = f(a + 2 * h);
		T fe = f(a + 3 * h);
		T sab = (b - a) * (7 * (fa + fb) + 32 * (fc + fe) + 12 * fd) / 90;

		T appInt, appErr;
		int nEval;
		std::tie(appInt, appErr, nEval) =
			adaptiveBoole(f, sab, fa, fc, fd, fe, fb, a, b, tol);
		return {appInt, appErr, 5 + nEval};
	}

	/**
	 * Approximation of the definite integral or an unary function
	 * with specified error tolerance using adaptive quadrature
	 * based on the three-point Gaussian quadrature rule.
	 */
	static std::tuple<T, T, unsigned int> AdaptiveGaussLegendre(
		FuncT f,
		const T & a,
		const T & b,
		const T & tol
		)
	{
		constexpr T w = std::sqrt(T(3) / T(5));
		T h2 = (b - a) / 2;
		T fm = f(a + h2 - w * h2);
		T fc = f(a + h2);
		T fp = f(a + h2 + w * h2);
		T sab = h2 * (5 * fm + 8 * fc + 5 * fp) / 9;

		T appInt, appErr;
		int nEval;
		std::tie(appInt, appErr, nEval) =
			adaptiveGaussLegendre(f, sab, a, b, tol);
		return {appInt, appErr, 3 + nEval};
	}
};


/**
 * Numerical integration (front-end).
 *
 * @param	f	function to integrate
 * @param	a	lower integration bound
 * @param	b	upper integration bound
 * @param	tol	tolerance
 * @param	method	Integration method to use (default: adaptive Boole)
 * @return	tuple (approximate integral, approximate error, number of evals)
 */
template<typename T,typename FuncT> std::tuple<T,T,unsigned int> Integrate(
	FuncT f,
	const T& a,
	const T& b,
	const T& tol,
	QuadratureMethod method = QuadratureMethod::AdaptiveBoole)
{
	return Integrator<T,FuncT>::Run(f,a,b,tol,method);
}

// vim: fenc=utf-8 noet:
