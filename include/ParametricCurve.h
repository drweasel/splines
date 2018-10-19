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
#include "Lerp.h"
#include <array>

template<typename T,unsigned int D> class ParametricCurve
{
private:
	/** polynomial curves for individual dimensions */
	Curve<T> curve_[D];

public:
	inline ParametricCurve()
	{
		static_assert(D>0,"number of dimensions has to be positive");
	}

	inline Curve<T>& operator[](unsigned int i) { return curve_[i]; }
	inline const Curve<T>& operator[](unsigned int i) const { return curve_[i]; }

	inline std::array<T,D> Eval(T t) const
	{
		std::array<T,D> values;
		for (unsigned int d=0; d<D; ++d)
			values[d] = curve_[d](t);
		return values;
	}

	inline static constexpr unsigned int GetDim() { return D; }

	inline T GetParam(unsigned int s) const
	{
		return curve_[0].GetAbscissa(s);
	}

	inline T GetMinParam() const
	{
		return curve_[0].GetAbscissa(0);
	}

	inline T GetMaxParam() const
	{
		return curve_[0].GetAbscissa(GetNumSegments());
	}

	inline unsigned int GetNumSegments() const
	{
		return curve_[0].GetNumSegments();
	}

	/**
	 * Orthogonal projection of a point to the spline.
	 *
	 * Bibliography:
	 *  H. Wang, J. Kearney, K. Atkinson, 2002: Robust and Efficient
	 *  Computation of the Closest Point on a Spline Curve; Curve and
	 *  Surface Design, Saint-Malo 2002
	 */
	T ProjectOn(const std::array<T,D>& p) const;
};

template<typename T,unsigned int D> T ParametricCurve<T,D>::ProjectOn(
	const std::array<T,D>& p
	) const
{
	// quadratic distance from curve(t) to point p
	auto dist2 = [&](const T& t)->T
	{
		T d2(0);
		for (unsigned int k=0; k<D; ++k)
		{
			T d = curve_[k].Eval(t) - p[k];
			d2 += d*d;
		}
		return d2;
	};

	// Newton iteration
	auto newton = [&](T t)->T
	{
		constexpr T tol = std::numeric_limits<T>::epsilon()*T(10);
		// up to 10 iterations (quadratic convergence)
		for (unsigned int i=0; i<10; ++i)
		{
			T Ddist2(0);
			T DDdist2(0);
			for (unsigned int k=0; k<D; ++k)
			{
				T d = curve_[k].Eval(t) - p[k];
				T Dd = curve_[k].Eval(t,1);
				Ddist2 += T(2) * Dd * d;
				DDdist2 += T(2) * (curve_[k].Eval(t,2)*d + Dd*Dd);
			}
			// Newton update
			T t_new = t - Ddist2/DDdist2;
			if (t_new != t_new)
				break;
			// convergence?
			if (std::abs(t_new - t) > tol)
				t = t_new;
			else
				break;
		}
		return t;
	};

	// quadratic approximation of parameter t at minimal distance
	auto quad_approx_min = [](
		const T& t1, const T& t2, const T& t3,
		const T& Dt1, const T& Dt2, const T& Dt3
		)->T
	{
		T y12 = t1*t1 - t2*t2, y23 = t2*t2 - t3*t3, y31 = t3*t3 - t1*t1;
		T t12 = t1 - t2, t23 = t2 - t3, t31 = t3 - t1;
		return T(0.5) * (y23*Dt1 + y31*Dt2 + y12*Dt3) /
			(t23*Dt1 + t31*Dt2 + t12*Dt3);
	};

	// approximate and minimize a quadratic polynomial (robust)
	Lerp<T> L(GetMinParam(),GetMaxParam(),4*GetNumSegments()+1);
	T t_min = GetMinParam();
	T d_min = std::numeric_limits<T>::infinity();
	for (unsigned int k=1; k<L.Size()-1; ++k)
	{
		T t[4] = { L[k-1],0,L[k+1],0 };
		T d[4] = { dist2(t[0]),0,dist2(t[2]),0 };
		for (unsigned int r=0; r<3; ++r)
		{
			t[1] = T(0.5)*(t[0]+t[2]);
			d[1] = dist2(t[1]);

			t[3] = quad_approx_min(t[0],t[1],t[2],d[0],d[1],d[2]);
			d[3] = t[3]==t[3] ? dist2(t[3])
				: std::numeric_limits<T>::infinity();

			unsigned int i_max = 0;
			for (unsigned int j=1; j<4; ++j)
				if (d[j] > d[i_max])
					i_max = j;
			if (i_max == 3)
				break;
			for (int j=2; j>=0; --j)
				if (t[j+1]<t[j])
				{
					std::swap(t[j],t[j+1]);
					std::swap(d[j],d[j+1]);
				} else break;
		}
		// result of quadratic approximation
		unsigned int i_min = 0;
		for (unsigned int j=1; j<4; ++j)
			if (d[j] < d[i_min])
				i_min = j;
		// Newton iteration for accuracy
		T t_min_k = newton(t[i_min]);
		T d_min_k = dist2(t_min_k);
		if (d_min_k < d_min)
		{
			t_min = t_min_k;
			d_min = d_min_k;
		}
	}
	if (t_min < GetMinParam())
		return GetMinParam();
	if (t_min > GetMaxParam())
		return GetMaxParam();
	return t_min;
}

// vim: fenc=utf-8 noet:
