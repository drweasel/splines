#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "ParametricCurve.h"
#include "BezierPolyLine.h"
#include <array>
#include <memory>

/**
 * Bézier-Curve
 *
 * A Bézier-Curve is simply a two dimensional ParametricCurve,
 * obtained from a two dimensional Bézier-Poly-Line.
 *
 * @author	Michael Weitzel <mich@elweitzel.de>
 */
template<typename T>
class BezierCurve : public ParametricCurve<T, 2>
{
public:
	using ParametricCurve<T, 2>::GetDim;

private:
	/** allocated memory for running parameter */
	std::unique_ptr<T[]> t_;
	/** allocated memory for polynomial coefficients */
	std::array<std::unique_ptr<T[][4]>, GetDim()> C_;

public:
	BezierCurve() = default;

	/** copy */
	BezierCurve(const BezierCurve & cpy);
	BezierCurve & operator=(const BezierCurve & cpy);
	/** move */
	BezierCurve(BezierCurve && cpy) = default;
	BezierCurve & operator=(BezierCurve && cpy) = default;

	BezierCurve(const BezierPolyLine<T> & bpl);
};

template<typename T>
BezierCurve<T>::BezierCurve(const BezierCurve<T> & cpy)
{
	t_ = std::make_unique<T[]>(cpy.GetNumSegments());
	for (unsigned int s = 0; s < cpy.GetNumSegments(); ++s)
		t_[s] = cpy.t_[s];

	for (unsigned int d = 0; d < GetDim(); ++d)
	{
		C_[d].reset(new T[cpy.GetNumSegments()][4]);
		for (unsigned int s = 0; s < cpy.GetNumSegments(); ++s)
			for (unsigned int i = 0; i < 4; ++i)
				C_[d][s][i] = cpy.C_[d][s][i];
		ParametricCurve<T, 2>::operator[](d) =
			Curve<T>(cpy.GetNumSegments(), t_.get(), C_[d].get());
	}
}

template<typename T>
BezierCurve<T> & BezierCurve<T>::operator=(const BezierCurve & cpy)
{
	return operator=(std::move(BezierCurve(cpy)));
}

template<typename T>
BezierCurve<T>::BezierCurve(const BezierPolyLine<T> & bpl)
{
	t_ = Lerp<T>(T(0), T(bpl.GetNumSegments()),
		bpl.GetNumNodes()).Render();

	for (unsigned int d = 0; d < GetDim(); ++d)
	{
		C_[d].reset(new T[bpl.GetNumSegments()][4]);
		for (int k = 0; k < int(bpl.GetNumSegments()); ++k)
		{
			T p0 = bpl.Get(k)[d];
			T p1 = bpl.GetRightCtrl(k)[d];
			T p2 = bpl.GetLeftCtrl(k + 1)[d];
			T p3 = bpl.Get(k + 1)[d];

			// convert Bézier control points into polynomial coefficients:
			C_[d][k][0] = p0;                                // a,t^0
			C_[d][k][1] = T(3) * (p1 - p0);                  // b,t^1
			C_[d][k][2] = T(3) * p2 - T(6) * p1 + T(3) * p0; // c,t^2
			C_[d][k][3] = p3 - T(3) * (p2 - p1) - p0;        // d,t^3
		}
		ParametricCurve<T, 2>::operator[](d) =
			Curve<T>(bpl.GetNumSegments(), t_.get(), C_[d].get());
	}
}

// vim: fenc=utf-8 noet:
