#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "Lerp.h"
#include <cstddef>
#include <cmath>
#include <algorithm>

class LinearPartition
	: private Lerp<double>
{
public:
	struct Interval
	{
		const int index;
		const int low;
		const int high;
		int pos;

		inline Interval() : index(-1), low(0), high(0), pos(-1) { }
		inline Interval(int index, int low, int high)
			: index(index), low(low), high(high), pos(-1) { }
		inline std::size_t size() const { return high-low+1; }
		inline int operator*() const { return pos; }
		inline Interval& operator++() { ++pos; return *this; }
		inline Interval begin() const { Interval b(*this); b.pos = low; return b; }
		inline Interval end() const { Interval e(*this); e.pos = high+1; return e; }
		inline bool operator!=(const Interval& r) const { return pos!=r.pos; }
	};

	LinearPartition() = default;

	inline LinearPartition(std::size_t dim, std::size_t n)
		: Lerp<double>(0,dim>0?int(dim):1,std::min(int(n+1),
			(dim>0?int(dim)-1:0)+2))
	{
		if (dim == 0 || n == 0)
			Lerp<double>();
	}

	inline LinearPartition(int low, int high, std::size_t n)
		: Lerp<double>(low,low<=high?high+1:low+1,
			std::min(int(n+1),(low<=high?high:low)-low+2))
	{
		if (high == low || n == 0)
			Lerp<double>();
	}

	inline Interval operator[](int b) const
	{
		return Interval(b,int(std::floor(Lerp<double>::operator[](b)+.5)),
			int(std::floor(Lerp<double>::operator[](b+1)+.5))-1);
	}

	inline Interval operator*() const { return operator[](Step()); }
	inline Interval operator->() const { return operator*(); }
	inline LinearPartition& operator++()
	{
		Lerp<double>::operator++();
		return *this;
	}
	inline LinearPartition& operator--()
	{
		Lerp<double>::operator--();
		return *this;
	}

	using Lerp<double>::operator bool;

	inline LinearPartition begin() const
	{
		LinearPartition rp;
		static_cast<Lerp<double>&>(rp) = Lerp<double>::begin();
		return rp;
	}

	inline LinearPartition end() const
	{
		LinearPartition rp;
		static_cast<Lerp<double>&>(rp) = Lerp<double>::end();
		if (rp.Step() > 0)
			--rp;
		return rp;
	}

	inline bool operator==(const LinearPartition& r) const
	{
		return Lerp<double>::operator==(r);
	}
	inline bool operator!=(const LinearPartition& r) const
	{
		return Lerp<double>::operator!=(r);
	}

	inline std::size_t Size() const
	{
		return Lerp<double>::Size()>0 ? Lerp<double>::Size() - 1 : 0;
	}

	using Lerp<double>::Step;

};

// vim: fenc=utf-8 noet:
