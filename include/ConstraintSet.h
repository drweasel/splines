#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "Vector.h"
#include <vector>
#include <array>
#include <limits>

template<unsigned int D,typename T> class ConstraintSet
{
public:
	struct Constraint
	{
		Vector<T,D> C;
		T d;

		bool CheckLE(const Vector<T,D>& v) const
		{
			return C*v <= d;
		}
	};

private:
	std::vector<Constraint> Cd_;

public:
	inline ConstraintSet(
		const std::initializer_list<Constraint>& values
		) { for (const auto& c : values) Add(c); }

	inline void Add(const Constraint& Cd) { Cd_.push_back(Cd); }
	inline void Add(const std::array<T,D>& C, T d) { Cd_.push_back({C,d}); }
	inline const T& operator()(std::size_t i, std::size_t j) const
	{
		return Cd_[i].C[j];
	}
	inline const T& operator[](std::size_t i) const
	{
		return Cd_[i].d;
	}

	inline constexpr std::size_t cols() const { return std::size_t(D); }
	inline std::size_t rows() const { return Cd_.size(); }

	bool CheckLE(const Vector<T,D>& v)
	{
		for (const Constraint& cn : Cd_)
			if (!cn.CheckLE(v))
				return false;
		return true;
	}

	Vector<T,D> CutBack(
		const Vector<T,D>& p,
		const Vector<T,D>& delta_p
		) const
	{
		// update direction
		Vector<T,D> r = delta_p.dir();
		std::size_t i_min = 0;
		T d_min = std::numeric_limits<T>::infinity();
		for (std::size_t i=0; i<rows(); ++i)
		{
			T rho = (Cd_[i].C * p - Cd_[i].d) / (Cd_[i].C * r);
			Vector<T,D> pc = p - r*rho; // intersection with constraint
			T d_pc = std::abs(rho); // distance to intersection
			if (r*(pc-p) >= T(0) && d_pc < d_min)
			{
				// pc in same direction as r (resp. delta_p) and lower distance to pc
				d_min = d_pc;
				i_min = i;
			}
		}
		if (d_min < delta_p.norm2())
			return r*d_min;
		return delta_p;
	}
};

// vim: fenc=utf-8 noet:
