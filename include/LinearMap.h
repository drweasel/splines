#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */

template<typename T> class LinearMap
{
private:
	T s_, t_;
public:
	inline LinearMap() : s_(1), t_(0) { }
	inline LinearMap(T s, T t) : s_(s), t_(t) { }
	inline LinearMap(T il, T ih, T ol, T oh)
		: s_((oh-ol)/(ih-il)), t_((ih*ol-il*oh)/(ih-il)) { }
	inline T operator()(T v) const { return v*s_+t_; }
	inline LinearMap operator~() const { return LinearMap(1/s_,-t_/s_); }

	inline T scale(T v) const { return v*s_; }
};

// vim: fenc=utf-8 noet:
