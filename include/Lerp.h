#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "LinearMap.h"
#include <cstddef>
#include <memory>

template<typename T> class Lerp
{
private:
	int k_;
	unsigned int n_;
	T s_, t_;

	inline Lerp(const Lerp& r, int k)
		: k_(k), n_(r.n_), s_(r.s_), t_(r.t_) { }
public:
	inline Lerp()
		: k_(0), n_(0), s_(0), t_(0) { }
	inline Lerp(T lo, T hi, unsigned int n)
		: k_(0), n_(n<2?2:n), s_((hi-lo)/(n_-1)), t_(lo) { }
	inline Lerp& operator++() { ++k_; return *this; }
	inline Lerp& operator--() { --k_; return *this; }
	inline T operator*() const { return k_*s_+t_; }
	inline Lerp begin() const { return Lerp(*this,0); }
	inline Lerp end() const { return Lerp(*this,int(n_)); }
	inline T operator[](int k) const { return k*s_+t_; }
	inline bool operator==(const Lerp& r) const { return k_==r.k_; }
	inline bool operator!=(const Lerp& r) const { return k_!=r.k_; }
	inline std::size_t Step() const { return k_; }
	inline std::size_t Size() const { return n_; }
	inline operator bool() const { return Step()>=0 && Step()<Size(); }
	inline operator LinearMap<T>() const { return LinearMap<T>(s_,t_); }
	inline void Render(T* p) const
	{
		for (T v : *this)
			*(p++) = v;
	}
	inline std::unique_ptr<T[]> Render() const
	{
		auto sp = std::make_unique<T[]>(n_);
		Render(sp.get());
		return sp;
	}
};

// vim: fenc=utf-8 noet:
