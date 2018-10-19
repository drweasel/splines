#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include <cstdint>

namespace fixed {

template< int > struct Int { typedef void type; };
template<> struct Int<8> { typedef int8_t type; };
template<> struct Int<16> { typedef int16_t type; };
template<> struct Int<32> { typedef int32_t type; };
template<> struct Int<64> { typedef int64_t type; };

/**
 * A class template implementing basic Q<I,F> fixed point arithmetic,
 * cf. https://en.wikipedia.org/wiki/Q_(number_format)
 *
 * @param	I	number of integral digits
 * @param	F	number of fractional digits
 * @author	Michael Weitzel <mich@elweitzel.de>
 */
template< int I,int F > class Q
{
	typedef typename Int<I+F>::type intB_t;
	typedef typename Int<2*(I+F)>::type int2B_t;

private:
	intB_t intVal_;

	static inline Q create(const intB_t& v)
	{
		Q f; f.intVal_ = v;
		return f;
	}

public:
	inline Q()
	{
		static_assert(I+F == 8 || I+F == 16 || I+F == 32 || I+F == 64,
			"error: I+F must be a power of 2");
	}
	inline /*explicit*/ constexpr Q(int8_t v) : intVal_(v << F) { }
	inline /*explicit*/ constexpr Q(int16_t v) : intVal_(v << F) { }
	inline /*explicit*/ constexpr Q(int32_t v) : intVal_(v << F) { }
	inline /*explicit*/ constexpr Q(float v) : intVal_(intB_t(v*(1<<F)+.5f)) { }
	inline /*explicit*/ constexpr Q(double v) : intVal_(intB_t(v*(1<<F)+.5)) { }

	inline explicit operator float() const
	{
		return float(intVal_)/float(1<<F);
	}
	inline explicit operator double() const
	{
		return double(intVal_)/double(1<<F);
	}
	inline explicit operator int8_t () const
	{
		return int8_t(intVal_ >> F);
	}
	inline explicit operator int16_t () const
	{
		return int16_t(intVal_ >> F);
	}
	inline explicit operator int32_t () const
	{
		return int32_t(intVal_ >> F);
	}
	inline explicit operator int64_t () const
	{
		return int64_t(intVal_ >> F);
	}

	// relational operators

	inline bool operator==(const Q& r) const { return intVal_==r.intVal_; }
	inline bool operator!=(const Q& r) const { return intVal_!=r.intVal_; }
	inline bool operator<=(const Q& r) const { return intVal_<=r.intVal_; }
	inline bool operator>=(const Q& r) const { return intVal_>=r.intVal_; }
	inline bool operator< (const Q& r) const { return intVal_< r.intVal_; }
	inline bool operator> (const Q& r) const { return intVal_> r.intVal_; }

	// arithmetic operators

	inline Q operator-() const { return create(-intVal_); }

	inline Q operator+(const Q& r) const { return create(intVal_+r.intVal_); }
	inline Q operator-(const Q& r) const { return create(intVal_-r.intVal_); }
	inline Q operator*(const Q& r) const
	{
		return create((intB_t)((int2B_t(intVal_) * r.intVal_) >> F));
	}
	inline Q operator/(const Q& r) const
	{
		return create((int2B_t)((int2B_t(intVal_) << F) / r.intVal_));
	}

	inline Q& operator+=(const Q& r) { intVal_ += r.intVal_; return *this; }
	inline Q& operator-=(const Q& r) { intVal_ -= r.intVal_; return *this; }
	inline Q& operator*=(const Q& r)
	{
		intVal_ = intB_t((int2B_t(intVal_)*r.intVal_) >> F);
		return *this;
	}
	inline Q& operator/=(const Q& r)
	{
		intVal_ = intB_t((int2B_t(intVal_) << F) / r.intVal_);
		return *this;
	}

	// raw access

	const intB_t& raw() const { return intVal_; }

	static inline constexpr int bitsI() { return I; }
	static inline constexpr int bitsF() { return F; }

}; // class Q<I,F>

} // namespace fixed

// vim: fenc=utf-8 noet:
