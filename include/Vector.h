#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include <array>
#include <cmath>

namespace detail
{

template<typename T, size_t I, size_t... J> struct ndarray_ctor
{
	using type = std::array<typename ndarray_ctor<T,J...>::type,I>;
};
template<typename T, size_t I> struct ndarray_ctor<T,I>
{
	using type = std::array<T,I>;
};
template<typename T, size_t I, size_t... J> using ndarray =
	typename ndarray_ctor<T,I,J...>::type;


template<typename T, size_t D, size_t... Ds>
class ndview
{
private:
	ndarray<T,D,Ds...> * ptr_;

	template<typename NDArray, typename... Args>
	inline T & unfold(NDArray & A, int i, Args... args)
	{
		return unfold(A[i], std::forward<Args>(args)...);
	}
	inline T & unfold(T & v) { return v; }

public:
	ndview(T * ptr) : ptr_((ndarray<T,D,Ds...>*)ptr) { }

	template< typename... Args > inline T& operator()(Args... args)
	{
		return unfold(*ptr_,std::forward<Args>(args)...);
	}

}; // class ndview<T,D...>

template<typename T, size_t D, size_t... Ds>
class ndmatrix : public ndarray<T,D,Ds...>
{
private:
	template<typename NDArray, typename... Args>
	inline T & unfold(NDArray & A, size_t i, Args... args)
	{
		return unfold(A[i], std::forward<Args>(args)...);
	}
	inline T & unfold(T & v) { return v; }

	template<typename NDArray, typename... Args>
	inline const T & unfold(NDArray & A, size_t i, Args... args) const
	{
		return unfold(A[i], std::forward<Args>(args)...);
	}
	inline const T & unfold(const T & v) const { return v; }

public:
	using ndarray<T,D,Ds...>::ndarray;

	template< typename... Args > inline T& operator()(Args... args)
	{
		return unfold(*this,std::forward<Args>(args)...);
	}
	template< typename... Args > inline const T& operator()(Args... args) const
	{
		return unfold(*this,std::forward<Args>(args)...);
	}

};

} // namespace detail

template<typename T, size_t I, size_t J> using MatrixView = detail::ndview<T,I,J>;
//template<typename T, size_t I, size_t J> using Matrix = detail::ndmatrix<T,I,J>;

template<typename T, size_t I, size_t J> class Matrix
	: public detail::ndmatrix<T,I,J>
{
public:
	Matrix() = default;

	Matrix(std::initializer_list<std::initializer_list<T>> values)
	{
		T * ptr = &(this->operator()(0,0));
		for (const auto & row : values)
			for (T value : row)
				*ptr++ = value;
	}

	Matrix(std::initializer_list<T> values)
	{
		T * ptr = &(this->operator()(0,0));
		for (T value : values)
			*ptr++ = value;
	}

	template< size_t K, size_t L >
	Matrix<T,I,L> operator*(const Matrix<T,K,L>& B) const
	{
		static_assert(J==K,"invalid arguments for matrix multiplication");
		Matrix<T,I,L> C;
		for (size_t i=0; i<I; ++i)
			for (size_t j=0; j<L; ++j)
			{
				T c_ij(0);
				for (size_t k=0; k<J; ++k)
					c_ij += this->operator()(i,k) * B(k,j);
				C(i,j) = c_ij;
			}
		return C;
	}

	inline T& operator[](size_t i)
	{
		return *(&(this->operator()(0,0)) + i);
	}

};

template<typename T, unsigned int D>
class Vector
{
private:
	std::array<T,D> v_;

public:
	Vector() = default;

	Vector(const Vector&) = default;
	Vector& operator=(const Vector&) = default;

	Vector(Vector&&) = default;
	Vector& operator=(Vector&&) = default;

	inline Vector(const std::array<T,D>& v) : v_(v) { }
	inline Vector& operator=(const std::array<T,D>& v) { v_ = v; }
	inline operator std::array<T, D> &() { return v_; }
	inline operator const std::array<T, D> &() const { return v_; }

	inline Vector(const T & value)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] = value;
	}

	template < typename... Args > Vector( Args... args )
		: v_{ static_cast< T >(std::forward<Args>( args ) )... } { }

	inline T& operator[](unsigned int k) { return v_[k]; }
	inline const T& operator[](unsigned int k) const { return v_[k]; }
	inline T& operator()(unsigned int k) { return v_[k]; }
	inline const T& operator()(unsigned int k) const { return v_[k]; }

	inline bool operator==(const Vector& b) const { return v_==b.v_; }
	inline bool operator!=(const Vector& b) const { return v_!=b.v_; }

	inline T Cross(const Vector<T, 2> & b)
	{
		return v_[0] * b.v_[1] - v_[1] * b.v_[0];
	}

	inline Vector<T, 3> Cross(const Vector<T, 3> & b)
	{
		return Vector<T, 3>{
			v_[1] * b.v_[2] - v_[2] * b.v_[1],
			v_[2] * b.v_[0] - v_[0] * b.v_[2],
			v_[0] * b.v_[1] - v_[1] * b.v_[0] };
	}

	inline Vector operator+(const Vector & r) const
	{
		Vector v;
		for (unsigned int d = 0; d < D; ++d) v[d] = v_[d] + r[d];
		return v;
	}

	inline Vector & operator+=(const Vector & r)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] += r[d];
		return *this;
	}

	inline Vector operator-(const Vector & r) const
	{
		Vector v;
		for (unsigned int d = 0; d < D; ++d) v[d] = v_[d] - r[d];
		return v;
	}

	inline Vector & operator-=(const Vector & r)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] -= r[d];
		return *this;
	}

	inline Vector operator-() const
	{
		Vector u;
		for (unsigned int d = 0; d < D; ++d) u.v_[d] = -v_[d];
		return u;
	}

	inline Vector operator*(const T & f) const
	{
		Vector v;
		for (unsigned int d = 0; d < D; ++d) v[d] = v_[d] * f;
		return v;
	}

	friend inline Vector operator*(const T & f, const Vector & b)
	{
		return b * f;
	}

	inline Vector & operator*=(const T & f)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] *= f;
		return *this;
	}

	inline Vector operator/(const T & f) const
	{
		Vector v;
		for (unsigned int d = 0; d < D; ++d) v[d] = v_[d] / f;
		return v;
	}

	inline Vector & operator/=(const T & f)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] /= f;
		return *this;
	}

	inline T operator*(const Vector & r) const
	{
		T dp(v_[0] * r[0]);
		for (unsigned int d = 1; d < D; ++d) dp += v_[d] * r[d];
		return dp;
	}

	inline Vector Scale(const Vector & b) const
	{
		Vector s;
		for (unsigned int d = 0; d < D; ++d) s.v_[d] = v_[d] * b._[d];
		return s;
	}

	inline T Norm2() const { return (*this) * (*this); }

	inline T Norm() const { return std::sqrt(Norm2()); }

	inline T Dist2(const Vector & r) const { return (*this - r).Norm2(); }

	inline T Dist(const Vector & r) const { return std::sqrt(Dist2(r)); }

	inline Vector Normalize() const { return *this / Norm(); }

	inline T Sum() const
	{
		T s(0);
		for (unsigned int d = 0; d < D; ++d) s += v_[d];
		return s;
	}

	inline Vector Min(const Vector & b) const
	{
		Vector m;
		for (unsigned int d = 0; d < D; ++d) m.v_[d] = std::min(v_[d], b.v_[d]);
		return m;
	}

	inline Vector Max(const Vector & b) const
	{
		Vector m;
		for (unsigned int d = 0; d < D; ++d) m.v_[d] = std::max(v_[d], b.v_[d]);
		return m;
	}

	inline Vector & MinMutable(const Vector & b)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] = std::min(v_[d], b.v_[d]);
		return *this;
	}

	inline Vector & MaxMutable(const Vector & b)
	{
		for (unsigned int d = 0; d < D; ++d) v_[d] = std::max(v_[d], b.v_[d]);
		return *this;
	}

	inline T Max() const
	{
		T v(v_[0]);
		for (unsigned int d = 1; d < D; ++d)
			if (v_[d] > v) v = v_[d];
		return v;
	}

	inline T Min() const
	{
		T v(v_[0]);
		for (int i = 1; i < D; ++i)
			if (v_[i] < v) v = v_[i];
		return v;
	}

	inline Vector Abs() const
	{
		Vector a;
		for (unsigned int d = 0; d < D; ++d) a.v_[d] = std::abs(v_[d]);
		return a;
	}

	inline Vector Floor() const
	{
		Vector f;
		for (unsigned int d = 0; d < D; ++d) f.v_[d] = std::floor(v_[d]);
		return f;
	}

	inline Vector Ceil() const
	{
		Vector c;
		for (unsigned int d = 0; d < D; ++d) c.v_[d] = std::ceil(v_[d]);
		return c;
	}

	inline Vector Round() const
	{
		Vector r;
		for (unsigned int d = 0; d < D; ++d) r.v_[d] = std::round(v_[d]);
		return r;
	}

	inline Vector Dir(const Vector & r) const
	{
		return (r - *this).Normalize();
	}

	inline Vector ProjectOn(const Vector & b) const
	{
		// a_prj = (a*b)*b/(|b|^2).
		return ((*this * b) / b.Norm2()) * b;
	}

	/** Projects the vector to line p+v.t */
	inline Vector ProjectToLine(const Vector & p, const Vector & v) const
	{
		return p + (*this - p).ProjectOn(v);
	}

	/** Projects the vector to line segment (p,q) */
	inline Vector ProjectToSegment(const Vector & p, const Vector & q) const
	{
		return ProjectToLine(p, q - p);
	}

	/** Projects the vector into a plane (p,n); assumes ||n||=1 */
	inline Vector ProjectToPlane(const Vector & p, const Vector & n) const
	{
		return (*this) - n * ((*this - p) * n);
	}

	/** Linear interpolation with parameter t to point b */
	inline Vector Lerp(const Vector & b, const T & t) const
	{
		return *this + (b - *this) * t;
	}

	/** Project to line p+v.t and compute the lerp parameter */
	inline T LerpOnLine(const Vector & p, const Vector & v) const
	{
		Vector u(ProjectToLine(p, v) - p);
		T lerp = (v * u) / (v * v);
		if (lerp != lerp) return T(.5); // repair NaN
		return lerp;
	}

	/** Project to line segment (p,q) and compute the lerp parameter */
	inline T LerpOnSegment(const Vector & p, const Vector & q) const
	{
		return LerpOnLine(p, q - p);
	}

	/** Returns true, if distance to line (p+v.t) is <= d_max */
	inline bool IsOnLine(const Vector & p, const Vector & v, T d_max) const
	{
		return Dist2(ProjectToLine(p, v)) <= d_max * d_max;
	}

	/** Returns true, if vector is somewhere between p and q */
	inline bool IsBetween(const Vector & p, const Vector & q, T tol) const
	{
		T t = LerpOnSegment(p, q);
		return t >= -tol && t <= (T(1) + tol);
	}

	/** Returns true, if vector is on a line segment (p,q) */
	inline bool IsOnSegment(const Vector & p, const Vector & q, T tol) const
	{
		return IsOnLine(p, q - p, tol) && IsBetween(p, q, tol);
	}

	inline bool IsAnyNaN() const
	{
		bool bAnyNaN = false;
		for (int j = 0; j < D && !bAnyNaN; ++j) bAnyNaN |= std::isnan(v_[j]);
		return bAnyNaN;
	}

	inline bool IsAnyInf() const
	{
		bool i = false;
		for (unsigned int d = 0; d < D && !i; ++d) i |= std::isinf(v_[d]);
		return i;
	}

	inline bool IsFinite() const
	{
		bool f = true;
		for (unsigned int d = 0; d < D && f; ++d) f = f && std::isfinite(v_[d]);
		return f;
	}

	inline bool IsZero() const
	{
		bool z = true;
		for (unsigned int d = 0; d < D && z; ++d) z = z && (v_[d] == T(0));
		return z;
	}

};

// vim: fenc=utf-8 noet:
