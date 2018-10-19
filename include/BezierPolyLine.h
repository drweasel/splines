#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "ConstraintSet.h"
#include "Vector.h"
#include "Curve.h"
#include "ParametricCurve.h"
#include <vector>
#include <stdexcept>

enum class BezierNodeType { Invalid,LeftCtrl,Node,RightCtrl };

/**
 * Editable Bézier-Spline
 *
 * @author Michael Weitzel <mich@elweitzel.de>
 */
template<typename T>
class BezierPolyLine
{
private:
	/** Interpolated (I) and approximated (A) nodes; pattern (IAA)*I */
	std::vector<Vector<T, 2>> nodes_;
	/** Smoothness flags for all interpolated nodes  */
	std::vector<bool> is_smooth_;
	/** Node constraints */
	//std::vector<ConstraintSet<2,T>> constraints_;

	void segment_to_path(
		std::vector<Vector<T, 2>> & path,
		Vector<T, 2> P0,
		Vector<T, 2> P0r,
		Vector<T, 2> P1l,
		Vector<T, 2> P1,
		const T & collinearity_eps,
		int level) const;

public:
	BezierPolyLine() = default;

	BezierPolyLine(const BezierPolyLine&) = default;
	BezierPolyLine& operator=(const BezierPolyLine&) = default;

	BezierPolyLine(BezierPolyLine&&) = default;
	BezierPolyLine& operator=(BezierPolyLine&&) = default;

	/** Create 2d-Bézier-Poly-Line from an 1d curve */
	BezierPolyLine(const Curve<T>& spline);

	/** Create 2d-Bézier-Poly-Line from a 2d parameteric curve */
	BezierPolyLine(const ParametricCurve<T,2>& spline);

	std::size_t GetNumSegments() const;
	std::size_t GetNumNodes() const;

	void Append(const Vector<T, 2> & xy, bool smooth);

	void Insert(std::size_t k, const Vector<T, 2> & xy, bool smooth);

	/** Subdivide at parameter t and return index of inserted node
		(if successful) */
	std::size_t Subdivide(T t);

	/** Smoothly subdivide segment s by creating a new node s+1.
		Parameter t must be in (0,1) */
	bool Subdivide(std::size_t s, const T & t);

	void Delete(std::size_t k);

	/** Returns the postition of the left control point of node k */
	const Vector<T, 2> & GetLeftCtrl(std::size_t k) const;

	/** Returns the postition of the right control point of node k */
	const Vector<T, 2> & GetRightCtrl(std::size_t k) const;

	/** Returns the postition of node k */
	const Vector<T, 2> & Get(std::size_t k) const;

	const Vector<T, 2> & MoveNodeTo(std::size_t k, const Vector<T, 2> & xy);

	/** Drag the point associated with parameter t to position (x,y) */
	void MoveSegmentTo(T t, const Vector<T, 2> & xy);

	const Vector<T, 2> & MoveLeftCtrlTo(std::size_t k,
		const Vector<T, 2> & xy);

	const Vector<T, 2> & MoveRightCtrlTo(std::size_t k,
		const Vector<T, 2> & xy);

	bool IsSmoothNode(std::size_t k) const;

	bool IsStraightSegment(std::size_t s) const;

	void MakeNodeSmooth(std::size_t k);

	/** Node k */
	void MakeNodeSharp(std::size_t k);

	/** Segment (k,k+1) */
	void MakeSegmentStraight(std::size_t k);

	void MakeSegmentSmooth(std::size_t k);

	void MakeNodeSymmetric(std::size_t k);

	//void AddConstraint(std::size_t k);

	std::vector<Vector<T, 2>> ToPath(
		const T & collinearity_eps
		) const;

	std::pair<int, BezierNodeType> Find(const Vector<T, 2> & p,
			const Vector<T, 2> & max_dist) const;

	void Dump() const;

};

// Douglas-Peucker-like recursion to linearize the Bézier curve
template<typename T>
void BezierPolyLine<T>::segment_to_path(
	std::vector<Vector<T, 2>> & path,
	Vector<T, 2> P1,
	Vector<T, 2> P2,
	Vector<T, 2> P3,
	Vector<T, 2> P4,
	const T & collinearity_eps,
	int level
	) const
{
	constexpr int max_recursion_level = 8;
	if (level > max_recursion_level)
		return;
	if (P1.Dist2(P4) <= collinearity_eps*collinearity_eps)
		return;

	// de Casteljau
	constexpr T half(0.5);
	Vector<T, 2> P12 = (P1 + P2) * half;
	Vector<T, 2> P23 = (P2 + P3) * half;
	Vector<T, 2> P34 = (P3 + P4) * half;
	Vector<T, 2> P123 = (P12 + P23) * half;
	Vector<T, 2> P234 = (P23 + P34) * half;
	Vector<T, 2> P1234 = (P123 + P234) * half;

	// can the curve be linearized with sufficiently accuracy?
	auto distance_point_to_line = [](
			const Vector<T,2>& p1,
			const Vector<T,2>& p2,
			const Vector<T,2>& q
			)->T
	{
		Vector<T,2> d(p2-p1);
		return std::abs( d[1]*q[0] - d[0]*q[1] + p2[0]*p1[1] - p2[1]*p1[0] ) /
			std::sqrt( d[0]*d[0] + d[1]*d[1] );
	};

	if (distance_point_to_line(P1,P4,P1234) <= collinearity_eps)
		return;

	// recursion
	segment_to_path(path,P1,P12,P123,P1234,collinearity_eps,level+1);
	path.push_back(P1234);
	segment_to_path(path,P1234,P234,P34,P4,collinearity_eps,level+1);
}

template<typename T>
BezierPolyLine<T>::BezierPolyLine(
	const Curve<T> & spline
	)
{
	using Vec2 = Vector<T, 2>;
	using Vec3 = Vector<T, 3>;

	for (unsigned int s = 0; s < spline.GetNumSegments(); ++s)
	{
		T sx = spline.GetAbscissa(s+1) - spline.GetAbscissa(s);
		// map polynomial coefficients to [0,1]
		T C[] = { spline.GetCoeff(s,0),spline.GetCoeff(s,1)*sx,
			spline.GetCoeff(s,2)*sx*sx,spline.GetCoeff(s,3)*sx*sx*sx };

		// divide coeffs c1,c2 by binomial coeffs (3 1), (3 2)
		// x(t) = t, y(t) = ..., w(t) = 1

		Vec3 c00(T(0), C[0], T(1));
		Vec3 c10(T(1) / T(3), C[1] / T(3), T(0));
		Vec3 c20(T(0), C[2] / T(3), T(0));
		Vec3 c30(T(0), C[3], T(0));

		Vec3 c01 = c00 + c10;
		Vec3 c11 = c10 + c20;
		Vec3 c21 = c20 + c30;

		Vec3 c02 = c01 + c11;
		Vec3 c12 = c11 + c21;

		Vec3 c03 = c02 + c12;

		// map Bézier control points back to original interval
		LinearMap<T> mp_back(T(0),T(1),
			spline.GetAbscissa(s),spline.GetAbscissa(s+1));
		if (s == 0) {
			c00[0] = mp_back(c00[0] / c00[2]);
			c00[1] /= c00[2];
			nodes_.push_back(Vector<T,2>(c00[0], c00[1]));
		}
		c01[0] = mp_back(c01[0] / c01[2]);
		c01[1] /= c01[2];
		nodes_.push_back(Vector<T,2>(c01[0], c01[1]));
		c02[0] = mp_back(c02[0] / c02[2]);
		c02[1] /= c02[2];
		nodes_.push_back(Vector<T,2>(c02[0], c02[1]));
		c03[0] = mp_back(c03[0] / c03[2]);
		c03[1] /= c03[2];
		nodes_.push_back(Vector<T,2>(c03[0], c03[1]));

		bool C0 = spline.IsContinuous(s, 0, T(1e-5));
		bool C1 = spline.IsContinuous(s, 1, T(1e-5));
		is_smooth_.push_back(C0 && C1);
	}
	is_smooth_.push_back(false);
}

template<typename T>
BezierPolyLine<T>::BezierPolyLine(
	const ParametricCurve<T, 2> & spline
	)
{
	using Vec3 = Vector<T, 3>;
	for (unsigned int s=0; s<spline[0].GetNumSegments(); ++s)
	{
		// divide coeffs c1,c2 by binomial coeffs (3 1), (3 2)
		// x(t) = t, y(t) = ..., w(t) = 1
		Vec3 c00(spline[0].GetCoeff(s,0)     , spline[1].GetCoeff(s,0)     , T(1));
		Vec3 c10(spline[0].GetCoeff(s,1)/T(3), spline[1].GetCoeff(s,1)/T(3), T(0));
		Vec3 c20(spline[0].GetCoeff(s,2)/T(3), spline[1].GetCoeff(s,2)/T(3), T(0));
		Vec3 c30(spline[0].GetCoeff(s,3)     , spline[1].GetCoeff(s,3)     , T(0));

		Vec3 c01 = c00 + c10;
		Vec3 c11 = c10 + c20;
		Vec3 c21 = c20 + c30;

		Vec3 c02 = c01 + c11;
		Vec3 c12 = c11 + c21;

		Vec3 c03 = c02 + c12;

		if (s == 0)
		{
			c00[0] /= c00[2];
			c00[1] /= c00[2];
			nodes_.push_back(Vector<T,2>(c00[0],c00[1]));
		}
		c01[0] /= c01[2];
		c01[1] /= c01[2];
		nodes_.push_back(Vector<T,2>(c01[0],c01[1]));
		c02[0] /= c02[2];
		c02[1] /= c02[2];
		nodes_.push_back(Vector<T,2>(c02[0],c02[1]));
		c03[0] /= c03[2];
		c03[1] /= c03[2];
		nodes_.push_back(Vector<T,2>(c03[0],c03[1]));

		bool C0 = spline[0].IsContinuous(s, 0, T(1e-5));
		bool C1 = spline[0].IsContinuous(s, 1, T(1e-5));
		is_smooth_.push_back(C0 && C1);
	}
	is_smooth_.push_back(false);
}

template<typename T>
std::size_t BezierPolyLine<T>::GetNumSegments() const
{
	return (is_smooth_.size() <= 1) ? 0 : is_smooth_.size()-1;
}

template<typename T>
std::size_t BezierPolyLine<T>::GetNumNodes() const
{
	return is_smooth_.size();
}

template<typename T>
void BezierPolyLine<T>::Append(
	const Vector<T,2> & xy,
	bool smooth
	)
{
	if (nodes_.empty())
		nodes_.emplace_back(xy);
	else
	{
		constexpr T oneThird = T(1)/T(3);
		constexpr T twoThird = T(2)/T(3);
		if (!smooth || GetNumNodes() == 1)
		{
			// non-smooth or just one point => no smoothing
			const Vector<T,2>& p0 = nodes_.back();
			Vector<T,2> p3(xy);
			Vector<T,2> p1 = p0.Lerp(p3,oneThird);
			Vector<T,2> p2 = p0.Lerp(p3,twoThird);

			nodes_.push_back(p1);
			nodes_.push_back(p2);
			nodes_.push_back(p3);
		}
		else
		{
			// smooth and GetNumNodes() > 1
			const Vector<T,2>& p1m = nodes_[nodes_.size()-2];
			const Vector<T,2>& p0 = nodes_.back();
			Vector<T,2> p3(xy);

			// p1 is obtained by mirroring on p0
			Vector<T,2> p1 = p0 + (p0-p1m);
			Vector<T,2> p2 = p0.Lerp(p3,twoThird);

			nodes_.push_back(p1);
			nodes_.push_back(p2);
			nodes_.push_back(p3);
		}
	}
	is_smooth_.emplace_back(smooth);
}

template<typename T>
void BezierPolyLine<T>::Insert(
	std::size_t k,
	const Vector<T, 2> & xy,
	bool smooth
	)
{
	if (k>GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	if (k == GetNumNodes())
	{
		Append(xy,smooth);
		return;
	}

	constexpr T oneThird = T(1)/T(3);
	constexpr T twoThird = T(2)/T(3);
	if (k > 0)
	{
		Vector<T,2> pp = Get(k-1);
		Vector<T,2> pq(xy);
		Vector<T,2> pr = Get(k);

		nodes_.insert(nodes_.begin()+3*k-1,pq.Lerp(pr,oneThird));
		nodes_.insert(nodes_.begin()+3*k-1,pq);
		nodes_.insert(nodes_.begin()+3*k-1,pp.Lerp(pq,twoThird));
	}
	else
	{
		Vector<T,2> pp(xy);
		Vector<T,2> pq = nodes_[0];

		nodes_.insert(nodes_.begin(),pp.Lerp(pq,twoThird));
		nodes_.insert(nodes_.begin(),pp.Lerp(pq,oneThird));
		nodes_.insert(nodes_.begin(),pp);
	}
	is_smooth_.insert(is_smooth_.begin()+k, smooth);
	if (smooth)
		MakeNodeSmooth(k);
}

template<typename T>
std::size_t BezierPolyLine<T>::Subdivide(T t)
{
	std::size_t s = (std::size_t)t;
	if (Subdivide(s,t-s))
		return s+1;
	return s;
}

template<typename T>
bool BezierPolyLine<T>::Subdivide(
	std::size_t s,
	const T & t
	)
{
	if (s >= GetNumSegments() || t <= T(0) || t >= T(1))
		return false;

	T it(1-t);
	Matrix<T, 4, 4> Dl{
		{1, 0, 0, 0},
		{it, t, 0, 0},
		{it * it, 2 * t * it, t * t, 0},
		{it * (it * it), 3 * t * (it * it), 3 * (t * t) * it, t * (t * t)}
		};
	Matrix<T, 4, 4> Dr{
		{it * (it * it), 3 * (it * it) * t, 3 * it * (t * t), t * (t * t)},
		{0, it * it, 2 * it * t, t * t},
		{0, 0, it, t},
		{0, 0, 0, 1}
		};

	const Vector<T,2>& P0 = nodes_[3*s]; // Get(s)
	const Vector<T,2>& P1 = nodes_[3*s+1]; // GetRightCtrl(s)
	const Vector<T,2>& P2 = nodes_[3*(s+1)-1]; //GetLeftCtrl(s+1)
	const Vector<T,2>& P3 = nodes_[3*(s+1)]; //Get(s+1)

	Matrix<T,4,2> P { {P0[0],P0[1]},{P1[0],P1[1]},{P2[0],P2[1]},{P3[0],P3[1]} };

	// Pl = (P0,P11,P12,P13)
	Matrix<T,4,2> Pl = Dl*P;
	// Pr = (P13,P22,P31,P3)
	Matrix<T,4,2> Pr = Dr*P;

	nodes_[3*s+1] = Vector<T,2>{ Pl(1,0),Pl(1,1) }; // adapt GetRightCtrl(s) of P11
	nodes_[3*(s+1)-1] = Vector<T,2>{ Pr(2,0),Pr(2,1) }; // adapt GetLeftCtrl(s+1) of P31
	nodes_.insert(nodes_.begin()+3*s+2,{
		Vector<T,2>{ Pl(2,0),Pl(2,1) }, // P12 (leftCtrl of P13)
		Vector<T,2>{ Pl(3,0),Pl(3,1) }, // new vertex P13
		Vector<T,2>{ Pr(1,0),Pr(1,1) }  // P22 (rightCtrl of P13)
		});
	// the newly created vertex is always smooth
	is_smooth_.insert(is_smooth_.begin()+s+1,true);
	return true;
}

template<typename T>
void BezierPolyLine<T>::Delete(std::size_t k)
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	if (GetNumNodes() == 1)
	{
		is_smooth_.clear();
		nodes_.clear();
	}
	else if (k == 0)
	{
		is_smooth_.erase(is_smooth_.begin());
		nodes_.erase(nodes_.begin());
		nodes_.erase(nodes_.begin());
		nodes_.erase(nodes_.begin());
	}
	else if (k + 1 == GetNumNodes())
	{
		is_smooth_.erase(is_smooth_.end()-1);
		nodes_.erase(nodes_.end()-1);
		nodes_.erase(nodes_.end()-1);
		nodes_.erase(nodes_.end()-1);
	}
	else
	{
		bool was_smooth = is_smooth_[k];
		is_smooth_.erase(is_smooth_.begin() + k);
		nodes_.erase(nodes_.begin() + (3*k-1));
		nodes_.erase(nodes_.begin() + (3*k-1));
		nodes_.erase(nodes_.begin() + (3*k-1));
		if (was_smooth)
			MakeSegmentSmooth(k-1);
		else
			MakeSegmentStraight(k-1);
	}
}

template<typename T>
const Vector<T, 2> & BezierPolyLine<T>::GetLeftCtrl(std::size_t k) const
{
	if (k < 1 || k + 1 > GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	return nodes_[3*k-1];
}

template<typename T>
const Vector<T, 2> & BezierPolyLine<T>::GetRightCtrl(std::size_t k) const
{
	if (k + 1 >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	return nodes_[3*k+1];
}

template<typename T>
const Vector<T, 2> & BezierPolyLine<T>::Get(std::size_t k) const
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	return nodes_[3*k];
}

template<typename T>
const Vector<T, 2> & BezierPolyLine<T>::MoveNodeTo(
	std::size_t k,
	const Vector<T, 2> & xy
	)
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	// parallel shift of the control points
	Vector<T,2> p(xy);
	Vector<T,2> d = p - nodes_[3*k];
	nodes_[3*k] = p;
	if (k>0)
		nodes_[3*k-1] += d;
	if (k+1 < GetNumNodes())
		nodes_[3*k+1] += d;
	return nodes_[3*k];
}

template<typename T>
void BezierPolyLine<T>::MoveSegmentTo(
	T t,
	const Vector<T, 2> & xy
	)
{
	std::size_t k = (std::size_t)t;
	t -= T(k);

	// Compute u (ratio between start- and end-node of segment s)
	T t2 = t * t;
	T t3 = t2 * t;
	T it = T(1) - t;
	T it2 = it * it;
	T it3 = it2 * it;
	T u = it3 / (t3 + it3);

	//LOG_INFO(L"k=%zi t=%g, u=%g x=%g y=%g", k, t, u, xy[0], xy[1]);

	const Vector<T, 2> & P0 = nodes_[3 * k];
	Vector<T, 2> & P1 = nodes_[3 * k + 1];
	Vector<T, 2> & P2 = nodes_[3 * (k + 1) - 1];
	const Vector<T, 2> & P3 = nodes_[3 * (k + 1)];

	// Compute C, the point at ratio u between start- and end-node
	Vector<T, 2> C = P0 * u + P3 * (T(1) - u);

	// Compute B, the point on the curve associated with t
	Vector<T, 2> B =
		P0 * it3 + P1 * (T(3) * t * it2) + P2 * (T(3) * t2 * it) + P3 * t3;

	//LOG_WARN(L"B: %g %g", B[0],B[1]);

	// Compute A, the point on the "hat" of B
	T v = std::abs((t3 + it3 - T(1)) / (t3 + it3));
	Vector<T, 2> A = B + (B - C) / v;

	// Destination point
	Vector<T, 2> Bnew(xy);
	// Corresponding A
	Vector<T, 2> Anew = Bnew + (B - C) / v;

	// Compute new control points P1, P2
	Vector<T, 2> e1 = (P0 * it + P1 * t) * it + (P1 * it + P2 * t) * t;
	Vector<T, 2> e2 = (P1 * it + P2 * t) * it + (P2 * it + P3 * t) * t;
	Vector<T, 2> v1 = Anew + (e1 - Anew) / t;
	Vector<T, 2> v2 = Anew + (e2 - Anew) / it;
	P1 = v1 + (v1 - P0) / t;
	P2 = v2 + (v2 - P3) / it;
}

template<typename T>
const Vector<T, 2> & BezierPolyLine<T>::MoveLeftCtrlTo(
	std::size_t k,
	const Vector<T, 2> & xy
	)
{
	if (k<1 || k+1>GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	if (!is_smooth_[k] || k+1==GetNumNodes())
		return nodes_[3*k-1] = xy;
	T lenR = nodes_[3*k].Dist(nodes_[3*k+1]);
	nodes_[3*k+1] = nodes_[3*k] - (nodes_[3*k].Dir(nodes_[3*k-1]))*lenR;

	return nodes_[3*k-1] = xy;
}

template<typename T>
const Vector<T, 2> & BezierPolyLine<T>::MoveRightCtrlTo(std::size_t k,
	const Vector<T, 2> & xy
	)
{
	if (k+1 >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	if (!is_smooth_[k] || k==0)
		return nodes_[3*k+1] = xy;
	T lenL = nodes_[3*k].Dist(nodes_[3*k-1]);
	nodes_[3*k-1] = nodes_[3*k] - (nodes_[3*k].Dir(nodes_[3*k+1]))*lenL;

	return nodes_[3*k+1] = xy;
}

template<typename T>
bool BezierPolyLine<T>::IsSmoothNode(std::size_t k) const
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	return is_smooth_[k];
}

template<typename T>
bool BezierPolyLine<T>::IsStraightSegment(std::size_t s) const
{
	if (s >= GetNumSegments())
		throw std::invalid_argument("index out of bounds");
	const Vector<T, 2> & p = Get(s);
	const Vector<T, 2> & q = Get(s+1);
	const Vector<T, 2> & l = GetRightCtrl(s);
	const Vector<T, 2> & r = GetLeftCtrl(s + 1);
	Vector<T, 2> pq = q-p;
	return l.IsOnLine(p, pq, T(1e-5)) && l.LerpOnLine(p, pq) >= T(0) &&
		r.IsOnLine(q, pq, T(1e-5)) && r.LerpOnLine(q, pq) <= T(0);
}

template<typename T>
void BezierPolyLine<T>::MakeNodeSmooth(std::size_t k)
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	if (k==0 || k+1==GetNumNodes())
		return;
	is_smooth_[k] = true;
	Vector<T,2> dir = nodes_[3*k-1].Dir(nodes_[3*k+1]);
	T lenL = nodes_[3*k-1].Dist(nodes_[3*k]);
	T lenR = nodes_[3*k+1].Dist(nodes_[3*k]);
	nodes_[3*k-1] = nodes_[3*k] - dir*lenL;
	nodes_[3*k+1] = nodes_[3*k] + dir*lenR;
}

template<typename T>
void BezierPolyLine<T>::MakeNodeSharp(std::size_t k)
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	is_smooth_[k] = false;
}

/** Segment (s,s+1) */
template<typename T>
void BezierPolyLine<T>::MakeSegmentStraight(std::size_t s)
{
	if (s >= GetNumSegments())
		throw std::invalid_argument("index out of bounds");
	is_smooth_[s] = is_smooth_[s+1] = false;
	constexpr T oneThird = T(1)/T(3);
	constexpr T twoThird = T(2)/T(3);
	nodes_[3*s+1] = nodes_[3*s].Lerp(nodes_[3*(s+1)],oneThird);
	nodes_[3*s+2] = nodes_[3*s].Lerp(nodes_[3*(s+1)],twoThird);
}

/** Segment (s,s+1) */
template<typename T>
void BezierPolyLine<T>::MakeSegmentSmooth(std::size_t s)
{
	if (s >= GetNumSegments())
		throw std::invalid_argument("index out of bounds");
	T len = nodes_[3*s].Dist(nodes_[3*(s+1)])/T(3);
	if (s == 0)
		nodes_[3*s+1] = nodes_[3*s].Lerp(nodes_[3*(s+1)],T(1)/T(3));
	else
		nodes_[3*s+1] = nodes_[3*s] + nodes_[3*s-1].Dir(nodes_[3*s])*len;

	if (s+1 == GetNumSegments())
		nodes_[3*s+2] = nodes_[3*(s+1)].Lerp(nodes_[3*s],T(1)/T(3));
	else
		nodes_[3*s+2] = nodes_[3*(s+1)] + nodes_[3*(s+1)+1].Dir(nodes_[3*(s+1)])*len;
}

template<typename T>
void BezierPolyLine<T>::MakeNodeSymmetric(std::size_t k)
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");
	if (k==0 || k+1==GetNumNodes())
		return;
	Vector<T,2> dir = nodes_[3*k-1].Dir(nodes_[3*k+1]);
	T lenL = nodes_[3*k-1].Dist(nodes_[3*k]);
	T lenR = nodes_[3*k+1].Dist(nodes_[3*k]);
	T len = (lenL + lenR)/T(2);
	nodes_[3*k-1] = nodes_[3*k] - dir*len;
	nodes_[3*k+1] = nodes_[3*k] + dir*len;

}

/*
template<typename T>
void BezierPolyLine<T>::AddConstraint(std::size_t k)
{
	if (k >= GetNumNodes())
		throw std::invalid_argument("index out of bounds");

	if (constraints_.size()
}
*/

template<typename T>
std::vector<Vector<T,2>> BezierPolyLine<T>::ToPath(
	const T& collinearity_eps
	) const
{
	std::vector<Vector<T,2>> path;
	if (GetNumSegments() == 0)
		return path;

	for (std::size_t s=0; s<GetNumSegments(); ++s)
	{
		path.push_back(Get(s));
		segment_to_path(path,
			Get(s), GetRightCtrl(s), GetLeftCtrl(s + 1), Get(s + 1),
			collinearity_eps, 0);
	}
	path.push_back(Get(GetNumSegments()));
	return path;
}

template<typename T>
std::pair<int, BezierNodeType> BezierPolyLine<T>::Find(
	const Vector<T, 2> & p,
	const Vector<T, 2> & max_dist
	) const
{
	T d2_min = std::nextafter<T>(T(1), T(2));
	auto sqr = [](const T & v) -> T { return v * v; };
	T r2_x = sqr(max_dist[0]);
	T r2_y = sqr(max_dist[0]);

	int min_idx = -1;
	int idx = 0;
	for (const Vector<T, 2> & q : nodes_)
	{
		// ellipsoid distance
		T d2 = sqr(p[0] - q[0]) / r2_x + sqr(p[1] - q[1]) / r2_y;
		if (d2 < d2_min) {
			d2_min = d2;
			min_idx = idx;
		}
		idx++;
	}

	BezierNodeType type = BezierNodeType::Invalid;
	if (min_idx >= 0)
	{
		switch (min_idx % 3)
		{
		case 0:
			type = BezierNodeType::Node;
			min_idx /= 3;
			break;
		case 1:
			type = BezierNodeType::RightCtrl;
			min_idx = (min_idx - 1) / 3;
			break;
		case 2:
			type = BezierNodeType::LeftCtrl;
			min_idx = (min_idx + 1) / 3;
			break;
		default: break;
		}
	}
	return {type == BezierNodeType::Invalid ? -1 : min_idx, type};
}

template<typename T>
void BezierPolyLine<T>::Dump() const
{
	for (std::size_t k=0; k<GetNumNodes(); ++k)
	{
		if (k>0)
			printf("(%g,%g)--", GetLeftCtrl(k)[0], GetLeftCtrl(k)[2]);

		printf("[%g,%g;%c]", Get(k)[0], Get(k)[1], IsSmoothNode(k)?'S':'^');

		if (k+1<GetNumNodes())
			printf("--(%g,%g)--", GetRightCtrl(k)[0], GetRightCtrl(k)[1]);
	}
	printf("\n");
}

// vim: fenc=utf-8 noet:
