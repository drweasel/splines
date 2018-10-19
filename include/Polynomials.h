#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include <complex>
#include <limits>

/**
 * Francis' implicitly shifted QR algorithm.
 * Computes all (complex) eigenvalues of an upper Hessenberg matrix H.
 *
 * @param	H	[in] upper n x n Hessenberg matrix
 * @param	n	[in] dimension
 * @param	W	[out] computed eigenvalues
 */
template<typename T> bool HessenbergEig(
	T* H,
	int n,
	std::complex<T>* W
	)
{
	auto sign = [](const T& a, const T& b)->T
	{
		return b >= T(0) ? std::abs(a) : -std::abs(a);
	};

	T Hnorm(0);
	// Compute matrix norm for possible use in locating single small
	// subdiagonal element.
	for (int i=0; i<n; ++i)
	{
		W[i] = std::complex<T>(T(0));
		for (int j=std::max(i-1,0); j<n; ++j)
			Hnorm += std::abs(H[i*n+j]);
	}

	int nn(n);
	T t(0); // Gets changed only by an exceptional shift.
	while (nn >= 1)
	{
		// Begin search for next eigenvalue.
		int iters = 0;
		int l;
		do
		{
			for (l=nn; l>=2; --l)
			{
				// Begin iteration: look for single small
				// subdiagonal element.
				T s = std::abs(H[(l-2)*n+(l-2)]) + std::abs(H[(l-1)*n+(l-1)]);
				if (s == T(0))
					s = Hnorm;
				if (std::abs(H[(l-1)*n+(l-2)]) + s == s)
					break;
			}
			T x = H[(nn-1)*n+(nn-1)];

			if (l == nn)
			{
				// One root found.
				W[nn-1] = std::complex<T>(x+t,T(0));
				nn--;
				continue;
			}

			T z,y,w,r,q,p;
			y = H[(nn-2)*n+(nn-2)];
			w = H[(nn-1)*n+(nn-2)] * H[(nn-2)*n+(nn-1)];
			if (l == nn-1)
			{
				// Two roots found...
				p = T(0.5)*(y-x);
				q = p*p + w;
				z = std::sqrt(std::abs(q));
				x += t;
				if (q >= T(0))
				{
					// ...a real pair.
					z = p + sign(z,p);
					W[nn-2] = std::complex<T>(x+z,T(0));
					W[nn-1] = z != T(0) ? std::complex<T>(x-w/z,T(0)) : W[nn-2];
				}
				else
				{
					// ...a complex pair.
					W[nn-1] = std::complex<T>(x+p,z);
					W[nn-2] = std::conj<T>(W[nn-1]);
				}
				nn -= 2;
			}
			else
			{
				// No roots found. Continue iteration.
				if (iters == 30)
					return false; // too many interations
				if (iters == 10 || iters == 20)
				{
					// Form exceptional shift.
					t += x;
					for (int i=0; i<nn; ++i)
						H[i*n+i] -= x;
					T s = std::abs(H[(nn-1)*n+(nn-2)]) + std::abs(H[(nn-2)*n+(nn-3)]);
					y = x = T(0.75)*s;
					w = T(-0.4375)*s*s;
				}
				++iters;

				int m;
				for (m=nn-2; m>=l; --m)
				{
					// Form shift and then look for 2 consecutive
					// small subdiagonal elements.
					z = H[(m-1)*n+(m-1)];
					r = x-z;
					T s = y-z;
					p = (r*s-w)/H[m*n+(m-1)] + H[(m-1)*n+m]; // Equation (11.6.23).
					q = H[m*n+m] - z - r - s;
					r = H[(m+1)*n+m];
					s = std::abs(p) + std::abs(q) + std::abs(r);
					// Scale to prevent overflow or underflow.
					p /= s;
					q /= s;
					r /= s;

					if (m == l)
						break;

					T u = std::abs(H[(m-1)*n+(m-2)]) * (std::abs(q)+std::abs(r));
					T v = std::abs(p) * (std::abs(H[(m-2)*n+(m-2)]) + std::abs(z)+std::abs(H[m*n+m]));

					if (u+v == v)
						break; // Equation (11.6.26).
				}

				for (int i=m+1; i<nn; ++i)
				{
					H[i*n+(i-2)] = T(0);
					if (i != m+1)
						H[i*n+(i-3)] = T(0);
				}

				for (int k=m; k<=nn-1; ++k)
				{
					// Double QR step on rows l to nn and columns m to nn.
					if (k != m)
					{
						p = H[(k-1)*n+(k-2)];
						// Begin setup of Householder vector.
						q = H[k*n+(k-2)];
						r = T(0);
						if (k != nn-1)
							r = H[(k+1)*n+(k-2)];
						if ((x=std::abs(p)+std::abs(q)+std::abs(r)) != T(0))
						{
							// Scale to prevent overflow or underflow.
							p /= x;
							q /= x;
							r /= x;
						}
					}

					T s(sign(std::sqrt(p*p + q*q + r*r),p));
					if (s != T(0))
					{
						if (k == m)
						{
							if (l != m)
								H[(k-1)*n+(k-2)] = -H[(k-1)*n+(k-2)];
						}
						else
							H[(k-1)*n+(k-2)] = -s*x;
						p += s;
						// Equations (11.6.24).
						x = p/s;
						y = q/s;
						z = r/s;
						q /= p;
						r /= p;
						for (int j=k-1; j<nn; ++j)
						{
							// Row modification.
							p = H[(k-1)*n+j] + q * H[k*n+j];
							if (k != nn-1)
							{
								p += r*H[(k+1)*n+j];
								H[(k+1)*n+j] -= p*z;
							}
							H[k*n+j] -= p*y;
							H[(k-1)*n+j] -= p*x;
						}

						int mmin = nn<k+3 ? nn : k+3;
						for (int i=l-1; i<mmin; ++i)
						{
							// Column modification.
							p = x * H[i*n+(k-1)] + y * H[i*n+k];
							if (k != nn-1)
							{
								p += z * H[i*n+(k+1)];
								H[i*n+(k+1)] -= p*r;
							}
							H[i*n+k] -= p*q;
							H[i*n+(k-1)] -= p;
						}
					}
				}
			}

		}
		while (l < nn-1);
	}
	return true;
}

/**
 * Determines all N-1 complex roots of a polynomial of degree N:
 *   p(x) = a(0).x^0 + a(1).x^1 + ... + a(N-1).x^(N-1) + a(N).x^N
 *
 * The algortihm computes the eigenvalues of the (N-1)x(N-1)
 * Companion-Matrix C, an upper Hessenberg matrix:
 *
 *     | 0 0 ... 0   -a(0)/a(N) |
 *     | 1 0 ... 0   -a(1)/a(N) |
 * C = | 0 1 ... 0   -a(2)/a(N) |
 *     | ... ... ...     ...    |
 *     | 0 0 ... 1 -a(N-1)/a(N) |
 *
 * @param	p	[in] N coefficients of the polynomial
 * @param	N	[in] number of polynomial coefficients (degree +1)
 * @param	W	[out] N-1 complex roots
 * @param	work	room for (N-1)*(N-1) temporary values
 * @return	true, on success
 */
template<typename T> bool PolynomialRoots(
	const T* p,
	int N,
	std::complex<T>* W,
	T* work
	)
{
	T* C = work;
	for (int k=0, M=(N-1)*(N-1); k<M; ++k)
		C[k] = T(0);
	--N;
	// Companion-Matrix:
	for (int i=1; i<N; ++i)
		C[i*N+(i-1)] = T(1);
	for (int i=0; i<N; ++i)
		C[i*N+(N-1)] = -p[i]/p[N];
	bool result = HessenbergEig(C,N,W);
	return result;
}

// vim: fenc=utf-8 noet:
