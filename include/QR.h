#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include <cmath>
#include <algorithm>
#include <type_traits>
#include <limits>

/**
 * constexpr computation of the machine epsilon for floating-point type T.
 * The machine epsilon is the smallest number eps which can be discriminated
 * from 0 (in particular: the smallest number eps with 1.0 + eps != 1.0).
 */
template< typename T > constexpr T macheps(const T& x = T(0))
{
#if 0
	static_assert(std::is_floating_point<T>::value,"T is not a fp type!");

	T o_eps(1), eps = T(0.5);
	while (x + eps != x)
	{
		o_eps = eps;
		eps *= T(0.5);
	}
	return o_eps;
#else
	return std::nextafter(x,std::numeric_limits<T>::max());
#endif
}

/**
 * Computes the Householder vector v for a given column vector x
 * and the associated scaling factor beta = 2/(v^T.v).
 *
 * @param	x		[in] input vector (dim M)
 * @param	v		[out] computed Householder vector (dim M)
 * @param	beta	[out] scaling factor (beta) for Householder matrix
 * @param	M		[in] dimension of x and v
 */
template< typename T > void householder(
	const T * x,
	T * v,
	T * beta,
	int M
	)
{
	if (M == 1)
	{
		*v = T(1);
		*beta = T(0);
		return;
	}
	T sigma = x[1]*x[1];
	v[0] = T(1);
	v[1] = x[1];
	for (int k=2; k<M; ++k)
	{
		sigma += x[k]*x[k];
		v[k] = x[k];
	}
	if (sigma == T(0))
		*beta = T(0);
	else
	{
		T mu(std::sqrt(x[0]*x[0] + sigma));
		if (x[0] <= T(0))
			v[0] = x[0] - mu;
		else
			v[0] = -sigma / (x[0] + mu);

		*beta = T(2)*(v[0]*v[0]) / (sigma + (v[0]*v[0]));
		T v0(v[0]);
		for (int k=0; k<M; ++k)
			v[k] /= v0;
	}
}

/**
 * Householder QR Factorization.
 *
 * Factorizes MxN matrix (with M>=N) into orthonormal MxM matrix Q and the
 * upper/right MxN triangular matrix R (with up to N*(N+1)/2 non-zeros).
 *
 * After factorization A contains:
 *  1. the triangular matrix R on and above the diagonal
 *  2. the Householder vectors v below the diagonal
 *
 * @param	A		[in] dim MxN matrix, [out] QR factorization of A
 * @param	beta	[out] dim M vector; scaling factors of the outer products
 * 		of the Householder vectors
 * @param	M		[in] number of rows of A (M>=N)
 * @param	N		[in] number of columns of A (N<=M)
 * @param	work	[in,out] room for M*(N+2) temp. values
 */
template< typename T > void QR(
	T * A,
	T * beta,
	int M,
	int N,
	T * work
	)
{
	// Platz für M*(N+2) temp. Werte
	T * v = work;
	T * x = work + M;
	T * C = work + 2*M;

	for (int j=0; j<N; ++j)
	{
		for (int i=j; i<M; ++i)
			x[i-j] = A[i*N+j];
		householder(x,v,&beta[j],M-j);
		// A(j:m,j:n) = (eye(m-j+1)-beta(j)*v*v') * A(j:m,j:n);
		for (int ii=j; ii<M; ++ii)
		{
			for (int jj=j; jj<N; ++jj)
			{
				T & c_ij = C[ii*N+jj];
				c_ij = T(0);
				for (int kk=j; kk<M; ++kk)
					c_ij += (T(ii==kk?1:0) - beta[j]*v[ii-j]*v[kk-j]) * A[kk*N+jj];
			}
		}
		for (int ii=j; ii<M; ++ii)
			for (int jj=j; jj<N; ++jj)
				A[ii*N+jj] = C[ii*N+jj];

		// A(j+1:m,j) = v(2:m-j+1)
		for (int i=j+1; i<M; ++i)
			A[i*N+j] = v[i-j];
	}
}

/**
 * Computes the orthonormal matrix Q for a given QR factorization.
 *
 * @param	A		[in] QR factorization of a matrix A (dim MxN)
 * @param	beta	[in] scaling factors of the outer products of
 * 		Householder vectors (dim N)
 * @param	Q		[out] orthonormal matrix Q (dim MxM)
 * @param	M		[in] number of rows (dim of Q)
 * @param	N		[in] number of columns
 */
template< typename T > void QfromQR(
	const T* A,
	const T* beta,
	T* Q,
	int M,
	int N
	)
{
	// Q = eye(m);
	for (int i=0; i<M; ++i)
		for (int j=0; j<M; ++j)
			Q[i*M+j] = (i==j)*T(1); // indeffizient!

	// Platz für M*(M+1) temp. Werte
	T* v = new T[M];
	T* C = new T[M*M];

	v[0] = T(1);
	for (int j=N-1; j>=0; --j)
	{
		// v = [ 1 ; A(j+1:m,j) ];
		for (int i=j+1; i<M; ++i)
			v[i-j] = A[i*N+j];

		// Q(j:m,j:m) = (eye(m-j+1)-beta(j)*v*v')*Q(j:m,j:m);
		for (int ii=j; ii<M; ++ii)
		{
			for (int jj=j; jj<M; ++jj)
			{
				T & c_ij = C[ii*M+jj];
				c_ij = T(0);
				for (int kk=j; kk<M; ++kk)
					c_ij += (T(ii==kk?1:0) - beta[j]*v[ii-j]*v[kk-j]) * Q[kk*M+jj];
			}
		}
		for (int ii=j; ii<M; ++ii)
			for (int jj=j; jj<M; ++jj)
				Q[ii*M+jj] = C[ii*M+jj];
	}

	delete[] v;
	delete[] C;
}

/**
 * R ist eine MxN rechte obere Dreiecksmatrix. Weil M>=N, hat R nicht mehr als
 * N*(N+1)/2 non-zeros (=> NxN rechte obere Dreiecksmatrix).
 *
 * @param	A	[in] QR-Zerlegung einer Matrix A
 * @param	R	[out] Speicherplatz für R (N*(N+1)/2 Werte)
 * @param	N	[in] Anzahl der Spalten von A, R
 */
template< typename T > void RfromQR(
	const T* A,
	T* R,
	int N
	)
{
	for (int i=0,k=0; i<N; ++i)
		for (int j=i; j<N; ++j)
			R[k++] = A[i*N+j];
}

/** Allgemeiner Fall, d.h. auch M<N mit MxN-Matrix R */
template< typename T > void RfromQR(
	const T* A,
	T* R,
	int M,
	int N
	)
{
	for (int i=0; i<M; ++i)
	{
		for (int j=0; j<i; ++j)
			R[i*N+j] = T(0);
		for (int j=i; j<N; ++j)
			R[i*N+j] = A[i*N+j];
	}
}

/**
 * Berechnet die orthonormale Matrix Q aus der QR(P)-Zerlegung von A,
 * A.P = Q.R (QR-Zerlegung mit Spaltenpivotisierung).
 *
 * @param	A		[in,out] QR-Zerlegung einer Matrix A (dim MxN)
 * @param	beta	[in,out] Skalierungsfaktoren des Outer Products der
 *		Householder-Vektoren (dim N)
 * @param	P		[out]
 * @param	M		[in] Anzahl der Zeilen
 * @param	N		[in] Anzahl der Spalten
 * @param	work	Platz für M*(N+2)+N temp. Werte
 */
template< typename T > void QRP(
	T * A,
	T * beta,
	int * P,
	int M,
	int N,
	T * work
	)
{
	// Platz für M*(N+2) temp. Werte
	T * v = work; // M
	T * x = work + M; // M
	T * C = work + 2*M; // M*N
	T * c = work + M*(N+2); // N

	// P = 1:N, [tau,k] = max(c)
	int k = -1;
	T tau = T(-1);
	for (int j=0; j<N; ++j)
	{
		P[j] = j;
		c[j] = A[0*N+j]*A[0*N+j];
		for (int i=1; i<M; ++i)
			c[j] += A[i*N+j]*A[i*N+j];

		if (c[j] > tau)
		{
			k = j;
			tau = c[j];
		}
	}

	// j:=r
	for (int j=0; j<N && tau>macheps<T>(); ++j)
	{
		std::swap(P[j],P[k]);
		std::swap(c[j],c[k]);
		for (int i=0; i<M; ++i)
			std::swap(A[i*N+j],A[i*N+k]);

		for (int i=j; i<M; ++i)
			x[i-j] = A[i*N+j];
		householder(x,v,&beta[j],M-j);
		// A(j:m,j:n) = (eye(m-j+1)-beta(j)*v*v') * A(j:m,j:n);
		for (int ii=j; ii<M; ++ii)
		{
			for (int jj=j; jj<N; ++jj)
			{
				T & c_ij = C[ii*N+jj];
				c_ij = T(0);
				for (int kk=j; kk<M; ++kk)
					c_ij += (T(ii==kk?1:0) - beta[j]*v[ii-j]*v[kk-j]) * A[kk*N+jj];
			}
		}
		for (int ii=j; ii<M; ++ii)
			for (int jj=j; jj<N; ++jj)
				A[ii*N+jj] = C[ii*N+jj];

		// A(j+1:m,j) = v(2:m-j+1)
		for (int i=j+1; i<M; ++i)
			A[i*N+j] = v[i-j];

		if (j<M)
		{
			for (int i=j+1; i<N; ++i)
				c[i] -= A[j*N+i]*A[j*N+i];
		}
		if (j+1<N)
		{
			// tau = max{ c(j+1),...,c(n) }
			// finde kleinstes k mit j+1<=k<=n und c(j)=tau
			//[tau,k] = max(c(j+1:n)); k = k+j;
			k = j+1;
			tau = c[j+1];
			for (int i=j+2; i<N; ++i)
			{
				if (c[i] > tau)
				{
					tau = c[i];
					k = i;
				}
			}
		}
		//else tau = T(0); // hinfällig
	}
}

/**
 * Rangbestimmung auf Basis der QRP-Zerlegung von A.
 *
 * @param	A	QR(Matrix)
 * @param	M	Anzahl der Zeilen von A
 * @param	N	Anzahl der Spalten von A
 * @return	Rang der Matrix A
 */
template< typename T > int rankFromQR(
	const T* A,
	int M,
	int N
	)
{
	const T tol((M>N?M:N) * T(10) * macheps<T>());
	for (int rank=M<N?M:N; rank>0; --rank)
		for (int j=rank-1; j<N; ++j)
			if (std::abs(A[(rank-1)*N+j]) > tol)
				return rank;
	return 0;
}

/**
 * Berechnet Q^T.b auf Basis der QR-Zerlegung von A.
 *
 * @param	A		[in] QR(Modellmatrix) (dim MxN, M Samples, N Unbekannte)
 * @param	beta	[in] Skalierungsfaktoren für das Outer-Product der
 *		Householder-Vektoren (dim M, M Samples)
 * @param	b		[in] rechte Seite, [out] Q'.b
 * @param	x		[out] Least-Squares-Lösung (dim N, N Unbekannte)
 * @param	M		[in] Anzahl der Messungen (Zeilen der Modellmatrix A)
 * @param	N		[in] Anzahl der Unbekannten (Spalten von Modellmatrix A)
 * @param	work	[in,out] Platz für 2*M temporäre Werte
 */
template< typename T > void QTtimes(
	const T* A,
	const T* beta,
	T* b,
	int M, int N,
	T* work
	)
{
	T* v = work;
	T* b_ = work + M;

	// berechne Q^T.b
	for (int j=0; j<N; ++j)
	{
		// v(j)=1; v(j+1:m) = A(j+1:m,j)
		v[j] = T(1);
		for (int i=j+1; i<M; ++i)
			v[i] = A[i*N+j];
		// b(j:m) = (I_{m-j+1} - beta_j.v.v^T).b(j:m)
		for (int ii=j; ii<M; ++ii)
		{
			b_[ii] = T(0);
			for (int jj=j; jj<M; ++jj)
				b_[ii] += (T(ii==jj?1:0) - beta[j] * v[ii]*v[jj]) * b[jj];
		}
		// von b_ nach b kopieren
		for (int i=j; i<M; ++i)
			b[i] = b_[i];
	}
}

/**
 * Berechnet Q.b auf Basis der QR-Zerlegung von A (sonst wie oben).
 */
template< typename T > void Qtimes(
	const T* A,
	const T* beta,
	T* b,
	int M, int N,
	T* work
	)
{
	T* v = work;
	T* b_ = work + M;

	// berechne Q.b
	for (int j=N-1; j>=0; --j)
	{
		// v(j)=1; v(j+1:m) = A(j+1:m,j)
		v[j] = T(1);
		for (int i=j+1; i<M; ++i)
			v[i] = A[i*N+j];
		// b(j:m) = (I_{m-j+1} - beta_j.v.v^T).b(j:m)
		for (int ii=j; ii<M; ++ii)
		{
			b_[ii] = T(0);
			for (int jj=j; jj<M; ++jj)
				b_[ii] += (T(ii==jj?1:0) - beta[j] * v[ii]*v[jj]) * b[jj];
		}
		// von b_ nach b kopieren
		for (int i=j; i<M; ++i)
			b[i] = b_[i];
	}
}

/**
 * Least-Squares-Lösung auf Basis der QR-Zerlegung.
 *
 * @param	A		[in] QR(Modellmatrix) (dim MxN, M Samples, N Unbekannte)
 * @param	beta	[in] Skalierungsfaktoren für das Outer-Product der
 *		Householder-Vektoren (dim M, M Samples)
 * @param	b		[in] rechte Seite, [out] Q'.b
 * @param	x		[out] Least-Squares-Lösung (dim N, N Unbekannte)
 * @param	M		[in] Anzahl der Messungen (Zeilen der Modellmatrix A)
 * @param	N		[in] Anzahl der Unbekannten (Spalten von Modellmatrix A)
 * @param	work	[in,out] Platz für 2*M temporäre Werte
 */
template< typename T > void QRLSSolve(
	const T* A,
	const T* beta,
	T* b,
	T* x,
	int M,
	int N,
	T* work
	)
{
	// berechne Q^T.b
	QTtimes(A,beta,b,M,N,work);
	// Backward; löse R(1:n,1:n).x = b(1:n)
	for (int i=N-1; i>=0; --i)
	{
		x[i] = b[i];
		for (int j=i+1; j<N; ++j)
			x[i] -= A[i*N+j]*x[j];
		x[i] /= A[i*N+i];
	}
}

/**
 * Least-Squares-Löser (Householder-QR, Frontend).
 *
 * Seiteneffekte:
 *   (1) in A wird die QR-Zerlegung von A gespeichert
 *   (2) in b wird Q'.b gespeichert
 *
 * @param	A	MxN-Modellmatrix
 * @param	b	M-Vektor (rechte Seite)
 * @param	x	N-Vektor (Least-Squares-Lösung)
 * @param	M	Anzahl der Zeilen ("Messungen")
 * @param	N	Anzahl der Spalten ("Parameter")
 * @param	work	Platz für M*(N+3) temp. Werte
 */
template< typename T > void QRLSSolve(
	T *A,
	T *b,
	T *x,
	int M,
	int N,
	T *work
	)
{
	// Lösung mit Rückwärtseinsetzen (suboptimal, siehe mqrsolve.m)
	T* beta = work;
	work += M;
	// QR-Zerlegung von A (wird zu QR)
	QR(A,beta,M,N,work);
	QRLSSolve(A,beta,b,x,M,N,work);
}

// vim: fenc=utf-8 noet:
