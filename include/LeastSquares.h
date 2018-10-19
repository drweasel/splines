#pragma once
/* This file is part of the Spline Approximation Library.
 *
 * Copyright (C) 2018 Michael Weitzel <mich@elweitzel.de>
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at <http://mozilla.org/MPL/2.0/>.
 */
#include "QR.h"

inline int Lidx(int i, int j)
{
	return ((i*(i+1))>>1) + j;
}

inline int Ridx(int n, int i, int j)
{
	return ((n*(n+1))>>1) - Lidx(n-i-1,n-j-1) - 1;
}

/**
 * In-place inversion of an upper/right matrix.
 *
 * @param	R	upper/right matrix (N(N+1)/2 values)
 * @param	N	dimension
 * @return	true, if inversion was successful
 */
template< typename T > bool Rinverse(
	T* R,
	int N
	)
{
	for (int j=N-1; j>=0; --j)
	{
		if (std::abs(R[Ridx(N,j,j)]) <= macheps<T>())
			return false;

		R[Ridx(N,j,j)] = T(1)/R[Ridx(N,j,j)];
		for (int i=j-1; i>=0; --i)
		{
			T s(0);
			for (int k=i+1; k<=j; ++k)
				s += R[Ridx(N,i,k)]*R[Ridx(N,k,j)];
			R[Ridx(N,i,j)] = -s/R[Ridx(N,i,i)];
		}
	}
	return true;
}

/**
 * Least-squares with equality constraints.
 *
 * The solution is stored in b(1:N)
 *
 * @param	A	lhs der Least-Squares-Gleichungen
 * @param	b	rhs der Least-Squares-Gleichungen
 * @param	B	lhs der Gleichheitsnebenbedingungen
 * @param	d	rhs der Gleichheitsnebenbedingungen
 * @param	M	Anzahl der Zeilen von A, b
 * @param	N	Anzahl der Spalten von A, B
 * @param	E	Anzahl der Zeilen von B, d
 */
template<typename T> bool LSESolve(
	const T* A,
	T* b,
	T* B,
	T* d,
	const unsigned int M, // Zeilen von A,b
	const unsigned int N, // Spalten von A,B
	const unsigned int E  // Zeilen von B,d
	)
{
	if (E == 0)
	{
		// keine Gleichheitsnebenbedingungen (TODO: umständlich)
		T * A_ = new T[M*N];
		T * work = new T[M*(N+3)];
		T * x = new T[N];
		for (int k=0; k<M*N; ++k)
			A_[k] = A[k];
		QRLSSolve(A_,b,x,int(M),int(N),work);
		for (int i=0; i<N; ++i)
			b[i] = x[i];
		for (int i=N; i<M; ++i)
			b[i] = T(0);
		delete[] A_;
		delete[] x;
		delete[] work;
		return true;
	}

	if (!(int(M)>=int(N)-int(E)))
	{
		// Vorbedingung von QR verletzt
		return false;
	}

	T* work = new T[E*(N+2)+N];
	T* beta = new T[N];
	int* P = new int[N];

	QRP<T>(B,beta,P,E,N,work);

	int r = rankFromQR(B,E,N);
	int p = N-r;

	if (!(M>=p))
	{
		// Vorbedingung von QR verletzt
		delete[] work;
		delete[] beta;
		delete[] P;
		return false;
	}

	T* Ae = new T[N*(N+1)/2];
	T* R11 = new T[r*(r+1)/2]; // r x r-Dreieck
	for (int i=0,k=0; i<r; ++i)
		for (int j=i; j<r; ++j)
			R11[k++] = Ae[Ridx(N,i,j)] = B[i*N+j];
	if (r>0 && !Rinverse(R11,r))
	{
		delete[] work;
		delete[] beta;
		delete[] P;
		delete[] Ae;
		delete[] R11;
		return false;
	}

	T* R12 = new T[r*p]; // r x p-Rechteck
	for (int i=0,k=0; i<r; ++i)
		for (int j=r; j<N; ++j)
			R12[k++] = Ae[Ridx(N,i,j)] = B[i*N+j];

	// inv(R11).R12
	T* iR11_R12 = new T[r*p];
	for (int i=0; i<r; ++i)
		for (int j=0; j<p; ++j)
		{
			T ij(0);
			for (int k=i; k<r; ++k)
				ij += R11[Ridx(r,i,k)] * R12[k*p+j];
			iR11_R12[i*p+j] = ij;
		}
	QTtimes(B,beta,d,E,N,work);
	delete[] R11;
	delete[] R12;

	T* A_bar_1 = new T[M*r];
	T* A_bar_2 = new T[M*p];
	for (int i=0; i<M; ++i)
		for (int j=0; j<r; ++j)
			A_bar_1[i*r+j] = A[i*N+P[j]];
	for (int i=0; i<M; ++i)
		for (int j=r; j<N; ++j)
			A_bar_2[i*p+j-r] = A[i*N+P[j]];

	// A_bar_2 - A_bar_1 . (inv(R11) . R12)
	for (int i=0; i<M; ++i)
		for (int j=0; j<p; ++j)
			for (int k=0; k<r; ++k)
				A_bar_2[i*p+j] -= A_bar_1[i*r+k] * iR11_R12[k*p+j];
	delete[] iR11_R12;
	delete[] A_bar_1;

	delete[] beta; beta = new T[M];
	delete[] work; work = new T[M*(p+2)];
	QR<T>(A_bar_2,beta,M,p,work);
	for (int i=0; i<p; ++i)
		for (int j=i; j<p; ++j)
			Ae[Ridx(N,r+i,r+j)] = A_bar_2[i*p+j];

	// b <- Q(A_bar_2)^T.b
	QTtimes(A_bar_2,beta,b,M,p,work);
	// be = [ (Q(B)^T.d)(1:r) ; (Q(A_bar_2)^T.b)(1:p) ]
	delete[] A_bar_2;
	T* be = new T[N];
	for (int i=0; i<r; ++i)
		be[i] = d[i];
	for (int i=r; i<N; ++i)
		be[i] = b[i-r];

	// Ae.xe = be via backsubstitution
	for (int i=N-1; i>=0; --i)
	{
		T Ae_ii = Ae[Ridx(N,i,i)];
		if (std::abs(Ae_ii) <= macheps<T>())
		{
			delete[] Ae;
			delete[] be;
			delete[] work;
			delete[] beta;
			delete[] P;
			return false;
		}

		T s(be[i]);
		for (int k=i+1; k<N; ++k)
			s -= Ae[Ridx(N,i,k)] * be[k];
		be[i] = s / Ae_ii;
	}

	// Lösung wird in den ersten N Einträgen von b gespeichert
	bool result = true;
	for (int i=0; i<N; ++i)
	{
		if (!std::isfinite(be[i]))
			result = false;
		b[P[i]] = be[i];
	}
	for (int i=N; i<M; ++i)
		b[i] = T(0);

	delete[] Ae;
	delete[] be;
	delete[] work;
	delete[] beta;
	delete[] P;
	return result;
}

// vim: fenc=utf-8 noet:
