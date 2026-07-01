/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#pragma once
#include <exception>
#include "Matrix.hpp"
#include "Vector.hpp"

template<typename T, int N, int M, int O>
FixedMatrix<T, N, O>
operator*(const FixedMatrix<T, N, M>& a, const FixedMatrix<T, M, O>& b) {
	FixedMatrix<T, N, O> c;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < O; j++) {
			c(i, j) = 0.0;
			for (int k = 0; k < M; k++) {
				c(i, j) += a(i, k) * b(k, j);
			}
		}
	}
	return c;
}

template<typename T, int N, int M, typename E>
typename boost::enable_if<
	boost::type_traits::ice_and<
		is_vector<E>::value,
		E::SIZE == M
	>,
	Vector<T, N>
>::type
operator*(const FixedMatrix<T, N, M> A, const E& x) {
	Vector<T, N> c;
	for (int i = 0; i < N; i++) {
		c[i] = 0.0;
		for (int k = 0; k < M; k++) {
			c[i] += A(i, k) * x[k];
		}
	}
	return c;
}

// super-special overload for using homogeneous coordinates with 2D vectors
template<typename T, typename E>
typename boost::enable_if<
	boost::type_traits::ice_and<
		is_vector<E>::value,
		E::SIZE == 2
	>,
	Vector<T, 2>
>::type
operator*(const FixedMatrix<T, 3, 3> A, const E& x) {
	Vector<T, 2> c;
	for (int i = 0; i < 2; i++) {
		c[i] = 0.0;
		for (int k = 0; k < 2; k++) {
			c[i] += A(i, k) * x[k];
		}
		c[i] += A(i, 2);
	}
	return c;
}
