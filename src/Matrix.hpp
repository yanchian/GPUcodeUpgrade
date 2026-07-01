/*
  Copyright © Cambridge Numerical Solutions Ltd 2013
*/
#pragma once
#include <iostream>
#include <cmath>

template<typename T, int N, int M>
class FixedMatrix {
private:
	T data[N][M];
public:
	T typedef valueType;

	valueType operator()(const int i, const int j) const { return data[i][j]; }
	valueType& operator()(const int i, const int j) { return data[i][j]; }

	template<typename S>
	FixedMatrix<T, N, M>& operator=(const S c) {
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				(*this)(i, j) = c;
			}
		}
		return *this;
	}
};

