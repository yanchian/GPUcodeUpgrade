#pragma once
#include <iostream>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <cmath>

template<typename T> struct is_vector { static const bool value = false; };


template<typename T, int N>
class Vector {
private:
	T data[N];
public:
	T typedef valueType;
	static const int SIZE = N;

	__host__ __device__ __forceinline__ Vector() {
	}

	template<typename E>
	__host__ __device__ __forceinline__ Vector(const E& vec) {
		for (int i = 0; i < vec.size(); i++) {
			data[i] = vec[i];
		}
	}

	__host__ __device__ __forceinline__ Vector(const T v) {
		for (int i = 0; i < size(); i++) {
			data[i] = v;
		}
	}
	__host__ __device__ __forceinline__ Vector(const T v[N]) {
		for (int i = 0; i < size(); i++) {
			data[i] = v[i];
		}
	}
	__host__ __device__ __forceinline__ Vector(const T a, const T b) {
		data[0] = a;
		data[1] = b;
	}
	__host__ __device__ __forceinline__ Vector(const T a, const T b, const T c) {
		data[0] = a;
		data[1] = b;
		data[2] = c;
	}

	__host__ __device__ __forceinline__ Vector<T, N>& operator=(const T v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			data[i] = v;
		}
		return *this;
	}

	template<typename E>
	__host__ __device__ __forceinline__ Vector<T, N>& operator+=(const E& v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			data[i] += v[i];
		}
		return *this;
	}
	template<typename E>
	__host__ __device__ __forceinline__ Vector<T, N>& operator-=(const E& v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			data[i] -= v[i];
		}
		return *this;
	}

	__host__ __device__ __forceinline__ Vector<T, N>& operator*=(const T v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			data[i] *= v;
		}
		return *this;
	}
	__host__ __device__ __forceinline__ Vector<T, N>& operator/=(const T v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			data[i] /= v;
		}
		return *this;
	}

	__host__ __device__ __forceinline__ T& operator[](const int i) { return data[i]; }
	__host__ __device__ __forceinline__ T  operator[](const int i) const { return data[i]; }
	__host__ __device__ __forceinline__ int size() const { return N; }
};

template<typename T, int N>
struct is_vector<Vector<T, N> > { static const bool value = true; };

template<typename T, int N>
std::ostream& operator<<(std::ostream& out, const Vector<T, N>& v) {
	out << "{";
	for (int i = 0; i < N; i++) {
		out << v[i];
		if (i != N - 1) {
			out << ", ";
		}
	}
	out << "}";
	return out;
}

template<typename T, int N>
class StridedVector {
private:
	T* const data;
	const int stride;
public:
	typedef T valueType;
	static const int SIZE = N;

	__host__ __device__ __forceinline__ StridedVector(T* const data_, const int stride_) : data(data_), stride(stride_) {
	}

	template<typename E>
	__host__ __device__ __forceinline__ StridedVector<T, N>& operator=(const E& vec) {
		#pragma unroll
		for (int i = 0; i < vec.size(); i++) {
			(*this)[i] = vec[i];
		}
		return *this;
	}

	template<typename S>
	__host__ __device__ __forceinline__ StridedVector<T, N>& operator*=(const S v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			(*this)[i] *= v;
		}
		return *this;
	}

	template<typename S>
	__host__ __device__ __forceinline__ StridedVector<T, N>& operator/=(const S v) {
		#pragma unroll
		for (int i = 0; i < size(); i++) {
			*(data + i * stride) /= v;
		}
		return *this;
	}

	__host__ __device__ __forceinline__ T& operator[](const int i) { return *(data + i * stride); }
	__host__ __device__ __forceinline__ T  operator[](const int i) const { return *(data + i * stride); }
	__host__ __device__ __forceinline__ int size() const { return N; }
};

template<typename T, int N>
struct is_vector<StridedVector<T, N> > { static const bool value = true; };

// VECTOR ADDITION
template<typename E1, typename E2>
class VectorAddition {
	const E1& u;
	const E2& v;
public:
	typedef typename E1::valueType valueType;
	static const int SIZE = E1::SIZE;

	__host__ __device__ __forceinline__ VectorAddition(const E1& u_, const E2& v_) : u(u_), v(v_) { }

	__host__ __device__ __forceinline__ int size() const { return u.size(); }
	__host__ __device__ __forceinline__ typename E1::valueType operator[](const int i) const { return u[i] + v[i]; }
};

template<typename E1, typename E2>
struct is_vector<VectorAddition<E1, E2> > { static const bool value = true; };

template<typename E1, typename E2>
__host__ __device__ __forceinline__ const typename boost::enable_if<is_vector<E1>, VectorAddition<E1, E2> >::type
operator+(const E1& u, const E2& v) {
	return VectorAddition<E1, E2>(u, v);
}

// VECTOR DIFFERENCE
template<typename E1, typename E2>
class VectorDifference  {
	const E1& u;
	const E2& v;
public:
	typedef typename E1::valueType valueType;
	static const int SIZE = E1::SIZE;

	__host__ __device__ __forceinline__ VectorDifference(const E1& u_, const E2& v_) : u(u_), v(v_) { }

	__host__ __device__ __forceinline__ int size() const { return u.size(); }
	__host__ __device__ __forceinline__ typename E1::valueType operator[](const int i) const { return u[i] - v[i]; }
};

template<typename E1, typename E2>
struct is_vector<VectorDifference<E1, E2> > { static const bool value = true; };

template<typename E1, typename E2>
__host__ __device__ __forceinline__ const typename boost::enable_if<is_vector<E1>, VectorDifference<E1, E2> >::type
operator-(const E1& u, const E2& v) {
	return VectorDifference<E1, E2>(u, v);
}

// VECTOR CROSS PRODUCT
template<typename E1, typename E2>
class VectorCross {
	const E1& u;
	const E2& v;
public:
	typedef typename E1::valueType valueType;
	static const int SIZE = E1::SIZE;

	__host__ __device__ __forceinline__ VectorCross(const E1& u_, const E2& v_) : u(u_), v(v_) { }

	__host__ __device__ __forceinline__ int size() const { return u.size(); }
	__host__ __device__ __forceinline__ typename E1::valueType operator[](const int i) const {
		switch (i) {
			case 0: return u[1] * v[2] - u[2] * v[1];
			case 1: return u[2] * v[0] - u[0] * v[2];
			case 2: return u[0] * v[1] - u[1] * v[0];
		}
		return 0; // avoid error
	}
};

template<typename E1, typename E2>
struct is_vector<VectorCross<E1, E2> > { static const bool value = true; };

template<typename E1, typename E2>
__host__ __device__ __forceinline__
const typename boost::enable_if<boost::type_traits::ice_and<is_vector<E1>::value, is_vector<E2>::value >, VectorCross<E1, E2> >::type
cross(const E1& u, const E2& v) {
	return VectorCross<E1, E2>(u, v);
}

// SCALAR MULTIPLICATION
template<typename E1, typename E2>
class VectorScale {
	const E1 u;
	const E2& v;
public:
	typedef typename E2::valueType valueType;
	static const int SIZE = E2::SIZE;

	__host__ __device__ __forceinline__ VectorScale(const E1 u_, const E2& v_) : u(u_), v(v_) { }

	__host__ __device__ __forceinline__ int size() const { return v.size(); }
	__host__ __device__ __forceinline__ typename E2::valueType operator[](const int i) const { return u * v[i]; }
};

template<typename E1, typename E2>
struct is_vector<VectorScale<E1, E2> > { static const bool value = true; };

template<typename E1, typename E2>
__host__ __device__ __forceinline__ const typename boost::enable_if<is_vector<E2>, VectorScale<E1, E2> >::type
operator*(const E1 u, const E2& v) {
	return VectorScale<E1, E2>(u, v);
}
template<typename E1, typename E2>
__host__ __device__ __forceinline__ const typename boost::enable_if<is_vector<E2>, VectorScale<E1, E2> >::type
operator*(const E2& v, const E1 u) {
	return VectorScale<E1, E2>(u, v);
}
template<typename E1, typename E2>
__host__ __device__ __forceinline__ const typename boost::enable_if<is_vector<E2>, VectorScale<E1, E2> >::type
operator/(const E2& v, const E1 u) {
	return VectorScale<E1, E2>(1.0/u, v);
}

template<typename E1, typename E2>
__host__ __device__ __forceinline__ typename E1::valueType dot(const E1& v, const E2& u) {
	typename E1::valueType sum = 0;
	#pragma unroll
	for (int i = 0; i < v.size(); i++) {
		sum += v[i] * u[i];
	}
	return sum;
}
template<typename E>
__host__ __device__ __forceinline__ typename E::valueType abs(const E& v) {
	return sqrt(dot(v, v));
}
template<typename E>
__host__ __device__ __forceinline__ typename E::valueType norm(const E& v) {
	return sqrt(dot(v, v));
}
template<typename E>
__host__ __device__ __forceinline__ typename E::valueType norm2(const E& v) {
	return dot(v, v);
}
template<typename E>
__host__ __device__ __forceinline__ void normalise(E& v) {
	v /= abs(v);
}
template<typename E>
__host__ __device__ __forceinline__ typename E::valueType vectorAngle(const E& v, const E& u) {
	return std::acos(dot(v, u) / (abs(v) * abs(u)));
}

template<typename E1, typename E2>
__host__ __device__ __forceinline__ typename E1::valueType distance(const E1& u, const E2& v) {
	return norm(u - v);
}

template<typename E1, typename E2>
__host__ __device__ __forceinline__
const typename boost::enable_if<boost::type_traits::ice_and<is_vector<E1>::value, is_vector<E2>::value >, bool>::type
operator==(const E1& u, const E2& v) {
	for (int i = 0; i < E1::SIZE; i++) {
		if (u[i] != v[i]) return false;
	}
	return true;
}

template<typename E1, typename E2>
__host__ __device__ __forceinline__
const typename boost::enable_if<boost::type_traits::ice_and<is_vector<E1>::value, is_vector<E2>::value >, bool>::type
operator!=(const E1& u, const E2& v) {
	return !(u == v);
}

