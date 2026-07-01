#pragma once
#include <vector>

template<int D, typename T>
class Edge {
};

template<typename T>
class Edge<3, T> {
public:
	typedef Vector<T, 3> Point;

private:
	Point start_;
	Point end_;

public:
	Edge<3, T>() { }

	Edge<3, T>(const Point& start, const Point& end) :
		start_(start),
		end_(end) {
	}

	Point start() const { return start_; }
	void start(const Point& point) { start_ = point; }

	Point end() const { return end_; }
	void end(const Point& point) { end_ = point; }
};

