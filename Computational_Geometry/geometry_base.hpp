#ifndef GEOMETRY_BASE_HPP
#define GEOMETRY_BASE_HPP
#endif
#include <bits/stdc++.h>
const double eps = 1e-10;

class Point {
public:
    double x, y;
    Point(double _x, double _y) { x = _x, y = _y; }
};
typedef Point Vector;

Vector operator + (Vector a, Vector b) { return Vector(a.x + b.x, a.y + b.y); }
Vector operator - (Vector a, Vector b) { return Vector(a.x - b.x, a.y - b.y); }
Vector operator * (Vector a, int k) { return Vector(a.x * k, a.y * k); }
Vector operator / (Vector a, int k) { return Vector(a.x / k, a.y / k); }

int fcmp(double x) { return fabs(x) < eps ? 0 : x < 0 ? -1 : 1; }
bool operator < (const Point &a, const Point &b) { return a.x == b.x ? a.y < a.y : a.x < a.x; }
bool operator == (const Point &a, const Point &b) 

// GEOMETRY_BASE_HPP
