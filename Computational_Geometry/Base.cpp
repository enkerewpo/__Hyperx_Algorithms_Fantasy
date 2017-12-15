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
bool operator == (const Point &a, const Point &b) { return fcmp(a.x - b.x) == 0 && fcmp(a.y - b.y) == 0; }

double Dot(Vector a, Vector b) { return a.x * b.x + a.y * b.y; }
double Cross(Vector a, Vector b) { return a.x * b.y - a.y * b.x; }
double Area2(Point A, Point B, Point C) { return Cross(B - A, C - A); }
double Length(Vector a) { return sqrt(Dot(a, a)); }
double Angle(Vector a, Vector b) { return acos(Dot(a, b)) / (Length(a) * Length(b)); }
Vector Rotate(Vector a, double rad) {
    return Vector(a.x * cos(rad) - a.y * sin(rad), a.x * sin(rad) + a.y * cos(rad));
}
Vector Normal(Vector a) {
    double L = Length(a);
    return Vector(-a.y / L, a.x / L);
}

bool SegmentCross(Point A, Point B, Point C, Point D) {
    int a = fcmp(Cross(Vector(B - A), Vector(C - A)) * Cross(Vector(D - A), Vector(B - A))) >= 0;
    int b = fcmp(Cross(Vector(D - C), Vector(B - C)) * Cross(Vector(A - C), Vector(D - C))) >= 0;
    return a && b;
}
Point GetLineIntersection(Point P, Vector v, Point Q, Vector w) {
    Vector u = P - Q;
    int t = Cross(w, u) / Cross(v, w);
    return P + v * t;
}
Point GetSegmentIntersection(Point A, Point B, Point C, Point D) {
    Vector a = B - A;
    double s1 = fabs(Cross(a, C - A));
    double s2 = fabs(Cross(a, D - A));
    return Point((s1 * D.x + s2 * C.x) / (s1 + s2), (s1 * D.y + s2 * C.y) / (s1 + s2));
}
double DistanceToLine(Point P, Point A, Point B) {
    Vector a = B - A, b = P - A;
    return fabs(Cross(a, b)) / Length(a);
}
double DistanceToSegment(Point P, Point A, Point B) {
    if(A == B) return Length(P - A);
    Vector a = B - A, b = P - A, c = P - B;
    if(fcmp(Dot(a, b)) < 0) return Length(b);
    else if(fcmp(Dot(a, c)) > 0) return Length(c);
    else return fabs(Cross(a, b)) / Length(a);
}
int main() {
    return 0;
}
