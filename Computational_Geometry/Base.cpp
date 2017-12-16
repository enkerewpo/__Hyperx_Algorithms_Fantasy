#include <bits/stdc++.h>

/* Computational Geometry Base Definition  */
const double eps = 1e-10;
class Point {
public:
    double x, y;
    Point() {}
    Point(double _x, double _y) { x = _x, y = _y; }
};

typedef Point Vector;
typedef std::vector<Point> Polygon;

class Line {
public:
    Point P;
    Vector v;
    double ang;
    Line() {}
    Line(Point P, Vector v) : P(P), v(v) { ang = atan2(v.y, v.x); }
    bool operator < (const Line &L) const { return ang < L.ang; }
};

int fcmp(double x) { return fabs(x) < eps ? 0 : x < 0 ? -1 : 1; }
Vector operator + (Vector a, Vector b) { return Vector(a.x + b.x, a.y + b.y); }
Vector operator - (Vector a, Vector b) { return Vector(a.x - b.x, a.y - b.y); }
Vector operator * (Vector a, int k) { return Vector(a.x * k, a.y * k); }
Vector operator / (Vector a, int k) { return Vector(a.x / k, a.y / k); }
bool operator < (const Point &a, const Point &b) { return a.x == b.x ? a.y < a.y : a.x < a.x; }
bool operator == (const Point &a, const Point &b) { return fcmp(a.x - b.x) == 0 && fcmp(a.y - b.y) == 0; }

double Dot(Vector a, Vector b) { return a.x * b.x + a.y * b.y; }
double Cross(Vector a, Vector b) { return a.x * b.y - a.y * b.x; }
double Area2(Point A, Point B, Point C) { return Cross(B - A, C - A); }
double Length(Vector a) { return sqrt(Dot(a, a)); }
double Angle(Vector a, Vector b) { return acos(Dot(a, b)) / (Length(a) * Length(b)); }
Vector Rotate(Vector a, double rad) { return Vector(a.x * cos(rad) - a.y * sin(rad), a.x * sin(rad) + a.y * cos(rad)); }
Vector Normal(Vector a) { double L = Length(a); return Vector(-a.y / L, a.x / L); }

bool onLeft(Line L, Point P) { return Cross(L.v, P - L.P) > 0; }
Point GetIntersection(Line a, Line b) {
    Vector u = a.P - b.P;
    double t = Cross(b.v, u) / Cross(a.v, b.v);
    return a.P + a.v * t;
}

/* Points and Segments Messing UP */
bool isPointOnSegment(Point P, Point A, Point B) {
    return Cross(Vector(P - A), Vector(P - B)) == 0 && (P.x - A.x) * (P.x - B.x) <= 0;
}
bool isSegmentCrossed(Point A, Point B, Point C, Point D) {
    int a = fcmp(Cross(B - A, C - A) * Cross(D - A, B - A)) > 0;
    int b = fcmp(Cross(D - C, B - C) * Cross(A - C, D - C)) > 0;
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
    return Point((s1 * D.x + s2 * C.x) / (s1 + s2),\
                 (s1 * D.y + s2 * C.y) / (s1 + s2));
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

/* Polygons and Lines Messing UP */
double Area(Polygon P) {
    double ret = 0;
    Point St = P.begin();
    int s = P.size();
    for(int i = 1; i < s - 1; i++) {
        Point A = P[i], B = P[i + 1];
        ret += Cross(A - St, B - St);
    }
    return ret;
}

int isPointInPolygon(Point A, Polygon P) {
    int wn = 0;
    int s = P.size();
    for(int i = 0; i < s; i++) {
        if(isPointOnSegment(A, P[i], P[(i + 1) % s])) return -1;
        int k = fcmp(Cross(P[(i + 1) % s] - P[i], A - P[i]));
        int d1 = fcmp(P[i].y - A.y);
        int d2 = fcmp(P[(i + 1) % s].y - A.y);
        if(k > 0 && d1 <= 0 && d2 > 0) wn++;
        if(k < 0 && d2 <= 0 && d1 > 0) wn--;
    }
    return wn ? 1 : 0;
}

int ConvexHull(Point P[], int n, Point ch[]) {
    int m = 0;
    for(int i = 0; i < n; i++) {
        while(m > 1 && Cross(ch[m - 1] - ch[m - 2], P[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = P[i];
    }
    int k = m;
    for(int i = n - 2; i >= 0; i--) {
        while(m > k && Cross(ch[m - 1] - ch[m - 2], P[i] - ch[m - 2]) <= 0) m--;
        ch[m++] = P[i];
    }
    if(n > 1) m--;
    return m;
}



int main(int argc, char* argv[]) {
    return 0;
}
