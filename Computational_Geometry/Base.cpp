 /**
 *   base.cpp   Copyright(C) 2017 Kvar_ispw17 All rights reserved.
 */

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

class Circle {
public:
    Point c;
    double r;
    Circle(Point c, double r) : c(c), r(r) {}
    Point point(double a) { return Point(c.x * cos(a) * r, c.y * sin(a) * r); }
};

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
double angle(Vector a) { return atan2(a.y, a.x); }
double Angle(Vector a, Vector b) { return acos(Dot(a, b)) / (Length(a) * Length(b)); }
Vector Rotate(Vector a, double rad) { return Vector(a.x * cos(rad) - a.y * sin(rad), a.x * sin(rad) + a.y * cos(rad)); }
Vector Normal(Vector a) { double L = Length(a); return Vector(-a.y / L, a.x / L); }
double Dist(Point A, Point B) { return Length(B - A); }

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
    Point St = *P.begin();
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

int HalfplainIntersection(Line L[], int n, Point Poly[]) {
    std::sort(L, L + n);
    int hd, tl;
    Point *P = new Point[n];
    Line *q = new Line[n];
    q[hd = tl = 0] = L[0];
    for(int i = 1; i < n; i++) {
        while(hd < tl && !onLeft(L[i], P[tl - 1])) tl--;
        while(hd < tl && !onLeft(L[i], P[hd])) hd++;
        q[++tl] = L[i];
        if(fabs(Cross(q[tl].v, q[tl - 1].v) < eps)) {
            tl--;
            if(onLeft(q[tl], L[i].P)) q[tl] = L[i];
        }
        if(hd < tl) P[tl - 1] = GetIntersection(q[tl - 1], q[tl]);
    }
    while(hd < tl && !onLeft(q[hd], P[tl - 1])) tl--;
    if(tl - hd <= 1) return 0;
    P[tl] = GetIntersection(q[tl], q[hd]);
    int m = 0;
    for(int i = hd; i <= tl; i++) Poly[m++] = P[i];
    return m;
}

double rotating_calipers(Point P[], int n) {
    int x = 1;
    double ans = 0;
    P[n] = P[0];
    for(int i = 0; i < n; i++) {
        while(Cross(P[i + 1] - P[i], P[x + 1] - P[i]) > Cross(P[i + 1] - P[i], P[x] - P[i])) x = (x + 1) % n;
        ans = std::max(ans, Dist(P[x], P[i]));
        ans = std::max(ans, Dist(P[x + 1], P[i + 1]));
    }
    return ans;
}

using namespace std;
/* Circle and &^%^%#@ messing up */
int getLineCircleIntersection(Line L, Circle C, double &t1, double &t2. vector<Point> &sol) {
    double a = L.v.x, b = L.p.x - C.c.x, c = L.v.y, d = L.p.y - C.c.y;
    double e = a * a + c * c, f = 2 * (a * b + c * d), g = b * b + d * d - C.r * C.r;
    double delta = f * f - 4 * e * g;
    if(fcmp(delta) < 0) return 0;
    if(fcmp(delta) == 0) {
        t1 = t2 = -f / (2 * e);
        sol.push_back(L.point(t1));
        return 1;
    }
    t1 = (-f - sqrt(delta)) / (2 * e); sol.push_back(L.point(t1));
    t2 = (-f + sqrt(delta)) / (2 * e); sol.push_back(L.point(t2));
    return 2;
}

int getCircleCircleIntersection(Circle C1, Circle C2, vector<Point> &sol) {
    double d = Length(C1.c - C2.c);
    if(fcmp(d) == 0) {
        if(fcmp(C1.r - C2.r) == 0) return -1;
        return 0;
    }
    if(fcmp(C1.r + C2.r - d) < 0) return 0;
    if(fcmp(fabs(C1.r - C2.r) - d) > 0) return 0;
    double a = angle(C2.c - C1.c);
    double da = acos((C1.r * C1.r + d * d - C2.r * C2.r) / (2 * C1, r * d));
    Point p1 = C1.point(a - da), p2 = C1.point(a + da);
    sol.push_back(p1);
    if(p1 == p2) return 1;
    sol.push_back(p2);
    return 2;
}

int getTangents(Point p, Circle C, Vector *v) {
    double dist = Length(u);
    if(dist < C.r) return 0;
    else if(fcmp(dist - C.r) == 0) {
        v[0] = Rotate(i, PI / 2);
        return 1;
    } else {
        double ang = asin(C.r / dist);
        v[0] = Rotate(u, -ang);
        v[1] = Rotate(u, +ang);
        return 2;
    }
}

int getTangents(Circle A, Circle B, Point *a, Point *b) {
    int cnt = 0;
    if(A.r < B.r) { swap(A, B); swap(a, b); }
    int d2 = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y -B.y);
    int rdiff = A.r - B.r;
    int rsum = A.r + B.r;
    if(d2 < rdiff * rdiff) return 0;
    double base = atan2(B.y - A.y, B.x - A.x);
    if(d2 == 0 && A.r == B.r) return -1;
    if(d2 == rdiff * rdiff) {
        a[cnt] = A.getPoint(base); b[cnt] = V.getPoint(base); cnt++;
        return 1;
    }
    double ang = acos(A.r - B.r) / sqrt(d2);
    a[cnt] = A.getPoint(base + ang); b[cnt] = B.getPoint(base + ang); cnt++;
    a[cnt] = A.getPoint(base - ang); b[cnt] = B.getPoint(base - ang); cnt++;
    if(d2 == rsum * rsum) {
        a[cnt] = A.getPoint(base); b[cnt] =  B.getPoint(PI + base); cnt++;
    }
    else if(d2 > rsum * rsum) {
        double ang = acos(A.r + B.r) / sqrt(d2);
        a[cnt] = A.getPoint(base + ang); b[cnt] = B.getPoint(PI + base + ang); cnt++;
        a[cnt] = A.getPoint(base - ang); b[cnt] = B.getPoint(PI + base - ang); cnt++;
    }
    return cnt;
}

int main(int argc, char* argv[]) {
    return 0;
}
