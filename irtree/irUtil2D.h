#pragma once
#include <cstdlib>
#include <limits>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <list>
#include <cmath>

#define INF_N -numeric_limits<double>::max()
#define INF_P numeric_limits<double>::max()

namespace irstar {

using namespace std;

const double EPS = 1e-9;

typedef double Coord;

typedef struct Line
{
    Line()
        : a(0), b(0), vertical(false), c(0)
        { }

    // <a,b> representing y = ax + b
    Coord a;
    Coord b;

    // vertical line is represented by x = c
    bool vertical;
    Coord c;

    void print()
    {
        if(vertical)
            cout << "x = " << c << endl;
        else
            cout << "y = " << a << "x " << (b > 0 ? "+ " : "- ") << abs(b) << endl;
    }
}Line2D;

typedef vector<Line> Line_V;

typedef struct Vertex
{
    Vertex()
        : x(0), y(0)
        { }

    Coord x;
    Coord y;

    bool operator< (const Vertex& p) const;
}Vertex;

inline bool Vertex::operator< (const Vertex& v) const
{
    if(x == v.x)
        return y < v.y;
    return x < v.x;
}

typedef vector<Vertex> Vertex_V;
typedef set<Vertex> Vertex_S;

typedef double Radians;

class Util2D
{
    public:
        Util2D() { }

        static bool lineCross(Line& line1, Line& line2, Vertex& cross); // intersection of two lines
        static Coord y(Line& line, Coord x);
        static void linePass(Vertex& vertex1, Vertex& vertex2, Line& line); // line passing through two vertex

        static const int STRICT_ABOVE = 1;
        static const int STRICT_BELOW = -1;
        static const int AROUND = 0;
        static const int STRICT_LEFT = -2;
        static const int STRICT_RIGHT = 2;
        static int isAboveLine(Line& line, Vertex& p);

        static const int LEFT_TURN = 1;
        static const int RIGHT_TURN = -1;
        static const int COLLINEAR = 0;
        static int turn(Vertex& p1, Vertex& p2, Vertex& p3);

        static Coord distance2(Vertex& p1, Vertex& p2);
        static Coord perpendicularDistance2(Line& line, Vertex& vertex);
        static Coord distance2(Line& line, Vertex& from, Vertex& to, Vertex& vertex);

        constexpr static const Radians H_PI = 0.5 * M_PI; // half pi
        constexpr static const Radians O_H_PI = 1.5 * M_PI; // one and half pi
        constexpr static const Radians T_PI = 2 * M_PI; // two pi
        static Radians angle(Coord x0, Coord y0, Coord x, Coord y); // [0:T_PI)

        static double getRandom(); // [0, 1]
};

inline bool Util2D::lineCross(Line& line1, Line& line2, Vertex& cross)
{
    // if both lines are vertical, no intersection
    if(line1.vertical && line2.vertical)
        return false;
    // if they are parallel, no intersection
    if(line1.a == line2.a)
        return false;

    if(line1.vertical)
    {
        cross.x = line1.c;
        cross.y = line2.a * cross.x + line2.b;
        return true;
    }
    else
    {
        if(line2.vertical)
        {
            cross.x = line2.c;
            cross.y = line1.a * cross.x + line1.b;
            return true;
        }
        else
        {
            cross.x = (line2.b - line1.b) / (line1.a - line2.a);
            cross.y = line1.a * cross.x + line1.b;
            return true;
        }
    }
}

inline Coord Util2D::y(Line& line, Coord x)
{
    //assert(!line.vertical);
    return line.a * x + line.b;
}

inline void Util2D::linePass(Vertex& vertex1, Vertex& vertex2, Line& line)
{
    if(vertex1.x == vertex2.x)
    {
        line.vertical = true;
        line.c = vertex1.x;
    }
    else
    {
        line.a = (vertex1.y - vertex2.y) / (vertex1.x - vertex2.x);
        line.b = vertex1.y - line.a * vertex1.x;
    }
}

inline int Util2D::isAboveLine(Line& line, Vertex& p)
{
    if(line.vertical)
    {
        if(p.x > line.c + EPS)
            return STRICT_RIGHT;
        if(p.x + EPS < line.c)
            return STRICT_LEFT;
    }
    else
    {
        Coord y = line.a * p.x + line.b;
        if(p.y > y + EPS)
            return STRICT_ABOVE;
        if(p.y + EPS < y)
            return STRICT_BELOW;
    }
    return AROUND;
}

inline int Util2D::turn(Vertex& p1, Vertex& p2, Vertex& p3)
{
    //ofstream debug;
    //debug.open ("data/debug/debug.txt", ios::out | ios::app );
    //assert(debug.is_open());

    double num1 = (p2.x - p1.x) * (p3.y - p1.y);
    double num2 = (p2.y - p1.y) * (p3.x - p1.x);
    double res = num1 - num2;
    //debug << "[" << res << "]";
    if(abs(res) < EPS)
        return COLLINEAR;
    if(res > 0)
        return LEFT_TURN;
    return RIGHT_TURN;

    //debug.close();
}

inline Coord Util2D::distance2(Vertex& p1, Vertex& p2)
{
    Coord diffx = p1.x - p2.x, diffy = p1.y - p2.y;
    return diffx * diffx + diffy * diffy;
}

inline Coord Util2D::perpendicularDistance2(Line& line, Vertex& vertex)
{
    Coord num = line.a * vertex.x - vertex.y + line.b;
    return num * num / (line.a * line.a + 1);
}

inline Coord Util2D::distance2(Line& line, Vertex& from, Vertex& to, Vertex& vertex)
{
    Coord minDis = min(distance2(vertex, from), distance2(vertex, to));
    minDis = min(minDis, perpendicularDistance2(line, vertex));
    return minDis;
}

inline Radians Util2D::angle(Coord x0, Coord y0, Coord x, Coord y)
{
    if(x0 == x)
    {
        if(y > y0)
            return H_PI;
        else if(y < y0)
            return O_H_PI;
        else
            assert(false);
    }
    if(y0 == y)
    {
        if(x > x0)
            return 0;
        else if(x < x0)
            return M_PI;
        else
            assert(false);
    }
    double value = (y - y0) / (x - x0);
    double radians = atan(value);
    if(x > x0)
    {
        if(y > y0)
            return radians;
        else
            return radians + T_PI;
    }
    else
        return radians + M_PI;
}

inline double Util2D::getRandom()
{
    double num = rand();
    return num / RAND_MAX;
}

typedef struct HalfSpace
{
    HalfSpace() { }

    Line line;
    Vertex from;
    Vertex to;
    int isAbove; // is q above line

    static void halfSpace(Vertex& p, Vertex& q, HalfSpace& hSpace);
}HalfSpace;

inline void HalfSpace::halfSpace(Vertex& p, Vertex& q, HalfSpace& hSpace)
{
    if(p.y == q.y)
    {
        hSpace.line.vertical = true;
        Coord mid = (p.x + q.x) / 2;
        hSpace.line.c = mid;
        if(q.x > mid + EPS)
            hSpace.isAbove = Util2D::STRICT_RIGHT;
        else if(q.x + EPS < mid)
            hSpace.isAbove = Util2D::STRICT_LEFT;
        else
            assert(false);
    }
    else
    {
        hSpace.line.a = (q.x - p.x) / (p.y - q.y);
        Coord midx = (p.x + q.x) / 2, midy = (p.y + q.y) / 2;
        hSpace.line.b = midy - hSpace.line.a * midx;
        hSpace.isAbove = Util2D::isAboveLine(hSpace.line, q);
        assert(hSpace.isAbove != Util2D::AROUND);
    }
}

typedef vector<HalfSpace> HalfSpace_V;

}
