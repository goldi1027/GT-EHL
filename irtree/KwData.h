#pragma once
#include <sys/time.h>
#include <stdio.h>
#include <unistd.h>
#include <cmath>
#include "IRStarTree.h"

namespace irstar {
const double NORMALIZE_SCALE = 100000;

typedef struct Point
{
    Coord coord[DIM];

    friend istream& operator>>(istream& in, Point& p);

    bool operator< (const Point& p) const;

    Point();
    Point(Coord c[DIM]);
    Point(Coord x, Coord y);

    void print();

    double distance2(Point& p);
}Point;

inline double Point::distance2(Point& p)
{
    double dis = 0;
    for(size_t dim = 0; dim < DIM; dim++)
    {
        double diff = coord[dim] - p.coord[dim];
        dis += diff * diff;
    }
    return dis;
}

typedef vector<Point> Point_V;
typedef vector<LeafNodeEntry> Entry_V;
}
