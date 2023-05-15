#include "KwData.h"

/*******************************************************************************
 * Point
 ******************************************************************************/
namespace irstar{
  Point::Point() { }

  Point::Point(Coord c[DIM])
  {
      for(size_t dim = 0; dim < DIM; dim++)
          coord[dim] = c[dim];
  }

  Point::Point(Coord x, Coord y)
  {
      coord[0] = x;
      coord[1] = y;
  }

  istream& operator>>(istream& in, Point& p)
  {
      for(size_t dim = 0; dim < DIM; dim++)
          in >> p.coord[dim];
      return in;
  }

  bool Point::operator< (const Point& p) const
  {
      for(size_t dim = 0; dim < DIM; dim++)
          if(coord[dim] != p.coord[dim])
              return coord[dim] < p.coord[dim];
      return true;
  }

  void Point::print()
  {
      cout << "<" << coord[0];
      for(size_t dim = 0; dim < DIM; dim++)
          cout << "," << coord[dim];
      cout << ">" << endl;
  }
}
