#include "meshdefs.h"
#include "dog_math.h"

//  Signed distance function:
//
//      SignedDistance(x,y) < 0 inside the region
//      SignedDistance(x,y) = 0 on the boundary
//      SignedDistance(x,y) > 0 outside of the region
//
double SignedDistance(point pt)
{
  double xin = pt.x;
  double yin = pt.y;

  double d1 = 0.95e0 - yin; //top edge
  double d2 = -0.3e0 + xin; //-1*left edge
  double d3 = 0.6e0 - xin; //right edge
  double d4 = -0.65e0 + yin; //-1*bottom edge

  double dist = -Min(Min(d1,d2),Min(d3,d4));

  return dist;
}
