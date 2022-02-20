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

  double xct = 0.5e0;
  double yct = 0.7e0;

  double rad = sqrt(pow(xin-xct,2)+pow(yin-yct,2));
  double dist = rad - 0.15e0;

  return dist;
}
