#include "gda2020.h"
#include <math.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  projection proj;
  proj.a = 6378137.0;
  proj.invf = 298.257222101;
  proj.m0 = 0.9996;
  proj.false_easting = 500000;
  proj.false_northing = 10000000;

  geopoint point;
  point.easting = 232681.853;
  point.northing = 5867898.032;

  grid_to_geo(&point, 147, proj);
  printf("%f\n", point.latitude);
  printf("%f\n", point.longitude);

  geopoint point2;
  point2.latitude = point.latitude;
  point2.longitude = point.longitude;

  geo_to_grid(&point2, 147, proj);
  printf("%f\n%f\n", point2.easting, point2.northing);

  return 0;
}