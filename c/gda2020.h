#include <math.h>
#define radians(angle_degrees) ((angle_degrees) * M_PI / 180.0)
#define degrees(angle_radians) ((angle_radians) * 180.0 / M_PI)

/*
 a; ellipsoid semi - major axis
 invf: inverse of flattening 1/f
 m0: central scale factor
 false_easting: false easting for projection
 false_northing: false northing for projection 
 */
typedef struct {
  double a;
  double invf;
  double m0;
  double false_easting;
  double false_northing;
} projection;

/*
 Geographic Point Structure
 */
typedef struct {
  double easting;
  double northing;
  double latitude;
  double longitude;
  double m;
  double grid_conv;
} geopoint;

/*
 Convert decimal degrees to dd.mmss
 */
double dec_to_dms(double decdeg) {
  int sign;
  if (decdeg < 0) {
    sign = -1;
    decdeg = decdeg * sign;
  } else {
    sign = 1;
  }

  double deg = (int)decdeg;
  double min = (int)((decdeg - deg) * 60);
  double sec = (((decdeg - deg) * 60) - min) * 60;

  return (deg + min * .01 + sec * .0001) * sign;
}

/* 
 Convert dd.mmss to decimal degrees
 */
double dms_to_dec(double dmsdeg) {

  int sign;
  double dec;
  int deg;

  if (dmsdeg < 0) {
    sign = -1;
    dmsdeg = dmsdeg * sign;
  } else {
    sign = 1;
  }

  deg = (int)dmsdeg;
  int min = floor((dmsdeg - deg) * 100.0); //,14-size);
  double sec = (((dmsdeg - deg) * 100.0) - min) * 100.0;

  dec = deg + min * 1.0 / 60 + sec * 1.0 / 3600;
  return dec * sign;
}

/*
Convert grid coordinates to geographic using Kruger n-series equations

        See Sec 4.1.1 (pg 36) of Geocentric Datum of Australia 2020
        Technical Manual 1.2

        Parameters
        ----------
        easting: float
            easting of the grid coordinate
        northing: float
            northing of the grid coordinate
        central_meridian: float
            central meridian in degrees
        projection: Proj class
            a class which contains the projection information

        Returns
        -------
        latitude: float
            the latitude in dd.mmssss
        longitude: float
            the longitude in dd.mmssss
 */
void grid_to_geo(geopoint *point, double central_meridian, projection proj) {
  double a = proj.a;
  double invf = proj.invf;
  double m0 = proj.m0;
  double false_easting = proj.false_easting;
  double false_northing = proj.false_northing;

  /* Calculate geometrical constants */
  double f = 1.0 / invf;
  double b = a * (1 - f); /* semi - minor axis */

  /* Eccentricity (eq 18) */
  double e2 = 2 * f - f * f; // = f*(2-f) = (a^2-b^2)/a^2
  double e = sqrt(e2);

  // Compute 3rd flattening and powers (eq 19)
  double n = (a - b) / (a + b);
  double n2 = n * n;
  double n3 = n * n2;
  double n4 = n2 * n2;
  double n5 = n3 * n2;
  double n6 = n2 * n4;
  double n7 = n4 * n3;
  double n8 = n4 * n4;

  // Rectifying Radius A (eq 20)
  double A = (a / (1 + n)) * (1 + (1.0 / 4.0) * n2 + (1.0 / 64) * n4 +
                              (1.0 / 256) * n6 + (25.0 / 16384) * n8);

  // calculate beta values (eq 37)
  double b2 = (-1.0 / 2) * n + (2.0 / 3) * n2 - (37.0 / 96) * n3 +
              (1.0 / 360) * n4 + (81.0 / 512) * n5 - (96199.0 / 604800) * n6 +
              (5406467.0 / 38707200) * n7 - (7944359.0 / 67737600) * n8;

  double b4 = (-1.0 / 48) * n2 - (1.0 / 15) * n3 + (437.0 / 1440) * n4 -
              (46.0 / 105) * n5 + (1118711.0 / 3870720) * n6 -
              (51841.0 / 1209600) * n7 - (24749483.0 / 348364800) * n8;

  double b6 = (-17.0 / 480) * n3 + (37.0 / 840) * n4 + (209.0 / 4480) * n5 -
              (5569.0 / 90720) * n6 - (9261899.0 / 58060800) * n7 +
              (6457463.0 / 17740800) * n8;

  double b8 = (-4397.0 / 161280) * n4 + (11.0 / 504) * n5 +
              (830251.0 / 7257600) * n6 - (466511.0 / 2494800) * n7 -
              (324154477.0 / 7664025600) * n8;

  double b10 = (-4583.0 / 161280) * n5 + (108847.0 / 3991680) * n6 +
               (8005831.0 / 63866880) * n7 - (22894433.0 / 124540416) * n8;

  double b12 = (-20648693.0 / 638668800) * n6 + (16363163.0 / 518918400) * n7 +
               (2204645983.0 / 12915302400) * n8;

  double b14 =
      (-219941297.0 / 5535129600) * n7 + (497323811.0 / 12454041600) * n8;

  double b16 = (-191773887257.0 / 3719607091200) * n8;

  // transverse Mercator coordinates X,Y (eq 38, 39)
  double X = (point->easting - false_easting) / m0;
  double Y = (point->northing - false_northing) / m0;

  // Calculate partial solutions for q and p to make calcs easier
  double x1 = cos(2 * (Y / A)) * sinh(2 * (X / A));
  double x2 = cos(4 * (Y / A)) * sinh(4 * (X / A));
  double x3 = cos(6 * (Y / A)) * sinh(6 * (X / A));
  double x4 = cos(8 * (Y / A)) * sinh(8 * (X / A));
  double x5 = cos(10 * (Y / A)) * sinh(10 * (X / A));
  double x6 = cos(12 * (Y / A)) * sinh(12 * (X / A));
  double x7 = cos(14 * (Y / A)) * sinh(14 * (X / A));
  double x8 = cos(16 * (Y / A)) * sinh(16 * (X / A));

  double y1 = sin(2 * (Y / A)) * cosh(2 * (X / A));
  double y2 = sin(4 * (Y / A)) * cosh(4 * (X / A));
  double y3 = sin(6 * (Y / A)) * cosh(6 * (X / A));
  double y4 = sin(8 * (Y / A)) * cosh(8 * (X / A));
  double y5 = sin(10 * (Y / A)) * cosh(10 * (X / A));
  double y6 = sin(12 * (Y / A)) * cosh(12 * (X / A));
  double y7 = sin(14 * (Y / A)) * cosh(14 * (X / A));
  double y8 = sin(16 * (Y / A)) * cosh(16 * (X / A));

  // gauss-schreiber ratios (eq 42)
  double nn = (X / A) + (b2 * x1) + (b4 * x2) + (b6 * x3) + (b8 * x4) +
              (b10 * x5) + (b12 * x6) + (b14 * x7) + (b16 * x8);

  // gauss-schreiber ratios (eq 43)
  double ee = (Y / A) + b2 * y1 + b4 * y2 + b6 * y3 + b8 * y4 + b10 * y5 +
              b12 * y6 + b14 * y7 + b16 * y8;

  // (eq 44)
  double t_dash = sin(ee) / (sqrt(sinh(nn) * sinh(nn) + cos(ee) * cos(ee)));

  double t = t_dash;

  // solve t = tan(phi) by Newton-Raphson iteration
  // initial value for t = t'
  for (int i = 0; i < 5; i++) {
    // (eq 46)
    double sigma = sinh(e * atanh(e * t / sqrt(1 + t * t)));
    // (eq 48)
    double ft = t * sqrt(1 + sigma * sigma) - sigma * sqrt(1 + t * t) - t_dash;
    // (eq 49)
    double ft_dash = (sqrt(1 + sigma * sigma) * sqrt(1 + t * t) - sigma * t) *
                     ((1 - e2) * sqrt(1 + t * t)) / (1 + (1 + e2) * t * t);
    // (eq 47)
    t = t - (ft / ft_dash);
  }

  // longitude difference omega (eq 51)
  double omega = atan(sinh(nn) / cos(ee));

  // longitude calculation (eq 52)
  double calculated_longitude = degrees(radians(central_meridian) + omega);

  // latitude calculation (eq 50)
  double calculated_latitude = degrees(atan(t));

  point->latitude = dec_to_dms(calculated_latitude);
  point->longitude = dec_to_dms(calculated_longitude);
}

/*  
Convert geographic to grid coordinates using Kruger n-series equations

    See Sec 4.1.1 (pg 36) of Geocentric Datum of Australia
    2020 Technical Manual 1.2

    Parameters
    ----------
    latitude: float
        latitude in dd.mmssss
    longitude: float
        longitude in dd.mmssss
    central_meridian: float
        central meridian in degrees
    projection: Proj class
        a class which contains the projection information

    Returns
    -------
    easting: float
        the easting grid coordinate
    northing: float
        the northing grid coordinate
    m: float
        the point scale factor
    grid_conv: float
        the grid convergence
*/
void geo_to_grid(geopoint *point, double central_meridian, projection proj) {

  /* Change lat/long to decimal degrees and convert to radians */
  double rlat = radians(dms_to_dec(point->latitude));
  double rlong = radians(dms_to_dec(point->longitude));
  double rcentral_meridian = radians(central_meridian);

  double a = proj.a;
  double invf = proj.invf;
  double m0 = proj.m0;
  double false_easting = proj.false_easting;
  double false_northing = proj.false_northing;

  /* Calculate geometrical constants */
  double f = 1.0 / invf;
  /* semi - minor axis */
  double b = a * (1 - f);

  /* Eccentricity (eq 18) */
  /* = f*(2-f) = (a^2-b^2)/a^2 */
  double e2 = 2 * f - f * f;
  double e = sqrt(e2);

  /* Compute 3rd flattening and powers (eq 19) */
  double n = (a - b) / (a + b);
  double n2 = n * n;
  double n3 = n * n2;
  double n4 = n2 * n2;
  double n5 = n3 * n2;
  double n6 = n2 * n4;
  double n7 = n4 * n3;
  double n8 = n4 * n4;

  /*  Rectifying Radius A (eq 20) */
  double A = (a / (1 + n)) * (1 + (1.0 / 4.0) * n2 + (1.0 / 64) * n4 +
                              (1.0 / 256) * n6 + (25.0 / 16384) * n8);

  /* Calculate conformal latitude (eq 22, 23) */
  double sigma =
      sinh(e * atanh((e * tan(rlat)) / (sqrt(1 + tan(rlat) * tan(rlat)))));

  double conformal_lat = tan(rlat) * sqrt(1 + sigma * sigma) -
                         sigma * sqrt(1 + tan(rlat) * tan(rlat));

  /* Compute the coefficients (eq 21) */
  double a2 = (1.0 / 2) * n - (2.0 / 3) * n2 + (5.0 / 16) * n3 +
              (41.0 / 180) * n4 - (127.0 / 288) * n5 + (7891.0 / 37800) * n6 +
              (72161.0 / 387072) * n7 - (18975107.0 / 50803200) * n8;

  double a4 = (13.0 / 48) * n2 - (3.0 / 5) * n3 + (557.0 / 1440) * n4 +
              (281.0 / 630) * n5 - (1983433.0 / 1935360) * n6 +
              (13769.0 / 28800) * n7 + (148003883.0 / 174182400) * n8;

  double a6 = (61.0 / 240) * n3 - (103.0 / 140) * n4 + (15061.0 / 26880) * n5 +
              (167603.0 / 181440) * n6 - (67102379.0 / 29030400) * n7 +
              (79682431.0 / 79833600) * n8;

  double a8 = (49561.0 / 161280) * n4 - (179.0 / 168) * n5 +
              (6601661.0 / 7257600) * n6 + (97445.0 / 49896) * n7 -
              (40176129013.0 / 7664025600) * n8;

  double a10 = (34729.0 / 80640) * n5 - (3418889.0 / 1995840) * n6 +
               (14644087.0 / 9123840) * n7 + (2605413599.0 / 622702080) * n8;

  double a12 = (212378941.0 / 319334400) * n6 - (30705481.0 / 10378368) * n7 +
               (175214326799.0 / 58118860800) * n8;

  double a14 =
      (1522256789.0 / 1383782400) * n7 - (16759934899.0 / 3113510400) * n8;

  double a16 = (1424729850961.0 / 743921418240) * n8;

  /* Find w by subtracting central meridian from longitude (eq 24) */
  double w = rlong - rcentral_meridian;

  /* Compute the gauss-Schreiber coordinates (eq 25, 26) */
  double u = a * atan(conformal_lat / cos(w));
  double v = a * asinh(sin(w) /
                       (sqrt(conformal_lat * conformal_lat + cos(w) * cos(w))));

  /* Calculate partial solutions for x and y to
  make calcs easier (eq 27, 28) */
  double x1 = cos(2 * (u / a)) * sinh(2 * (v / a));
  double x2 = cos(4 * (u / a)) * sinh(4 * (v / a));
  double x3 = cos(6 * (u / a)) * sinh(6 * (v / a));
  double x4 = cos(8 * (u / a)) * sinh(8 * (v / a));
  double x5 = cos(10 * (u / a)) * sinh(10 * (v / a));
  double x6 = cos(12 * (u / a)) * sinh(12 * (v / a));
  double x7 = cos(14 * (u / a)) * sinh(14 * (v / a));
  double x8 = cos(16 * (u / a)) * sinh(16 * (v / a));

  double y1 = sin(2 * (u / a)) * cosh(2 * (v / a));
  double y2 = sin(4 * (u / a)) * cosh(4 * (v / a));
  double y3 = sin(6 * (u / a)) * cosh(6 * (v / a));
  double y4 = sin(8 * (u / a)) * cosh(8 * (v / a));
  double y5 = sin(10 * (u / a)) * cosh(10 * (v / a));
  double y6 = sin(12 * (u / a)) * cosh(12 * (v / a));
  double y7 = sin(14 * (u / a)) * cosh(14 * (v / a));
  double y8 = sin(16 * (u / a)) * cosh(16 * (v / a));

  /* Calculate partial solutions for q and p to make calcs easier */
  double q1 = sin(2 * (u / a)) * sinh(2 * (v / a));
  double q2 = sin(4 * (u / a)) * sinh(4 * (v / a));
  double q3 = sin(6 * (u / a)) * sinh(6 * (v / a));
  double q4 = sin(8 * (u / a)) * sinh(8 * (v / a));
  double q5 = sin(10 * (u / a)) * sinh(10 * (v / a));
  double q6 = sin(12 * (u / a)) * sinh(12 * (v / a));
  double q7 = sin(14 * (u / a)) * sinh(14 * (v / a));
  double q8 = sin(16 * (u / a)) * sinh(16 * (v / a));

  double p1 = cos(2 * (u / a)) * cosh(2 * (v / a));
  double p2 = cos(4 * (u / a)) * cosh(4 * (v / a));
  double p3 = cos(6 * (u / a)) * cosh(6 * (v / a));
  double p4 = cos(8 * (u / a)) * cosh(8 * (v / a));
  double p5 = cos(10 * (u / a)) * cosh(10 * (v / a));
  double p6 = cos(12 * (u / a)) * cosh(12 * (v / a));
  double p7 = cos(14 * (u / a)) * cosh(14 * (v / a));
  double p8 = cos(16 * (u / a)) * cosh(16 * (v / a));

  /* Calculate q and p for calculating point scale factor (eq 33, 34) */
  double q = -(2 * a2 * q1 + 4 * a4 * q2 + 6 * a6 * q3 + 8 * a8 * q4 +
               10 * a10 * q5 + 12 * a12 * q6 + 14 * a14 * q7 + 16 * a16 * q8);
  double p =
      1 + (2 * a2 * p1 + 4 * a4 * p2 + 6 * a6 * p3 + 8 * a8 * p4 +
           10 * a10 * p5 + 12 * a12 * p6 + 14 * a14 * p7 + 16 * a16 * p8);

  /* Calculate point scale factor m (eq 35) */
  double m = m0 * (A / a) * sqrt(q * q + p * p) *
             (sqrt(1 + (tan(rlat) * tan(rlat))) *
              sqrt(1 - e2 * (sin(rlat) * sin(rlat)))) /
             sqrt(conformal_lat * conformal_lat + cos(w) * cos(w));

  /* compute grid convergence (eq 36) */
  double grid_conv =
      atan(q / p) +
      atan((conformal_lat * tan(w)) / sqrt(1 + conformal_lat * conformal_lat)) *
          180 / M_PI;

  /* Calculate transverse Mercator coordinates (eq 29, 30) */
  double X = A * ((v / a) + a2 * x1 + a4 * x2 + a6 * x3 + a8 * x4 + a10 * x5 +
                  a12 * x6 + a14 * x7 + a16 * x8);

  double Y = A * ((u / a) + a2 * y1 + a4 * y2 + a6 * y3 + a8 * y4 + a10 * y5 +
                  a12 * y6 + a14 * y7 + a16 * y8);

  /* Calculate scaled coordinates with false easting and northing
      (eq 31, 32)
  */
  double easting = false_easting + m0 * X;
  double northing = false_northing + m0 * Y;

  point->easting = easting;
  point->northing = northing;
  point->m = m;
  point->grid_conv = grid_conv;
  // return easting, northing, m, grid_conv
}