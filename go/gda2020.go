package gda2020

import (
    "math"
)

type Projection struct {
	A float64
	Invf float64
	M0 float64
	FalseEasting float64
	FalseNorthing float64
}

// Convert degrees to radians.
func Radians(degrees float64) float64 {
	return (degrees * math.Pi) / 180
}

// Convert radians to degrees.
func Degrees(Radians float64) float64 {
	return (Radians * 180) / math.Pi
}

// Convert degrees, minutes, seconds (dd.mmss) to decimal degrees
func Dms2dec(dms float64) float64 {
	sign := 1.0
	if dms < 0 {
		sign = -1.0
		dms = dms * sign
	}

	degrees := float64(int(dms))
	minutes := float64(int((dms * 100))) - degrees*100.0
	seconds := (((dms - degrees) * 100) - minutes) * 100
	decdeg := degrees + float64(minutes)/60 + float64(seconds)/3600
	return decdeg * sign
}

// Convert decimal degrees to degrees, minutes, seconds (dd.mmss)
func Dec2dms(decdeg float64) float64 {
	sign := 1.0
	if decdeg < 0 {
		sign = -1.0
		decdeg = decdeg * sign
	}

	degrees := float64(int(decdeg))
	minutes := float64(int(decdeg*60)) - degrees*60.0
	seconds := (((decdeg - degrees) * 3600) - minutes*60.0)
	dms := degrees + minutes/100 + seconds/10000
	return dms * sign
}

// Convert geodetic to cartesian
func GeodeticToCartesian(latitude, longitude, h float64, proj Projection)(float64, float64, float64){
    f := 1/proj.Invf 
    e2 := 2*f-f*f // eccentricity squared of reference ellipsoid (eq 12)

    // convert dd.mmss to radians
    latitude_radians := Radians(Dms2dec(latitude))
    longitude_radians := Radians(Dms2dec(longitude))

    // radius of curvature at a point on ellipsoid with respect to prime vertical through
    // that point (eq 11)
    part := math.Pow(math.Sin(latitude_radians),2)
    v := proj.A/(math.Sqrt(1- e2*part))

    // cartesian coordinates
    X := (v+h)*math.Cos(latitude_radians)*math.Cos(longitude_radians) // (eq 8)
    Y := (v+h)*math.Cos(latitude_radians)*math.Sin(longitude_radians) // (eq 9)
    Z := ((1-e2)*v + h)*math.Sin(latitude_radians) // (eq 10)

    return X, Y, Z
}

// Convert cartesian coordinates to geodetic
func CartesianToGeodetic(X, Y, Z float64, proj Projection)(float64, float64, float64){
	a := proj.A // semi-major axis
    invf := proj.Invf
    f := 1/invf // flattening
    e2 := 2*f-f*f // eccentricity squared of reference ellipsoid
    subf := 1 - f

    p := math.Sqrt(X*X + Y*Y) // radius of curvature at a point on ellipsoid with respect to meridian (eq 4)
    r := math.Sqrt(p*p + Z*Z) // (eq 6)
    u := math.Atan((Z/p)*(subf + (e2 * a/r))) // parametric latitude (eq 5)

    longitude_radians := math.Atan2(Y, X) // (eq 1)

    // (eq 2)
    top := Z*subf+e2*a*math.Pow(math.Sin(u),3)
    bot := subf*(p-e2*a*math.Pow(math.Cos(u), 3))
    latitude_radians := math.Atan(top/bot)

    // (eq 3)
    part1 := p * math.Cos(latitude_radians)
    part2 := Z * math.Sin(latitude_radians)
    part3 := a * math.Sqrt(1 - e2 * math.Pow(math.Sin(latitude_radians), 2))
    h := part1 + part2 - part3

    return Dec2dms(Degrees(latitude_radians)), Dec2dms(Degrees(longitude_radians)), h
}


func GeoToGrid(latitude, longitude, central_meridian float64, proj Projection)(float64, float64, float64, float64){


    // Change lat/long to decimal degrees and convert to radians
    rlat := Radians(Dms2dec(latitude))
    rlong := Radians(Dms2dec(longitude))
    rcentral_meridian := Radians(central_meridian)

    a := proj.A
    invf := proj.Invf
    m0 := proj.M0
    false_easting := proj.FalseEasting
    false_northing := proj.FalseNorthing

    // Calculate geometrical consTants
    f := 1.0/invf
    b := a * (1-f)  // semi - minor axis

    // Eccentricity (eq 18)
    e2 := 2 * f - f * f  // = f*(2-f) = (a^2-b^2)/a^2
    e := math.Sqrt(e2)

    // Compute 3rd flattening and powers (eq 19)
    n := (a - b)/(a + b)
    n2 := n * n
    n3 := n * n2
    n4 := n2 * n2
    n5 := n3 * n2
    n6 := n2 * n4
    n7 := n4 * n3
    n8 := n4 * n4

    // Rectifying Radius A (eq 20)
    A := (a/(1+n)) * (1+(1.0/4.0) * n2 + (1.0/64) * n4 + (1.0/256) * n6 + (25.0/16384) * n8)

    // Calculate conformal latitude (eq 22, 23)
    sigma := math.Sinh( e*math.Atanh(( e * math.Tan(rlat)) / (math.Sqrt( 1 + math.Tan(rlat) * math.Tan(rlat)))))
    conformal_lat := math.Tan(rlat) * math.Sqrt(1 + sigma * sigma) - sigma * 
                                                                   math.Sqrt(1 + math.Tan(rlat) * math.Tan(rlat))

    // Compute the coefficients (eq 21)
    a2 := (1.0/2) * n - (2.0/3) * n2 + (5.0/16) * n3 + (41.0/180) * n4 - (127.0/288) * n5 + (7891.0/37800) * n6 +
         (72161.0/387072) * n7 - (18975107.0/50803200) * n8
    a4 := (13.0/48) * n2 - (3.0/5) * n3 + (557.0/1440) * n4 + (281.0/630) * n5 - (1983433.0/1935360) * n6 + 
         (13769.0/28800) * n7 + (148003883.0/174182400) * n8
    a6 := (61.0/240) * n3 - (103.0/140) * n4 + (15061.0/26880) * n5 + (167603.0/181440) * n6 - 
         (67102379.0/29030400) * n7 + (79682431.0/79833600) * n8
    a8 := (49561.0/161280) * n4 - (179.0/168) * n5 + (6601661.0/7257600) * n6 + (97445.0/49896) * 
                                                                               n7 - (40176129013.0/7664025600) * n8
    a10 := (34729.0/80640)* n5-(3418889.0/1995840)* n6+(14644087.0/9123840)* n7+(2605413599.0/622702080) * n8
    a12 := (212378941.0/319334400) * n6 - (30705481.0/10378368) * n7 + (175214326799.0/58118860800) * n8
    a14 := (1522256789.0/1383782400) * n7 - (16759934899.0/3113510400) * n8
    a16 := (1424729850961.0/743921418240) * n8

    // Find w by subtracting central meridian from longitude (eq 24)
    w := rlong - rcentral_meridian // Converted to radians

    // Compute the gauss-Schreiber coordinates (eq 25, 26)
    u := a * math.Atan( conformal_lat/math.Cos(w))
    v := a * math.Asinh(math.Sin(w) / (math.Sqrt( conformal_lat * conformal_lat + math.Cos(w) * math.Cos( w))))

    // Calculate partial solutions for x and y to make calcs easier (eq 27, 28)
    x1 := math.Cos(2 * (u / a)) * math.Sinh(2 * (v / a))
    x2 := math.Cos(4 * (u / a)) * math.Sinh(4 * (v / a))
    x3 := math.Cos(6 * (u / a)) * math.Sinh(6 * (v / a))
    x4 := math.Cos(8 * (u / a)) * math.Sinh(8 * (v / a))
    x5 := math.Cos(10 * (u / a)) * math.Sinh(10 * (v / a))
    x6 := math.Cos(12 * (u / a)) * math.Sinh(12 * (v / a))
    x7 := math.Cos(14 * (u / a)) * math.Sinh(14 * (v / a))
    x8 := math.Cos(16 * (u / a)) * math.Sinh(16 * (v / a))

    y1 := math.Sin(2 * (u / a)) * math.Cosh(2 * (v / a))
    y2 := math.Sin(4 * (u / a)) * math.Cosh(4 * (v / a))
    y3 := math.Sin(6 * (u / a)) * math.Cosh(6 * (v / a))
    y4 := math.Sin(8 * (u / a)) * math.Cosh(8 * (v / a))
    y5 := math.Sin(10 * (u / a)) * math.Cosh(10 * (v / a))
    y6 := math.Sin(12 * (u / a)) * math.Cosh(12 * (v / a))
    y7 := math.Sin(14 * (u / a)) * math.Cosh(14 * (v / a))
    y8 := math.Sin(16 * (u / a)) * math.Cosh(16 * (v / a))

    // Calculate partial solutions for q and p to make calcs easier
    q1 := math.Sin(2 * (u / a)) * math.Sinh(2 * (v / a))
    q2 := math.Sin(4 * (u / a)) * math.Sinh(4 * (v / a))
    q3 := math.Sin(6 * (u / a)) * math.Sinh(6 * (v / a))
    q4 := math.Sin(8 * (u / a)) * math.Sinh(8 * (v / a))
    q5 := math.Sin(10 * (u / a)) * math.Sinh(10 * (v / a))
    q6 := math.Sin(12 * (u / a)) * math.Sinh(12 * (v / a))
    q7 := math.Sin(14 * (u / a)) * math.Sinh(14 * (v / a))
    q8 := math.Sin(16 * (u / a)) * math.Sinh(16 * (v / a))

    p1 := math.Cos(2 * (u / a)) * math.Cosh(2 * (v / a))
    p2 := math.Cos(4 * (u / a)) * math.Cosh(4 * (v / a))
    p3 := math.Cos(6 * (u / a)) * math.Cosh(6 * (v / a))
    p4 := math.Cos(8 * (u / a)) * math.Cosh(8 * (v / a))
    p5 := math.Cos(10 * (u / a)) * math.Cosh(10 * (v / a))
    p6 := math.Cos(12 * (u / a)) * math.Cosh(12 * (v / a))
    p7 := math.Cos(14 * (u / a)) * math.Cosh(14 * (v / a))
    p8 := math.Cos(16 * (u / a)) * math.Cosh(16 * (v / a))

    // Calculate q and p for calculating point scale factor (eq 33, 34)
    q := - (2 * a2 * q1 + 4 * a4 * q2 + 6* a6* q3 + 8 * a8 * q4 +
           10 * a10 * q5 + 12 * a12 * q6 + 14 * a14 * q7 + 16 * a16 * q8)
    p := 1 + (2 * a2 * p1 + 4 * a4 * p2 + 6 * a6 * p3 + 8 * a8 *
             p4 + 10 * a10 * p5 + 12 * a12 * p6 + 14 * a14 * p7 + 16 * a16 * p8)


    // Calculate point scale factor m (eq 35)
    m := m0 * (A / a)*math.Sqrt(q * q + p * p) * (math.Sqrt(1+(math.Tan(rlat)*math.Tan(rlat))) * math.Sqrt(1 - e2*(math.Sin(rlat)*math.Sin(rlat))))/ 
        math.Sqrt(conformal_lat * conformal_lat+math.Cos(w)*math.Cos(w))

    // compute grid convergence (eq 36)
    grid_conv := math.Atan(q / p)+math.Atan((conformal_lat*math.Tan(w))/math.Sqrt(1 + conformal_lat * conformal_lat))*180/math.Pi

    // Calculate transverse Mercator coordinates (eq 29, 30)
    X := A*((v / a) + a2 * x1 + a4 * x2 + a6 * x3 + a8 * x4 + a10 * x5 + a12 * x6 + a14 * x7 + a16 * x8)
    Y := A*((u / a) + a2 * y1 + a4 * y2 + a6 * y3 + a8 * y4 + a10 * y5 + a12 * y6 + a14 * y7 + a16 * y8)

    // Calculate scaled coordinates with false easting and northing
    // (eq 31, 32)
    easting := false_easting + m0 * X
    northing := false_northing + m0 * Y

	return easting, northing, m, grid_conv
}