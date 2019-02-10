import numpy as np
import math

class Projection:
    """ simple class to hold a geographic projection
    """
    def __init__(self, a, invf, m0, false_easting, false_northing):
        self.a = a # ellipsoid semi-major axis
        self.invf = invf # 1/f
        self.m0 = m0 # central scale factor
        self.false_easting = false_easting # false easting
        self.false_northing = false_northing # false northing

def dms2dec(dms):
    """ Convert dd.mmss to decimal degrees
    """
    if dms < 0:
        sign = -1
        dms = dms * sign
    else:
        sign =1

    degrees = int(dms)
    minutes = int((dms*100)-degrees*100)
    seconds = (((dms-degrees)*100) - minutes) * 100
    decdeg = degrees + float(minutes)/60 + float(seconds)/3600
    return decdeg*sign


def dec2dms(decdeg):
    """ Convert decimal degrees to dd.mmss
    """
    if decdeg < 0:
        sign = -1
        decdeg = decdeg * sign
    else:
        sign = 1

    degrees = int(decdeg)
    minutes = int((decdeg*60) - degrees*60)
    seconds = (((decdeg - degrees) * 3600) - minutes * 60)
    dms = degrees + float(minutes)/100 + seconds/10000
    return dms*sign

def geo_to_grid(latitude, longitude, central_meridian, proj):
    """ Convert geographic to grid coordinates using Kruger n-series equations

        See Sec 4.1.1 (pg 36) of Geocentric Datum of Australia 2020 Technical Manual 1.2

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
    """
    # Change lat/long to decimal degrees and convert to radians
    rlat = np.radians(dms2dec(latitude))
    rlong = np.radians(dms2dec(longitude))
    rcentral_meridian = np.radians(central_meridian)

    a = proj.a
    invf = proj.invf
    m0 = proj.m0
    false_easting = proj.false_easting
    false_northing = proj.false_northing

    # Calculate geometrical constants
    f = 1.0/invf
    b = a * (1-f)  # semi - minor axis

    # Eccentricity (eq 18)
    e2 = 2 * f - f * f  # = f*(2-f) = (a^2-b^2)/a^2
    e = math.sqrt(e2)

    # Compute 3rd flattening and powers (eq 19)
    n = (a - b)/(a + b)
    n2 = n * n
    n3 = n * n2
    n4 = n2 * n2
    n5 = n3 * n2
    n6 = n2 * n4
    n7 = n4 * n3
    n8 = n4 * n4

    # Rectifying Radius A (eq 20)
    A = (a/(1+n)) * (1+(1.0/4.0) * n2 + (1.0/64) * n4 + (1.0/256) * n6 + (25.0/16384) * n8)

    # Calculate conformal latitude (eq 22, 23)
    sigma = math.sinh( e*math.atanh(( e * math.tan(rlat)) / (math.sqrt( 1 + math.tan(rlat) * math.tan(rlat)))))
    conformal_lat = math.tan(rlat) * math.sqrt(1 + sigma * sigma) - sigma * \
                                                                   math.sqrt(1 + math.tan(rlat) * math.tan(rlat))

    # Compute the coefficients (eq 21)
    a2 = (1.0/2) * n - (2.0/3) * n2 + (5.0/16) * n3 + (41.0/180) * n4 - (127.0/288) * n5 + (7891.0/37800) * n6 +\
         (72161.0/387072) * n7 - (18975107.0/50803200) * n8
    a4 = (13.0/48) * n2 - (3.0/5) * n3 + (557.0/1440) * n4 + (281.0/630) * n5 - (1983433.0/1935360) * n6 + \
         (13769.0/28800) * n7 + (148003883.0/174182400) * n8
    a6 = (61.0/240) * n3 - (103.0/140) * n4 + (15061.0/26880) * n5 + (167603.0/181440) * n6 - \
         (67102379.0/29030400) * n7 + (79682431.0/79833600) * n8
    a8 = (49561.0/161280) * n4 - (179.0/168) * n5 + (6601661.0/7257600) * n6 + (97445.0/49896) * \
                                                                               n7 - (40176129013.0/7664025600) * n8
    a10 = (34729.0/80640)* n5-(3418889.0/1995840)* n6+(14644087.0/9123840)* n7+(2605413599.0/622702080) * n8
    a12 = (212378941.0/319334400) * n6 - (30705481.0/10378368) * n7 + (175214326799.0/58118860800) * n8
    a14 = (1522256789.0/1383782400) * n7 - (16759934899.0/3113510400) * n8
    a16 = (1424729850961.0/743921418240) * n8

    # Find w by subtracting central meridian from longitude (eq 24)
    w = rlong - rcentral_meridian # Converted to radians

    # Compute the gauss-Schreiber coordinates (eq 25, 26)
    u = a * math.atan( conformal_lat/math.cos(w))
    v = a * math.asinh(math.sin(w) / (math.sqrt( conformal_lat * conformal_lat + math.cos(w) * math.cos( w))))

    # Calculate partial solutions for x and y to make calcs easier (eq 27, 28)
    x1 = math.cos(2 * (u / a)) * math.sinh(2 * (v / a))
    x2 = math.cos(4 * (u / a)) * math.sinh(4 * (v / a))
    x3 = math.cos(6 * (u / a)) * math.sinh(6 * (v / a))
    x4 = math.cos(8 * (u / a)) * math.sinh(8 * (v / a))
    x5 = math.cos(10 * (u / a)) * math.sinh(10 * (v / a))
    x6 = math.cos(12 * (u / a)) * math.sinh(12 * (v / a))
    x7 = math.cos(14 * (u / a)) * math.sinh(14 * (v / a))
    x8 = math.cos(16 * (u / a)) * math.sinh(16 * (v / a))

    y1 = math.sin(2 * (u / a)) * math.cosh(2 * (v / a))
    y2 = math.sin(4 * (u / a)) * math.cosh(4 * (v / a))
    y3 = math.sin(6 * (u / a)) * math.cosh(6 * (v / a))
    y4 = math.sin(8 * (u / a)) * math.cosh(8 * (v / a))
    y5 = math.sin(10 * (u / a)) * math.cosh(10 * (v / a))
    y6 = math.sin(12 * (u / a)) * math.cosh(12 * (v / a))
    y7 = math.sin(14 * (u / a)) * math.cosh(14 * (v / a))
    y8 = math.sin(16 * (u / a)) * math.cosh(16 * (v / a))

    # Calculate partial solutions for q and p to make calcs easier
    q1 = math.sin(2 * (u / a)) * math.sinh(2 * (v / a))
    q2 = math.sin(4 * (u / a)) * math.sinh(4 * (v / a))
    q3 = math.sin(6 * (u / a)) * math.sinh(6 * (v / a))
    q4 = math.sin(8 * (u / a)) * math.sinh(8 * (v / a))
    q5 = math.sin(10 * (u / a)) * math.sinh(10 * (v / a))
    q6 = math.sin(12 * (u / a)) * math.sinh(12 * (v / a))
    q7 = math.sin(14 * (u / a)) * math.sinh(14 * (v / a))
    q8 = math.sin(16 * (u / a)) * math.sinh(16 * (v / a))

    p1 = math.cos(2 * (u / a)) * math.cosh(2 * (v / a))
    p2 = math.cos(4 * (u / a)) * math.cosh(4 * (v / a))
    p3 = math.cos(6 * (u / a)) * math.cosh(6 * (v / a))
    p4 = math.cos(8 * (u / a)) * math.cosh(8 * (v / a))
    p5 = math.cos(10 * (u / a)) * math.cosh(10 * (v / a))
    p6 = math.cos(12 * (u / a)) * math.cosh(12 * (v / a))
    p7 = math.cos(14 * (u / a)) * math.cosh(14 * (v / a))
    p8 = math.cos(16 * (u / a)) * math.cosh(16 * (v / a))

    # Calculate q and p for calculating point scale factor (eq 33, 34)
    q = - (2 * a2 * q1 + 4 * a4 * q2 + 6* a6* q3 + 8 * a8 * q4 +
           10 * a10 * q5 + 12 * a12 * q6 + 14 * a14 * q7 + 16 * a16 * q8)
    p = 1 + (2 * a2 * p1 + 4 * a4 * p2 + 6 * a6 * p3 + 8 * a8 *
             p4 + 10 * a10 * p5 + 12 * a12 * p6 + 14 * a14 * p7 + 16 * a16 * p8)


    # Calculate point scale factor m (eq 35)
    m = m0 * (A / a)*math.sqrt(q * q + p * p) * (math.sqrt(1+(math.tan(rlat)*math.tan(rlat))) *
        math.sqrt(1 - e2*(math.sin(rlat)*math.sin(rlat))))/\
        math.sqrt(conformal_lat * conformal_lat+math.cos(w)*math.cos(w))

    # compute grid convergence (eq 36)
    grid_conv = math.atan(q / p)+math.atan((conformal_lat*math.tan(w))/math.sqrt(1 + conformal_lat * conformal_lat))*\
                180/math.pi

    # Calculate transverse Mercator coordinates (eq 29, 30)
    X = A*((v / a) + a2 * x1 + a4 * x2 + a6 * x3 + a8 * x4 + a10 * x5 + a12 * x6 + a14 * x7 + a16 * x8)
    Y = A*((u / a) + a2 * y1 + a4 * y2 + a6 * y3 + a8 * y4 + a10 * y5 + a12 * y6 + a14 * y7 + a16 * y8)

    # Calculate scaled coordinates with false easting and northing
    # (eq 31, 32)
    easting = false_easting + m0 * X
    northing = false_northing + m0 * Y

    return easting, northing, m, grid_conv

def grid_to_geo(easting, northing, central_meridian, projection):
    """ Convert grid coordinates to geographic using Kruger n-series equations

        See Sec 4.1.1 (pg 36) of Geocentric Datum of Australia 2020 Technical Manual 1.2

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
    """
    a = projection.a
    invf = projection.invf
    m0 = projection.m0
    false_easting = projection.false_easting
    false_northing = projection.false_northing

    # Calculate geometrical constants
    f = 1.0/invf
    b = a * (1-f)  # semi - minor axis

    # Eccentricity (eq 18)
    e2 = 2 * f - f * f  # = f*(2-f) = (a^2-b^2)/a^2
    e = math.sqrt(e2)

    # Compute 3rd flattening and powers (eq 19)
    n = (a - b)/(a + b)
    n2 = n * n
    n3 = n * n2
    n4 = n2 * n2
    n5 = n3 * n2
    n6 = n2 * n4
    n7 = n4 * n3
    n8 = n4 * n4

    # Rectifying Radius A (eq 20)
    A = (a/(1+n)) * (1+(1.0/4.0) * n2 + (1.0/64) * n4 + (1.0/256) * n6 + (25.0/16384) * n8)

    # calculate beta values (eq 37)
    b2 = (-1.0/2) * n + (2.0/3) * n2 - (37.0/96) * n3 + (1.0/360) * n4 + (81.0/512) * n5 - (96199.0/604800) * n6 \
        + (5406467.0/38707200) * n7 - (7944359.0/67737600) * n8
    b4 = (-1.0/48) * n2 - (1.0/15) * n3 + (437.0/1440) * n4 - (46.0/105) * n5 + (1118711.0/3870720) * n6 \
        - (51841.0/1209600) * n7 - (24749483.0/348364800) * n8

    b6 = (-17.0/480) * n3 + (37.0/840) * n4 + (209.0/4480) * n5 - (5569.0/90720) * n6 - (9261899.0/58060800) * n7 \
        + (6457463.0/17740800) * n8

    b8 = (-4397.0/161280) * n4 + (11.0/504) * n5 + (830251.0/7257600) * n6 - (466511.0/2494800) * n7 \
        - (324154477.0/7664025600) * n8

    b10 = (-4583.0/161280) * n5 + (108847.0/3991680) * n6 + (8005831.0/63866880) * n7 - (22894433.0/124540416) * n8

    b12 = (-20648693.0/638668800) * n6 + (16363163.0/518918400) * n7 + (2204645983.0/12915302400) * n8

    b14 = (-219941297.0/5535129600) * n7 + (497323811.0/12454041600) * n8

    b16 = (-191773887257.0/3719607091200) * n8

    # transverse Mercator coordinates X,Y (eq 38, 39)
    X = (easting - false_easting)/m0
    Y = (northing - false_northing)/m0

    # Calculate partial solutions for q and p to make calcs easier
    x1 = math.cos(2 * (Y / A)) * math.sinh(2 * (X / A))
    x2 = math.cos(4 * (Y / A)) * math.sinh(4 * (X / A))
    x3 = math.cos(6 * (Y / A)) * math.sinh(6 * (X / A))
    x4 = math.cos(8 * (Y / A)) * math.sinh(8 * (X / A))
    x5 = math.cos(10 * (Y / A)) * math.sinh(10 * (X / A))
    x6 = math.cos(12 * (Y / A)) * math.sinh(12 * (X / A))
    x7 = math.cos(14 * (Y / A)) * math.sinh(14 * (X / A))
    x8 = math.cos(16 * (Y / A)) * math.sinh(16 * (X / A))

    y1 = math.sin(2 * (Y / A)) * math.cosh(2 * (X / A))
    y2 = math.sin(4 * (Y / A)) * math.cosh(4 * (X / A))
    y3 = math.sin(6 * (Y / A)) * math.cosh(6 * (X / A))
    y4 = math.sin(8 * (Y / A)) * math.cosh(8 * (X / A))
    y5 = math.sin(10 * (Y / A)) * math.cosh(10 * (X / A))
    y6 = math.sin(12 * (Y / A)) * math.cosh(12 * (X / A))
    y7 = math.sin(14 * (Y / A)) * math.cosh(14 * (X / A))
    y8 = math.sin(16 * (Y / A)) * math.cosh(16 * (X / A))

    # gauss-schreiber ratios (eq 42)
    nn = (X/A) + b2 * x1 + b4 * x2 + b6 * x3 + b8 * x4 + b10 * x5 + b12 * x6 + b14 * x7 + b16 * x8

    # gauss-schreiber ratios (eq 43)
    ee = (Y/A) + b2 * y1 + b4 * y2 + b6 * y3 + b8 * y4 + b10 * y5 + b12 * y6 + b14 * y7 + b16 * y8

    # (eq 44)
    t_dash = math.sin(ee)/(math.sqrt(math.sinh(nn) * math.sinh(nn) + math.cos(ee) * math.cos(ee)))

    t = t_dash

    # solve t = tan(phi) by Newton-Raphson iteration
    # initial value for t = t'
    for i in range(0,5):
        # (eq 46)
        sigma = math.sinh(e*math.atanh(e*t/math.sqrt(1+t**2)))
        # (eq 48)
        ft = t*math.sqrt(1+sigma**2)-sigma*math.sqrt(1+t**2)-t_dash
        # (eq 49)
        ft_dash = (math.sqrt(1+sigma**2)*math.sqrt(1+t**2)-sigma*t)*((1-e2)*math.sqrt(1+t**2))/(1+(1+e2)*t*t)
        # (eq 47)
        t = t - (ft/ft_dash)

    # longitude difference omega (eq 51)
    omega = math.atan(math.sinh(nn)/math.cos(ee))

    # longitude calculation (eq 52)
    calculated_longitude = np.degrees(np.radians(cm) + omega)

    # latitude calculation (eq 50)
    calculated_latitude = np.degrees(math.atan(t))

    return dec2dms(calculated_latitude), dec2dms(calculated_longitude)

def geodetic_to_cartesian(latitude, longitude, h, projection):
    """ convert geodetic coordinates to cartesian

        Parameters
        ----------
        latitude: float
            latitude in dd.mmss hp format
        longitude: float
            longitude in dd.mmss hp format
        h: float
            ellipsoidal height
        projection: Projection class
            the projection details (semi minor axis, flattening etc)

        Returns
        -------
        X: float
            X location in cartesian format
        Y: float
            Y location in cartesian format
        Z: float
            Z location in cartesian format
    """
    a = projection.a # semi-major axis
    invf = projection.invf
    f = 1/invf # flattening
    e2 = 2*f-f*f # eccentricity squared of reference ellipsoid (eq 12)

    # convert dd.mmss to radians
    latitude_radians = math.radians(dms2dec(latitude))
    longitude_radians = math.radians(dms2dec(longitude))

    # radius of curvature at a point on ellipsoid with respect to prime vertical through
    # that point (eq 11)
    v = a/(math.sqrt(1-e2*math.sin(latitude_radians)**2))

    # cartesian coordinates
    X = (v+h)*math.cos(latitude_radians)*math.cos(longitude_radians) # (eq 8)
    Y = (v+h)*math.cos(latitude_radians)*math.sin(longitude_radians) # (eq 9)
    Z = ((1-e2)*v + h)*math.sin(latitude_radians) # (eq 10)

    return X, Y, Z


def cartesian_to_geodetic(X, Y, Z, projection):
    """ convert cartesian coordinates to geodetic

        Parameters
        -------
        X: float
            X location in cartesian format
        Y: float
            Y location in cartesian format
        Z: float
            Z location in cartesian format
        projection: Projection class
            the projection details (semi minor axis, flattening etc)

        Returns
        -------
        latitude: float
            latitude in dd.mmss format
        longitude: float
            longitude in dd.mmss format
        h: float
            ellipsoid height
    """
    a = projection.a # semi-major axis
    invf = projection.invf
    f = 1/invf # flattening
    e2 = 2*f-f*f # eccentricity squared of reference ellipsoid
    subf = 1 - f

    p = math.sqrt(X**2 + Y**2) # radius of curvature at a point on ellipsoid with respect to meridian (eq 4)
    r = math.sqrt(p**2 + Z**2) # (eq 6)
    u = math.atan((Z/p)*(subf + (e2 * a/r))) # parametric latitude (eq 5)

    longitude_radians = math.atan2(Y, X) # (eq 1)

    # (eq 2)
    top = (Z*subf+e2*a*math.sin(u)**3)
    bot = (subf*(p-e2*a*math.cos(u)**3))
    latitude_radians = math.atan(top/bot)

    # (eq 3)
    part1 = p * math.cos(latitude_radians)
    part2 = Z * math.sin(latitude_radians)
    part3 = a * math.sqrt(1 - e2 * math.sin(latitude_radians)**2)
    h = part1 + part2 - part3

    return dec2dms(math.degrees(latitude_radians)), dec2dms(math.degrees(longitude_radians)), h

if __name__=="__main__":
    # Alice Springs GDA94 (ALIC)
    # Example 3.1.1 pg 26 Geocentric Datum of Australia 2020 Technical Manual Version 1.2

    a = 6378137.0 # ellipsoid semi-major axis
    invf = 298.257222101 # inverse flattening
    m0 = 0.9996 # central scale factor
    false_easting = 500000 # grid false easting
    false_northing = 10000000 # grid false northing
    projection = Projection(a, invf, m0, false_easting, false_northing)

    latitude = -23.4012446019
    longitude = 133.5307847844
    h = 603.3466 # ellipsoidal height

    # convert geodetic to cartesian
    x,y,z = geodetic_to_cartesian(latitude, longitude, h, projection)

    print('Calculated cartesian')
    print('X: {0}'.format(x))
    print('Y: {0}'.format(y))
    print('Z: {0}'.format(z))

    # convert cartesian to geodetic
    lat_calc, lon_calc, h_calc = cartesian_to_geodetic(x, y, z, projection)

    print('\nCalculated geodetic')
    print('latitude: {0}'.format(lat_calc))
    print('longitude: {0}'.format(lon_calc))
    print('h: {0}'.format(h))

    # smeaton
    # zone 55 gda 2020
    cm = 147 # zone 55
    east = 232681.853
    north = 5867898.032
    lat = -37.174973133
    lon = 143.590316715

    # convert from lat/long to mga2020
    calc_east, calc_north, psf, grid_conv = geo_to_grid(lat, lon, cm, projection)
    print(psf, grid_conv)
    print("calculated easting: {0} difference: {1:.5f}".format(calc_east, calc_east-east))
    print("calculated northing: {0} difference: {1:.5f}".format(calc_north, calc_north-north))

    # convert from mga2020 to lat/long
    calc_lat, calc_lon = grid_to_geo(calc_east, calc_north, cm, projection)

    print("calculated latitude: {0} difference: {1:.5f}".format(calc_lat, calc_lat-lat))
    print("calculated longitude: {0} difference: {1:.5f}".format(calc_lon, calc_lon-lon))
