package gda2020

import (
    "testing"
    "math"
    "fmt"
)

// Test dms2dec
func TestDms2Dec1(t *testing.T){
    dms := 20.30
    expected := 20.50
    result := Dms2dec(dms)
    tolerance := 0.001
    if math.Abs(result - expected) > tolerance{
        t.Error("Expected 20.50 got: ", result)
    }
}

// Test dms2dec
func TestDms2Dec2(t *testing.T){
    dms := 20.3124
    expected := 20.5233
    result := Dms2dec(dms)
    tolerance := 0.001
    if math.Abs(result - expected) > tolerance{
        t.Error("Expected 20.5233 got: ", result)
    }
}

// Test radians conversion
func TestRadians(t *testing.T){
    dec := 30.55892356
    expected := 0.5333538320
    result := Radians(dec)
    tolerance := 0.000000001
    if math.Abs(result - expected) > tolerance{
        t.Error("Expected  0.5333538320 got: ", result)
    }
}

// Alice Springs GDA94 (ALIC)
// Example 3.1.1 pg 26 Geocentric Datum of Australia 2020 Technical Manual Version 1.2
// test GeodeticToCartesian and CartesianToGeodetic conversions
func TestAlice(t *testing.T){
    a := 6378137.0 // ellipsoid semi-major axis
    invf := 298.257222101 // inverse flattening
    m0 := 0.9996 // central scale factor
    var  false_easting float64 = 500000 // grid false easting
    var false_northing float64 = 10000000 // grid false northing
    projection := Projection{A: a, Invf:invf, M0:m0, FalseEasting: false_easting, FalseNorthing: false_northing}

    latitude := -23.4012446019
    longitude := 133.5307847844
    h := 603.3466 // ellipsoidal height

    tolerance := 0.0000001

    // convert geodetic to cartesian
    x, y, z := GeodeticToCartesian(latitude, longitude, h, projection)

    // convert cartesian to geodetic
    lat_calc, lon_calc, h_calc := CartesianToGeodetic(x, y, z, projection)

    if math.Abs(lat_calc - latitude) > tolerance{
        t.Error("Latitude calc incorrect", lat_calc)
    }

    if math.Abs(lon_calc - longitude) > tolerance{
        t.Error("Longitude calc incorrect", lon_calc)
    }

    if math.Abs(h_calc - h) > tolerance{
        t.Error("Ellipsoid calc incorrect", h_calc)
    }
}

func TestSmeaton(t *testing.T){
    a := 6378137.0 // ellipsoid semi-major axis
    invf := 298.257222101 // inverse flattening
    m0 := 0.9996 // central scale factor
    var  false_easting float64 = 500000 // grid false easting
    var false_northing float64 = 10000000 // grid false northing
    
    projection := Projection{A: a, Invf:invf, M0:m0, FalseEasting: false_easting, FalseNorthing: false_northing}
    // smeaton
    // zone 55 gda 2020
    var cm float64 = 147 // zone 55
    east := 232681.853
    north := 5867898.032
    lat := -37.174973133
    lon := 143.590316715

    tolerance := 0.001

    // convert from lat/long to mga2020
    east_calc, north_calc, psf, grid_conv := GeoToGrid(lat, lon, cm, projection)

    if math.Abs(east_calc - east) > tolerance{
        t.Error("Easting calc incorrect", east_calc)
    }

    if math.Abs(north_calc - north) > tolerance{
        t.Error("Northing calc incorrect", north_calc)
    }

    fmt.Printf("%.6f\n", psf)
    fmt.Printf("%.6f\n", grid_conv)
}

func TestFlinders(t *testing.T){
    a := 6378137.0 // ellipsoid semi-major axis
    invf := 298.257222101 // inverse flattening
    m0 := 0.9996 // central scale factor
    var  false_easting float64 = 500000 // grid false easting
    var false_northing float64 = 10000000 // grid false northing
    
    projection := Projection{A: a, Invf:invf, M0:m0, FalseEasting: false_easting, FalseNorthing: false_northing}
    // smeaton
    // zone 55 gda 2020
    var cm float64 = 141 // zone 55
    east := 758173.797
    north := 5828674.337
    lat := -37.391015619
    lon := 143.553538390

    tolerance := 0.001

    // convert from lat/long to mga2020
    east_calc, north_calc, psf, grid_conv := GeoToGrid(lat, lon, cm, projection)

    if math.Abs(east_calc - east) > tolerance{
        t.Error("Easting calc incorrect", east_calc)
    }

    if math.Abs(north_calc - north) > tolerance{
        t.Error("Northing calc incorrect", north_calc)
    }

    fmt.Printf("%.6f\n", psf)
    fmt.Printf("%.6f\n", grid_conv)
}