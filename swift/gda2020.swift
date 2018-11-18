/**
 Convert dd.mmss to decimal degrees
 
 - Parameter dms: the input in dd.mmss
 
 - Returns: The decimal degrees for the input
 */
func dms2dec(dms: Float64)->Float64{
    var sign: Float64
    var dms = dms
    
    if dms < 0{
        sign = -1
        dms = dms * sign
    }
    else{
        sign = 1
    }
    
    let degrees = Int(dms)
    let minutes = Int((dms*100) - Float64(degrees*100))
    let seconds = (((dms-Float64(degrees))*100) - Float64(minutes)) * 100
    let decDeg = Float64(degrees) + Float64(minutes)/60 + Float64(seconds)/3600
    return decDeg*sign
}

/**
 Convert dd.mmss to decimal degrees
 
 - Parameter decdeg: the input in decimal degrees
 
 - Returns: The dd.mmss for the input
 */
func dec2dms(decdeg: Float64)->Float64{
    var sign: Float64
    var decdeg = decdeg
    
    if decdeg < 0{
        sign = -1
        decdeg = decdeg * sign
    }
    else{
        sign = 1
    }
    
    let degrees = Int(decdeg)
    let minutes = Int((decdeg*60) - Float64(degrees)*60)
    let seconds = (((decdeg - Float64(degrees)) * 3600) - Float64(minutes) * 60)
    let dms = Float64(degrees) + Float64(minutes)/100 + Float64(seconds)/10000
    return dms*sign
}