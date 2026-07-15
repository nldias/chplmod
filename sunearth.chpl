use Math only sin,cos,pi,tan,acos,asin;
const TwoPi = 2*pi;
// -----------------------------------------------------------------------------
// --> ddse: ddse: declination and distance sun-earth a funcion of (year, month,
// day)
//
// based on
//
// Van Flandern, T. C. and Pulkkinen, K. F. (1979) "Low Precision Formulae for
// Planetary Positions" - The Astronomical Journal Supplement Series,
// 41,391:411.
//
// 2019-05-04T13:46:15 fumbling around
//
// 2020-03-30T15:38:03 reverting to a faithful reproduction of my dear ddstx in
//                     weather.c
//
// 2021-04-06T13:56:56 translating from Python to Chapel
// -----------------------------------------------------------------------------
proc ddse(yea: int, mon: int, day: int): (real,real) {
// -----------------------------------------------------------------------------
// At GMT noon: this is done with purely integer arithmetic
// -----------------------------------------------------------------------------
   var JD = 367 * yea - ( 7 * (yea + (mon + 9) / 12 ) / 4 ) 
          + ( 275 * mon / 9 + day ) + 1721014;
   var tee = (JD - 2451545): real; // tee == thousands of Julian years from 2000
   var TC = tee/36525.0 + 1.0;     // TC == hundreds of Julian years from 1900
   var LS = 0.779072 + 0.00273790931 * tee;  // mean longitude, Sun
   var GS = 0.993126 + 0.00273777850 * tee;  // mean anomaly, Sun
   var G5 = 0.056531 + 0.00023080893 * tee;  // mean anomaly, Jupiter
   var OM = 0.347343 - 0.00014709391 * tee;  // long of lunar ascending mode
// -----------------------------------------------------------------------------
// extracts fractional part
// -----------------------------------------------------------------------------
   LS = LS - trunc(LS);
   GS = GS - trunc(GS);
   G5 = G5 - trunc(G5);
   OM = OM - trunc(OM);
// -----------------------------------------------------------------------------
// converts to radians
// -----------------------------------------------------------------------------
   LS = LS * TwoPi;
   GS = GS * TwoPi;
   G5 = G5 * TwoPi;
   OM = OM * TwoPi;
// -----------------------------------------------------------------------------
// obtains VS
// -----------------------------------------------------------------------------
   var VS = + 0.39785 * sin( LS )
        - 0.01000 * sin( LS - GS )
        + 0.00333 * sin( LS + GS )
        - 0.00021 * TC * sin( LS ) 
        + 0.00004 * sin( LS + 2.0 * GS ) 
        - 0.00004 * cos( LS ) 
        - 0.00004 * sin( OM - LS ) 
        + 0.00003 * TC * sin( LS - GS );
// -----------------------------------------------------------------------------
// obtains US
// -----------------------------------------------------------------------------
   var US = + 1.0 
        - 0.03349 * cos( GS ) 
        - 0.00014 * cos( 2.0 * GS ) 
        + 0.00008 * TC * cos( GS ) 
        - 0.00003 * sin( GS - G5 );
// -----------------------------------------------------------------------------
// Sun's declination
// -----------------------------------------------------------------------------
   var delta = asin( VS / sqrt(US) ) ;
// -----------------------------------------------------------------------------
// distance Sun-Earth in the form (r/a) where a is the length of the largest
// semi-axis of the Earth's orbit, i.e.: the equivalent to one astronomical
// unit, and rr is the Sun-Earth distance
// -----------------------------------------------------------------------------
   var rr = 1.00021 * sqrt( US );
   return (delta,rr);
}

// -----------------------------------------------------------------------------
// -->sunman: only the sun mean anomaly a funcion of (year, month,day)
//
// based on
//
// Van Flandern, T. C. and Pulkkinen, K. F. (1979) "Low Precision Formulae for
// Planetary Positions" - The Astronomical Journal Supplement Series,
// 41,391:411.
//
// 2021-04-06T14:19:53 Translating to Chahpel
// -----------------------------------------------------------------------------
proc sunman(yea,mon,day,sec=0.0): real {
// -----------------------------------------------------------------------------
// At GMT noon: this is done with purely integer arithmetic
// -----------------------------------------------------------------------------
   var IJD = 367 * yea 
      - 7*(yea + (mon + 9)/12)/4 
      - 3*((yea + (mon - 9)/7)/100 + 1)/4 
      + 275*mon/9 + day + 1721029;
// -----------------------------------------------------------------------------
// trying to get more accuracy by calculating seconds
// -----------------------------------------------------------------------------
   var JD = IJD + sec/86400.0;
// -----------------------------------------------------------------------------
// Obtains tee,
// -----------------------------------------------------------------------------
   var tee = JD - 2451545.0;
// -----------------------------------------------------------------------------
// other variables
// -----------------------------------------------------------------------------
   var GS = 0.993126 + 0.00273777850 * tee;  // mean anomaly, Sun
   GS = GS - trunc(GS);                  // fractional part
// -----------------------------------------------------------------------------
// converts to radians
// -----------------------------------------------------------------------------
   GS = GS * TwoPi ;
   return(GS);
}

// -----------------------------------------------------------------------------
// --> rsds: extra-atmospheric solar radiation and maximum number of hours of
// sunshine as a function of latitude, the Sun-Earth distance and the Sun's
// declination, where:
//
// lat   -- latitude (radians)
// rr    -- Sun-Earth distance in A.U.  
// delta -- Sun's declination (radians)
//
// returns
//
// rsea  -- extra-atmospheric solar radiation ( W/m2 ) 
// dsmax -- maximum sunshine duration (hours)
//
// References: Sellers, W.D. (1965) " Physical Climatology " . The University of
// Chicago Press / Chicago & London
//
// 2020-03-30T16:15:43 since my MSc thesis the solar constant has changed! Now
// it is 1361.5 (according to Wikipedia)
// -----------------------------------------------------------------------------
proc rsds(
   const in lat:real,
   const in rr: real,
   const in delta: real
   ): (real,real) {
   const rs0 = 1361.5;    // solar constant 
// -----------------------------------------------------------------------------
// calculates H = half-duration of a day in radians 
// -----------------------------------------------------------------------------
   var H = acos( - tan(lat) * tan(delta) );
// ----------------------------------------------------------------------------- 
// calculates dsmax in hours
// -----------------------------------------------------------------------------
   var dsmax = 24.0 * H / pi;
// -----------------------------------------------------------------------------
// calculates rsea 
// -----------------------------------------------------------------------------
   var rsea = ( rs0 / pi ) * (1.0/rr**2) *
      ( H * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(H) );
   return (rsea,dsmax);
}

// -----------------------------------------------------------------------------
// --> rsdsZ: extra-atmospheric solar radiation and maximum number of hours of
// sunshine as a function of latitude, the Sun-Earth distance and the Sun's
// declination, where:
//
// lat   -- latitude (radians)
// rr    -- Sun-Earth distance in A.U.  
// delta -- Sun's declination (radians)
//
// returns
//
// rsea  -- extra-atmospheric solar radiation ( W/m2 ) 
// dsmax -- maximum sunshine duration (hours)
// Z     -- zenith angle of the sun at solar noon
//
// References: Sellers, W.D. (1965) " Physical Climatology " . The University of
// Chicago Press / Chicago & London
//
// 2020-03-30T16:15:43 since my MSc thesis the solar constant has changed! Now
// it is 1361.5 (according to Wikipedia)
// -----------------------------------------------------------------------------
proc rsdsZ(lat,rr,delta): (real,real,real) {
   const rs0 = 1361.5;    // solar constant 
// -----------------------------------------------------------------------------
// calculates H = half-duration of a day in radians 
// -----------------------------------------------------------------------------
   var H = acos( - tan(lat) * tan(delta) );
// ----------------------------------------------------------------------------- 
// calculates dsmax in hours
// -----------------------------------------------------------------------------
   var dsmax = 24.0 * H / pi;
// -----------------------------------------------------------------------------
// calculates rsea 
// -----------------------------------------------------------------------------
   var rsea = ( rs0 / pi ) * (1.0/rr**2) *
      ( H * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(H) );
// -----------------------------------------------------------------------------
// Zenith angle at solar noon
// -----------------------------------------------------------------------------
   var Z = lat - delta;
//   if lat < 0 then Z = -Z;
   return (rsea,dsmax,Z);
}

// -----------------------------------------------------------------------------
// --> rseat: instantaneous extra-atmospheric solar radiation and maximum number
// of hours of sunshine as a function of latitude, the Sun-Earth distance and
// the Sun's declination, where:
//
// tnoon -- time of solar noon, in decimal hours
// t     -- time, in decimal hours (9.47, etc.)
// lat   -- latitude (radians)
// rr    -- Sun-Earth distance in A.U.  
// delta -- Sun's declination (radians)
//
// returns
//
//    rsea  -- extra-atmospheric solar radiation at time t ( W/m2 ) 
//    cosZ  -- cosine of Zenith angle of the Sun
//
// References: Sellers, W.D. (1965) " Physical Climatology " . The University of
// Chicago Press / Chicago & London 2020-08-13T15:01:51 since my MSc thesis the
// solar constant has changed! Now it is 1361.5 (according to Wikipedia)
//
// returns instantaneous values, adequate for hourly or sub-hourly time
// intervals
// -----------------------------------------------------------------------------
proc rseat(
   tnoon: real,
   t: real,
   lat: real,
   rr: real,
   delta: real ): (real,real) {
   const rs0 = 1361.5;        // solar constant
// -----------------------------------------------------------------------------
// calculates H = half-duration of a day in radians 
// -----------------------------------------------------------------------------
   var H = acos( - tan(lat) * tan(delta) );
// -----------------------------------------------------------------------------
// the first thing we need to calculate is the hour angle!
// -----------------------------------------------------------------------------
   var h = (pi/12.0)*(tnoon - t);
   if ( h < -H)  || (h > H) then {   // control nightime periods
      return(0.0,0.0);
   }                               // daytime!
// -----------------------------------------------------------------------------
// cosine of Zenith angle
// -----------------------------------------------------------------------------
   var cosZ = (sin(delta)*sin(lat)+cos(delta)*cos(lat)*cos(h));
// -----------------------------------------------------------------------------
// extra-terrestrial radiation at time t
// -----------------------------------------------------------------------------
   var rsea = rs0*cosZ/(rr**2); 
   return (rsea,cosZ);
}
// -----------------------------------------------------------------------------
// --> rseac: calculates the instantaneous clear-sky solar radiation, following
//
// - Crawford, T. M. & Duchon, C. E. An improved parameterization for estimating
//   effective atmospheric emissivity for use in calculating daytime downwelling
//   longwave radiation J Appl Meteorol, 1999, 38, 474-480
//
// - Meyers, T. & Dale, R. Predicting daily insolation with hourly cloud height
//   and coverage Journal of climate and applied meteorology, 1983, 22, 537-545
//
// - Smith, W. L. Note on the relationship between total precipitable water and
//   surface dew point Journal of Applied Meteorology, 1966, 5, 726-727
//
// as a function of year, month, day, latitude (radians), local solar
// noon time and local time, Smith's lambda (Table 1), atmospheric
// pressure (Pa) and Dew Point temperature (K)
//
// Nelson Luís Dias
// 2020-08-13T16:05:54
// 2020-08-13T16:05:57
// -----------------------------------------------------------------------------
proc rseac(year,mon,day,lat,tnoon,t,slambda,pa,Td): real {
   var (delta,rr) = ddse(year,mon,day); // sun's decl and sun-earth distance
// -----------------------------------------------------------------------------   
// instantaneous extra-terrestrial solar radiation and Zenith angle
// -----------------------------------------------------------------------------
   var (rsea,cosZ) = rseat(tnoon,t,lat,rr,delta);
   if ( rsea == 0.0 ) then {  // nighttime? return
      return 0.0;
   }
// -----------------------------------------------------------------------------
// optical air mass at 101325 Pa
// -----------------------------------------------------------------------------
   var m = 35*(1224*cosZ**2 + 1)**(-0.5);
   pa /= 1000.0;               // use pressure in kPa in formulas
   var TrTg = 1.021 - 0.084*(m*(949*pa*1.0e-5+0.051))**(0.5);
// -----------------------------------------------------------------------------
// Dew point temperature
// -----------------------------------------------------------------------------
   Td = Td - 273.15;           // from Kelvin to Celsius
   Td = Td*9.0/5.0 + 32.0;     // from Celsius to Fahrenheit
   var u = exp(0.1133 - log(slambda + 1.0) + 0.0393*Td);
   var Tw = 1 - 0.077*(u*m)**0.3;
   var Ta = 0.935**m;
   return rsea*(TrTg*Tw*Ta);
}
// -----------------------------------------------------------------------------
// --> EoT: equation of time
// -----------------------------------------------------------------------------
inline proc EoT(
   const in jd: int           // the julian day 1 <= jd <= 365
   ): real {
   var B = (360.0/365)*(jd - 81);
   var aux = 9.87*sin(2*B) - 7.53*cos(2*B) - 1.5*sin(B);
   return aux;
}
// -----------------------------------------------------------------------------
// --> Loc2GMT: convert local solar time hh.mm.sec to GMT (include
// fractions of seconds)
// -----------------------------------------------------------------------------
inline proc Loc2GMT(
   const in local: time,
   const in long: real
   ): time {
   return;
}
   


   
