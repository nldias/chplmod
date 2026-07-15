const TwoPi = 2*pi;
// -----------------------------------------------------------------------------
// --> ddse: ddse: declination and distance sun-earth a funcion of (year, month,
// day)
//
// based on Stull, 1995, Meteorology Today for Scientists and Engineers.
//
// 2021-05-12T10:25:54 a new star is born
// -----------------------------------------------------------------------------

private const
   NNDias: [0..1][1..12] int
   = [ [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
       [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31] ];
inline proc Bissexto(Ano: int ): bool {
   return ( ( ( Ano % 4 == 0 ) && ( Ano % 100 != 0 ) ) || ( Ano % 400 == 0 ) ) ;
}

inline proc DataCor(
   Bis: bool,
   Mes: int,
   Dia: int): int {
   var DiAno = Dia ;
   var i = 1 ;
   while i < Mes do {
      DiAno += NNDias[Bis][i];
      i += 1;
   }
   return DiAno;
}



proc ddse(yea: int, mon: int, day: int): (real,real) {
// -----------------------------------------------------------------------------
// O Ano é bissexto?
// -----------------------------------------------------------------------------
   var da,dr,dp: int;
   const e = 0.0167;
   if Bissexto(yea) then {
      da = 366;
      dr = 174;
   }
   else {
      da = 365;
      dr = 173;
   }
   dp = 3;
// -----------------------------------------------------------------------------
// O dia juliano é ...
// -----------------------------------------------------------------------------
   var d = DataCor(Bissexto(yea),mon,day);
// -----------------------------------------------------------------------------
// Sun's declination
// -----------------------------------------------------------------------------
   var delta = 0.408*cos(2*pi*(d-dr)/da);
// -----------------------------------------------------------------------------
// distance Sun-Earth in the form (r/a) where a is the length of the largest
// semi-axis of the Earth's orbit, i.e.: the equivalent to one astronomical
// unit, and rr is the Sun-Earth distance
// -----------------------------------------------------------------------------
   var rr = (1 - e**2)/(1 + e*cos(2*pi*(d-dp)/da));
   return (delta,rr);
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
proc rsds(lat,rr,delta): (real,real) {
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
      
