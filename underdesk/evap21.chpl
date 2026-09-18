// -----------------------------------------------------------------------------
// evap21 is the legacy code of evap: increase in Ra due to cloudiness
// is hard-wired
//
// 2022-07-07T09:22:28
// -----------------------------------------------------------------------------
const sigma = 5.670374419e-8;  // Stefan-Boltzmann constant, W/m2/K4
const P0 = 101325.0;           // atm pressure at sea level
use atmgas, ssr;
use IO;
// -----------------------------------------------------------------------------
// static local variables
// -----------------------------------------------------------------------------
private var
   aP = 0.25,      // constant, Prescott's equation
   bP = 0.50,      // constant, Prescott's equation
   aB = 0.625,     // constant, Brutsaert's longwave radiation equation
   bB = 0.131,     // constant, Brutsaert's longwave radiation equation
   aV,             // constante, coef. transf. vapor   
   bV,             // constante, coef. transf. vapor   
   eps             // emissividade                     
   : real;         
var Ph: real;      // atmospheric pressure at altitude of a standard atmosphere
var Th: real;      // temperature at altitude of a standard atmosphere
// -----------------------------------------------------------------------------
// --> Inipar: inicializa o procedimento, pegando a latitude e a altitude da
// estacao e calculando a pressao atmosferica em uma atmosfera padrao
// -----------------------------------------------------------------------------
proc IniPar(
   Altitude: real,        // altitude 
   Emissividade: real     // emissividade
) {
 eps = Emissividade ;
 Th    = 288.15 - 0.0065*Altitude;
 Ph    = P0*(Th/288.15)**5.256;
}
// -----------------------------------------------------------------------------
// --> Prescott: define os parametros a,b da equacao de Prescott:
//
// Rs = Rsea(a + b S)
// -----------------------------------------------------------------------------
proc Prescott(
   a: real,
   b: real
) {
 aP = a;
 bP = b;
}
// -----------------------------------------------------------------------------
// --> Brutsaert: define the parameters of Brutsaert's equation
//
// As a note, Duarte et al. use epsac = (0.625)*((ea/Ta)**(0.131), so a = 0.625
// -----------------------------------------------------------------------------
proc Brutsaert(
   a: real,
   b: real
) {
   aB = a;
   bB = b;
}
   
// -----------------------------------------------------------------------------
// --> Wind: define os parametros a,b do coeficiente de transferencia de massa:
//   
// f(u) = a + b u
// -----------------------------------------------------------------------------
proc Wind(
   a: real,
   b: real
)  {
 aV = a;
 bV = b;
}
// -----------------------------------------------------------------------------
// --> Bissexto: Diz se um ano e' bissexto ou nao
// -----------------------------------------------------------------------------
proc Bissexto(Ano: int ): bool {
   return ( ( Ano % 4 == 0 ) && ( Ano % 100 != 0 ) || ( Ano % 400 == 0 ) ) ;
}
// -----------------------------------------------------------------------------
// --> LakeAlbedo: calculates the albedo of a lake using Morton's 1983 CRLE
// equations. Note: Z = lat - delta, where latitude is the latitude in radians
// and delta is the Sun's declination
// -----------------------------------------------------------------------------
// proc LakeAlbedo(
//    in Z: real,          // zenith angle at solar noon
//    const in S: real           // sunshine duration
// ): real {
//    const b = 0.6875;
//    const az = 0.05;          // snow-free clear-sky albedo
//    var num, den;              // Arnfields' generalized integrals
//    Z = abs(Z);
//    var aux = sin(Z) + b*cos(Z);
//    num = (exp(pi*b/2) - exp(b*Z)*aux)/(1+b**2);
//    den = 1 - sin(Z);
//    var a0 = az*num/den;
//    var alb = a0*(S + (1.0 - S)*(1 - Z/330.0));
//    return alb;
// }
// -----------------------------------------------------------------------------
// --> WaterAlbedo: Implementation of table 5 of Cogley, J. G. The
// albedo of water as a function of latitude Monthly Weather Review,
// 1979, 107, 775-781
//
// Data are interpolated between latitudes and the 15th of each month
// -----------------------------------------------------------------------------
use DateTime;
private const Lats = [0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0]*(pi/180);
private var Galb: [0..9,1..12] real;    // albedos by latitudes and months
Galb[9,1..12] = [ NAN, NAN,30.1,29.3,17.1,14.8,16.0,24.6,34.2, NAN, NAN, NAN]/100.0;
Galb[8,1..12] = [ NAN,30.1,31.9,22.5,16.0,13.1,14.5,20.6,29.4,30.5, NAN, NAN]/100.0;
Galb[7,1..12] = [30.1,33.8,22.9,14.8,11.6,11.2,11.4,13.4,20.2,31.3,30.1, NAN]/100.0;
Galb[6,1..12] = [33.9,24.0,15.5,10.5, 8.8, 8.4, 8.6, 9.8,13.6,21.6,32.1,35.5]/100.0;
Galb[5,1..12] = [22.0,16.1,10.8, 8.4, 7.5, 7.3, 7.4, 8.0, 9.9,14.4,21.0,24.1]/100.0;
Galb[4,1..12] = [14.5,11.1, 8.5, 7.3, 6.8, 6.7, 6.8, 7.1, 8.0,10.3,13.8,16.1]/100.0;
Galb[3,1..12] = [10.3, 8.6, 7.3, 6.7, 6.5, 6.4, 6.4, 6.6, 7.1, 8.2,10.0,11.1]/100.0;
Galb[2,1..12] = [ 8.3, 7.4, 6.7, 6.4, 6.3, 6.3, 6.3, 6.4, 6.6, 7.2, 8.1, 8.7]/100.0;      
Galb[1,1..12] = [ 7.2, 6.7, 6.4, 6.3, 6.4, 6.4, 6.4, 6.3, 6.3, 6.6, 7.1, 7.4]/100.0;      
Galb[0,1..12] = [ 6.6, 6.4, 6.3, 6.4, 6.6, 6.8, 6.7, 6.4, 6.3, 6.4, 6.6, 6.8]/100.0;      
proc WaterAlbedo(
   const in year: int,
   const in month: int,
   const in day: int,
   in lat: real
   ): real {
   var today = new date(year,month,day);
   if lat < 0.0 then {
      var dt = new timedelta(183);
      today += dt;                 // shift by 6 months approx 183 days
   }
// -----------------------------------------------------------------------------
// how do I interpolate the latitude "vertically"? Can Chapel do the magic?
// -----------------------------------------------------------------------------
   lat = abs(lat);
   var (il,iu) = indx_interp(lat,Lats);
   var deltal = (lat-Lats[il])/(Lats[iu] - Lats[il]);
// -----------------------------------------------------------------------------
// find months to interpolate from
// -----------------------------------------------------------------------------
   var amon,bmon,aux: int;
   var ayear=today.year;
   var byear=today.year;
   if today.day <= 15 then {
      amon = today.month-1;
      bmon = today.month;
      if amon < 1 then {
         amon = 12;
         ayear -= 1;
      } 
   }
   else {
      amon = today.month;
      bmon = today.month+1;
      if bmon > 12 then {
         bmon = 1;
         byear += 1;
      }
   }
// -----------------------------------------------------------------------------
// now find the distance between current day and 15th of amon
// -----------------------------------------------------------------------------
   var lastd = new date(ayear,amon,15);
   var nextd = new date(byear,bmon,15);
   var diffl = today - lastd;
   var difmo = nextd - lastd;
   var deltad = (diffl.days):real/(difmo.days):real;
// -----------------------------------------------------------------------------
// now interpolate vertically (deltal) and horizontally (deltad) in that order
// -----------------------------------------------------------------------------
   var alba = (1-deltal)*Galb[il,amon] + deltal*Galb[iu,amon];
   var albb = (1-deltal)*Galb[il,bmon] + deltal*Galb[iu,bmon];
   var alb =  (1-deltad)*alba + deltad*albb;
   return alb;
}
// -----------------------------------------------------------------------------
// --> Sprescott: S from Prescott parameters; returns S = n/N
// -----------------------------------------------------------------------------
inline proc SPrescott(
   const in Rsea: real,       // extra-terrestrial solar radiation (W/m2)
   const in Rs: real          // solar radiation (W/m^2)
   ): real {
   assert(Rs > 0.0);
   var S = (Rs/Rsea - aP)/bP;
   S = max(0.0,S);
   S = min(1.0,S);
   return S;
}

// -----------------------------------------------------------------------------
// --> Radiation: calcula a radiacao solar, a radiacao atmosferica e a radiacao
// emitida pela superficie a partir dos dados sobre dia, mes e ano, numero de
// horas de brilho de sol, pressao de vapor, temperatura do ar e temperatura da
// superficie
//
// The atmospheric irradiance components are calculated after Duarte,
// H. F.; Dias, N. L. & Maggiotto, S. R. Assessing daytime downward
// longwave radiation estimates for clear and cloudy skies in Southern
// Brazil Agric For Meteorol, 2006, 139, 171-181
// -----------------------------------------------------------------------------
proc Radiation(
   const in alb: real,        // albedo
   const in ea: real,         // vapor pressure of the air
   const in Ta: real,         // air termperature (Kelvin)
   const in T0: real,         // surface temperature (Kelvin)
   const in S: real,          // relative sunshine duration
   const in Rs: real,         // solar radiation (W/m^2)
   out Ra: real,              // atmospheric radiation (W/m^2)
   out Re: real,              // emitted radiation     (W/m^2)
   out Rl0: real              // net radiation     (W/m^2)
) {
   var epsAC = aB*((ea/Ta)**bB);   // clear-sky atmospheric emissivity
   Ra = (epsAC)*(sigma)*Ta**4;     // clear-sky atmospheric radiation
// -----------------------------------------------------------------------------
// aumento da radiacao atmosferica devido aa presenca de nuvens segundo a
// equacao (22) de Duarte et al (2006)
//
//  checking going back to Bolz
// -----------------------------------------------------------------------------
   var c = 1.0 - S;
   Ra = Ra*(1.0 - c**0.671) + 0.990*(c**0.671)*sigma*Ta**4;
//   Ra = Ra*(1 + 0.22*c**2);
   Re = eps*sigma*T0**4;                // radiacao emitida
   Rl0 = Rs*(1.0-alb) + eps*Ra - Re;    // radiação líquida
}
// -----------------------------------------------------------------------------
// --> Epen: evaporacao segundo Penman
// -----------------------------------------------------------------------------
proc EPen(
   const in Ta: real,    // temperatura do ar Kelvin 
   const in Rl: real,    // net available energy
   const in u: real,     // velocidade do vento (m/s)
   const in ea: real,    // pressao de vapor (Pa)    
   out Ep,               // evaporacao kg/m^2/s      
   out Tp                // temperatura de equilibrio
) {
// -----------------------------------------------------------------------------
// calcula gamma
// -----------------------------------------------------------------------------
   var q = sphum(ea,Ph);
   var cp = spheat(q);
   var L = latent(Ta);
   var gamma = cp*Ph/(0.622*L);
// -----------------------------------------------------------------------------
// pressao de saturacao de vapor d'agua aa temp. ar e sua inclinacao
// -----------------------------------------------------------------------------
   var  
      esa,                  // pressao saturacao vapor temp. ar
      da,                   // inclinacao da curva es x T aa temp. ar
      ft:                   // coef. transf. vapor
      real;
   (esa,da) = svpd(Ta);
// -----------------------------------------------------------------------------
// coeficiente de transferencia de vapor
// -----------------------------------------------------------------------------
   ft = aV  +  bV*u;
// -----------------------------------------------------------------------------
// evaporacao
// -----------------------------------------------------------------------------
   Ep = da/(da + gamma) * Rl / L  
         + gamma/(da + gamma) * ft * ( esa - ea ) ;
// -----------------------------------------------------------------------------
// temperatura de equilibrio
// -----------------------------------------------------------------------------
   Tp = Ta + (Rl - L*ft*(esa - ea))/(L*ft*(da + gamma)) ;
}
// -----------------------------------------------------------------------------
// --> EKeP: evaporacao com a equacao de Kohler e Parmele
// -----------------------------------------------------------------------------
proc EKeP(
   const in Ta: real,    // temp. ar (K)
   const in Rla: real,   // radiacao liquida temperatura do ar (W/m^2)
   const in u: real,     // velocidade do vento (m/s)
   const in ea: real,    // pressao de vapor do ar (Pa)
   out Ek: real,         // evaporacao (kg/m^2/s)
   out Tk: real          // temperatura de equilibrio
) {
// -----------------------------------------------------------------------------
// calcula gamma
// -----------------------------------------------------------------------------
   var q = sphum(ea,Ph);
   var cp = spheat(q);
   var L = latent(Ta);
   var gamma = cp*Ph/(0.622*L);
// -----------------------------------------------------------------------------
// pressao de saturacao de vapor aa temperatura do ar e sua inclinacao
// -----------------------------------------------------------------------------
   var esa: real;          // pressao saturacao vapor temp. ar
   var da: real;           // inclinacao da curva es x T aa temp. ar
   var ft: real;           // coef. transf. vapor
   var lamb: real;         // coeficiente de transferencia de calor
   (esa,da) = svpd(Ta);
// -----------------------------------------------------------------------------
// coeficiente de transferencia de massa
// -----------------------------------------------------------------------------
   ft = aV  +  bV*u;
// -----------------------------------------------------------------------------
// coeficiente de transferencia de calor
// -----------------------------------------------------------------------------
   lamb = gamma + 
      4.0 * eps * sigma * ( Ta * Ta * Ta ) / ( L * ft ) ;
// -----------------------------------------------------------------------------
// evaporacao
// -----------------------------------------------------------------------------
   Ek = da/(da + lamb) * Rla / L  
      + lamb/(da + lamb) * ft * ( esa - ea ) ;
// -----------------------------------------------------------------------------
// temperatura de equilibrio
// -----------------------------------------------------------------------------
   Tk = Ta + (Rla - L*ft*(esa - ea))/(L*ft*(da + lamb)) ;
}
// -----------------------------------------------------------------------------
// --> EPaT: evaporacao pela formula de Priestley-Taylor
// -----------------------------------------------------------------------------
proc EPT(
   const in alfa: real,       // cte de Priestley-Taylor, default = 1.26
   const in T: real,          // temperatura, K
   const in Rl: real,         // radiacao liquida, W/m^2
   out Es: real               // evaporacao, kg/m^2/s
) {
/* ------------------------------------------------------------------------
   censura o parametro alfa 
   ------------------------------------------------------------------------ */
   var (e,d) = svpd(T);
// -----------------------------------------------------------------------------
// calcula gamma
// -----------------------------------------------------------------------------
   var q = sphum(e,Ph);
   var cp = spheat(q);
   var L = latent(T);
   var gamma = cp*Ph/(0.622*L);
// -----------------------------------------------------------------------------
// evaporacao de Priestley-Taylor 
// -----------------------------------------------------------------------------
   Es = alfa * ( d / ( d + gamma ) ) * (Rl/L);
}
