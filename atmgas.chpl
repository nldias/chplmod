// =============================================================================
// ==> atmgas: various properties of atmospheric gases
//
// Nelson Luís Dias (nelsonluisdias@gmail.com)
// 2024-03-09T10:51:58
// =============================================================================
use Math only exp;
// -----------------------------------------------------------------------------
// a table of molar masses in kg/mol
// -----------------------------------------------------------------------------
const Mmol = [
   "DRY" =>  28.966413e-3,
   "N2"  =>  28.013400e-3,
   "O2"  =>  31.99880e-3,
   "AR"  =>  39.94800e-3,   
   "H2O" =>  18.0153e-3,
   "CO2" =>  44.00950e-3,
   "NE"  =>  20.1797e-3,
   "HE"  =>   4.002602e-3,
   "KR"  =>  83.798e-3,
   "XE"  => 131.293e-3,
   "CH4" =>  16.04250e-3,
   "H2"  =>   2.01588e-3,
   "O3"  =>  47.99820e-3,
   "N2O" =>  44.01280e-3
   ];
// -----------------------------------------------------------------------------
// the universal gas constant
// -----------------------------------------------------------------------------
const Ru  = 8.31446261815;    // J/K/mol universal gas constant
// -----------------------------------------------------------------------------
// the individual gas constants are handy
// -----------------------------------------------------------------------------
const Rgas = Ru/Mmol;
// -----------------------------------------------------------------------------
// a table of useful abbreviations
// ------------------------------------------------------------------------------
const Rd  = Rgas["DRY"];      // gas constant for dry air
const Rv  = Rgas["H2O"];      // gas constant for water vapor
const Rc  = Rgas["CO2"];      // gas constant for CO2
// ------------------------------------------------------------------------------
// --> x2rho: converts x in molar fraction (volumetric concentration) (parts per
// thousand, million, billion, etc.) to density (kg/m^3).
//
// Example: if = 12 ppm of CO2 at an atmospheric pressure of 101325.0 Pa, and a
// temperature of 300.0 K, do:
// rho = x2rho(12e-6,101325.0,300.0,'CO2') 
//
// gas needs to be one of the strings in the Mmol associative array
//  ------------------------------------------------------------------------------
inline proc x2rho(
   const in x: real,          // the molar fraction
   const in p: real,          // atm pressure (Pa)
   const in T: real,          // temperature (K)
   const in gas: string       // the gas name in the Mmol associative array
   ): real {                  // the gas density in kg/m^3
// ------------------------------------------------------------------------------
// obtains the gas constant
// ------------------------------------------------------------------------------
      var R = Rgas[gas];
      return (x*p)/(R*T);
}
// ------------------------------------------------------------------------------
// these quantities are needed for viscosity and diffusivity calculations
// ------------------------------------------------------------------------------
private const mu0 = 18.18e-6;   // dynamical viscosity of air at 20o C (Pa * s)
private const T0m = 293.15;     // 20o C = 293.15 K
private const C0m = 120.00;     // C = 120 K
// ------------------------------------------------------------------------------
// --> viscair: returns the viscosity of air mu as a function of temperature T
// (in K). mu in Pa s. From: Montgomery, R. B. Viscosity and thermal
// conductivity of air and diffusivity of water vapor in air J of Meteorology,
// 1947, 4, 193-196, and Wikipedia (somewhere...)
// ------------------------------------------------------------------------------
inline proc viscair(
   T: real                    // temperature (K)
   ): real {                  // viscosity  (Pa s)
   return mu0*((T0m + C0m)/(T + C0m))*((T/T0m)**1.5);
}
// -----------------------------------------------------------------------------
// --> difmom: retorna a viscosidade cinemática do ar em função da temperatura T
// T em K, da pressão atmosférica p em Pa e da umidade específica q em
// kg/kg. nu_u em m^2/s
// -----------------------------------------------------------------------------
inline proc difmom(
   const in p: real,          // pressure (Pa)                   
   const in T: real,          // temperature (K)                 
   const in q: real           // specific humidity (kg/kg)       
   ): real {                  // kinematic viscosity of air (m^2/s)
   var rho = rho_air(p,T,q);
   return viscair(T)/rho;
}

// -----------------------------------------------------------------------------
// --> difheat: thermal diffusivity of air as a function of pressure,
// temperature and specific humidity.  See Montgomery, R. B. Viscosity and
// thermal conductivity of air and diffusivity of water vapor in air J of
// Meteorology, 1947, 4, 193-196.
// -----------------------------------------------------------------------------
inline proc difheat(
   const in p: real,          // pressure (Pa)                   
   const in T: real,          // temperature (K)                 
   const in q: real           // specific humidity (kg/kg)       
   ): real {                  // thermal diffusivity of air (m^2/s)
   return (difmom(p,T,q)/0.711);
}

// -----------------------------------------------------------------------------
// --> difvap: molecular diffusivity of water vapor in air as a function of
// pressure, temperature and specific humidity.  See Montgomery, R. B. Viscosity
// and thermal conductivity of air and diffusivity of water vapor in air J of
// Meteorology, 1947, 4, 193-196.
// -----------------------------------------------------------------------------
inline proc difvap(
   const in p: real,          // pressure (Pa)                   
   const in T: real,          // temperature (K)                 
   const in q: real           // specific humidity (kg/kg)       
   ): real {                  // diffusivity of water vapor in air m^2/s
   return (difmom(p,T,q)/0.596);
}

// -----------------------------------------------------------------------------
// --> rho2x: converts rhoi (kg/m^3) to molar fraction (ppm, ppb, etc.)
// 
// Example: if rhoi = 16 g/m^3 de H2O at p =101325 Pa, and a uma temperatura de
// 300 K,
//
// x = rho2x(16.0e-3,101325.0,300.0,'H2O');
// -----------------------------------------------------------------------------
inline proc rho2x(
   const in rhoi: real,       // gas density (kg/m^3)
   const in p: real,          // pressure (Pa)
   const in T: real,          // temperature (K)
   const in gas: string       // the gas name in the Mmol associative array
   ): real {                  // molar fraction
   return (rhoi*Rgas[gas]*T)/p;
}

// -----------------------------------------------------------------------------
// --> pressalt: atmospheric pressure (Pa) as a function of altitude (m) in a
// standard atmosphere
// -----------------------------------------------------------------------------
inline proc pressalt(
   const in alt: real         // altitude (m)
   ): real {                  // pressure (Pa)
   return 101325.0 * ((288 - 0.0065*alt)/288)**5.256;
}

// -----------------------------------------------------------------------------
// --> eexact: the exact eqn for e from T, y (and p)
// -----------------------------------------------------------------------------
proc eexact(T: real,y: real,p: real): real {
   var es = svp(T);
   return (y * es)/(1 + (y-1)*(es/p));
}
// -----------------------------------------------------------------------------
// --> sphum(e,p): specific humidity (in kg/kg) as a function of water vapor
// pressure and atmospheric pressure (both in Pa)
// -----------------------------------------------------------------------------
inline proc sphum(e: real,p: real): real {
   return e / ( 1.608*(p-e) + e );
}

// -----------------------------------------------------------------------------
// --> spheat: specific heat of humid air
// -----------------------------------------------------------------------------
inline proc spheat(q: real): real {
   return 1005.0 + 845.0*(q);
}
// -----------------------------------------------------------------------------
// R_air: gas constant for moist air
// -----------------------------------------------------------------------------
inline proc R_air(q: real): real {
   return q* Rv + ( 1.0 - q )*Rd;
}
// -----------------------------------------------------------------------------
// humid air density, kg/m^3
// -----------------------------------------------------------------------------
inline proc rho_air(p: real,T: real,q: real): real {
   return p/(R_air(q)*T);
}
   
// -----------------------------------------------------------------------------
// --> tempot: the potential temperature
// -----------------------------------------------------------------------------
inline proc tempot(p0: real,p: real,T: real,q: real): real {
   var R = R_air(q);
   var cp = spheat(q);
   return T * (p0/p)**(R/cp);
}

// -----------------------------------------------------------------------------
// --> latent: Latent heat of evaporation (in J kg^{-1}) as a function of 
// thermodynamic temperature (in K)
//
// From Dake 1972 "Evaporative Cooling of a body of water", Water Resour Res,
// (8)1087--1091.
// -----------------------------------------------------------------------------
inline proc latent(T: real): real {
   return (3142689.0 - 2356.01 * T);
}



// -------------------------------------------------------------------
// --> svp(T): saturation vapor pressure of water and its derivative
// as a function of the absolute temperature T (K)
//
// input: T 
// output (es)
//
// es == saturation vapor pressure (Pa)
//
// References:
//
// Brutsaert, W. (1982) "Evaporation Into the Atmosphere". D. Reidel 
// Publishing Company, Dordrecht, Holland.
//
// Richards, J. M. (1971) "Simple expression for the saturation vapour
// pressure of water in the range -50o to 140o", Journal of Physics D:
// Applied Physics, v. 4, n. 4, L15--L18
//
// Richards, J. M. (1971) "Simple expression for the saturation vapour
// pressure of water in the range -50o to 140o", Journal of Physics D:
// Applied Physics, v. 4, n. 6, 876
//
// Murray, F. M. (1966) "On the computation of saturation vapor
// pressure", Journal of Applied Meteorolgy v. 6, 203---204.
// -------------------------------------------------------------------
proc svp(
   const in T: real,
   param eqn: string="Richards"
   ): real {
   if eqn == "Richards" then {
      const a1 = 13.3185;
      const a2 = -1.9760;
      const a3 = -0.6445;
      const a4 = -0.1229;
      var tr = 1.0 - 373.15 / T ;
      var aux1 = ((( a4*tr + a3)*tr + a2)*tr + a1)*tr ;
      return (101325.0 * exp(aux1));
   }
   else if eqn == "Tetens" then {
      if T >= 273.15 then {
         return 610.78*exp(17.2693882*(T - 273.16)/(T - 35.86));
      }                                                      
      else {
         return 610.78*exp(21.8745584*(T - 273.16)/(T - 7.66));
      }
   }
   else {
      halt ("svp(T,eqn): eqn must be: omit == 'Richards', 'Richards', or 'Tetens'");
   }
}

// -----------------------------------------------------------------------------
// --> svpd(T): saturation vapor pressure of water and its derivative as a
// function of the absolute temperature T (K)
//
// input: T 
// output (es,ds)
//
// es == saturation vapor pressure (Pa)
// ds == derivative of es with respect to T (Pa/K)
// 
// References:
// 
// Brutsaert, W. (1982) "Evaporation Into the Atmosphere". D. Reidel Publishing
// Company, Dordrecht, Holland.
// -----------------------------------------------------------------------------
proc svpd(T: real, eqn: string="Richards"): (real, real) {
   if eqn == "Richards" then {
      const a1 = 13.3185;
      const a2 = -1.9760;
      const a3 = -0.6445;
      const a4 = -0.1229;
      const b0 = 13.3185;    // b0 = a1
      const b1 = -3.9520;    // b1 = 2.0 * a2
      const b2 = -1.9335;    // b2 = 3.0 * a3
      const b3 = -0.4916;    // b3 = 4.0 * a4
      var tr = 1.0 - 373.15 / T ;
      var aux1 = ((( a4*tr + a3)*tr + a2)*tr + a1)*tr;
      var aux2 = (( b3*tr + b2)*tr + b1)*tr + b0;
      var es = 101325.0 * exp(aux1);
      var ds = 373.15 * (es) * aux2 / (T*T);
      return (es,ds);
   }
   else if eqn == "Tetens" then {
      const b = 17.2693882;
      const b0 = 21.8745584;
      const T1 = 273.16;
      const T2 = 35.86;
      const T20 = 7.66;
      if T >= 273.15 then {
         var es = 610.78*exp(b*(T - T1)/(T - T2));
         var ds = b*(T1-T2)*es/(T - T2)**2;
         return (es,ds);
      }                                                      
      else {
         var es = 610.78*exp(b0*(T - T1)/(T - T20));
         var ds = b0*(T1-T2)*es/(T - T20)**2;
         return (es,ds);
      }
   }
   else {
      halt ("svp(T,eqn): eqn must be: 'Richards', or 'Tetens'");
   }
}

// =============================================================================
// examples
// =============================================================================
config const examples: bool ;
if examples then {
   writef('Mmol["N2"] = %12.8er\n',Mmol["N2"]);
   writef('Rgas["N2"] = %12.8er\n',Rgas["N2"]);
}
