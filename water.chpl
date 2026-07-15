config const describe = false;
const doc = "\
=================================================================================\
==> water: thermodynamic properties of water                                     \
                                                                                 \
based on Batchelor's tables of rhow and cp: only interpolates                    \
viscosities from:                                                                \
Kestin and Wakeham,                                                              \
Viscosity of Liquid Water in the Range -8°C to 150°C,                            \
Journal of Physical and Chemical Reference v. 7 pp 941--948                      \
                                                                                 \
=================================================================================\
";
if describe then {
   writeln(doc);
   exit(0);
}
// =============================================================================
// 2021-04-16T14:41:23 a new star is born
// 2024-11-21T10:57:00 viscosity table
// =============================================================================
use ssr only interp;
private const Tw =
   [ 0.0,      5.0,   10.0,   15.0,   20.0,   25.0,   30.0,   35.0,   40.0,   50.0,   60.0,   70.0,   80.0,   90.0,  100.0];
private const rhow =
   [ 999.9, 1000.0,  999.7,  999.1,  998.2,  997.1,  995.7,  994.1,  992.3,  988.1,  983.2,  977.8,  971.8,  965.3,  958.4];
private const cpw =
   [4217.0, 4202.0, 4192.0, 4186.0, 4182.0, 4179.0, 4178.0, 4178.0, 4178.0, 4180.0, 4184.0, 4189.0, 4196.0, 4205.0, 4216.0];
private const muw =
   [1791.5, 1519.3, 1307.0, 1138.3, 1002.0, 890.2, 797.3, 719.1, 652.7, 547.1, 467.0, 404.6, 355.1, 315.0, 282.1]*1.0e-6;
// -----------------------------------------------------------------------------
// --> rhow: density of water as a function of temperature and pressure.
// -----------------------------------------------------------------------------
proc rho_water(
   in T: real       // temperature, K (for 0--100\degree{C})
   ): real  {
   T -= 273.15;     // formulae are in Celsius
   var rho = interp(T,Tw,rhow);
   return rho;
}
// -----------------------------------------------------------------------------
// --> cp_water: specific heat at constant pressure of water as a
// function of temperature
// -----------------------------------------------------------------------------
proc cp_water(
   in T: real       // temperature, K (for 0--100\degree{C})
   ): real  {
   T -= 273.15;     // formulae are in Celsius
   var cp = interp(T,Tw,cpw);
   return cp;
}
// -----------------------------------------------------------------------------
// --> mu_water: viscosity of water at constant pressure as a function
// of temperature
// -----------------------------------------------------------------------------
proc mu_water(
   in T: real
   ): real {
   T -= 273.15;
   var mu = interp(T,Tw,muw);
   return mu;
}
      