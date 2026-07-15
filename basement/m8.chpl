// =============================================================================
// ==> eleets-eb: good old energy-budget method
//
// 2021-04-26T10:02:59
//
// 2021-04-28T14:06:16 daily energy budget !
// =============================================================================
use IO, DateTime;                  // system modules
use rdayeleets;                    // reads all daily data
assert(firstmonth == 1  &&  firstday == 1);
assert(lastmonth == 12  &&  lastday == 31);
use solterra;
use atmgas;
use angles;
use evap;     
use nstat;
var dlat = todec(lat);             // convert latitude to decimal
var rlat = dlat*(pi/2)/90.0;       // convert latitude to radians
IniPar(alt,0.97);
Prescott(a=0.25,b=0.50);           // so says FAO
// -----------------------------------------------------------------------------
// Output a series of monthly evaporation now
// -----------------------------------------------------------------------------
// use numint only trapepsilon;
use water;
// -----------------------------------------------------------------------------
// for simplicity, rho, cp, gamma and L0 will all be held constant
// -----------------------------------------------------------------------------
// const (nbT,ST0) = butsum(NAN,T0);
// const T0mean = ST0/nbT;
// const e0 = svp(T0mean+273.15);
// const q0 = 0.622*e0/P;
// const cp = spheat(q0);
// const L0 = latent(T0mean+273.15);
// const Tv = (1 + 0.61*q0)*T0mean;
// const rho = rho_air(P,T0mean,q0);
const L0 = 2.464e6 ;
const cp = 1005;
const gamma = cp*P/(0.622*L0);
const rho = rho_air(P,288.15,0.0);
// -----------------------------------------------------------------------------
// now loop over all days of the record
// -----------------------------------------------------------------------------
var firstdate = new date(firstyear,firstmonth,firstday);
var lastdate =  new date(lastyear,lastmonth,lastday);
var d15days = new timedelta(15);
for iy in 1..ny do {
   var year = iy + firstyear - 1;
   for mo in 1..12 do {
      for day in 1..daysInMonth(year,mo) do {
         var (delta,rr) = ddse(year,mo,day);
         var (Rsea,dsmax,Z) = rsdsZ(rlat,rr,delta);
         var S = min(1.0,nn[iy,mo,day]/dsmax);
         var alb = LakeAlbedo(Z,S);
         Radiation(alb,ea[iy,mo,day],Ta[iy,mo,day]+273.15,
                   T0[iy,mo,day]+273.15,Rsea,S,Rs[iy,mo,day],
                   Ra[iy,mo,day],Re[iy,mo,day],Rl0[iy,mo,day]);
         e0[iy,mo,day] = svp(T0[iy,mo,day]+273.15); 
// -----------------------------------------------------------------------------
// the calculation of D: looks for 15 days ahead and behind
// -----------------------------------------------------------------------------
         var today = new date(year,mo,day);
         var yedel = today - d15days;
         if yedel < firstdate then {
            yedel = firstdate;
         }
         var todel = today + d15days;
         if todel > lastdate then {
            todel = lastdate;
         }
         var Vlow = min(Vl[yedel.year-firstyear+1,yedel.month,yedel.day],
                        Vl[todel.year-firstyear+1,todel.month,todel.day]);
         // if fio then {
         //    var fC: Acte;
         //    (Vlow,eps) = trapepsilon(1.0,zfundo,zlow,fC);
         // }
         // else {
         //    var fA: Atot;
         //    (Vlow,eps) = trapepsilon(1.0,zfundo,zlow,fA);
         // }
         var rhow = rho_water(T0[iy,mo,day]+273.15);
         var cpw = cp_water(T0[iy,mo,day]+273.15);
         DD[iy,mo,day] = ((rhow*cpw)/(Al[iy,mo,day]*Deltat_d*31))*
                 (T0[todel.year-firstyear+1,todel.month,todel.day]-
                  T0[yedel.year-firstyear+1,yedel.month,yedel.day])*Vlow;
// -----------------------------------------------------------------------------
// the energy-budget method
// -----------------------------------------------------------------------------
         var Bo = gamma*(T0[iy,mo,day] - Ta[iy,mo,day])/(e0[iy,mo,day] - ea[iy,mo,day]);
         LE[iy,mo,day] = (1.0/(1.0+Bo))*(Rl0[iy,mo,day] - DD[iy,mo,day]);
         HH[iy,mo,day] = Bo*LE[iy,mo,day];
         El[iy,mo,day] = LE[iy,mo,day]*Deltat_d/L0;
      }
   }
}
// -----------------------------------------------------------------------------
// there are 46 periods of 8 days (the last has fewer) in a year
// -----------------------------------------------------------------------------
const n8 = 46;
var El8: [1..ny,0..#46] real = 0.0;
var Elm: [1..ny,1..12] real = 0.0;      // monthly lake evaporation
var Erm: [1..ny,1..12] real = 0.0;      // monthly land evapotranspiration
// -----------------------------------------------------------------------------
// cumulative evaporation in n8 periods
// -----------------------------------------------------------------------------
for iy in 1..ny do {
   var jday = 0;
   for mo in 1..12 do {
      for day in 1..daysInMonth(firstyear+iy-1,mo) do {
         El8[iy,jday/8] += El[iy,mo,day];
         jday += 1;
      }
   }
}
for iy in 1..ny do {
   for mo in 1..12 do {
      for day in 1..daysInMonth(firstyear+iy-1,mo) do {
         Elm[iy,mo] += El[iy,mo,day];
         Erm[iy,mo] += Er[iy,mo,day];
      }
   }
}
// -----------------------------------------------------------------------------
// mean monthly evaporation and evapotranspiration
// -----------------------------------------------------------------------------
var Elmon: [1..12] real = 0.0;
var Ermon: [1..12] real = 0.0;
for mo in 1..12 do {
   for iy in 1..ny do {
      Elmon[mo] += Elm[iy,mo];
      Ermon[mo] += Erm[iy,mo];
   }
   Elmon[mo] /= ny;
   Ermon[mo] /= ny;
}
use ssr;
var EL = sum(Elmon);     // mean yearly lake evaporation
var ER = sum(Ermon);     // mean yearly land evapotranspiratoin
// -----------------------------------------------------------------------------
// print results for daily data
// -----------------------------------------------------------------------------
var cod = openwriter(lake+"_D-eb.out");
cod.writef("%s\n","# "+lake);
cod.writef("# EL  = %+6.2dr\n",EL);
cod.writef("# ER  = %+6.2dr\n",ER);
cod.writef("# EN  = %+6.2dr\n",EL-ER);
cod.writef("#year mo dd  Rs(W/m2)  Ra(W/m2)  Re(W/m2) Rlo(W/m2)  D (W/m2)    H(W/m2)   LE(W/m2)  T0(oC)    Ta(oC)    e0(Pa)    ea(Pa)    El(mm)\n");
for iy in 1..ny do {
   var year = firstyear + iy - 1;
   for mo in 1..12 do {
      for day in 1..daysInMonth(year,mo) do {
         cod.writef("%5i %02i %02i"+"%10.2dr"*11+"\n",
                    year, mo, day,
                    Rs[iy,mo,day], Ra[iy,mo,day], Re[iy,mo,day], Rl0[iy,mo,day], DD[iy,mo,day],
                    HH[iy,mo,day], LE[iy,mo,day], T0[iy,mo,day], Ta[iy,mo,day],
                    e0[iy,mo,day], ea[iy,mo,day], El[iy,mo,day]);
      }
   }
}
cod.close();
// -----------------------------------------------------------------------------
// now monthly results
// -----------------------------------------------------------------------------
var dm = {1..ny,1..12};
var Rsm: [dm] real = 0.0;
var Ram: [dm] real = 0.0;
var Rem: [dm] real = 0.0;
var Rlm: [dm] real = 0.0;
var DDm: [dm] real = 0.0;
var HHm: [dm] real = 0.0;
var LEm: [dm] real = 0.0;
var T0m: [dm] real = 0.0;
var Tam: [dm] real = 0.0;
var e0m: [dm] real = 0.0;
var eam: [dm] real = 0.0;
for iy in 1..ny do {
   var year = firstyear + iy - 1;
   for mo in 1..12 do {
      var nd = daysInMonth(year,mo);
      for day in 1..nd do {
         Rsm[iy,mo] += Rs[iy,mo,day];
         Ram[iy,mo] += Ra[iy,mo,day];
         Rem[iy,mo] += Re[iy,mo,day];
         Rlm[iy,mo] += Rl0[iy,mo,day];
         HHm[iy,mo] += HH[iy,mo,day];
         LEm[iy,mo] += LE[iy,mo,day];
         T0m[iy,mo] += T0[iy,mo,day];
         Tam[iy,mo] += Ta[iy,mo,day];
         e0m[iy,mo] += e0[iy,mo,day];
         eam[iy,mo] += ea[iy,mo,day];
      }
      Rsm[iy,mo] /= nd;
      Ram[iy,mo] /= nd;
      Rem[iy,mo] /= nd;
      Rlm[iy,mo] /= nd;
      HHm[iy,mo] /= nd;
      LEm[iy,mo] /= nd;
      T0m[iy,mo] /= nd;
      Tam[iy,mo] /= nd;
      e0m[iy,mo] /= nd;
      eam[iy,mo] /= nd;
// -----------------------------------------------------------------------------
// the rate of change of stored enthalpy is different!
// -----------------------------------------------------------------------------      
      DDm[iy,mo] = Rlm[iy,mo] - HHm[iy,mo] - LEm[iy,mo];
   }
}
// -----------------------------------------------------------------------------
// print results for monthly data
// -----------------------------------------------------------------------------
var com = openwriter(lake+"_M-eb.out");
com.writef("%s\n","# "+lake);
com.writef("# EL  = %+6.2dr\n",EL);
com.writef("#year mo   Rs(W/m2)  Ra(W/m2)  Re(W/m2) Rlo(W/m2)  D (W/m2)    H(W/m2)   LE(W/m2)  T0(oC)    Ta(oC)    ea(Pa)    El(mm)\n");
for iy in 1..ny do {
   var year = firstyear + iy - 1;
   for mo in 1..12 do {
      com.writef("%5i %02i "+"%10.2dr"*11+"\n",
                 year, mo, 
                 Rsm[iy,mo], Ram[iy,mo], Rem[iy,mo], Rlm[iy,mo], DDm[iy,mo],
                 HHm[iy,mo], LEm[iy,mo], T0m[iy,mo], Tam[iy,mo], eam[iy,mo],Elm[iy,mo]);
   }
}
com.close();
// -----------------------------------------------------------------------------
// now 8-day results
// -----------------------------------------------------------------------------
var d8 = {1..ny,0..#46};
var Rs8: [d8] real = 0.0;
var Ra8: [d8] real = 0.0;
var Re8: [d8] real = 0.0;
var Rl8: [d8] real = 0.0;
var DD8: [d8] real = 0.0;
var HH8: [d8] real = 0.0;
var LE8: [d8] real = 0.0;
var T08: [d8] real = 0.0;
var Ta8: [d8] real = 0.0;
var ea8: [d8] real = 0.0;
var e08: [d8] real = 0.0;
var da8: [d8] date;
for iy in 1..ny do {
   var year = firstyear + iy - 1;
   var jday = 0;
   for mo in 1..12 do {
      var nd = daysInMonth(year,mo);
      for day in 1..nd do {
         Rs8[iy,jday/8] += Rs[iy,mo,day];
         Ra8[iy,jday/8] += Ra[iy,mo,day];
         Re8[iy,jday/8] += Re[iy,mo,day];
         Rl8[iy,jday/8] += Rl0[iy,mo,day];
         HH8[iy,jday/8] += HH[iy,mo,day];
         LE8[iy,jday/8] += LE[iy,mo,day];
         T08[iy,jday/8] += T0[iy,mo,day];
         Ta8[iy,jday/8] += Ta[iy,mo,day];
         e08[iy,jday/8] += ea[iy,mo,day];
         ea8[iy,jday/8] += ea[iy,mo,day];
         jday += 1;
      }
   }
}
for iy in 1..ny do {
   var year = firstyear + iy - 1;
   var jday = 0;
   var leap: bool = isLeapYear(year);
   for mo in 1..12 do {
      var nd = daysInMonth(year,mo);
      for day in 1..nd do {
         if jday % 8 == 3 then {
            da8[iy,jday/8] = new date(year,mo,day);
         }
         var full8 = (jday % 8 == 7);
         var endy = (mo == 12 && day == 31) ;
         var ndays: int;
         if full8 then {
            ndays = 8 ;
         }
         else if endy then {
            if leap then {
               ndays = 6;
            }
            else {
               ndays = 5;
            }
         }
         if full8 || endy then {
            Rs8[iy,jday/8] /= ndays;
            Ra8[iy,jday/8] /= ndays;
            Re8[iy,jday/8] /= ndays; 
            Rl8[iy,jday/8] /= ndays;
            HH8[iy,jday/8] /= ndays;
            LE8[iy,jday/8] /= ndays;
            T08[iy,jday/8] /= ndays;
            Ta8[iy,jday/8] /= ndays;
            e08[iy,jday/8] /= ndays;
            ea8[iy,jday/8] /= ndays;
// -----------------------------------------------------------------------------
// the rate of change of stored enthalpy is different!
// -----------------------------------------------------------------------------      
            DD8[iy,jday/8] = Rl8[iy,jday/8] - HH8[iy,jday/8] - LE8[iy,jday/8];
         }
         jday += 1;
      }
   }
}
// -----------------------------------------------------------------------------
// print results for 8-day data
// -----------------------------------------------------------------------------
var co8 = openwriter(lake+"_8-eb.out");
co8.writef("%s\n","# "+lake);
co8.writef("# EL  = %+6.2dr\n",EL);
co8.writef("#year mo dd   Rs(W/m2)  Ra(W/m2)  Re(W/m2) Rlo(W/m2)  D (W/m2)    H(W/m2)   LE(W/m2)  T0(oC)    Ta(oC)    e0(Pa)    ea(Pa)    El(mm)\n");
for iy in 1..ny do {
   var year = firstyear + iy - 1;
   for j8 in 0..#46 do {
      co8.writef("%5i %02i %02i "+"%10.2dr"*11+"\n",
                 da8[iy,j8].year, da8[iy,j8].month, da8[iy,j8].day,
                 Rs8[iy,j8], Ra8[iy,j8], Re8[iy,j8], Rl8[iy,j8], DD8[iy,j8],
                 HH8[iy,j8], LE8[iy,j8], T08[iy,j8], Ta8[iy,j8],
                 e08[iy,j8], ea8[iy,j8], El8[iy,j8]);
   }
}
co8.close();
writeln((+ reduce DD8)/(46*17));
