//!/home/nldias/miniconda3/bin/python3
// -*- coding: iso-8859-1 -*-
// ====================================================================
// ==> smog: spectral estimates from smooth ogives
//
// 2016-10-08T09:46:10
// 2016-12-17T14:13:44 cosmetic changes
// 2016-12-27T13:39:05 the omega fix
// 2016-12-29T10:48:08 finally, the variances
// 2017-01-06T10:25:39 freezing the gamma factor correction
// 2017-01-10T13:57:58 completely revamped the classical correction; no
// omega yet.
// 2017-01-20T11:15:05 fixing bug: Og[k-ow[k]] - Og[k+ow[k]+1] in smog
// 2017-01-27T16:04:19 numba does not like variables created inside
// ifs and loops. Moreover, numba is uncomfortable with *args kinds of
// things. Simplified things for numba: tilogvar uses numpy only, while
// _otilg now does the numba optimization. Seems to be working.
// 2021-03-22T15:25:24 translating from Python to Chapel
//
// 2021-03-23T09:24:38 all private procs assume 0-based arrays with
// the correct lengths. The public procs must check consistency
//
// fix for k in unbias: 3 --> 2
// ====================================================================
// --------------------------------------------------------------------
// License
//
// Copyright Nelson Luís da Costa Dias, 2017
//
// This file is part of pub-smog.
// pub-smog is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// pub-smog is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// You should have received a copy of the GNU General Public License
// along with pub-smog.  If not, see <http://www.gnu.org/licenses/>.
// -------------------------------------------------------------------



// -------------------------------------------------------------------
// --> ogive(deltaf,aG,aO): use very simple integration to calculate
// the ogive from a spectrum G with data sampled at frequency deltaf.
//
// 2017-01-10T09:40:31 going back to a single frequency for Os and Gs
//
// 2016-10-08T09:47:12 re-created with numpy
//
// 2021-03-22T15:59:57 re-programmed (it's simpler!) in Chapel   
// -------------------------------------------------------------------
proc ogive(
   const in deltaf: real,
   const ref aG: [] real,
   ref aO: [] real
) {
   assert(aG.rank == 1);
   assert(aO.rank == 1);
   const M = aG.size-1;
   assert(aO.size == M+1);
   ref G = aG.reindex(0..M);
   ref Og = aO.reindex(0..M);
   Og[M] = G[M];
   for k in 0..M-1 by -1 do {
      Og[k] = G[k] + Og[k+1];
   }
   Og *= deltaf;
}

// -------------------------------------------------------------------
// --> wog returns the widths of the smoothing over array a:
// 0,0,1,2, ...,p-2,p-1,p, p, p, ..., p,p-1,p-2, ...,2,1,0,0
//
// you must run wog *before* smog
// -------------------------------------------------------------------
proc wog(
   const in p: int,           // the size of the half-window
   const ref aa: [] real,     // the array to be smoothed
   ref aw: [] int             // the widths of the smoothing array
) {
   assert(aa.rank == 1);
   assert(aw.rank == 1);
   const M = aa.size-1;
   assert(aw.size == M+1);
   ref a = aa.reindex(0..M);
   ref w = aw.reindex(0..M);
   w[0] = 0;
   w[1..p] = 0..p-1;
   w[p+1..M-p-1] = p;
   w[M-p..M-1] = 0..p-1 by -1;
   w[M] = 0;
}

// -------------------------------------------------------------------
// corrections
// -------------------------------------------------------------------
record correction {
   param M: int;
   var me: [0..M] real = 0.0;
   var gamma_O: [0..M] real = 1.0;
   var gamma_G: [0..M] real = 1.0;
   var omega = 0.0;
}
// -------------------------------------------------------------------
// -->zmog: smoothes a raw ogive Or with a windows of half-width ow
// and sampling frequency deltaf and returns a smoothed ogive Os, and
// a smoothed spectrum Gs.
// -------------------------------------------------------------------
proc zmog(
   const in p: int,
   const in deltaf: real,     // frequency interval
   const ref aGr: [] real,    // the raw spectrum
   const ref aOr: [] real,    // the raw ogive
   ref aGs: [] real,          // the smoothed spectrum
   ref aOs: [] real           // the smoothed ogive
) {
   assert(aGr.rank == 1);
   assert(aOr.rank == 1);
   assert(aGs.rank == 1);
   assert(aOs.rank == 1);
   const M = aOr.size - 1;    // check their size
   assert(aGr.size == M+1);
   assert(aGs.size == M+1);
   assert(aOs.size == M+1);
   ref Gr = aGr.reindex(0..M);
   ref Or = aOr.reindex(0..M);
   ref Os = aOs.reindex(0..M);
   ref Gs = aGs.reindex(0..M);
   var sign: real;
   if Or[0] == 0 then {
      halt("--> smog: something terribly wrong happened");
   }
   else if Or[0] > 0.0 then { // prevailing sign of the cospectrum
      sign = 1.0;
   }
   else {
      sign = -1.0;
   }
   Os[0] = Or[0];             // endpoints unaltered
   Os[1] = Or[1];
   for k in 2..M-2 do {       // smoothing loop
      var q = min(p,k);
      var r = min(p,M-k);
      Os[k] = (+ reduce Or[k-q..k+r])/(q+r+1);
   }
   Os[M-1] = Or[M-1];
   Os[M] = Or[M];           // endpoints unaltered
// --------------------------------------------------------------------
// og-derived smooth spectrum
// --------------------------------------------------------------------
   Gs[0] = 0.0;               // endpoints unaltered
   Gs[1] = (Or[1]-Or[2])/deltaf;
   for k in 2..M-2 do {
      var q = min(p,k);
      var r = min(p,M-k-1);
//      writeln("k, q, r ",k, "  ",q,"  ",r);
      Gs[k] = (Os[k-q] - Os[k+r+1])/((q+r+1)*deltaf);
   }
   Gs[M-1] = (Or[M-1]-Or[M])/deltaf;
   Gs[M] = Or[M]/deltaf;     // endpoints unaltered
}

// -------------------------------------------------------------------
// --> unbias: new bias correction algorithm based on exact arithmetic
// means for both Os and Gs:
//
// 2017-01-10T13:50:18
//
// 2021-03-22T18:23:34 this is a private proc; all arrays are assumed
// to be 1D, of the same size, and zero-based.
// -------------------------------------------------------------------
private proc unbias(
   const in deltaf: real,
   const in sign: real,
   const ref ow: [] int,
   ref Gs: [] real,
   ref Os: [] real,
   ref cor: correction
) {
   const M = Os.size - 1;          // == NN//2
   var pOs = sign*Os;              // positive copies of the ogive
   var pGs = sign*Gs;              // and the spectrum
   for k in 2..M-2 by -1 do {      // only k's for which ow[k] >= 2
      const fa = (k-ow[k])*deltaf; // left frequency
      const fb = (k+ow[k])*deltaf; // right frequency
      const fc = k*deltaf;         // center (arith mean) frequency
      const alpha = pOs[k];
      if ( alpha <= 0.0 ) then {
         continue;
      }
      const beta = pGs[k];
      if ( beta < 0.0 ) then {
         continue;
      }
      var m = mroot(fa,fb,beta/alpha);
      if ( m <= 0.0 ) then {            // sometimes m misbehaves
         m = 0.01;
      }
      else if ( m > 5.0 ) then {
         m = 5.0;
      }
// --------------------------------------------------------------------
// the bias correction factors for the ogive
// --------------------------------------------------------------------
      const fg = gmeanf(fa,fb);         // the geometric-mean frequency
      var gm = gmeano(fa,fb,1.0,m);
      var am = ameano(fa,fb,1.0,m);
      cor.me[k] = m;
      cor.gamma_O[k] = (gm/am)*(fg/fc)**m;
      Os[k] *= cor.gamma_O[k];          // bias correction
// --------------------------------------------------------------------
// the bias correction for the spectrum
// --------------------------------------------------------------------
      gm = gmeang(fa,fb,1.0,m);
      am = ameang(fa,fb,1.0,m);
      cor.gamma_G[k] = (gm/am)*(fg/fc)**(m+1);
   }
// --------------------------------------------------------------------
// the omega correction
// --------------------------------------------------------------------
   var omega_a = -10.0;
   var omega_b = +10.0;
   var omega = (omega_a + omega_b)/2.0;
   assert(omegazero(sign,deltaf,omega_a,Os[0],Gs,cor.gamma_G) > 0.0);
   assert(omegazero(sign,deltaf,omega_b,Os[0],Gs,cor.gamma_G) < 0.0);
   const eps = 1.0e-4;
   while (abs(omega_b - omega_a) > eps) do {
      if omegazero(sign,deltaf,omega,Os[0],Gs,cor.gamma_G) >= 0.0 then {
         omega_a = omega;
      }
      else {
         omega_b = omega;
      }
      omega = (omega_a + omega_b)/2.0;
   }
   cor.omega = omega;
   Gs *= (cor.gamma_G**omega);
}


// -------------------------------------------------------------------
// --> omegazero: returns sum (gamma_G[k])^omega Gs[k] - Os[0]
// omegazero is used in unbias
// 
// 2017-01-17T16:27:44
//
// 2021-03-23T09:27:29 assumes 0-based arrays
// -------------------------------------------------------------------
private proc omegazero(
   const in sign: real,
   const in deltaf: real,
   const in omega: real,
   const in Oszero: real,
   const ref Gs: [] real,
   const ref fGcor: [] real
): real {
   const M = Gs.size - 1;
   var locor = fGcor**omega;
   var Gcor = Gs*deltaf*locor;
   var ha = sign*(+ reduce Gcor[1..M]);
   var hb = sign*Oszero;
   return ha - hb;
}

// --------------------------------------------------------------------
// --> tilogvar_gxx: calculates the variance of the smoothed spectrum
// gtilvar and of smoothed ogive otilvar
// --------------------------------------------------------------------
proc tilogvar_gxx(
   const in deltaf: real,
   const ref aow: [] int,
   const ref Gxx: [] real,
   ref agtilvar: [] real,
   ref aotilvar: [] real
) {
   assert (aow.rank == 1);
   assert (Gxx.rank == 1);
   assert (agtilvar.rank == 1);
   assert (aotilvar.rank == 1);
   const M = aow.size - 1;
   assert (Gxx.size == M+1);
   assert (agtilvar.size == M+1);
   assert (aotilvar.size == M+1);
   ref ow = aow.reindex(0..M);
   ref gtilvar = agtilvar.reindex(0..M);
   ref otilvar = aotilvar.reindex(0..M);
   var sig2m: [0..M] real;
// --------------------------------------------------------------------
// Which case are we talking about?
// --------------------------------------------------------------------
   sig2m = Gxx**2;         // over the whole array
   otilg(deltaf,ow,sig2m,gtilvar,otilvar);
}

// --------------------------------------------------------------------
// --> tilogvar_gxx: calculates the variance of the cospectrum:
// --------------------------------------------------------------------
proc tilogvar_cxy(
   const in deltaf: real,
   const ref aow: [] int,
   const ref Gxx: [] real,
   const ref Gyy: [] real,
   const ref Coxy: [] real,
   const ref Quxy: [] real,
   ref agtilvar: [] real,
   ref aotilvar: [] real
) {
   assert (aow.rank == 1);
   assert (Gxx.rank == 1);
   assert (Gyy.rank == 1);
   assert (Coxy.rank == 1);
   assert (Quxy.rank == 1);
   assert (agtilvar.rank == 1);
   assert (aotilvar.rank == 1);
   const M = aow.size - 1;
   assert (Gxx.size == M+1);
   assert (Gyy.size == M+1);
   assert (Coxy.size == M+1);
   assert (Quxy.size == M+1);
   assert (agtilvar.size == M+1);
   assert (aotilvar.size == M+1);
   ref ow = aow.reindex(0..M);
   ref gtilvar = agtilvar.reindex(0..M);
   ref otilvar = aotilvar.reindex(0..M);
   var sig2m: [0..M] real;
// --------------------------------------------------------------------
// Which case are we talking about?
// --------------------------------------------------------------------
   sig2m = (Gxx*Gyy + Coxy**2 - Quxy**2)/2.0;
   otilg(deltaf,ow,sig2m,gtilvar,otilvar);
}


// --------------------------------------------------------------------
// --> tilogvar_qxy: calculates the variance of the quadspectrum.
// --------------------------------------------------------------------
proc tilogvar_qxy(
   const in deltaf: real,
   const ref aow: [] int,
   const ref Gxx: [] real,
   const ref Gyy: [] real,
   const ref Coxy: [] real,
   const ref Quxy: [] real,
   ref agtilvar: [] real,
   ref aotilvar: [] real
) {
   assert (aow.rank == 1);
   assert (Gxx.rank == 1);
   assert (Gyy.rank == 1);
   assert (Coxy.rank == 1);
   assert (Quxy.rank == 1);
   assert (agtilvar.rank == 1);
   assert (aotilvar.rank == 1);
   const M = aow.size - 1;
   assert (Gxx.size == M+1);
   assert (Gyy.size == M+1);
   assert (Coxy.size == M+1);
   assert (Quxy.size == M+1);
   assert (agtilvar.size == M+1);
   assert (aotilvar.size == M+1);
   ref ow = aow.reindex(0..M);
   ref gtilvar = agtilvar.reindex(0..M);
   ref otilvar = aotilvar.reindex(0..M);
   var sig2m: [0..M] real;
// --------------------------------------------------------------------
// Which case are we talking about?
// --------------------------------------------------------------------
   sig2m = (Gxx*Gyy - Coxy**2 + Quxy**2)/2.0;
   otilg(deltaf,ow,sig2m,gtilvar,otilvar);
}


// -------------------------------------------------------------------
// --> otilg: what do I do?
// -------------------------------------------------------------------
private proc otilg(
   const in deltaf: real,
   const ref ow: [] int,
   const ref sig2m: [] real,
   ref gtilvar: [] real,
   ref otilvar: [] real
) {
   const M = sig2m.size-1;
   gtilvar = 0.0;  // will return this
   otilvar = 0.0;  // will return this too
// --------------------------------------------------------------------
// the spectrum is simpler
// --------------------------------------------------------------------
   for k in 1..M do {    // calculating loop
      gtilvar[k] = (+ reduce sig2m[k-ow[k]..k+ow[k]])/(2*ow[k]+1)**2;
   }
// --------------------------------------------------------------------
// the ogive is more complicated
// --------------------------------------------------------------------
   for k in 0..M do {
      var ovar = 0.0;
      for l in k-ow[k]+1..k+ow[k] do {
         ovar += ((l - k + ow[k]+1)**2)*sig2m[l];
      }
      ovar /= (2*ow[k]+1)**2;
      for l in k+ow[k]..M do {
         ovar += sig2m[l];
      }
      otilvar[k] = ovar*(deltaf**2);
   }
}





// -------------------------------------------------------------------
// --> ameano: arithmetic mean, power-law ogive
// -------------------------------------------------------------------
private inline proc ameano(
   const in fa: real,
   const in fb: real,
   const in C: real,
   const in m: real): real {
   if m == 1.0 then {
      return (C/(fb-fa))*log(fb/fa);
   } 
   else {
      return (C/((fb-fa)*(-m+1.0)))*(fb**(-m+1.0) - fa**(-m+1.0));
   }
}

// -------------------------------------------------------------------
// --> ameang: arithmetic mean, power-law spectrum
// -------------------------------------------------------------------
private proc ameang(
   const in fa: real,
   const in fb: real,
   const in C: real,
   const in m: real): real {
   return (C/(fb-fa))*(fa**(-m) - fb**(-m));
}

// -------------------------------------------------------------------
// --> gmeano geometric mean, power-law ogive
// -------------------------------------------------------------------
private proc gmeano(
   const in fa: real,
   const in fb: real,
   const in C: real,
   const in m: real): real {
   var aux = log(C);
   aux += m;
   aux -= (m/(fb-fa))*(fb*log(fb) - fa*log(fa));
   return exp(aux);
}

// -------------------------------------------------------------------
// --> gmeang: geometric mean, power-law spectrum
// -------------------------------------------------------------------
private proc gmeang(
   const in fa: real,
   const in fb: real,
   const in C: real,
   const in m: real
): real {
   var aux = log(m*C);
   aux += (m+1.0);
   aux -= ((m+1.0)/(fb-fa))*(fb*log(fb) - fa*log(fa));
   return exp(aux);
}

// -------------------------------------------------------------------
// --> gmeanf: geometric-mean frequency
// -------------------------------------------------------------------
private proc gmeanf(
   const in fa: real,
   const in fb: real
): real {
   var aux = fb*log(fb) - fa*log(fa) + fa - fb;
   aux /= (fb - fa);
   return exp(aux);
}




// ---------------------------------------------------------
// --> mroot: calculates the numerical solution of
//
// goa = ameang(fa,fb,1.0,m)/ameano(fa,fb,1.0,m)
//
// for the unknown m, using bisection.
// ---------------------------------------------------------
private proc mroot(
   const in fa: real,
   const in fb: real,
   const in goa: real
) {
   var ma = -1.0;
   var mb = 5.0;
   var m = (ma + mb)/2.0;
   const epsm = 1.0e-4;
   var kontrol = 0;
   while abs(mb - ma) > epsm do {
      var res = ameang(fa,fb,1.0,m)/ameano(fa,fb,1.0,m) - goa;
      if res >= 0.0 then {
         mb = m;
      }
      else {
         ma = m;
      }
      kontrol += 1;
      if kontrol > 1000 then {
         halt('too many (1001) iterations');
      }
      m = (ma + mb)/2.0;
   }
   return m;
}
