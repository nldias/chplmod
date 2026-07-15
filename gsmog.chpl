// =============================================================================
// ==> gsmog: spectral estimates from smooth ogives
//
// 2021-03-30T14:46:38 total change of direction: now using geometric
// means
// 2021-08-26T10:34:47 cosmetic changes on gsmog
// 2022-12-15T08:44:56 cleaning up geometric running means
// =============================================================================
// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------

use nstat only gmom;
// -----------------------------------------------------------------------------
// --> ogive(deltaf,aG,aO): use very simple integration to calculate
// the ogive from a spectrum G with data sampled at frequency deltaf.
//
// 2017-01-10T09:40:31 going back to a single frequency for Os and Gs
//
// 2016-10-08T09:47:12 re-created with numpy
//
// 2021-03-22T15:59:57 re-programmed (it's simpler!) in Chapel   
// -----------------------------------------------------------------------------
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
// you must run wog *before* grm
// -------------------------------------------------------------------
proc wog(
   const in p: int,           // the size of the half-window
   const ref aa: [] real,     // the array to be smoothed
   ref aw: [] int             // the widths of the smoothing array
) {
   assert(aa.rank == 1);
   assert(aw.rank == 1);
   const M = aa.size-1;
   assert(aw.size == aa.size);
   ref w = aw.reindex(0..M);
   w[0] = 0;
   w[1..p] = 0..p-1;
   w[p+1..M-p-1] = p;
   w[M-p..M-1] = 0..p-1 by -1;
   w[M] = 0;
}

// -----------------------------------------------------------------------------
// --> grm: geometric running mean of x
// -----------------------------------------------------------------------------
proc grm(
   const ref aow: [] int,               // smoothing window indices
   const ref axr: [] real,              // input xr
   ref axs: [] real                     // output xs
   ) where (aow.rank == 1 && axr.rank == 1 && axs.rank == 1)  { 
   const M = aow.size - 1;              // check their size
   assert(axr.size == M+1);
   assert(axs.size == M+1);
   ref ow = aow.reindex(0..M);          // re-index
   ref xr = axr.reindex(0..M);
   ref xs = axs.reindex(0..M);
   for k in 0..M do {
      xs[k] = gmom(xr[k-ow[k]..k+ow[k]]);
   }
}

// -----------------------------------------------------------------------------
// --> gsmof: geometric mean of frequencies
// -----------------------------------------------------------------------------
proc gsmof(
   const ref aow: [] int,               // smoothing window indices
   ref afs: [] real                     // the smoothed frequencies
   ) where (aow.rank == 1 && afs.rank == 1)  { 
   const M = aow.size - 1;             // check their size
   assert(afs.size == M+1);
   ref ow = aow.reindex(0..M);         // re-index
   ref fs = afs.reindex(0..M);
   var fr = fs;                        // copy fs to fr
   for k in 0..M do {
      fs[k] = gmom(fr[k-ow[k]..k+ow[k]]);
   }
}
// -----------------------------------------------------------------------------
// --> gsmog: smoothes an array Gs with a windows of half-width ow
// -----------------------------------------------------------------------------
proc gsmog(
   const ref aow: [] int,          // smoothing window indices
   ref aGs: [] real                // the smoothed spectrum
   ) where (aow.rank == 1 && aGs.rank == 1 ) {
   assert(aow.size == aGs.size);
   const M = aow.size - 1;         // check their size
   ref ow = aow.reindex(0..M);     // re-index
   ref Gs = aGs.reindex(0..M);
   var Gr = Gs;                    // copy Gs to Gr
   for k in 0..M do {
      Gs[k] = gmom(Gr[k-ow[k]..k+ow[k]]);
   }
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





