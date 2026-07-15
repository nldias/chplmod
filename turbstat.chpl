// =============================================================================
// ==> turbstat: specific processing of turbulence statistics
//
// 2020-12-11T16:26:58 a new star is born
//
// 2021-03-19T09:14:19 including mdfast and kfastcontrol translated
// from the Python versions
//
// 2021-03-21T12:45:31 finally fixed ibeg and iend and now everything
// is working fine
// =============================================================================
use nstat, ssr, mdot;

// -----------------------------------------------------------------------------
// --> pmstat1: conditional positive and negative means
// -----------------------------------------------------------------------------
proc pmxystat1(
   const in xmin: real,       // the threshold
   const ref x: [] real,      // usually the vertical velocity for REA
   const ref y: [] real       // the data for the conditional means
   ) : (real,real) where x.rank == 1 && y.rank == 1 {
   assert (xmin >= 0.0);         // xmin is a threshold
   var nplus = 0;
   var nminus = 0;
   var yplus = 0.0;
   var yminus = 0.0;
   for (ax,ay) in zip(x,y) do {
      if ax > xmin then {
         yplus += ay;
         nplus += 1;
      }
      else if ax < -xmin then {
         yminus += ay;
         nminus += 1;
      }
   }
   yplus /= nplus ;
   yminus /= nminus;
   return (yplus,yminus);
}


// -------------------------------------------------------------------
// --> lindetrend: linear detrending: extract the linearly varying
// "mean" inplace
// -------------------------------------------------------------------
proc lindetrend(ref x: [] real) where x.rank == 1 {
   var n = x.size;
// --------------------------------------------------------------------
// be careful with empty arrays
// --------------------------------------------------------------------
   assert (n > 0);
// --------------------------------------------------------------------
// make an array of times
// --------------------------------------------------------------------
   var t = linspace(0.0,(n-1):real,n);
// --------------------------------------------------------------------
// first the mean
// --------------------------------------------------------------------
   var (a,b,r) = reglin(t,x);
   var mx = a*t + b;          // the amazing vectorized operaton
   x = x - mx;                // another amazing vectorized operation
}


// -------------------------------------------------------------------
// --> uvwrot: rotate 3 arrays of sonic velocities
// -------------------------------------------------------------------
use Math only atan2;
proc uvwrot(ref au: [] real, ref av: [] real, ref aw: [] real) {
   assert(au.rank == 1);
   assert(av.rank == 1);
   assert(aw.rank == 1);
   const n = au.size;
   assert (av.size == n);
   assert (aw.size == n);
   ref u = au.reindex(1..n);
   ref v = av.reindex(1..n);
   ref w = aw.reindex(1..n);
   const ubar = stat1(u);
   const vbar = stat1(v);
   const wbar = stat1(w);
   const alphax = atan2(vbar,ubar);
   const alphaz = atan2(wbar,sqrt(ubar**2 + vbar**2));
   var CC: [1..3,1..3] real;
   CC[1,..] = [ cos(alphax)*cos(alphaz),  sin(alphax)*cos(alphaz), sin(alphaz)];
   CC[2,..] = [-sin(alphax),              cos(alphax),             0.0        ];
   CC[3,..] = [-cos(alphax)*sin(alphaz), -sin(alphax)*sin(alphaz), cos(alphaz)];
   foreach k in 1..n do {
      var ax = [u[k],v[k],w[k]];
      var ay: [1..3] real;
      dot_mv(CC,ax,ay);
      u[k] = ay[1];
      v[k] = ay[2];
      w[k] = ay[3];
   }
}

// -------------------------------------------------------------------
// --> mdfast: calculates deviations around a running median AND
// flags spikes
//
//   m    - filter width
//   ndev - | x - xmedian | > ndev*mad flags a spike
//   x    - array of data of length n (x[1] .. x[n])
//
//   Nelson Luís Dias
//   2019-11-26T08:38:24 (I was wrong: simpler is always better!)
//   2021-03-19T09:20:05 finishing translation from Python
// -------------------------------------------------------------------
proc mdfast(
   const in m: int,           // size of each block (filter width)
   const in ndev:real,        // deviation around median defines spike
   ref ax: [] real,           // the data
   out xctrl: (real,real,real,real),    // see maxxmdn, ..., below
   out nspikes: int = 0       // number of spikes found
) {
   var n = ax.size;           // the size of ax
   if ! ( n % m == 0 ) then { // n must be multiple of m
      halt("array length must be a multiple of m");
   }
// -------------------------------------------------------------------
// painful reindexing
// -------------------------------------------------------------------
   assert (ax.rank == 1);
   ref x = ax.reindex(1..n);
// -------------------------------------------------------------------
// the elements of xctrl
// -------------------------------------------------------------------
   var maxxmdn: real;         // max median over blocks of size m
   var minxmdn: real;         // min median over blocks of size m
   var maxxmad: real;         // max mad over blocks of size m
   var minxmad: real;         // min mad over blocks of size m
// -------------------------------------------------------------------
// the obvious limits to calculate maxima and minima
// -------------------------------------------------------------------
   maxxmdn = min(real);
   minxmdn = max(real);
   maxxmad = min(real);
   minxmad = max(real);
// -------------------------------------------------------------------
// all passes
// -------------------------------------------------------------------
   for k in 1..n by m do {
      ref xtest = x[k..#m];        // a sub-block
      var xmedian = nanmedian(xtest); // and its median
      assert(!isnan(xmedian));
// -------------------------------------------------------------------
// now calculates the mean absolute deviation around the median for
// the sub-block, which is xmad
// -------------------------------------------------------------------
      var adev = abs(xtest - xmedian);
      var (mval, xmad)  = nanstat1(adev);
// -------------------------------------------------------------------
// the filter's maxima and minima may be useful
// -------------------------------------------------------------------
      if xmedian > maxxmdn then {
         maxxmdn = xmedian;
      }
      if xmedian < minxmdn then {
         minxmdn = xmedian;
      }
      if xmad > maxxmad then {
         maxxmad = xmad;
      }
      if xmad < minxmad then {
         minxmad = xmad;
      }
// -------------------------------------------------------------------
// find and censor with nan all spikes in the sub-block
// -----------------------------------------------------------------------------
      for i in k..#m do {
         if abs(x[i] - xmedian) > ndev*xmad then {
            x[i] = nan;
            nspikes += 1;
         }
      }
   }
   xctrl = (maxxmdn,minxmdn,maxxmad,minxmad);
   return;
}

// -----------------------------------------------------------------------------
// kfastcontrol: fast control for korea
//
// the ideia is to implement something simple, along the lines of
//
// E. Zahn, T. L. Chor, N. L. Dias. A Simple Methodology for Quality
// Control of Micrometeorological Datasets. American Journal of
// Environmental Engineering 2016, 6(4A): 135-142 DOI:
// 10.5923/s.ajee.201601.20,
//
// return codes (rc)
//
// 0 => success (all nans have been interpolated)
// 1 => more than 1% nans in x
// 2 => more than 1% nans in x after despiking
// 3 => largest mad is too small
// 4 => difference between largest and smallest block medians is too 
//      large 
// -----------------------------------------------------------------------------
proc kfastcontrol(
   const in m: int,           // size of each block   
   const in ndev: real,       // deviation around median defines spike
   const in madmin: real,     // minimum acceptable mad
   const in maxdiff: real,    // maximum difference between block medians
   ref ax: [] real,           // the data block
   out rc: int,               // the return code
   out nspikes: int,          // how many spikes?
   out xctrl: (real,real,real,real)     // directly from mdfast
) {
// -----------------------------------------------------------------------------
// take a look at x
// -----------------------------------------------------------------------------   
   assert (ax.rank == 1);               // must be 1D
   var nx = ax.size;                    // count elements
   ref x = ax.reindex(1..nx);           // reindex
   var ninv = countval(nan,x);          // number of nans (invalids) in x
   var fninv = ninv:real;               // same, as real
   if fninv/nx > 0.01 then {  // if there are more than 1% nans, fail
      x = nan;
      rc = 1;
      nspikes = 0;
      xctrl = (nan,nan,nan,nan);
      return;
   }
// -----------------------------------------------------------------------------
// de-spiking and its consequences
// -----------------------------------------------------------------------------   
   var                        // returns from mdfast: 
      maxxmdn,                // maximum median
      minxmdn,                // minimum median
      maxxmad,                // maximum mad
      minxmad:                // minimum mad
      real;
   mdfast(m,ndev,x,xctrl,nspikes); // finds spikes over sub-block medians
   (maxxmdn,minxmdn,maxxmad,minxmad) = xctrl;
   ninv += nspikes;           // add # of spikes to ninv (# of invalids)
   fninv = ninv:real;         // convert again
   if fninv/nx > 0.01 then {  // if over 1%, it is game over
      x = nan;
      rc = 2;
      return;
   }
   if maxxmad < madmin then { // test for xmads too small
      x = nan;
      rc = 3;
      return;
   }
   if ( maxxmdn - minxmdn ) > maxdiff then { // maxdiff test
      x = nan;
      rc = 4;
      return;
   }
// -------------------------------------------------------------------
// if passed all tests and there are no gaps, returns 0
// -------------------------------------------------------------------
   if ninv == 0 then {
      rc = 0;
      return;
   }
// -------------------------------------------------------------------
// if passed everything but there are still gaps, fills each gap with
// linear interpolation. first, find each gap:
// -------------------------------------------------------------------
   var db = {0..nx+1};
   var binv: [db] bool;       // true when nan in x 
   binv[0] = false;           // place sentinels around binv
   binv[1..nx] = isnan(x);    // and find nans inside
   binv[nx+1] = false;        // place sentinels around binv
   var dinv = diff(binv);     // differentiate binv
// -------------------------------------------------------------------
// now +1 marks the beginning of gaps, and -1 their end, in dinv. the
// number of -1s is equal to the number of +1s. ibeg and iend are the
// indices of the last valid datum before each gap and first valid
// datum after each gap in 1-based array x
// -------------------------------------------------------------------
   var ibeg = whereval(1,dinv);         // note that ibeg is 0-based
   var iend = whereval(-1,dinv) + 1;    // note that iend is 0-based
   var nruns = ibeg.size;               // how many gaps?
   assert (nruns == iend.size);
// -------------------------------------------------------------------
// finally, linear interpolation of gaps: we place two sentinels with
// the mean at both ends of x. Since there are nans in x, we need to
// use nanmean
// -------------------------------------------------------------------
   var (nnans, xmean): 2*real = nanstat1(x);
// -------------------------------------------------------------------
// xcat acts like x with two sentinels at the extremeties
// -------------------------------------------------------------------
   var xcat: [db] real;
   xcat[0] = xmean;
   xcat[1..nx] = x;
   xcat[nx+1] = xmean;
// -------------------------------------------------------------------
// (parallel!) loop over gaps: linear interpolation with linspace
// -------------------------------------------------------------------
   forall ir in 0..#nruns do {
      var irb = ibeg[ir];          // last valid position before
      var ire = iend[ir];          // first valid position after
      var xstart = xcat[irb];      // valid datum before
      var xstop  = xcat[ire];      // valid datum after
      var gaplen = ire - irb + 1;  // gap size
      var xfill = linspace(xstart,xstop,gaplen);  // linear interp
//      writeln(irb, " ",ire, " ",xstart," ",xstop," ",gaplen,"  ",x[irb+1..ire-1]," ",xfill);
      x[irb+1..ire-1] = xfill[1..gaplen-2];       // fill gaps
   }
// -----------------------------------------------------------------------------
// all is fine: passed quality control
// -----------------------------------------------------------------------------
   rc = 0;
   return;
}
// -----------------------------------------------------------------------------
// --> csat3shadow: correction of CSAT3 measurements for transducer shadowing
// -----------------------------------------------------------------------------
use Math only acos, pi, cos, sin;
proc csat3shadow(ref au: [] real, ref av: [] real, ref aw: [] real, const in report = false)
   where (au.rank == 1 && av.rank == 1 && aw.rank == 1) {
   const n = au.size;
   assert (av.size == n);
   assert (aw.size == n);
   ref u = au.reindex(1..n);
   ref v = av.reindex(1..n);
   ref w = aw.reindex(1..n);
   const phi = pi/3;
   const cphi = cos(phi);
   const sphi = sin(phi);
   const sqr3 = sqrt(3.0);
   foreach i in 1..n do {
// -----------------------------------------------------------------------------
// corrected cartesian velocities
// -----------------------------------------------------------------------------      
      var (uold, vold, wold) = (u[i],v[i],w[i]);            
      var (unew, vnew, wnew) = (uold,vold,wold);            
// -----------------------------------------------------------------------------
// measured sonic path velocities
// -----------------------------------------------------------------------------      
      const uam  = -uold*cphi + wold*sphi;                    
      const ubm  = ((uold + sqr3*vold)/2.0)*cphi + wold*sphi; 
      const ucm  = ((uold - sqr3*vold)/2.0)*cphi + wold*sphi; 
// -----------------------------------------------------------------------------
// corrected sonic path velocities
// -----------------------------------------------------------------------------      
      var (ua,ub,uc) = (uam,ubm,ucm);                  // corrected sonic
// -----------------------------------------------------------------------------
// absolute differences cartesian
// -----------------------------------------------------------------------------      
      var du, dv, dw : real = 0.1;                     
      if report then writeln("-"*40);
      while max(du,dv,dw) >= 0.001 do {
// -----------------------------------------------------------------------------
// update wind speed
// -----------------------------------------------------------------------------         
         var S = sqrt(uold**2 + vold**2 + wold**2);
// -----------------------------------------------------------------------------
// update cosine angles
// -----------------------------------------------------------------------------
         var theta_a = acos(ua/S);
         var theta_b = acos(ub/S);
         var theta_c = acos(uc/S);
// -----------------------------------------------------------------------------
// update path velocities
// -----------------------------------------------------------------------------
         ua = uam/(0.84 + 0.16*sin(theta_a));
         ub = ubm/(0.84 + 0.16*sin(theta_b));
         uc = ucm/(0.84 + 0.16*sin(theta_c));
// -----------------------------------------------------------------------------
// update cartesian velocities
// -----------------------------------------------------------------------------
         unew = (-2*ua + ub + uc)/(3*cphi);
         vnew = (ub - uc)/(sqr3*cphi);
         wnew = (ua + ub + uc)/(3*sphi);
// -----------------------------------------------------------------------------
// absolute values of differences between old and new
// -----------------------------------------------------------------------------         
         du = abs(unew - uold);
         dv = abs(vnew - vold);
         dw = abs(wnew - wold);
         if report then {
            writef(" %7.4dr"*3,uold,vold,wold);
            writef(" --- ");
            writef(" %7.4dr"*3,unew,vnew,wnew);
            writef(" --- ");
            writef(" %7.4dr"*2+"\n",max(du,dv,dw),S);
         }
         (uold,vold,wold) = (unew,vnew,wnew);
      }
      u[i] = unew;
      v[i] = vnew;
      w[i] = wnew;      
   }
}

// =====================================================================
// instat: instationarity test (Th. Foken & B.Wichura, 1995)
// =====================================================================

proc instat(
   const in nb: int,
   ref ax: [] real,           // the data block
   ref ay: [] real,
   out stst: real
) {
   
   if ax.rank != ay.rank then {
      return;
   }
   if ax.size != ay.size then {
      return;
   }
   
   var n = ax.size;           // count elements
   ref x = ax.reindex(1..n);  // reindex
   ref y = ay.reindex(1..n);  // reindex
   
   var cxyb: [1..nb] real;
   var cxyglob: real;
   var cxybm: real;
   var blksize = n/nb;
   
   var xm = stat1(x);
   var ym = stat1(y);
   
   for k in 1..nb do {
      var bgn = (k-1)*blksize+1 :int;
      var end = k*blksize :int;
      var xbm = stat1(x[bgn..end]): real;
      var ybm = stat1(y[bgn..end]): real;
      
      cxyb[k] = covar(xbm, ybm, x[bgn..end], y[bgn..end]); // create a block covars array
   } 
   cxybm = stat1(cxyb);
   cxyglob = covar(xm,ym,x,y);            // global covar
   stst = (abs((cxyglob-cxybm)/cxyglob)); // steady-state test
   
   return;
}
// -----------------------------------------------------------------------------
// --> arm: calculates the autorecursive mean from an array x as
//      
// m[i] = (1/l)*x[i] + ( 1 - (1/l))*m[i-1]
//   
// where l is the window width of the running mean.  An initial running mean
// m[0] is calculated as equal to x[0].
// -----------------------------------------------------------------------------
proc arm (
   const in l: int,           // window for calculation of running means
   const ref ax: [] real,     // data array 
   ref am: [] real            // autorecursive mean 
   ) where (ax.rank == 1) && (am.rank == 1) {
   const n = ax.size;
   assert( am.size == n);
   const fl: real = l;
   var   sum: real = 0.0;
   ref x = ax.reindex(0..n-1);
   ref m = am.reindex(0..n-1);
// -----------------------------------------------------------------------------
// obtains the first autorecursive mean
// -----------------------------------------------------------------------------
   m[0] = x[0];
   sum = x[0] * fl ; 
// -----------------------------------------------------------------------------
// loop to calculate rest of the arm vector
// -----------------------------------------------------------------------------
   for i in 1..n-1 do {
      sum  = sum + x[i] - m[i-1];
      m[i] = sum/fl;
   }
}
