// -*- chapel_mode -*-
// =============================================================================
// ==> nstat: useful statistics
//
// 2021-06-10T11:24:23 parallelizing stat1
// 2021-08-27T14:46:09 introducing gmom & friends
// 2022-02-17T14:44:24 fixing absolute values in gmom
// 2022-05-07T11:55:01 adding performance
// 2022-07-09T15:41:21 adding lowess_fast_estimate
// 2022-09-16T11:30:34 changing input data to "numeric"
// 2023-03-08T16:00:30 including steep
// =============================================================================
use smatrix;
use ssr only allequalto, amin, amax, aminz, interp, heapsort, indxsort, indxquickselect;
proc stat1(x: [] real): real {
// -----------------------------------------------------------------------------
// be careful with empty arrays
// -----------------------------------------------------------------------------
   assert (x.rank == 1);
   var n = x.size;
   if n == 0 then {
      halt("nstat-->stat1: empty array");
   }
   var fn = n:real; // just in case
// -----------------------------------------------------------------------------
// first the mean
// -----------------------------------------------------------------------------
   var xm = 0.0;
   xm = (+ reduce x);
   xm /= fn;
   return xm;
}
// -----------------------------------------------------------------------------
// --> nanstat1: calculate 1st moment, ignoring NANs
//
// 2021-03-20T12:31:37 
// -----------------------------------------------------------------------------
proc nanstat1(ref x: [] real): (int,real) {
   assert (x.rank == 1);
   const n = x.size;
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->nanstat1: empty array");
   }
   var m = 0;                 // how many not NANs?
   var xm = 0.0;              // the mean over them
   for xi in x do {
      if !isnan(xi) then {
         m += 1;
         xm += xi;
      }
   }
   if m == 0 then {           // is x all NANs?
      return (0,NAN);
   }
   xm /= m;                   // calculate mean over valid data
   return (m,xm);
}

// -----------------------------------------------------------------------------
// --> nanstat2: calculate 1st and 2nd moments, ignoring NANs
//
// 2021-04-20T15:35:44
// -----------------------------------------------------------------------------
proc nanstat2(ref x: [] real): (int,real,real) {
   assert (x.rank == 1);
   const n = x.size;
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->nanstat1: empty array");
   }
   var m = 0;                 // how many not NANs?
   var xm = 0.0;              // the mean over them
   for xi in x do {
      if !isnan(xi) then {
         m += 1;
         xm += xi;
      }
   }
   if m == 0 then {           // is x all NANs?
      return (0,NAN,NAN);
   }
   xm /= m;                   // calculate mean over valid data
   var xv = 0.0;              // the variance
   for xi in x do {
      if !isnan(xi) then {
         xv += (xi - xm)**2;
      }
   }
   xv /= m;
   return (m,xm,xv);
}

// -----------------------------------------------------------------------------
// --> wstat1: calculate weighted 1st moment
//
// 2021-03-26T13:46:25 w can be of generic type!
// 
// 2021-04-15T15:08:47 eliminating the analysis of NANs
// -----------------------------------------------------------------------------
proc wstat1(const ref x: [] real, const ref w: [] ?tw): real where (x.rank == 1 && w.rank == 1) {
   const n = x.size;
   assert (w.size == n);
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->wstat1: empty array");
   }
   var sw: tw = 0;            // the sum of the weights
   var sx = 0.0;              // the sum of the elements
   for (wi,xi) in zip(w,x) do {
      sw += wi;
      sx += wi*xi;
   }
   var xm = sx/sw;            // weighted mean 
   return xm; 
}

// -----------------------------------------------------------------------------
// --> wstat2: calculate weighted mean and variance
// -----------------------------------------------------------------------------
proc wstat2(const ref x: [] real, const ref w: [] ?tw): (real,real)
   where (x.rank == 1 && w.rank == 1) {
   const n = x.size;
   assert (w.size == n);
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->wstat2: empty array");
   }
   var sw: tw = 0;            // the sum of the weights
   var sx = 0.0;              // the sum of the elements
   for (wi,xi) in zip(w,x) do {
      sw += wi;
      sx += wi*xi;
   }
   var xm = sx/sw;            // weighted mean
   var v2x = 0.0;   
   for (wi,xi) in zip(w,x) do {
      v2x += wi*(xi-xm)**2;
   }
   v2x /= sw;                 // weighted variance
   return (xm,v2x);
}


// -----------------------------------------------------------------------------
// --> nanstat1: calculate weighted 1st moment, ignoring NANs
//
// 2021-03-26T13:46:25 w can be of generic type!
// -----------------------------------------------------------------------------
proc nanwstat1(ref x: [] real, ref w: [] ?tw): (int,real) {
   assert (x.rank == 1);
   const n = x.size;
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->nanwstat1: empty array");
   }
   var m = 0;                 // how many not NANs?
   var sw: tw = 0;            // the sum of the weights
   var sx = 0.0;              // the sum of the elements
   for (wi,xi) in zip(w,x) do {
      if !isnan(xi) then {
         m += 1;
         sw += wi;
         sx += wi*xi;
      }
   }
   if m == 0 then {           // is x all NANs?
      return (0,NAN);
   }
   var xm = sx/sw;            // weighted mean over valid data
   return (m,xm);
}
// -----------------------------------------------------------------------------
// --> gmean: the geometric mean, in very simple terms
// -----------------------------------------------------------------------------
proc gmean(
   ref x: [] real                  // the array being averaged
) : real where x.rank == 1 {
   var lsum = 0.0;                 // the sum of ln(x[i])
   for xi in x do {
      if xi <= 0.0 then {          // only allow positives
         halt("array x must hold only positive values");
      }
      lsum += log(xi);
   }
   var n = x.size;
   lsum /= n;
   var gm  = exp(lsum);
   return gm;
}
// -----------------------------------------------------------------------------
// --> log-modulus function sgn(x)*log(|x| + 1)
// -----------------------------------------------------------------------------
proc logmod(
   const in x: real
   ): real {
   return sgn(x)*log(abs(x)+1);
}
// -----------------------------------------------------------------------------
// --> exp-modulus function sgn(x)*(exp(|x|) - 1)
// -----------------------------------------------------------------------------
proc expmod(
   const in x: real
   ): real {
   return sgn(x)*(exp(abs(x))-1.0);
}
// -----------------------------------------------------------------------------
// --> gmom: geometric-modulus mean: calculates
//
// expmod( (sum logmod(x[i]))/n )
// -----------------------------------------------------------------------------
proc gmom(
   ref x: [] real                  // the array being averaged
): real where x.rank == 1 {
   var ax = abs(x);                // absolute values!
   var delta = aminz(ax)/16;       // minimum non-zero absolute value
   ax = x/delta;                   // change scale to avoid values close to 1
   var lsum = 0.0;
   for i in ax.domain do {         // mod-geometric mean: sum of L(x)
      lsum += logmod(ax[i]);
   }
   var n = x.size;
   lsum /= n;                 // mean(LogMod(x))
   var gm = expmod(lsum);     // ExpMod(mean(LogMod(x))
   gm *= delta;
   return gm;
}

// -----------------------------------------------------------------------------
// --> gmeanpn: calculates the geometric means of positives and negatives
// (with |x_i|) and adds them (with sign): gp - gn. This calculation follows
//
// Habib, Elsayed AE, (2012) Geometric mean for negative and zero values,
// International Journal of Research and Reviews in Applied Sciences}, 11(3)419--432.
// journal={
// File = {:/home/nldias/Dropbox/epapers/habib2012geometric.pdf},
// volume={11},
// number={3},
// pages={419--32},
// year={2012}

// -----------------------------------------------------------------------------
// proc gmeanpn(
//    ref x: [] real                  // the array being averaged
// ): real where x.rank == 1 {
//    const n = x.size;
//    const fn = n:real;
//    var gp = 0.0;
//    var np = 0;
//    var gn = 0.0;
//    var nn = 0;
//    for xi in x do {
//       if xi == 0.0 then {
//          continue;
//       }
//       else if xi < 0.0 then {
//          gn += log(abs(xi));
//          nn += 1;
//       }
//       else {
//          gp += log(xi);
//          np += 1;
//       }
//    }
//    gp = (np/fn)*exp(gp/fn);
//    gn = (nn/fn)*exp(gn/fn);
//    return gp - gn;
// }


proc stat2(x: []): (real,real) {
// -----------------------------------------------------------------------------
// be careful with empty arrays
// -----------------------------------------------------------------------------
   assert (x.rank == 1); // 
   var n = x.size;
   if n == 0 then {
      halt("--> stat2: array is empty");
   }
   var fn = n:real; // just in case
// -----------------------------------------------------------------------------
// first the mean
// -----------------------------------------------------------------------------
   var xm = 0.0;
   for xi in x do {
      xm += xi;
   }
   xm /= fn;
// -----------------------------------------------------------------------------
// now the variance
// -----------------------------------------------------------------------------
   var xv = 0.0;
   for xi in x do {
      xv += (xi - xm)**2;
   }
   xv /= fn;
   return (xm,xv);
}
// -----------------------------------------------------------------------------
// --> covar: calculate the covariance
// -----------------------------------------------------------------------------
proc covar(xmed: real, ymed: real, ax: [] real, ay: [] real): real {
   assert (ax.rank == 1);
   assert (ay.rank == 1);
   var n = ax.size;
   if n == 0 then {
      return NAN;
   }
   assert (n == ay.size);
   ref x = ax.reindex(1..n);
   ref y = ay.reindex(1..n);
   var fn = n:real; // just in case
   var sxy = 0.0;
   for i in 1..n do {
      sxy += (x[i] - xmed)*(y[i] - ymed);
   }
   return sxy/fn;
}
// -----------------------------------------------------------------------------
// --> trivar: calculate the trivariance
// -----------------------------------------------------------------------------
proc trivar(xmed: real, ymed: real, zmed: real, ax: [] real, ay: [] real, az: [] real): real
   where (ax.rank == 1) && (ay.rank == 1) && (az.rank == 1) {
   var n = ax.size;
   if n == 0 then {
      return NAN;
   }
   assert (n == ay.size);
   assert (n == az.size);
   ref x = ax.reindex(1..n);
   ref y = ay.reindex(1..n);
   ref z = az.reindex(1..n);
   var fn = n:real; // just in case
   var sxyz = 0.0;
   for i in 1..n do {
      sxyz += (x[i] - xmed)*(y[i] - ymed)*(z[i]-zmed);
   }
   return sxyz/fn;
}
// -------------------------------------------------------------------
// --> median: returns the median of array x
//
// 2012-08-21T08:55:53 Python version
// 2021-03-19T09:37:33 Chapel version
// -------------------------------------------------------------------
proc median(ref ax: [] ?at): at {
   assert (ax.rank == 1);
   var n = ax.size;
   ref x = ax.reindex(0..n-1);
   if n == 0 then {
      halt("nstat-->median: empty array");
   }
   var indx: [0..n-1] int;
   indxsort(x,indx);          // sort x by index
   var xmedian: at;
// -------------------------------------------------------------------
//  is n even or odd ?
//  ------------------------------------------------------------------
   if n == 1 then {
      xmedian = x[0];
   }
   else if ( (n % 2) == 0 ) then {
      xmedian = (x[indx[n/2 - 1]] + x[indx[n/2]]) / 2;
   }
   else {
      xmedian = x[indx[n/2]];
   }
   return xmedian;
}
proc nanmedian(ref ax: [] ?at): at {
   assert (ax.rank == 1);
   var x = purgeval(NAN,ax);
   var n = x.size;
   if n == 0 then {
      halt("nstat-->median: empty array");
   }
   var indx: [0..n-1] int;
   indxsort(x,indx);          // sort x by index
   var xmedian: at;
// -------------------------------------------------------------------
//  is n even or odd ?
//  ------------------------------------------------------------------
   if n == 1 then {
      xmedian = x[0];
   }
   else if ( (n % 2) == 0 ) then {
      xmedian = (x[indx[n/2 - 1]] + x[indx[n/2]]) / 2;
   }
   else {
      xmedian = x[indx[n/2]];
   }
   return xmedian;
}
// -----------------------------------------------------------------------------
// --> quartiles: returns the quartiles (q1,q2,q3) of an array; q2 is the median
// -----------------------------------------------------------------------------
proc quartiles(
   in ax: [] real                   // we need to make a copy!
   ): (real,real,real) {
   assert (ax.rank == 1);
   var n = ax.size;
   if n < 3 then {
      halt("nstat--> quartiles: need at least 3 data points...");
   }
   ref x = ax.reindex(0..n-1);
   var fn = [i in 0..n-1] (i:real)/(n-1);
   heapsort(x);                    // sort x: this only affects the local copy
   var q1 = interp(0.25,fn,x);     // just interpolate the 3 quartiles
   var q2 = interp(0.50,fn,x);
   var q3 = interp(0.75,fn,x);
   return (q1,q2,q3);
}
// -----------------------------------------------------------------------------
// --> whiskers: returns the whiskers limits [w0,w1,w2,w3,w4] of an array; w0 is
// the minimum, w1 is the first quartile, w2 is the median, w3 is the third
// quartile and w4 is the maximum
// -----------------------------------------------------------------------------
proc whiskers(
   in x: [?Dx] real                   // we need to make a copy!
   ): [0..4] real where x.rank == 1 {
   var n = x.size;
   writeln("whiskers: Dx = ",Dx);
   if n < 3 then {
      halt("nstat-->whiskers: need at least 3 data points...");
   }
   var fn = [i in 0..n-1] (i:real)/(n-1);
   heapsort(x);                    // sort x: this only affects the local copy
   var whisk: [0..4] real;
   whisk[0] = x[Dx.low];        // the minimum
   whisk[1] = interp(0.25,fn,x);   // just interpolate the 3 quartiles
   whisk[2] = interp(0.50,fn,x);
   whisk[3] = interp(0.75,fn,x);
   whisk[4] = x[Dx.high];       // the maximum
   return whisk;
}
// -----------------------------------------------------------------------------
// type vla is a variable-length array that is needed for xbins
// -----------------------------------------------------------------------------
record vla {
   var db: domain(1);
   var a: [db] real;
   proc this(k: int) ref {
      assert (db.contains(k));
      return a[k];
   }
}
// -----------------------------------------------------------------------------
// --> xbins: given the arrays x,y of data points, and the bin limits xbmin, db,
// xbmax, such that xb[j] = xbmin + (j+0.5)*db, returns yb where yb[j] is a
// (ragged) array that contains all data points y[i] such that xb[j]-dx <= x[i]
// < xb[j]+dx
// -----------------------------------------------------------------------------
proc xbins(
   const in xbmin: real,
   const in dxb: real,
   const in xbmax: real,   
   ref x: [] real,
   ref y: [] real
   ): [] vla where (x.rank == 1) && (y.rank == 1) {
   const n = x.size;                              // the size of x
   assert(n == y.size);
   const m = ((xbmax - xbmin)/dxb):int;       // the number of bins
   writeln("# of xbins: m = ",m);
// -----------------------------------------------------------------------------
// it is complicated to construct a ragged array
// -----------------------------------------------------------------------------
//   const dxb2 = dxb/2;
   var nb: [0..m-1] int = 0;                 // nb[j] is the size of yb[j]
   var yb: [0..m-1] vla;                     // yb is a ragged array
// -----------------------------------------------------------------------------
// initialize the domains in yb
// -----------------------------------------------------------------------------
   for j in 0..m-1 do {
      yb[j] = new vla(db = {0..9});
   }
   for (xi,yi) in zip(x,y) do {
      var j = ((xi - xbmin)/dxb):int;    // index of desired bin
// -----------------------------------------------------------------------------      
// if j < 0 || j >= m then continue;
// -----------------------------------------------------------------------------      
      if j*(m-1-j) < 0 then {
         continue;
      }
      var k = nb[j];
      dgrow(k,yb[j].db);
      yb[j][k] = yi;                       // fill jth bin
      nb[j] += 1;
   }
// -----------------------------------------------------------------------------
// now resize the second dimension of yb
// -----------------------------------------------------------------------------
   for j in 0..m-1 do {
      yb[j].db = {0..nb[j]-1};
   }
   return yb;
}

// -----------------------------------------------------------------------------
// --> outlimits: define limits for outliers using the interquartile range
// -----------------------------------------------------------------------------
proc outlimits(
   ref x: [] ?tx,
   in delta: real = 1.5
   ): (real,real) {
   assert (delta >= 0.0);
   assert (x.rank == 1);
//   writeln("x.size = ",x.size);
   var (q1,q2,q3) = quartiles(x);
   var iqr = q3 - q1;
   var minx = q1 - delta*iqr;
   var maxx = q3 + delta*iqr;
   return (minx,maxx);
}
// -----------------------------------------------------------------------------
// --> force_bounds is a convenient macro
// -----------------------------------------------------------------------------
inline proc force_bounds(
   ref x: real,
   ref xmin: real,
   ref xmax: real) {
   if x < xmin then {
      x = xmin;
   }
   else if x > xmax then {
      x = xmax;
   }
}
// -----------------------------------------------------------------------------
// --> reglin: calculates the linear regression y = ax + b, plus the
// correlation coefficient r
//
// Nelson Luís Dias
// 2010-02-11T15:23:42
// 2010-02-11T15:23:45
// 2019-07-12T11:20:43 prevent division by zero
// translating to Python
// ----------------------------------------------------------------------------
proc reglin(ref x: [] real,ref y: [] real): 3*real {
   assert (x.rank == 1);
   assert (y.rank == 1);
   var n = x.size;
   if n == 0 then {
      halt("nstat-->reglin: array x is empty");
   }
   assert (y.size == n);
   var (xavg,xvar) = stat2(x);
   var (yavg,yvar) = stat2(y);
   if ( xvar == 0.0 ) || ( yvar == 0.0) then {
      halt("nstat-->reglin: zero x or y variance");
   }
   var coxy = covar(xavg,yavg,x,y);
   var a = coxy / xvar;
   var b = yavg - (a * xavg);
   var r = coxy / sqrt(xvar*yvar);
   return (a,b,r);
}
// -----------------------------------------------------------------------------
// --> reglina: calculates the linear regression y = a x , the correlation
// coefficient r, the standard error of estimate se, and the standard deviation
// sa the estimator of a
// -----------------------------------------------------------------------------
proc reglina(
   ref x: [] real,
   ref y: [] real,
   out a: real,
   out r2: real,              // only r2 makes sense here!
   out se: real,              // std error of estimate?
   out sa: real
) {
   assert(x.rank == 1);
   assert(y.rank == 1);
   var n = x.size;
   var fn = n:real;
// ---------------------------------------------------------------------------------------
// starts by calculating central moments
// ---------------------------------------------------------------------------------------
   var (xavg,xvar) = stat2(x);
   var (yavg,yvar) = stat2(y);
   var coxy = covar(xavg,yavg,x,y);
// ----------------------------------------------------------------------------------------
// translates to non-central moments
// ----------------------------------------------------------------------------------------
   var sx20 = fn * ( xvar + (xavg*xavg) ) ;
   var sy20 = fn * ( yvar + (yavg*yavg) ) ;
   var sxy0 = fn * ( coxy + (xavg*yavg) ) ;
// ----------------------------------------------------------------------------------------
// obtains the slope
// ----------------------------------------------------------------------------------------
   a = sxy0 / sx20 ;
// ----------------------------------------------------------------------------------------
// 2007-09-21T01:40 each time I look at this I get a different result; at least this time
// I have documented the equation in stat.tex
// ----------------------------------------------------------------------------------------
   var s2 = ( sy20 - 2.0*a*sxy0 + a*a*sx20 ) / (fn-1.0)  ;
   sa = sqrt( s2 / sx20 ) ;
   se = sqrt(s2) ;
   r2 = 1.0 - s2/yvar;
//   printf(" s2 = %lf  yvar = %lf\n", s2, yvar) ;
} 
// -----------------------------------------------------------------------------
// --> reglinpar: calculates the linear regression y = ax + b, the correlation
// coefficient r, the standard deviation s of the distribution of y given x, and
// the standard deviations sa, sb of the estimators(?) of a and b
// -----------------------------------------------------------------------------
proc reglinpar(
   ref x: [] real,
   ref y: [] real,
   out a: real,
   out b: real,
   out r: real,
   out s: real,
   out sa: real,
   out sb: real
) {
   assert(x.rank == 1);
   assert(y.rank == 1);
   var n = x.size;
   assert(y.size == n);
   var fn = n:real;
   var dn = fn / ( fn - 1.0);
   var en = fn / ( fn - 2.0);
   var (xavg,xvar) = stat2(x);
   var (yavg,yvar) = stat2(y);
   var coxy = covar(xavg,yavg,x,y);
   xvar *= dn ;
   yvar *= dn ;
   coxy *= dn ;
   a = coxy / xvar ;
   b = yavg - a*xavg; 
   r = coxy / sqrt(xvar*yvar);
   var s2 = en*(1 - r**2)*yvar   ;
   s = sqrt(s2) ;
   sa = sqrt( s2 / (fn * xvar ) ) ;
   sb = sqrt( s2 * ( 1 + (xavg*xavg) / xvar ) / fn ) ;
}
// -----------------------------------------------------------------------------
// --> stat_ols: ordinary least squares multivariate regression: quite 
//     ordinary, indeed: wants to find the parameters C[m,1] such that:
//    
//     y[n,1] = x[n,m] c[m,1] =>
//     x'y = (x'x) c = a c    =>
//     c = a^-1 x' y          =>
//     c = b y                //
//
//     the variance of the residuals is
//
//     s^2 = [y - ye]'[y - ye]/(n-m)
//
//     and the covariance matrix of c is
//
//     cov{c} = s^2 a^-1
//
//     the vector of estimated y's returns:
//    
//     ye = x c
// -----------------------------------------------------------------------------
proc ols(
   ref ax: [] real,       // independent data matrix
   ref ay: [] real,       // dependent data vector  
   out S2: real,          // variance of residuals
   out R2: real,          // coefficient of determination
   ref ac: [] real,       // array of coefficients
   ref aye: [] real,      // vector of estimated ys
   ref accov: [] real,    // covariance matrix of c
   ref axcor: [] real     // covariance matrix of x
) {
   var (n,m) = ax.shape;
   assert(ay.shape == (n,));
   assert(aye.shape == (n,));
   assert(ac.shape == (m,));
   assert(accov.shape == (m,m));
   assert(axcor.shape == (m,m));
   assert(n > m);
   ref x = ax.reindex({1..n,1..m});
   ref y = ay.reindex(1..n);
   ref c = ac.reindex(1..m);
   ref ccov = accov.reindex({1..m,1..m});
   ref ye = aye.reindex(1..n);
   ref xcor = axcor.reindex({1..m,1..m});
// -----------------------------------------------------------------------------
// auxiliary storage S,b, is needed
// -----------------------------------------------------------------------------
   var S: [1..m,1..m] real;
   var b: [1..m,1..n] real;
// -----------------------------------------------------------------------------
// calculates S = (x'x) 
// -----------------------------------------------------------------------------
   dot_mtm(x,x,S);
   for i in 1..m do {         // the correlation matrix
      for j in 1..m do {
         xcor[i,j] = S[i,j]/sqrt(S[i,i]*S[j,j]);
      }
   }
// -----------------------------------------------------------------------------
// calculates S^-1
// -----------------------------------------------------------------------------
   minvgj(S);
// -----------------------------------------------------------------------------
// calculates b = S^-1 x'
// -----------------------------------------------------------------------------
   dot_mmt(S,x,b);
// -----------------------------------------------------------------------------
// calculates c = b y
// -----------------------------------------------------------------------------
   dot_mv(b,y,c);
// -----------------------------------------------------------------------------
// calculates ye = x c
// -----------------------------------------------------------------------------
   dot_mv(x,c,ye) ;
// -----------------------------------------------------------------------------
// variance of residuals
// -----------------------------------------------------------------------------
   S2 = 0.0;
   for i in 1..n do {
      S2 += (y[i] - ye[i])**2;
   }
   S2 /= (n - m);        // there you go
   ccov = S2*S;          // covariance matrix of c
// -----------------------------------------------------------------------------
// variance of observations
// -----------------------------------------------------------------------------
   var (ymed,yvar) = stat2(y);
   yvar *= n:real/(n-1.0);
   R2 = 1.0 - S2/yvar;
}
// /============================================================================
// ==> nlowess: my implementation of the lowess algorithm
//
// with most of the inspiration from
// https://towardsdatascience.com/loess-373d43b03564
//
// 2022-04-22T12:11:44 a new star is born
// 2022-04-22T17:34:40 essentially done ... hopefully
// 2022-04-23T18:54:06 moved into nstat
// /============================================================================
// -----------------------------------------------------------------------------
// I use tricubic as a scalar function only
// -----------------------------------------------------------------------------
private proc tricubic(const in x: real): real {
   var ax = abs(x) ;
   return if ax < 1.0 then (1-ax**3)**3 else 0.0 ;
}
// -----------------------------------------------------------------------------
// --> normalize_array: gess what: normalizes data
// -----------------------------------------------------------------------------
private proc normalize_array(
   ref a: [] real,            // raw data
   ref n_a: [] real           // normalized data
   ): (real,real)
   where (a.rank == 1) {
   assert (a.shape == n_a.shape);
   var a_min = amin(a);
   var a_max = amax(a);
   n_a = (a - a_min) / (a_max - a_min);
   return (a_min, a_max);
}
// -----------------------------------------------------------------------------
// private variables!
// -----------------------------------------------------------------------------
private var
   degree = 1;      // degree of linear regression of this lowess
private var
   window: int;     // window of this lowess
private var
   xx_min,          // minimum of xx data
   xx_max,          // maximum of xx data
   yy_min,          // minimum of yy data
   yy_max:          // maximum of yy data
   real;
private var
   ldom = {0..0};   // data domain initially empty
private var
   n_xx,            // normalized x data
   n_yy:            // normalized y data
   [ldom] real;
// -----------------------------------------------------------------------------
// --> init_lowess: normalize data arrays making local copies to n_xx and n_yy
// via normalize_array. For clarity, nlowess will use 1-based arrays throughout
// -----------------------------------------------------------------------------
proc init_lowess(
   ref xx: [] real,           // the raw x data
   ref yy: [] real,           // the raw y data
   in win: int,               // the window, once and for all, of this loess
   in deg: int = 1            // the degree of this loess
   )  where (xx.rank == 1) && (yy.rank == 1) {
   assert (xx.shape == yy.shape);
// -----------------------------------------------------------------------------
// set these global variables once and for all
// -----------------------------------------------------------------------------
   window = win;
   degree = deg;
   var nd: int = xx.size;
   ldom = {1..nd};            // allocate mem for n_xx and n_yy
   (xx_min,xx_max) = normalize_array(xx,n_xx);
   (yy_min,yy_max) = normalize_array(yy,n_yy);
   return;
}
// -----------------------------------------------------------------------------
// --> get_min_range: the indices of the window smallest
// -----------------------------------------------------------------------------
private proc get_min_range(
   ref distances: [] real,   // the local distances
   ref min_range: [] int,
   out maxdist: real
   ) where distances.rank == 1 {      
   var dinds = distances.indices;
   var indx = [i in dinds] i ;
   indxquickselect(indx,distances,dinds.first,dinds.last,window);
   min_range =  [k in 1..window] indx[k];
   maxdist = distances[indx[window]];
}
// -----------------------------------------------------------------------------
// --> get_weights: I suspect that max_distance == distances[indx[wind]]
// because of partial sorting
// -----------------------------------------------------------------------------
private proc get_weights(
   in maxdist: real,
   ref distances: [] real,
   ref min_range: [] int,
   ref weights: [] real
   ) where distances.rank == 1 && min_range.rank == 1 {
   assert(min_range.shape == weights.shape);
   var wind = min_range.size;
// -----------------------------------------------------------------------------
// normalize with weights
// -----------------------------------------------------------------------------
   for k in 1..wind do {
      weights[k] = tricubic(distances[min_range[k]]/maxdist);
   }
}
// -----------------------------------------------------------------------------
// we will need to come back to normalize
// -----------------------------------------------------------------------------
private inline proc normalize_x(
   const in x: real        // the data array?
   ): real {
   return (x - xx_min)/(xx_max - xx_min);
}

private inline proc denormalize_y(
   const in y: real) : real {
   return y * (yy_max - yy_min) + yy_min;
}
// -----------------------------------------------------------------------------
// --> estimate is the core loess algorithm
// -----------------------------------------------------------------------------
proc lowess_estimate(
   in x: real                 // point around which to do linear regression
   ) {
   const m: int = window;     // size of neighborhood == window
   const n: int = degree;     // degree of linear regression == degree
   var n_x = normalize_x(x);  // normalize independent variable
// -----------------------------------------------------------------------------
// do I really need to calculate all distances for each point?
// -----------------------------------------------------------------------------
   var nd = n_xx.size;
   var maxdist: real;
   var distances: [1..nd] real = abs(n_x - n_xx);
   var min_range: [1..window] int;
   var weights: [1..window] real;
   get_min_range(distances, min_range,maxdist);
   get_weights(maxdist,distances, min_range, weights);
// -----------------------------------------------------------------------------
// here I use m,n instead of window, degree
// estimate beta = [Xtr(n+1,m)W(m,m)X(m,n+1)]^(-1)Xtr(n+1,m)W(m,m)n_yy(m,1)
// -----------------------------------------------------------------------------      
   var W: [1..m,1..m] real = 0.0;
   var X: [1..m,0..n] real;
   var Y: [1..m] real;
// -----------------------------------------------------------------------------
// fill diagonal of W with weights
// -----------------------------------------------------------------------------
   foreach i in 1..m do {
      W[i,i] = weights[i];
   }
// -----------------------------------------------------------------------------
// trickiest part is to build X
// -----------------------------------------------------------------------------
   for i in 1..m do {
      X[i,0] = 1.0 ;
      X[i,1] = n_xx[min_range[i]];
      for j in 2..n do {
         X[i,j] = (n_xx[min_range[i]])**j;
      }
   }
// -----------------------------------------------------------------------------
// Y also needs to be built
// -----------------------------------------------------------------------------
   foreach i in 1..m do {
      Y[i] = n_yy[min_range[i]];
   }
// -----------------------------------------------------------------------------
// ready for linear regression
// -----------------------------------------------------------------------------
   var A: [0..n,1..m] real;
   var B: [0..n,0..n] real;
   var C: [0..n,1..m] real;
   var beta: [0..n] real;
// -----------------------------------------------------------------------------
// lots of matrix multiplications!
// -----------------------------------------------------------------------------
   dot_mtm(X,W,A);
   dot_mm(A,X,B);
   minvgj(B);
   dot_mmt(B,X,A);            // re-using A !!!
   dot_mm(A,W,C);
   dot_mv(C,Y,beta);          // finally the parameters of the LLR
// -----------------------------------------------------------------------------
// now we can estimate!
// -----------------------------------------------------------------------------
   var n_xp: [0..n] real ;
   n_xp[0] = 1.0;
   n_xp[1] = n_x ;
   foreach j in 2..n do {
      n_xp[j] = n_x**j;
   }      
   var y = dot_vtv(beta,n_xp);  // and the local estimate!
   return denormalize_y(y);
}
// -----------------------------------------------------------------------------
// --> estimate is the core loess algorithm
// -----------------------------------------------------------------------------
proc lowess_fast_estimate(
   in x: real                 // point around which to do linear regression
   ) {
   const m: int = window;     // size of neighborhood == window
   var n_x = normalize_x(x);  // normalize independent variable
// -----------------------------------------------------------------------------
// do I really need to calculate all distances for each point?
// -----------------------------------------------------------------------------
   var nd = n_xx.size;
   var maxdist: real;
   var distances: [1..nd] real = abs(n_x - n_xx);
   var min_range: [1..window] int;
   var weights: [1..window] real;
   get_min_range(distances, min_range, maxdist);
   get_weights(maxdist,distances, min_range, weights);
// -----------------------------------------------------------------------------
// arrays for linear regression
// -----------------------------------------------------------------------------      
   var X: [1..m] real;
   var Y: [1..m] real;
// -----------------------------------------------------------------------------
// build X
// -----------------------------------------------------------------------------
   foreach i in 1..m do {
      X[i] = n_xx[min_range[i]];
   }
// -----------------------------------------------------------------------------
// Y also needs to be built
// -----------------------------------------------------------------------------
   foreach i in 1..m do {
      Y[i] = n_yy[min_range[i]];
   }
// -----------------------------------------------------------------------------
// ready for linear regression
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// now we can estimate!
// -----------------------------------------------------------------------------
   var a: real;
   var b: real;
   var r: real;
   (a,b,r) = reglin(X,Y);
   var y = a*n_x + b;
   return denormalize_y(y);
}
// -----------------------------------------------------------------------------
// --> performance: given an array O of observations and an array P of
// predictions, calculates several performance statistics: correlation
// coefficient, BIAS, RMSE, MAE, and Willmott's refined performance index
//
// Willmott, C. J.; Robeson, S. M. & Matsuura, K. A refined index of
// model performance International Journal of Climatology, 2012,
// v. 32, p. 2088-2094
// -----------------------------------------------------------------------------
proc performance(O: [] real, P: [] real): (real,real,real,real,real,real) {
   assert (O.rank == 1);
   assert (P.rank == 1);
   var n = O.size;
   assert (P.size == n);
// -----------------------------------------------------------------------------
// first the means
// -----------------------------------------------------------------------------
   var 
      Omean = 0.0,
      Pmean = 0.0;
   for (xo,xp) in zip(O,P) do {
      Omean += xo;
      Pmean += xp;
   }
   Omean /= n;
   Pmean /= n;
// -----------------------------------------------------------------------------
// now the stats
// -----------------------------------------------------------------------------
   var
      absdif = 0.0,
      absdom = 0.0,
      mse = 0.0,
      varO = 0.0,
      varP = 0.0,
      cvOP = 0.0;
   for (xo,xp) in zip(O,P) do {
      var del_o = xo - Omean;
      var del_p = xp - Pmean;
      var del_op = xp - xo;
      absdif += abs(del_op);
      absdom += abs(del_o);
      mse += del_op**2;
      varO += del_o**2;
      varP += del_p**2;
      cvOP += del_o*del_p;
   }
// -----------------------------------------------------------------------------
// now the summaries
// -----------------------------------------------------------------------------
   var r = cvOP/(sqrt(varO*varP));      // coefficient of correlation
   var r2 = 1 - mse/varO;               // coefficient of determination
   var bias = Pmean - Omean;            // bias
   var mae = absdif/n;                  // mean absolute error
   mse /= n;                            // mean square error
   var rmse = sqrt(mse);                // root mean square error
// -----------------------------------------------------------------------------
// the refined willmott performance index is more complicated
// -----------------------------------------------------------------------------
   var dbldom = 2*absdom;
   var dr: real;
   if absdif <= dbldom then {
      dr = 1.0 - absdif/dbldom;
   }
   else {
      dr = dbldom/absdif - 1.0;
   }
   return(r,r2,bias,mae,rmse,dr);
}

// ------------------------------------------------------------------------------
// --> steep: nonlinear least squares by curve fitting with the steepest descent
// method
//
// Hopefully, backtracking will come from:
// https://www.cs.cmu.edu/~ggordon/10725-F12/slides/05-gd-revisited.pdf
// ------------------------------------------------------------------------------
proc steep(
   const ref ax: [] real,          // ind variables (used as arg to func) (m x 1)
   const ref ay: [] real,          // data to be fit by func(x,p) (m x 1)
   const ref aw: [] real,          // array, *not matrix*, of weights
   ref ap: [] real,                // initial guess of parameter values  (n x 1)
                                   // returns the estimated parameters
   ref asigp: [] real,             // standard  errors of the parameters
   ref acp: [] real,               // parameter covariance matrix
   func,                           // in the simulated model: y_hat = func(x,p)
   const in epsilon = 1.0e-6,      // stop criterion
   const in broyden = false        // use Broyden update?
   ) : (real,real,real)
   where ( ax.rank == 1 && ay.rank == 1 && ap.rank == 1
        && asigp.rank == 1 && acp.rank == 2) {
   const maxiter = 10000;
   const m = ax.size;
   const n = ap.size;
   assert (ay.size == m);
   assert (aw.size == m);
   assert (asigp.size == n);
   assert (acp.shape == (n,n));
// -----------------------------------------------------------------------------
// reindexing
// -----------------------------------------------------------------------------   
   ref x = ax.reindex(1..m);
   ref y = ay.reindex(1..m);
   ref w = aw.reindex(1..m);
   ref p = ap.reindex(1..n);
   ref sigp = asigp.reindex(1..n);
   ref cp = acp.reindex({1..n,1..n});
// -----------------------------------------------------------------------------
// local scalar variables
// -----------------------------------------------------------------------------   
   var eps = epsilon;
   var iiter = 0;
   var chi2p = 0.0;                // the current figure of merit
   var chi2p_o = 0.0;              // the previous figure of merit
// -----------------------------------------------------------------------------
// local array variables
// -----------------------------------------------------------------------------   
   var J: [1..m,1..n] real;        // the jacobian matrix
   var dely: [1..m] real;
   var vaux_m: [1..m] real;        // aux m-vector
   var vaux_n: [1..n] real;        // aux n-vector
   var maux_mn: [1..m,1..n] real;  // aux (m,n)-matrix
   var gradchi2: [1..n] real;      // the gradient
   var hh: [1..n] real;            // the step in p
   var yhat: [1..m] real;          // function estimates
   var yhat_o: [1..m] real;        // the old ones
// -----------------------------------------------------------------------------
// if all w are equal to -1, then there is no estimate of w: set it to 1
// -----------------------------------------------------------------------------
   var noweights = false;
   if allequalto(w,-1.0) then {
      w = 1.0;
      noweights = true;
   }
// -----------------------------------------------------------------------------
// main loop
// -----------------------------------------------------------------------------   
   while eps >= epsilon do {
      if iiter > maxiter then {
         writef("nstat-->steep: I have exceeded %i iterations",maxiter);
         halt();
      }
// -----------------------------------------------------------------------------
// recalculate the Jacobian
// -----------------------------------------------------------------------------
      if !broyden then {
         func(x,p,yhat);           // estimate yi's
         dely = yhat - y;          // [y - hat{y}]
         simplejacob();
      } else {
         chi2p_o = chi2p;
         yhat_o = yhat ;
         chi2p = chi2(p);
         func(x,p,yhat);           // estimate yi's
         dely = yhat - y;          // [y - hat{y}]
         if mod(iiter,2*n) == 0 || chi2p > chi2p_o then {
            simplejacob();          // recalculate the jacobian
         }
         else {                     // Broyden rank-1 update formula
            Broyden();
         }
      }
      dot_diagm_v(w,dely,vaux_m);
      dot_mtv(J,vaux_m,gradchi2);
      gradchi2 *= 2;
// -----------------------------------------------------------------------------
// backtracking line search
// -----------------------------------------------------------------------------      
      var modgr2 = dot_vtv(gradchi2,gradchi2);
      var t = 1.0;
      var pa = p - t*gradchi2;
      var chi2a = chi2(pa);
      var chi2b = chi2(p) - (t/2)*modgr2;
//      writef("       t: %8.2er   chi2a: %8.2er   chi2b: %8.2er\n",t,chi2a,chi2b);
      while chi2a > chi2b do {
         t *= 0.1;
         pa = p - t*gradchi2;
         chi2a = chi2(pa);
         chi2b = chi2(p) - (t/2)*modgr2;
//         writef("       t: %8.2er   chi2a: %8.2er   chi2b: %8.2er\n",t,chi2a,chi2b);
      }             // end backtracking
      hh = t*gradchi2;
      p = p - hh ;
      eps = sqrt(dot_vtv(hh,hh));
//      writeln("eps = ",eps);
      iiter += 1;
   }
// -----------------------------------------------------------------------------
// error statistics
// -----------------------------------------------------------------------------
   var sumw = (+ reduce w);
   var sum2yhat = chi2(p);
   var sig2yhat = sum2yhat/sumw;
   var (ym,yvar) = wstat2(y,w);
   var r2 = 1.0 - sig2yhat/yvar;
   var redchi2: real;
   if noweights then {
      redchi2 = 1.0;
   }
   else {
      redchi2 = sum2yhat/(m - n);
   }
// -----------------------------------------------------------------------------
// straighforward calculation of the parameter covariance matrix and the
// asymptotic standard parameter errors
// -----------------------------------------------------------------------------
   dot_diagm_m(w,J,maux_mn);       // WJ
   dot_mtm(J,maux_mn,cp);          // J'WJ
   minvgj(cp);                     // [J'WJ]^{-1}
   vdiag(cp,sigp);                 // diagonal of [J'WJ]^{-1}
   sigp = sqrt(sigp);              // square root of each element
// -----------------------------------------------------------------------------
// square root below is for the standard error of estimate, which is more
// readily grasped
// -----------------------------------------------------------------------------
   return (redchi2,sqrt(sig2yhat),r2);
// -----------------------------------------------------------------------------
// proc chi2: the figure of merit
// -----------------------------------------------------------------------------
   proc chi2(ref pa: [] real): real where (pa.rank == 1)  {
      assert(pa.size == n);
      var ya,delya: [1..m] real;
      func(x,pa,ya);
      delya = ya - y;
      dot_diagm_v(w,delya,vaux_m);
      var merit = dot_vtv(vaux_m,delya);
      // writeln(" merit = ", merit);
      return merit;      
   }
// -----------------------------------------------------------------------------
// brute-force calculation of the jacobian
// -----------------------------------------------------------------------------   
   proc simplejacob() {
      const delp: [1..n] real = 1.0e-6;
      var forwp: [1..n] real;
      var backp: [1..n] real;
      var yplus: [1..m] real;
      var yminus: [1..m] real;
      for k in 1..n do {
         forwp = p;
         backp = p ;
         forwp[k] += delp[k];
         backp[k] -= delp[k];
         func(x,forwp,yplus);
         func(x,backp,yminus);
         J[1..m,k] = (yplus[1..m] - yminus[1..m])/(2*delp[k]);
      }
      // for i in 1..10 do {
      //    writeln(J[i,1..2]);
      // }
   }
// -----------------------------------------------------------------------------
// Broyden's rank-1 update formula
// -----------------------------------------------------------------------------
   proc Broyden() {
      dot_mv(J,hh,vaux_m);
      vaux_m = yhat - yhat_o + vaux_m;
      dot_vvt(vaux_m,hh,maux_mn);
      maux_mn /= (dot_vtv(hh,hh));
      J += maux_mn;
      // for i in 1..10 do {
      //    writeln(J[i,1..2]);
      // }
   }
}

