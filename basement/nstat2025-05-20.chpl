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
// 2023-03-10T20:25:20 already included newton
// 2023-03-10T20:25:11 including levmar
// 2023-03-11T22:19:45 removed Broyden completely
// 2024-03-12T20:40:46 replacing record procedures with procedures using vec's
// 2024-03-15T11:13:04 using mat from ada for x in steep, gnewton and levmar
// 2024-10-11T13:51:00 including a procedure for barnes's analysis
// 2025-02-21T07:52:55 golay-savitzky (starting) [but not today]
// 2025-05-20T13:28:38 producing a backup file nstat2025-05-20.chpl
// =============================================================================
use smatrix;
use ssr only allequal, allequalto, amin, amax, aminz, interp, heapsort, indxsort, indxquickselect;
use IO only readln, stderr;
use ada;
use Math;
// -----------------------------------------------------------------------------
// --> stat1: mean
// -----------------------------------------------------------------------------
proc stat1(
   const ref x: [] real       // the data
   ): real                       // mean
   where x.rank == 1 {
   // --------------------------------------------------------------------------
   // be careful with empty arrays
   // --------------------------------------------------------------------------
   var n = x.size;
   if n == 0 then {
      halt("nstat-->stat1: empty array");
   }
   var fn = n:real; // just in case
   // --------------------------------------------------------------------------
   // first the mean
   // --------------------------------------------------------------------------
   var xm = 0.0;
   xm = (+ reduce x);
   xm /= fn;
   return xm;
}
// -----------------------------------------------------------------------------
// --> stat2: mean and variance (always biased)
// -----------------------------------------------------------------------------
proc stat2(
   const ref x: [] real       // the data
   ): (real,real)             // mean, variance
   where x.rank == 1 {
   // --------------------------------------------------------------------------
   // be careful with empty arrays
   // --------------------------------------------------------------------------
   var n = x.size;
   if n == 0 then {
      halt("nstat-->stat2: empty array");
   }
   var fn = n:real; // just in case
   // --------------------------------------------------------------------------
   // first the mean
   // --------------------------------------------------------------------------
   var xm = (+ reduce x);
   xm /= fn;
   // --------------------------------------------------------------------------
   // now the variance
   // --------------------------------------------------------------------------
   var xv = (+ reduce ((x-xm)**2) );
   xv /= fn;
   return (xm,xv);
}

// -----------------------------------------------------------------------------
// --> stat3: mean, variance and 3rd central moment (all biased)
// -----------------------------------------------------------------------------
proc stat3(
   const ref x: [] real       // the data
   ): (real,real,real)        // mean, variance, 3rd central moment
   where x.rank == 1 {
   // --------------------------------------------------------------------------
   // be careful with empty arrays
   // --------------------------------------------------------------------------
   var n = x.size;
   if n == 0 then {
      halt("nstat-->stat2: empty array");
   }
   // --------------------------------------------------------------------------
   // first the mean
   // --------------------------------------------------------------------------
   var xm = (+ reduce x);
   xm /= n;
   // --------------------------------------------------------------------------
   // now the variance
   // --------------------------------------------------------------------------
   var xv = (+ reduce ((x-xm)**2) );
   xv /= n;
   // --------------------------------------------------------------------------
   // now the 3rd
   // --------------------------------------------------------------------------
   var x3 = (+ reduce ((x-xm)**3) );
   x3 /= n;
   return (xm,xv,x3);
}


// -----------------------------------------------------------------------------
// --> nanstat1: calculate 1st moment, ignoring nans
// -----------------------------------------------------------------------------
proc nanstat1(
   ref x: [] real             // the data
   ): (int,real)                 // number of valid data, mean
   where x.rank == 1 {
   const n = x.size;
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->nanstat1: empty array");
   }
   var m = 0;                 // how many not nans?
   var xm = 0.0;              // the mean over them
   for xi in x do {
      if !isNan(xi) then {
         m += 1;
         xm += xi;
      }
   }
   if m == 0 then {           // is x all nans?
      return (0,nan);
   }
   xm /= m;                   // calculate mean over valid data
   return (m,xm);
}

// -----------------------------------------------------------------------------
// --> nanstat2: calculate 1st and 2nd moments, ignoring nans
// -----------------------------------------------------------------------------
proc nanstat2(
   ref x: [] real             // the data
   ): (int,real,real)         // number of valid data, mean, variance
   where x.rank == 1 {
   const n = x.size;
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->nanstat2: empty array");
   }
   var m = 0;                 // how many not nans?
   var xm = 0.0;              // the mean over them
   for xi in x do {
      if !isNan(xi) then {
         m += 1;
         xm += xi;
      }
   }
   if m == 0 then {           // is x all nans?
      return (0,nan,nan);
   }
   xm /= m;                   // calculate mean over valid data
   var xv = 0.0;              // the variance
   for xi in x do {
      if !isNan(xi) then {
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
// 2021-04-15T15:08:47 eliminating the analysis of nans
// -----------------------------------------------------------------------------
proc wstat1(
   const ref x: [] real,      // the data
   const ref w: [] ?tw        // the weights
): real                       // the weighted mean
where (x.rank == 1 && w.rank == 1) {
   const n = x.size;          // the size of x
   assert (w.size == n);      // must be the same as that of w
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
proc wstat2(
   const ref x: [] real,      // the data
   const ref w: [] ?tw        // the weights
): (real,real)                // mean, variance
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
// --> nanwstat1: calculate weighted 1st moment, ignoring nans
//
// 2021-03-26T13:46:25 w can be of generic type!
// -----------------------------------------------------------------------------
proc nanwstat1(
   const ref x: [] real,      // the data
   const ref w: [] ?tw        // the weights
): (int,real)                 // number of valid data points, weighted mean
where (x.rank == 1 && w.rank == 1) {
   const n = x.size;
   if n == 0 then {           // be careful with empty arrays
      halt("nstat-->nanwstat1: empty array");
   }
   var m = 0;                 // how many not nans?
   var sw: tw = 0;            // the sum of the weights
   var sx = 0.0;              // the sum of the elements
   for (wi,xi) in zip(w,x) do {
      if !isNan(xi) then {
         m += 1;
         sw += wi;
         sx += wi*xi;
      }
   }
   if m == 0 then {           // is x all nans?
      return (0,nan);
   }
   var xm = sx/sw;            // weighted mean over valid data
   return (m,xm);
}
// -----------------------------------------------------------------------------
// --> gmean: the geometric mean
// -----------------------------------------------------------------------------
proc gmean(
   ref x: [] real                  // the data
): real                            // the geometric mean
where x.rank == 1 {
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
inline proc logmod(
   const in x: real
): real {
   return sgn(x)*log(abs(x)+1);
}
// -----------------------------------------------------------------------------
// --> exp-modulus function sgn(x)*(exp(|x|) - 1)
// -----------------------------------------------------------------------------
inline proc expmod(
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
   ): real                         // the geometric-modulus mean
   where x.rank == 1 {
   var ax = abs(x);                // absolute values!
   var delta = aminz(ax)/16;       // minimum non-zero absolute value
   ax = x/delta;                   // change scale to avoid values close to 1
   var lsum = 0.0;
   for i in ax.domain do {         // mod-geometric mean: sum of logmod(x)
      lsum += logmod(ax[i]);
   }
   var n = x.size;                 // the size of x, of course
   lsum /= n;                      // mean(logmod(x))
   var gm = expmod(lsum);          // expmod(mean(logmod(x))
   gm *= delta;                    // rescale back
   return gm;
}

// -----------------------------------------------------------------------------
// --> covar: calculates the covariance
// -----------------------------------------------------------------------------
proc covar(
   const in xm: real,         // the mean of x
   const in ym: real,         // the mean of y
   const ref x: [] real,      // the x data
   const ref y: [] real       // the y data
   ): real                    // the covariance
   where (x.rank == 1 && y.rank == 1) {
   var n = x.size;
   if n == 0 then {
      halt("nstat-->covar: empty array x");
   }
   assert (n == y.size);
   var sxy = (+ reduce ((x-xm)*(y-ym)) );
   return sxy/n;
}
// -----------------------------------------------------------------------------
// --> trivar: calculate the trivariance
// -----------------------------------------------------------------------------
proc trivar(
   const in xm: real,         // the mean of x
   const in ym: real,         // the mean of y
   const in zm: real,         // the mean of z
   const ref x: [] real,      // the x data
   const ref y: [] real,      // the y data
   const ref z: [] real       // the z data
   ): real                       // the trivariance
   where (x.rank == 1) && (y.rank == 1) && (z.rank == 1) {
   var n = x.size;
   if n == 0 then {
      halt("nstat-->trivar: empty array x");
   }
   assert (n == y.size);
   assert (n == z.size);
   var sxyz = (+ reduce ((x-xm)*(y-ym)*(z-zm)) );
   return sxyz/n;
}
// -----------------------------------------------------------------------------
// --> median: returns the median of array x
//
// 2012-08-21T08:55:53 Python version
// 2021-03-19T09:37:33 Chapel version
// -----------------------------------------------------------------------------
proc median(
   const ref ax: [] ?at       // the data
   ): at                         // the median
   where (isNumericType(at) && ax.rank == 1) {
   var n = ax.size;
   ref x = ax.reindex(0..n-1);
   if n == 0 then {
      halt("nstat-->median: empty array");
   }
   var indx: [0..n-1] int;
   indxsort(x,indx);          // sort x by index
   var xmedian: at;
   // --------------------------------------------------------------------------
   //  is n even or odd ?
   //  -------------------------------------------------------------------------
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
// --> nanmedian: returns the median of array x, ignoring nans
//
// 2012-08-21T08:55:53 Python version
// 2021-03-19T09:37:33 Chapel version
// -----------------------------------------------------------------------------
proc nanmedian(
   const ref ax: [] ?at       // the data (may include NaNs)
): at                         // the median
where ax.rank == 1 {
   var x = purgeval(nan,ax);
   var n = x.size;
   if n == 0 then {
      halt("nstat-->nanmedian: empty array");
   }
   var indx: [0..n-1] int;
   indxsort(x,indx);          // sort x by index
   var xmedian: at;
// -----------------------------------------------------------------------------
//  is n even or odd ?
//  ----------------------------------------------------------------------------
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
   const in ax: [] real       // we need to make a copy!
   ): (real,real,real)        // the 3 quartiles
   where ax.rank == 1 {
   var n = ax.size;
   if n < 3 then {
      halt("nstat-->quartiles: need at least 3 data points...");
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
   const in x: [?Dx] real          // we need to make a copy!
   ): [0..4] real                  // an array with the whiskers
   where x.rank == 1 {
   var n = x.size;
   writeln("whiskers: Dx = ",Dx);
   if n < 3 then {
      halt("nstat-->whiskers: need at least 3 data points...");
   }
   var fn = [i in 0..n-1] (i:real)/(n-1);
   heapsort(x);                    // sort x: this only affects the local copy
   var whisk: [0..4] real;
   whisk[0] = x[Dx.low];           // the minimum
   whisk[1] = interp(0.25,fn,x);   // just interpolate the 3 quartiles
   whisk[2] = interp(0.50,fn,x);
   whisk[3] = interp(0.75,fn,x);
   whisk[4] = x[Dx.high];          // the maximum
   return whisk;
}
// -----------------------------------------------------------------------------
// --> xbins: given the arrays x,y of data points, and the bin limits xbmin,
// dxb, xbmax, such that xb[j] = xbmin + (j+0.5)*db, returns yb where yb[j] is a
// (ragged) array that contains all data points y[i] such that xb[j]-dx <= x[i]
// < xb[j]+dx
// -----------------------------------------------------------------------------
proc xbins(
   const in xbmin: real,      // lower x limit
   const in dxb: real,        // size of x bin
   const in xbmax: real,      // uper x limit
   const ref x: [] real,      // xdata
   const ref y: [] real       // ydata
   ): [] vec                  // y data in each bin
   where (x.rank == 1) && (y.rank == 1) {
   const n = x.size;                         // the size of x
   assert(n == y.size);
   const m = ((xbmax - xbmin)/dxb):int;      // the number of bins
   // writeln("# of xbins: m = ",m);
   // --------------------------------------------------------------------------
   // it is complicated to construct a ragged array
   // --------------------------------------------------------------------------
   var nb: [0..m-1] int = 0;                 // nb[j] is the size of yb[j]
   var yb: [0..m-1] vec;                     // yb is a ragged array
   // --------------------------------------------------------------------------
   // initialize the domains in yb
   // --------------------------------------------------------------------------
   for j in 0..m-1 do {
      yb[j] = new vec({0..9});
   }
   for (xi,yi) in zip(x,y) do {
      var j = ((xi - xbmin)/dxb):int;    // index of desired bin
      // -----------------------------------------------------------------------
      // if j < 0 || j >= m then continue;
      // -----------------------------------------------------------------------
      if j*(m-1-j) < 0 then {
         continue;
      }
      var k = nb[j];
      dgrow(k,yb[j].dom);
      yb[j][k] = yi;                       // fill jth bin
      nb[j] += 1;
   }
   // --------------------------------------------------------------------------
   // now resize the second dimension of yb
   // --------------------------------------------------------------------------
   for j in 0..m-1 do {
      yb[j].dom = {0..nb[j]-1};
   }
   return yb;
}
// -----------------------------------------------------------------------------
// --> outlimits: define limits for outliers using the interquartile range
// -----------------------------------------------------------------------------
proc outlimits(
   ref x: [] real,                 // the data
   in delta: real = 1.5            // interquartile range multiplier
   ): (real,real)                  // (minx, maxx)
   where x.rank == 1 {
   assert (delta >= 0.0);
   var (q1,q2,q3) = quartiles(x);  // the quartiles
   var iqr = q3 - q1;              // the interquartile range
   var minx = q1 - delta*iqr;      // the proposed minimum
   var maxx = q3 + delta*iqr;      // the proposed maximum
   return (minx,maxx);             
}


// -----------------------------------------------------------------------------
// --> force_bounds is a convenient macro
// -----------------------------------------------------------------------------
inline proc force_bounds(
   ref x: real,
   const in xmin: real,
   const in xmax: real) {
   if x < xmin then {
      x = xmin;
   }
   else if x > xmax then {
      x = xmax;
   }
}
// -----------------------------------------------------------------------------
// --> reglin: calculates the linear regression y = ax + b, plus the correlation
// coefficient r
//
// Nelson Luís Dias
// 2010-02-11T15:23:42
// 2010-02-11T15:23:45
// 2019-07-12T11:20:43 prevent division by zero
// ----------------------------------------------------------------------------
proc reglin(
   const ref x: [] real,      // x data
   const ref y: [] real       // y data
   ): (real,real,real)        // (a,b,r) where y = ax + b
   where (x.rank == 1 && y.rank == 1) {
   var n = x.size;
   if n == 0 then {
      halt("nstat-->reglin: array x is empty");
   }
   assert (y.size == n);
   var (xm,xvar) = stat2(x);
   var (ym,yvar) = stat2(y);
   if ( xvar == 0.0 ) || ( yvar == 0.0) then {
      halt("nstat-->reglin: zero x or y variance");
   }
   var coxy = covar(xm,ym,x,y);
   var a = coxy / xvar;
   var b = ym - (a * xm);
   var r = coxy / sqrt(xvar*yvar);
   return (a,b,r);
}
// -----------------------------------------------------------------------------
// --> reglina: calculates the linear regression y = a x , the correlation
// coefficient r, the standard error of estimate se, and the standard deviation
// sa the estimator of a
// -----------------------------------------------------------------------------
proc reglina(
   const ref x: [] real,      // x data
   const ref y: [] real,      // y data
   out a: real,               // the slope
   out cd: real,              // coeff of det: only cd makes sense here!
   out se: real,              // std error of estimate
   out sa: real               // error of estimate of a
   )
   where (x.rank == 1 && y.rank == 1) {
   var n = x.size;
   // --------------------------------------------------------------------------
   // starts by calculating central moments
   // --------------------------------------------------------------------------
   var (xm,xvar) = stat2(x);
   var (ym,yvar) = stat2(y);
   var coxy = covar(xm,ym,x,y);
   // --------------------------------------------------------------------------
   // translates to non-central moments
   // --------------------------------------------------------------------------
   var sx20 = n * ( xvar + (xm*xm) ) ;
   var sy20 = n * ( yvar + (ym*ym) ) ;
   var sxy0 = n * ( coxy + (xm*ym) ) ;
   // --------------------------------------------------------------------------
   // obtains the slope and its std deviation
   // --------------------------------------------------------------------------
   a = sxy0 / sx20 ;
   sa = sqrt(yvar/sx20);
   // --------------------------------------------------------------------------
   // 2007-09-21T01:40 and 2024-05-26T10:49:02 each time I look at this I
   // get a different result; at least this time I have documented the
   // equation in nstat.tex
   // --------------------------------------------------------------------------
   var se2 = ( sy20 - 2.0*a*sxy0 + a*a*sx20 ) / n  ;
   se = sqrt(se2) ;
   cd = 1.0 - se2/yvar;
} 
// -----------------------------------------------------------------------------
// --> reglinpar: calculates the linear regression y = ax + b, the correlation
// coefficient r, the standard deviation se of the distribution of y given x,
// and the standard deviations sa, sb of the estimators(?) of a and b
// -----------------------------------------------------------------------------
proc reglinpar(
   const ref x: [] real,      // the x data
   const ref y: [] real,      // the y data
   out a: real,               // the a of y = ax + b
   out b: real,               // the b of y = ax + b
   out r: real,               // correlation coefficient
   out se: real,              // standard error of estimate
   out sa: real,              // standard error of a
   out sb: real               // standard error of b
) 
where (x.rank == 1 && y.rank == 1) {
   var n = x.size;
   assert(y.size == n);
   // var fn = n:real;
   // var dn = n / ( n - 1.0);
   var en = n / ( n - 2.0);
   var (xavg,xvar) = stat2(x);
   var (yavg,yvar) = stat2(y);
   var coxy = covar(xavg,yavg,x,y);
   // xvar *= dn ;
   // yvar *= dn ;
   // coxy *= dn ;
   a = coxy / xvar ;
   b = yavg - a*xavg; 
   r = coxy / sqrt(xvar*yvar);
   var se2 = en*(1 - r**2)*yvar   ;
   se = sqrt(se2) ;
   sa = sqrt( se2 / (n * xvar ) ) ;
   sb = sqrt( se2 * ( 1 + (xavg*xavg) / xvar ) / n ) ;
}
// -----------------------------------------------------------------------------
// --> ols: ordinary least squares multivariate regression: quite ordinary,
//     indeed: wants to find the parameters C[m,1] such that:
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
   ref x: [] real,            // independent data matrix
   ref y: [] real,            // dependent data vector  
   out S2: real,              // variance of residuals
   out CD: real,              // adjusted coeff. of multiple determination
   ref c: [] real,            // array of coefficients
   ref ye: [] real,           // vector of estimated ys
   ref ccov: [] real,         // covariance matrix of c
   ref axcor: [] real         // correlation matrix of x
   )
   where (x.rank == 2 && y.rank == 1 &&
          c.rank == 1 && ccov.rank == 2 && axcor.rank == 2) {
   var (n,m) = x.shape;
   assert(y.shape == (n,));
   assert(ye.shape == (n,));
   assert(c.shape == (m,));
   assert(ccov.shape == (m,m));
   assert(axcor.shape == (m,m));
   assert(n > m);
   ref xcor = axcor.reindex({1..m,1..m});
   // --------------------------------------------------------------------------
   // auxiliary storage A and b is needed
   // --------------------------------------------------------------------------
   var A: [1..m,1..m] real;
   var b: [1..m,1..n] real;
   // --------------------------------------------------------------------------
   // calculates A = (x'x) 
   // --------------------------------------------------------------------------
   dot_mtm(x,x,A);
   for i in 1..m do {         // the correlation matrix
      for j in 1..m do {
         xcor[i,j] = A[i,j]/sqrt(A[i,i]*A[j,j]);
      }
   }
   // --------------------------------------------------------------------------
   // calculates A^-1
   // --------------------------------------------------------------------------
   minvgj(A);
   // --------------------------------------------------------------------------
   // calculates b = A^-1 x'
   // --------------------------------------------------------------------------
   dot_mmt(A,x,b);
   // --------------------------------------------------------------------------
   // calculates c = b y
   // --------------------------------------------------------------------------
   dot_mv(b,y,c);
   // --------------------------------------------------------------------------
   // calculates ye = x c
   // --------------------------------------------------------------------------
   dot_mv(x,c,ye) ;
   // --------------------------------------------------------------------------
   // variance of residuals
   // --------------------------------------------------------------------------
   S2 = (+ reduce ((y- ye)**2));   // reindexing no longer needed!
   S2 /= (n - m);                  // there you go
   ccov = S2*A;                    // covariance matrix of c
   // --------------------------------------------------------------------------
   // variance of observations
   // --------------------------------------------------------------------------
   var (ymed,yvar) = stat2(y);
   yvar *= (n:real)/(n-1.0);
   CD = 1.0 - S2/yvar;             // coefficient of determination
}
// =============================================================================
// ==> nlowess: my implementation of the lowess algorithm
//
// with most of the inspiration from
// https://medium.com/data-science-collective/loess-373d43b03564
//
// 2022-04-22T12:11:44 a new star is born
// 2022-04-22T17:34:40 essentially done ... hopefully
// 2022-04-23T18:54:06 moved into nstat
// 2025-02-19T08:33:44 is n_xx a single variable for all lowess_estimate calls?
//    yes, because n_xx is only used by lowess_estimate to calculate distances
//    locally
//
// 2025-05-05T13:54:12 using only the diagonal of W???
// =============================================================================
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
   const ref a: [] real,            // raw data
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
   ldom = {0..-1};  // data domain initially empty
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
   const in win: int,         // the window, once and for all, of this loess
   const in deg: int = 1      // the degree of this loess
   )  where (xx.rank == 1) && (yy.rank == 1) {
   assert (xx.shape == yy.shape);
   // --------------------------------------------------------------------------
   // set these global variables once and for all
   // --------------------------------------------------------------------------
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
   ): real {
   const m: int = window;     // size of neighborhood == window
   const n: int = degree;     // degree of linear regression == degree
   var n_x = normalize_x(x);  // normalize independent variable
   // --------------------------------------------------------------------------
   // do I really need to calculate all distances for each point?
   // --------------------------------------------------------------------------
   var nd = n_xx.size;
   var maxdist: real;
   var distances: [1..nd] real = abs(n_x - n_xx);
   var min_range: [1..window] int;
   var weights: [1..window] real;
   get_min_range(distances, min_range, maxdist);
   get_weights(maxdist, distances, min_range, weights);
   // --------------------------------------------------------------------------
   // here I use m,n instead of window, degree
   // estimate beta = [Xtr(n+1,m)W(m,m)X(m,n+1)]^(-1)Xtr(n+1,m)W(m,m)n_yy(m,1)
   // --------------------------------------------------------------------------
   var W: [1..m] real = 0.0;
   var X: [1..m,0..n] real;
   var Y: [1..m] real;
   // --------------------------------------------------------------------------
   // fill diagonal of W with weights
   // --------------------------------------------------------------------------
   // foreach i in 1..m do {
   //    W[i,i] = weights[i];
   // }
   foreach i in 1..m do {
      W[i] = weights[i];
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
   dot_mt_diagm(X,W,A);
   dot_mm(A,X,B);
   minvgj(B);
   dot_mmt(B,X,A);            // re-using A !!!
   dot_m_diagm(A,W,C);
   dot_mv(C,Y,beta);          // finally the parameters of the LLR
   // --------------------------------------------------------------------------
   // now we can estimate!
   // --------------------------------------------------------------------------
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
// --> lowess_fast_estimate: uses linear regression; this version does not use
//     weights!
// -----------------------------------------------------------------------------
proc lowess_fast_estimate(
   in x: real                 // point around which to do linear regression
): real {
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
// coefficient r, coefficient of determinantion Cd, BIAS, MAE, RMSE, and
// Willmott's refined performance index dr: (r,Cd,BIAS,MAE,RMSE,dr)
//
// 2023-05-18T18:08:51: Whenever r and/or Cd do not make sense, they return with
// an (impossible) value of +100
//
// Willmott, C. J.; Robeson, S. M. & Matsuura, K. A refined index of
// model performance International Journal of Climatology, 2012,
// v. 32, p. 2088-2094
// -----------------------------------------------------------------------------
proc performance(
   const ref O: [] real,
   const ref P: [] real
   ): (real,real,real,real,real,real)
   where (O.rank == 1) && (P.rank == 1) {
   var n = O.size;
   assert (P.size == n);
   // --------------------------------------------------------------------------
   // first the means
   // --------------------------------------------------------------------------
   var 
      Omean = 0.0,
      Pmean = 0.0;
   for (xo,xp) in zip(O,P) do {
      Omean += xo;
      Pmean += xp;
   }
   Omean /= n;
   Pmean /= n;
   // --------------------------------------------------------------------------
   // now the stats
   // --------------------------------------------------------------------------
   var
      absdif = 0.0,
      absdom = 0.0,
      cvOP = 0.0,
      mse = 0.0,
      r = 0.0,
      Cd = 0.0,
      varO = 0.0,
      varP = 0.0;
// -----------------------------------------------------------------------------
//  because all Os and all Ps may be equal, this is not the most efficient way
//  to calculate things, but ...
//  -----------------------------------------------------------------------------
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
   var allO = allequal(O);
   if allO then {
      r = 100.0;                   // impossible to calculate r
      Cd = 100.0;                  // impossible to calculate Cd
   }
   else if allequal(P) then {
      r = 100.0;                   // impossible to calculate r
      Cd = 1 - mse/varO;
   }
   else {
      r = cvOP/(sqrt(varO*varP));  // coefficient of correlation
      Cd = 1 - mse/varO;           // coefficient of determination
   }
   var bias = Pmean - Omean;       // bias
   var mae = absdif/n;             // mean absolute error
   mse /= n;                       // mean square error
   var rmse = sqrt(mse);           // root mean square error
// -----------------------------------------------------------------------------
// the refined willmott performance index is more complicated
// -----------------------------------------------------------------------------
   var dbldom: real;
   if allO then {
      dbldom = 0.0;
   }
   else {
      dbldom = 2*absdom;
   }
   var dr: real;
   if absdif <= dbldom then {
      dr = 1.0 - absdif/dbldom;
   }
   else {
      dr = dbldom/absdif - 1.0;
   }
   return(r,Cd,bias,mae,rmse,dr);
}

// ------------------------------------------------------------------------------
// --> steep: nonlinear least squares by curve fitting with the steepest descent
// method
//
// Hopefully, backtracking will come from:
// https://www.cs.cmu.edu/~ggordon/10725-F12/slides/05-gd-revisited.pdf
// ------------------------------------------------------------------------------
proc steep(
   ref x: mat,          // ind variables (used as arg to func) (m x ell)
   ref y: vec,          // data to be fit by func(x,p,y) (m x 1)
   ref w: [] real,      // array, *not matrix*, of weights (m x 1)
   ref p: vec,          // initial guess of parameter values  (n x 1)
                        // returns the estimated parameters
   ref sigp: [] real,             // standard  errors of the parameters
   ref cp: [] real,               // parameter covariance matrix
   const func: proc(ref ax: mat,   // the independent variables
                    ref ap: vec,   // the parameters
                    ref ay: vec),  // in the sim model we call func(x,p,y)
   const in epsilon = 1.0e-6       // stop criterion
   ) : (real,real,real)
   where ( w.rank == 1 && sigp.rank == 1 && cp.rank == 2) {
   const maxiter = 100000;         // steep may take a looong time to converge
   const m = x.shape(0);           // the number of data points
   const n = p.size;               // the number of parameters
// -----------------------------------------------------------------------------
// check all shapes and sizes
// -----------------------------------------------------------------------------   
   assert (y.size == m);
   assert (w.size == m);
   assert (sigp.size == n);
   assert (cp.shape == (n,n));
// -----------------------------------------------------------------------------
// local scalar variables
// -----------------------------------------------------------------------------   
   var eps = epsilon;
   var iiter = 0;                  // number of iterations
// -----------------------------------------------------------------------------
// local array variables
// -----------------------------------------------------------------------------   
   var J: [1..m,1..n] real;        // the jacobian matrix
   var dely: [1..m] real;          // yhat - y
   var vaux_m: [1..m] real;        // aux m-vector
   var vaux_n: [1..n] real;        // aux n-vector
   var maux_mn: [1..m,1..n] real;  // aux (m,n)-matrix
   var gradchi2: [1..n] real;      // the gradient
   var hh: [1..n] real;            // the step in p
// -----------------------------------------------------------------------------
// arguments to functions must be vecs
// -----------------------------------------------------------------------------   
   var yhat = new vec({1..m});     // function estimates
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
         writef("nstat-->steep: I have exceeded %i iterations\n",maxiter);
         return (-1.0, -1.0, -1.0);
      }
// -----------------------------------------------------------------------------
// recalculate the Jacobian
// -----------------------------------------------------------------------------
      func(x,p,yhat);                   // estimate yi's
      dely = yhat.arr - y.arr;          // [hat{y}-y]: dely is an array
      simplejacob();                    // update the jacobian matrix
      dot_diagm_v(w,dely,vaux_m);       // W[hat{y} - y]
      dot_mtv(J,vaux_m,gradchi2);       // J'W[hat{y} - y]
      gradchi2 *= 2;                    // 2J'W[hat{y} - y]
// -----------------------------------------------------------------------------
// backtracking line search
// -----------------------------------------------------------------------------      
      var modgr2 = dot_vtv(gradchi2,gradchi2);    // |grad chi^2|
      var t = 1.0;                                // backtracing parameter
      var pa = p - tovec(t*gradchi2);             // check two chi^2s
                                                  // pa is a vec
      var chi2a = chi2(pa);                       // chi2 needs a vec
      var chi2b = chi2(p) - (t/2)*modgr2;         // chi2 needs a vec
// -----------------------------------------------------------------------------
// backtracing loop
// -----------------------------------------------------------------------------      
      while chi2a > chi2b do {
         t *= 0.1;
         pa = p - tovec(t*gradchi2);
         chi2a = chi2(pa);
         chi2b = chi2(p) - (t/2)*modgr2;
      }
// -----------------------------------------------------------------------------
// found the right scaling for the size of the step
// -----------------------------------------------------------------------------      
      hh = t*gradchi2;
      p = p - tovec(hh) ;
      eps = sqrt(dot_vtv(hh,hh));
      iiter += 1;
   }
// -----------------------------------------------------------------------------
// error statistics
// -----------------------------------------------------------------------------
   var sumw = (+ reduce w);             // sum of weights               
   var sum2yhat = chi2(p);              // sum of squares of deviations 
   var sig2yhat = sum2yhat/sumw;        // variance of deviations       
   var (ym,yvar) = wstat2(y.arr,w);     // weighted variance of data    
   var r2 = 1.0 - sig2yhat/yvar;        // coefficient of determination 
   var redchi2: real;                   // reduced chi-square           
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
   mvdiag(cp,sigp);                // diagonal of [J'WJ]^{-1}
   sigp = sqrt(sigp);              // square root of each element
// -----------------------------------------------------------------------------
// square root below is for the standard error of estimate, which is more
// readily grasped
// -----------------------------------------------------------------------------
   return (redchi2,sqrt(sig2yhat),r2);
// -----------------------------------------------------------------------------
// proc chi2: the figure of merit
// -----------------------------------------------------------------------------
   proc chi2(
      ref pa: vec
   ): real  {
      assert(pa.size == n);
      var ya = new vec({1..m});
      var delya: [1..m] real;
      func(x,pa,ya);
      delya = ya.arr - y.arr;
      dot_diagm_v(w,delya,vaux_m);
      var merit = dot_vtv(vaux_m,delya);
      // writeln(" merit = ", merit);
      return merit;      
   }
// -----------------------------------------------------------------------------
// brute-force calculation of the jacobian matrix
// -----------------------------------------------------------------------------   
   proc simplejacob() {
      const delp: [1..n] real = 1.0e-6;
      var forwp = new vec({1..n});
      var backp = new vec({1..n});
      var yplus = new vec({1..m});
      var yminus = new vec({1..m});
      for k in 1..n do {
         forwp = p;
         backp = p ;
         forwp[k] += delp[k];
         backp[k] -= delp[k];
         func(x,forwp,yplus);
         func(x,backp,yminus);
         J[1..m,k] = (yplus.arr[1..m] - yminus.arr[1..m])/(2*delp[k]);
      }
   }
}




// ------------------------------------------------------------------------------
// --> gnewton: nonlinear least squares by curve fitting with the Gauss-Newton
// method
//
// Hopefully, backtracking will come from:
// https://www.cs.cmu.edu/~ggordon/10725-F12/slides/05-gd-revisited.pdf
// ------------------------------------------------------------------------------
proc gnewton(
   ref x: mat,                // ind variables (used as arg to func) (m x ell)
   ref y: vec,                // data to be fit by func(x,p) (m x 1)
   ref w: [] real,            // array, *not matrix*, of weights
   ref p: vec,                     // initial guess of parameter values  (n x 1)
                                   // returns the estimated parameters
   ref sigp: [] real,              // standard  errors of the parameters
   ref cp: [] real,                // parameter covariance matrix
   const func: proc(ref ax: mat,   // the independent variables
                    ref ap: vec,   // the parameters
                    ref ay: vec),  // in the sim model we call func(x,p,y)
   const in epsilon = 1.0e-6       // stop criterion
   ) : (real,real,real)
   where ( w.rank == 1 && sigp.rank == 1 && cp.rank == 2) {
   const maxiter = 100000;         // maximum number of iterations
   const m = x.shape(0);           // number of data points
   const n = p.size;               // number of parameters
   assert (y.size == m);
   assert (w.size == m);
   assert (sigp.size == n);
   assert (cp.shape == (n,n));
// -----------------------------------------------------------------------------
// local scalar variables
// -----------------------------------------------------------------------------   
   var eps = epsilon;              // size of the update vector hh in each iter
   var iiter = 0;                  // number of iterations
// -----------------------------------------------------------------------------
// local array variables
// -----------------------------------------------------------------------------   
   var J: [1..m,1..n] real;        // the jacobian matrix
   var dely: [1..m] real;          // yhat - y
   var vaux_m: [1..m] real;        // aux m-vector
   var vaux_n: [1..n] real;        // aux n-vector
   var maux_mn: [1..m,1..n] real;  // aux (m,n)-matrix
   var maux_nm: [1..n,1..m] real;  // aux (n,m)-matrix
   var maux_nn: [1..n,1..n] real;  // aux (n,n)-matrix
   var hh: [1..n] real;            // the step in p
   var yhat = new vec({1..m});     // function estimates
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
         writef("nstat-->gnewton: I have exceeded %i iterations\n",maxiter);
         return(-1.0,-1.0,-1.0);
      }
// -----------------------------------------------------------------------------
// recalculate the Jacobian
// -----------------------------------------------------------------------------
      func(x,p,yhat);              // estimate yi's
      dely = yhat.arr - y.arr;     // [hat{y} - y]: dely is an array
      simplejacob();               // recalculate the jacobian matrix
// -----------------------------------------------------------------------------
// calculates the RHS of the linear system
// -----------------------------------------------------------------------------
      dot_mt_diagm(J,w,maux_nm);   // J'W
      dot_mv(maux_nm,dely,vaux_n); // J'W[hat{y} - y]
// -----------------------------------------------------------------------------
// calculates the LHS of the linear system
// -----------------------------------------------------------------------------
      dot_mm(maux_nm,J,maux_nn);   // J'WJ
      gauss(maux_nn,vaux_n,hh);    // [J'WJ]h = J'W[hat{y} - y] => h 
// -----------------------------------------------------------------------------
// update, etc.
// -----------------------------------------------------------------------------      
      p = p - tovec(hh) ;         // p is vec, hh is array
      eps = sqrt(dot_vtv(hh,hh));
      iiter += 1;
   }
// -----------------------------------------------------------------------------
// error statistics
// -----------------------------------------------------------------------------
   var sumw = (+ reduce w);             // sum of weights               
   var sum2yhat = chi2(p);              // sum of squares of deviations 
   var sig2yhat = sum2yhat/sumw;        // variance of deviations       
   var (ym,yvar) = wstat2(y.arr,w);     // weighted variance of data    
   var r2 = 1.0 - sig2yhat/yvar;        // coefficient of determination 
   var redchi2: real;                   // reduced chi-square                 
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
   mvdiag(cp,sigp);                // diagonal of [J'WJ]^{-1}
   sigp = sqrt(sigp);              // square root of each element
// -----------------------------------------------------------------------------
// square root below is for the standard error of estimate, which is more
// readily grasped
// -----------------------------------------------------------------------------
   return (redchi2,sqrt(sig2yhat),r2);
// -----------------------------------------------------------------------------
// proc chi2: the figure of merit 
// -----------------------------------------------------------------------------
   proc chi2(ref pa: vec): real {
      assert(pa.size == n);
      var ya = new vec({1..m});
      var delya: [1..m] real;
      func(x,pa,ya);
      delya = ya.arr - y.arr;
      dot_diagm_v(w,delya,vaux_m);
      var merit = dot_vtv(vaux_m,delya);
      // writeln(" merit = ", merit);
      return merit;      
   }
// -----------------------------------------------------------------------------
// --> simplejacob: brute-force calculation of the jacobian
// -----------------------------------------------------------------------------   
// -----------------------------------------------------------------------------
// brute-force calculation of the jacobian matrix
// -----------------------------------------------------------------------------   
   proc simplejacob() {
      const delp: [1..n] real = 1.0e-6;
      var forwp = new vec({1..n});
      var backp = new vec({1..n});
      var yplus = new vec({1..m});
      var yminus = new vec({1..m});
      for k in 1..n do {
         forwp = p;
         backp = p ;
         forwp[k] += delp[k];
         backp[k] -= delp[k];
         func(x,forwp,yplus);
         func(x,backp,yminus);
         J[1..m,k] = (yplus.arr[1..m] - yminus.arr[1..m])/(2*delp[k]);
      }
   }
}


// ------------------------------------------------------------------------------
// --> levmar: nonlinear least squares by curve fitting with the
// Levenberg-Marquard method. Here x is a mat.
//
// this code is a free translation (with several modifications) of the Matlab
// code of
//
// Gavin, H. P. (2022). The Levenberg-Marquardt algorithm for nonlinear least
// squares curve-fitting prob- lems. Available at
// https://people.duke.edu/~hpgavin/ExperimentalSystems/lm.pdf.
// ------------------------------------------------------------------------------
proc levmar(
   ref x: mat,           // independent variables (used as arg to func) (m x ell)
   ref y: vec,           // data to be fit by func(x,p) (m x 1)
   ref w: [] real,       // array, *not matrix*, of weights (m x 1)
   ref p: vec,           // initial guess of parameter values  (n x 1)
                         // returns the estimated parameters
   ref sigp: [] real,    // standard  errors of the parameters (n x 1)
   ref cp: [] real,      // parameter covariance matrix (n x n)
   const ref func: proc(ref ax: mat,    // the independent variables
                    ref ap: vec,        // the parameters
                    ref yhat: vec),     // in the sim model call func(ax,ap,yhat)
   const in epsilon_p = 1.0e-6          // stop criterion
   ) : (real,real,real)  // (red chi sq, st err of estimate, coeff det)
   where ( w.rank == 1 && sigp.rank == 1 && cp.rank == 2) {
// =============================================================================
// initializations
// 2024-03-18T18:52:10 needed to change epsilon_r drastically for polyalb to
// work
// =============================================================================   
//   const epsilon_r = 1.0e-3;     // rho criterion
   const epsilon_r = 1.0e-10;      // rho criterion changed drastically now!
   const maxiter = 100000;         // the maximum number of iterations
   const m = x.shape(0);           // the number of data points
   const ell = x.shape(1);         // dimension of the data points
   const n = p.size;               // the number of parameters
   assert (p.dom == {1..n});       // but only this will work
   assert (y.size == m);
   assert (w.size == m);
   assert (sigp.size == n);
   assert (cp.shape == (n,n));
// -----------------------------------------------------------------------------
// local scalar variables
// -----------------------------------------------------------------------------   
   var eps = epsilon_p;            // convergence criterion in p
   var iiter = 0;                  // iteration counter
   var chi2p_o = 0.0;              // the previous figure of merit
   var chi2p = 1.0;                // the current figure of merit
   var lamb = 1.0;                 // initial value of Levenberg-Marquardt's
                                   // lambda
   var rho: real;                  // how to update lambda and p?
// -----------------------------------------------------------------------------
// local array variables
// -----------------------------------------------------------------------------   
   var J: [1..m,1..n] real;        // the jacobian matrix
   var dely: [1..m] real;          // hat{y} - y
   var vaux_m: [1..m] real;        // aux m-vector
   var vaux_n: [1..n] real;        // aux n-vector
   var vaux_n1: [1..n] real;       // aux n-vector
   var pa = new vec({1..n});       // auxiliary parameter vector
   var maux_nm: [1..n,1..m] real;  // aux (n,m)-matrix: maux_nm == J'W
   var
      A_nn,                        // LM system matrix
      D,                           // diagonal of J'WJ
      lambD,                       // lamb*D
      maux_nn: [1..n,1..n] real;   // J'WJ      
   var hh: [1..n] real;            // the step in p
   var yhat = new vec({1..m});     // function estimates
// =============================================================================
// algol starts
// =============================================================================
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
   while eps >= epsilon_p do {
      if iiter > maxiter then {
         writef("nstat-->levar: I have exceeded %i iterations\n",maxiter);
         return(-1.0,-1.0,-1.0);
      }
// -----------------------------------------------------------------------------
// note that chi2p_o and chi2p have already been calculated in the previous
// iteration (or defined to be different at the outset), so:
// a) if p has changed (chi2p != chi2p_o), update p, yhat and J
// b) otherwise, just pick up where lamb has changed
// -----------------------------------------------------------------------------         
      if chi2p != chi2p_o then {
// -----------------------------------------------------------------------------
// update the Jacobian matrix?
// -----------------------------------------------------------------------------         
         func(x,p,yhat);           // if p has changed, update yhat
         // writeln("-"*40);
         // writeln("p =    ",p);
         // writeln("x =    ",x[1..10]);
         // writeln("y =    ",y[1..10]);
         // writeln("yhat = ",yhat[1..10]);
         dely = yhat.arr - y.arr;  // [hat{y}- y]: dely is array
         simplejacob();            // new p -> new Jacobian matrix
         // writeln("J = ",J);
         // exit(1);
// -----------------------------------------------------------------------------
// with J updated, calculate the RHS of the LM linear system
// -----------------------------------------------------------------------------
         dot_mt_diagm(J,w,maux_nm);     // maux_nm == J'W
         dot_mv(maux_nm,dely,vaux_n);   // vaux_n == J'W[hat{y} - y]
// -----------------------------------------------------------------------------
// calculate (almost all of) the LHS of the LM linear system
// -----------------------------------------------------------------------------
         dot_mm(maux_nm,J,maux_nn);     // maux_nn == J'WJ
         mmdiag(maux_nn,D);             // D == its diagonal as a matrix
      }
// -----------------------------------------------------------------------------
// lambda-dependency only starts here!
// -----------------------------------------------------------------------------      
      lambD = lamb*D;              // lamb always changes, but not D
      A_nn = maux_nn + lambD;      // J'WJ + lamb D (LM system matrix)
      gauss(A_nn,vaux_n,hh);       // solve for hh
//      writeln("hh = ",hh);
// -----------------------------------------------------------------------------
// calculate denominator of rho
// -----------------------------------------------------------------------------
      dot_mv(lambD,hh,vaux_n1);    // vaux_n1 == lamb D h
      vaux_n1 += vaux_n;           // lamb D h + J'W[hat{y}-y]
      var rhoden = abs(dot_vtv(hh,vaux_n1));     // | above |
// -----------------------------------------------------------------------------
// what should I do with the updates?
// -----------------------------------------------------------------------------      
      chi2p_o = chi2(p);           // the old figure of merit
      pa = p - tovec(hh) ;         // do not update p yet!
//      writeln("pa = ",pa);
      chi2p = chi2(pa);            // the new figure of merit
      rho = (chi2p_o - chi2p)/rhoden;
//      writeln("rho = ",rho);
// -----------------------------------------------------------------------------
// delayed gratification: 5 and 2 instead of 10
// -----------------------------------------------------------------------------
      if rho > epsilon_r then {
         p = pa ;                  // accept the update
         lamb /= 5.0;              // decrease lamb
      }
      else {                       // do not update p!
         lamb *= 2.0;              // increase lamb
         chi2p = chi2p_o;          // undo chi2p
      }
      eps = sqrt(dot_vtv(hh,hh));  // |hh|
      iiter += 1;
//      writeln("p = ",p);
      // writeln("eps = ",eps);
      // writeln("rho = ",rho);
      // writeln("chi2p = ",chi2p);
   }      
// =============================================================================
// error statistics
// =============================================================================
   var sumw = (+ reduce w);             // sum of weights              
   var sum2yhat = chi2(p);              // sum of squares of deviations
   var sig2yhat = sum2yhat/sumw;        // variance of deviations      
   var (ym,yvar) = wstat2(y.arr,w);     // weighted variance of data   
   var r2 = 1.0 - sig2yhat/yvar;        // coefficient of determination
   var redchi2: real;                   // reduced chi-square          
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
   cp = maux_nn;                   // J'WJ
   minvgj(cp);                     // [J'WJ]^{-1}
   mvdiag(cp,sigp);                // diagonal of [J'WJ]^{-1}
   sigp = sqrt(sigp);              // square root of each element
// -----------------------------------------------------------------------------
// square root below is for the standard error of estimate, which is more
// readily grasped
// -----------------------------------------------------------------------------
   return (redchi2,sqrt(sig2yhat),r2);
// =============================================================================
// local procs
// =============================================================================
// -----------------------------------------------------------------------------
// proc chi2: the figure of merit
// -----------------------------------------------------------------------------
   proc chi2(ref pa: vec): real {
      assert(pa.size == n);
      var ya = new vec({1..m});
      var delya: [1..m] real;
      func(x,pa,ya);
      delya = ya.arr - y.arr;
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
      var forwp = new vec({1..n});
      var backp = new vec({1..n});
      var yplus = new vec({1..m});
      var yminus = new vec({1..m});
      for k in 1..n do {
         forwp = p;
         backp = p;
         forwp[k] += delp[k];
         backp[k] -= delp[k];
         func(x,forwp,yplus);
         func(x,backp,yminus);
         J[1..m,k] = (yplus.arr[1..m] - yminus.arr[1..m])/(2*delp[k]);
      }
   }
}
// -----------------------------------------------------------------------------
// --> barnes: barnes's analysis
//
// from a set of grid points and their values (xin, yin, zin), returns the
// "analyzed" value z at (x,y)
// -----------------------------------------------------------------------------
proc barnes(
   const in x: real,          // the abscissa of the desired point value
   const in y: real,          // the ordinate of the desired point value
   const in r0: real,         // the radius of the analysis
   const ref axin: [] real,   // the abscissas of the existing data (1D)
   const ref ayin: [] real,   // the ordinates of the existing data (1D)
   const ref azin: [] real    // the values of the existing data
   ) : real                   // the value at the desired point
   where (axin.rank == 1) && (ayin.rank == 1) && (azin.rank == 1) {                          
   // --------------------------------------------------------------------------
   // number of available data points
   // --------------------------------------------------------------------------
   var n = axin.size;
   assert( (ayin.size == n) && (azin.size) == n );
   // --------------------------------------------------------------------------
   // reindexing is part of my game
   // --------------------------------------------------------------------------
   const ref xin = axin.reindex(1..n);
   const ref yin = ayin.reindex(1..n);
   const ref zin = azin.reindex(1..n);
   var npt = 0;               // number of points used in the analysis
   var zb: [1..n] real;       // auxiliary array to perform analysis
   var wb: [1..n] real;       // auxiliary array to perform analysis
   // --------------------------------------------------------------------------
   // go over all data points
   // --------------------------------------------------------------------------
   for i in 1..n do {
      // -----------------------------------------------------------------------
      // distance to each existing data point
      // -----------------------------------------------------------------------
      var r = sqrt( (x - xin[i])**2 + (y - yin[i])**2 );
   // --------------------------------------------------------------------------
   // only analyze within radius
   // --------------------------------------------------------------------------
      if ( r <= r0 ) then {
         npt += 1;            // I have found one more point
         var rho = 3.034854*r/r0 ;
         zb[npt] = zin[i];
         wb[npt] = exp(-(rho**2)) ;  
      }
   } 
   // --------------------------------------------------------------------------
   // number of valid points
   // --------------------------------------------------------------------------
   if ( npt == 0 ) then {
      writef(stderr,"--> barnes: no valid points\n") ;
      writef(stderr,"    around (%10.4r,%10.4r)\n",x,y) ;
      halt();
   }
   // --------------------------------------------------------------------------
   // completed loop over all points inside radius of influence; estimate z at
   // point
   // --------------------------------------------------------------------------
   var sum0 = 0.0 ;
   var sumz = 0.0 ;
   for k in 1..npt do {
      sum0 += wb[k] ;
      sumz += wb[k]*zb[k] ;      
   }
   var z = sumz/sum0 ;
   return z;   
}

// -----------------------------------------------------------------------------
// --> savitzky-golay: the savitzky-golay filter? one day!
//
// given xdata and ydata, and a window size m ...
// -----------------------------------------------------------------------------
