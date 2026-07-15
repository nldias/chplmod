// =============================================================================
// ==> geo: things that vary spatially
//
// Nelson Luís Dias
// 2026-02-04T11:01:06 a new star is born
// =============================================================================
use IO only stderr, writeln;
use Math only isClose, phase,pi,sqrt;
use dgrow only dgrow;
use nstat only stat2, covar;
// -----------------------------------------------------------------------------
// --> inxy: the standard way of finding inside points
// -----------------------------------------------------------------------------
proc inxy(
   const in xp: real,              // is this point in or out?
   const in yp: real,              // is this point in or out?
   const ref axboun: [] real,      // abscissas of perimeter
   const ref ayboun: [] real       // ordinates of perimeter
   ): bool where (axboun.rank == 1) && (ayboun.rank == 1) {
   const n = axboun.size - 1;
   assert(axboun.size == ayboun.size) ;
   const ref xboun = axboun.reindex(0..n);
   const ref yboun = ayboun.reindex(0..n);
   // --------------------------------------------------------------------------
   // Enforce closed perimeter.
   // --------------------------------------------------------------------------
   if ( xboun[0] != xboun[n] || yboun[0] != yboun[n] ) then {
      stderr.writeln("geo-->inxy:");
      stderr.writeln("      first and last points of perimeter must be the same");
      halt();
   }
   var theta = 0.0;
   var z1,z2: complex;
   z1 = (xboun[0] - xp,yboun[0] - yp):complex;
   for i in 1..n do {
      z2 = (xboun[i] - xp,yboun[i] - yp):complex;
      var dtheta = phase(z2/z1);
      theta += dtheta;
      z1 = z2;
   }
   if isClose(abs(theta),0.0,absTol=1.0e-5) then {
      return false;
   }
   else if isClose(abs(theta),2*pi,absTol=1.0e-5) then {
      return true;
   }
   else {
      writef("geo-->inxy:\n");
      writef("      theta = %r\n",theta);
      halt("something wrong");
   }
}
// -----------------------------------------------------------------------------
// --> AreaPts: the "area" covered by a bunch of "missing" points.
// -----------------------------------------------------------------------------
proc AreaPts(
   const in deltagrid: real,       // the gridsize
   const in C: int,                // the desired cluster id (>= 0; 0 is noise)
   const ref xmis: [] real,        // the array of abscissas of missing pts
   const ref ymis: [] real,        // the array of ordinates of missing pts
   const ref ptlabel: [] int,      // the cluster id of each point.
   out xbar: real,                 // abscissa of centroid
   out ybar: real,                 // ordinate of centroid
   out lambda1: real,              // first eigenvalue of covariance matrix
   out lambda2: real,              // second eigenvalue of covariance matrix
   out v1: (real,real),            // first eigenvector of covariance matrix
   out v2: (real,real)             // second eigenvector of covariance matrix
   ) where (xmis.rank == 1) && (ymis.rank == 1) {
   // --------------------------------------------------------------------------
   // The usual and invaluable reindexing.
   // --------------------------------------------------------------------------
   var npts = xmis.size;
   assert (npts == ymis.size && npts == ptlabel.size);
   ref xh = xmis.reindex(1..npts);
   ref yh = ymis.reindex(1..npts);
   ref lh = ptlabel.reindex(1..npts);
   // --------------------------------------------------------------------------
   // Build two arrays xclu, yclu with points belonging to cluster C.
   // --------------------------------------------------------------------------
   var dclu = {1..10};
   var nclu = 0;
   var xclu, yclu: [dclu] real;
   for i in 1..npts do {
      if lh[i] == C then {
         nclu += 1;
         dgrow(nclu,dclu);
         xclu[nclu] = xh[i];
         yclu[nclu] = yh[i];
      }
   }
   dclu = {1..nclu};
   if nclu == 0 then {
      xbar = 0.0;
      ybar = 0.0;
      lambda1 = 0.0;
      lambda2 = 0.0;
      v1 = (0.0,0.0);
      v2 = (0.0,0.0);
      return;
   }
   // --------------------------------------------------------------------------
   // The covariance matrix.
   // --------------------------------------------------------------------------
   var Cxx, Cyy: real;
   (xbar,Cxx) = stat2(xclu);
   (ybar,Cyy) = stat2(yclu);
   var Cxy = covar(xbar,ybar,xclu,yclu);
   // --------------------------------------------------------------------------
   // set the Cxy tolerance to (1/100)^2 of (minimum, if grid is
   // anisotropic) gridsize^2.
   // --------------------------------------------------------------------------
   if isClose(Cxy,0.0,absTol=deltagrid/10000) then {
      lambda1 = Cxx;
      lambda2 = Cyy;
      v1 = (1.0,0.0);
      v2 = (0.0,1.0);
      return;
   }
   else {
      var Cplus = Cxx + Cyy;
      var Cminus = Cyy - Cxx;
      var sqdel = sqrt(Cminus**2 + 4*Cxy**2);
      lambda1 = (Cplus + sqdel)/2;
      lambda2 = (Cplus - sqdel)/2;
      var v1x = 1.0;
      var v1y = (Cminus + sqdel)/(2*Cxy);
      var v2x = 1.0;
      var v2y = (Cminus - sqdel)/(2*Cxy);
      var v1m = sqrt(v1x**2 + v1y**2);
      var v2m = sqrt(v2x**2 + v2y**2);
      v1x /= v1m;
      v1y /= v1m;
      v2x /= v2m;
      v2y /= v2m;
      v1 = (v1x,v1y);
      v2 = (v2x,v2y);
      return;
   }
}
// -----------------------------------------------------------------------------
// --> AreaEllipseEV: the area of the ellipse as found by the eigenvalues (EV)
// obtained from AreaPts. Note that the EVs are variances hence the sqrt below
// to get the correct dimension.
// -----------------------------------------------------------------------------
inline proc AreaEllipseEV(
   const in lambda1: real,
   const in lambda2: real
   ): real {
   const kstd = 1.64485;
   return pi*(kstd**2)*sqrt(lambda1*lambda2);
}
   