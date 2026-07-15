// =============================================================================
// ==> geo: things that vary spatially
//
// Nelson Luís Dias
// 2026-02-04T11:01:06 a new star is born
// =============================================================================
use IO only stderr, writeln;
use Math only phase,pi,isClose;
// -----------------------------------------------------------------------------
// --> findquad: finds the quadrant to which the point/vector (delx,dely)
// belongs
// -----------------------------------------------------------------------------
private proc findquad(
   const in delx: real,       // abscissa of point
   const in dely: real        // ordinate of point
   ): int {                   // quadrant of point (0, 1, 2, 3)
   var idx = sgn(delx);
   var idy = sgn(dely);
   var qpoint: int ;          // quadrant of this point
   // --------------------------------------------------------------------------
   // find the quadrant of this point
   // --------------------------------------------------------------------------
   select (idx,idy) {
      // -----------------------------------------------------------------------
      // most frequent cases first
      // -----------------------------------------------------------------------
      when (1,1)   do qpoint = 0;
      when (-1,1)  do qpoint = 1;
      when (-1,-1) do qpoint = 2;
      when (1,-1)  do qpoint = 3;
      // -----------------------------------------------------------------------
      // cases with zero are much less frequent; therefore, test them lastly
      // -----------------------------------------------------------------------
      when (1,0)  do qpoint = 0;
      when (0,1)  do qpoint = 1;
      when (-1,0) do qpoint = 2;
      when (0,-1) do qpoint = 3;
      // -----------------------------------------------------------------------
      // the rarest case? (0,0) is arbitrarily assigned to quadrant 0
      // -----------------------------------------------------------------------
      when (0,0) do qpoint = 0;
      otherwise {
         halt("this cannot happen");
      }
   }
   return qpoint;
}
// -----------------------------------------------------------------------------
// --> inxy: is the point (xp,yp) inside perimeter (axboun,ayboun)?
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
      stderr.writeln("first and last points of perimeter must be the same");
      halt();
   }
   var qsum = 0;                             // sum of quadrants
   var oldquad = findquad(xboun[0]-xp,yboun[0]-yp);
                                             // quadrant of first/current point
   var newquad: int;                         // quadrant of next point
   for i in 1..n do {
      newquad = findquad(xboun[i]-xp,yboun[i]-yp);
      var delquad = newquad - oldquad ;
      // -----------------------------------------------------------------------
      // adjust for 0->3 and 3->0 transitions
      // -----------------------------------------------------------------------
      if delquad == 3 then {
         delquad = -1;
      }
      else if delquad == -3 then {
         delquad = +1;
      }
      // -----------------------------------------------------------------------
      // quadrant arithmetic!
      // -----------------------------------------------------------------------
      qsum += delquad;
      // writeln("old, new, del, qsum: ",oldquad,newquad,delquad,qsum,sep=" ");
      oldquad = newquad ;
   }
   // writeln("qsum = ",qsum);
   return (qsum != 0);
}
// -----------------------------------------------------------------------------
// --> zinxy: the standard way of finding inside points
// -----------------------------------------------------------------------------
proc zinxy(
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
      stderr.writeln("first and last points of perimeter must be the same");
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
   if isClose(theta,0.0,absTol=1.0e-5) then {
      return false;
   }
   else if isClose(theta,2*pi,absTol=1.0e-5) then {
      return true;
   }
   else {
      writef("theta = %r\n",theta);
      halt("something wrong");
   }
}
      