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
