use vec;
use smatrix;
// ------------------------------------------------------------------------------
// --> steep: nonlinear least squares by curve fitting with the steepest descent
// method
//
// Hopefully, backtracking will come from:
// https://www.cs.cmu.edu/~ggordon/10725-F12/slides/05-gd-revisited.pdf
// ------------------------------------------------------------------------------
proc steep(
   const ref x: vec,          // ind variables (used as arg to func) (m x 1)
   const ref y: vec,          // data to be fit by func(x,p,y) (m x 1)
   const ref w: vec,          // array, *not matrix*, of weights (m x 1)
   ref p: vec,                // initial guess of parameter values  (n x 1)
                              // returns the estimated parameters
   ref asigp: [] real,             // standard  errors of the parameters
   ref acp: [] real,               // parameter covariance matrix
   const func: proc(ref _: vec,    // the independent variables
                    ref _: vec,    // the parameters
                    ref _: vec),   // in the sim model we call func(x,p,y)
   const in epsilon = 1.0e-6       // stop criterion
   ) : (real,real,real)
   where ( asigp.rank == 1 && acp.rank == 2) {
   const maxiter = 100000;         // steep may take a looong time to converge
   const m = x.size;              // the number of data points
   const n = p.size;              // the number of parameters
// -----------------------------------------------------------------------------
// check all shapes and sizes
// -----------------------------------------------------------------------------   
   assert (y.size == m);
   assert (w.size == m);
   assert (asigp.size == n);
   assert (acp.shape == (n,n));
// -----------------------------------------------------------------------------
// reindexing
// -----------------------------------------------------------------------------   
   x.reindex(1..m);
   y.reindex(1..m);
   w.reindex(1..m);
   p.reindex(1..n);
   ref sigp = asigp.reindex(1..n);
   ref cp = acp.reindex({1..n,1..n});
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
   if allequalto(w.arr,-1.0) then {
      w.arr = 1.0;
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
      dely = yhat - y;                  // [hat{y}-y]
      simplejacob();                    // update the jacobian matrix
      dot_diagm_v(w.arr,dely,vaux_m);   // W[hat{y} - y]
      dot_mtv(J,vaux_m,gradchi2);       // J'W[hat{y} - y]
      gradchi2 *= 2;                    // 2J'W[hat{y} - y]
// -----------------------------------------------------------------------------
// backtracking line search
// -----------------------------------------------------------------------------      
      var modgr2 = dot_vtv(gradchi2,gradchi2);    // |grad chi^2|
      var t = 1.0;                                // backtracing parameter
      var pa = new vec({1..n});                   // yet another vec!
      pa.arr = p.arr - t*gradchi2;                // check two chi^2s
      var chi2a = chi2(pa);
      var chi2b = chi2(p) - (t/2)*modgr2;
// -----------------------------------------------------------------------------
// backtracing loop
// -----------------------------------------------------------------------------      
      while chi2a > chi2b do {
         t *= 0.1;
         pa.arr = p.arr - t*gradchi2;
         chi2a = chi2(pa);
         chi2b = chi2(p) - (t/2)*modgr2;
      }
// -----------------------------------------------------------------------------
// found the right scaling for the size of the step
// -----------------------------------------------------------------------------      
      hh = t*gradchi2;
      p = p - hh ;
      eps = sqrt(dot_vtv(hh,hh));
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
   mvdiag(cp,sigp);                 // diagonal of [J'WJ]^{-1}
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
      delya = ya - y;
      dot_diagm_v(w.arr,delya,vaux_m);
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
         J[1..m,k] = (yplus[1..m] - yminus[1..m])/(2*delp[k]);
      }
   }
}
