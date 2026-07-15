// ===================================================================
// ==> tdma: tri-diagonal matrix algorithm
// ===================================================================
use IO;
// -------------------------------------------------------------------
// --> it3d: iterates over the 3 "diagonals"
// -------------------------------------------------------------------
iter it3d(n: int) {      // return tridiagonal indices                  @\label{lin:tdma-it3d}@
   for i in 1..n do {
      for j in max(i-1,1)..min(i+1,n) do {
         yield (i,j);               
      }
   }
}
// -------------------------------------------------------------------
// --> threediag: returns the sparse subdomain
// -------------------------------------------------------------------
proc threediag(n: int) {
   const D = {1..n,1..n};
   const SD: sparse subdomain(D) = [(i,j) in it3d(n)] (i,j);
   return SD;
}
// -------------------------------------------------------------------
// --> tridiag: the tdma algorithm
//
//     a[j]x[j-1] +   b[j]x[j] +     c[j]x[j+1] = d[j]
// A[i,j-1]x[i-1] + A[i,j]x[j] + A[i,j+1]x[j+1] = y[j]
//
// to translate from Numerical Recipes:
// a[j] == A[j-1,j]
// b[j] == A[j,j]
// c[j] == A[j,j+1]
// d[j] == y[j]
// -------------------------------------------------------------------
proc tridiag(const A: [] real, const ay: [] real, ax: []real) {
   const n = A.shape[0];           // the problem size                  @\label{lin:tdma-size}@
   assert (A.domain == threediag(n));   // is A tridiag?                @\label{lin:tdma-A-domain}@
   assert (ax.shape == (n,));           // check x shape                @\label{lin:tdma-x-shape}@
   assert (ay.shape == (n,));      // check y shape                     @\label{lin:tdma-y-shape}@
   ref x = ax.reindex({1..n});     // reindex ax->x                     @\label{lin:tdma-x-reind}@
   const ref y = ay.reindex({1..n});    // reindex ay->y                @\label{lin:tdma-y-reind}@
   var gam: [2..n] real;           // work vector                       @\label{lin:tdma-gam}@
   if A[1,1] == 0.0 then {         // algorithm can fail #1
      halt("tdma->tridiag, error 1");
   }
   var bet = A[1,1];
   x[1] = y[1]/bet;
   for j in 2..n do {              // decomposition & forward subst
      gam[j] = A[j-1,j]/bet;
      bet = A[j,j] - A[j-1,j]*gam[j];
      if bet == 0 then {           // algorithm can fail #2
         halt("tdma->tridiag, error 2");
      }
      x[j] = (y[j] - A[j-1,j]*x[j-1])/bet;
   }
   for j in 1..n-1 by -1 do {      // back substitution
      x[j] -= gam[j+1]*x[j+1];
   }
}

