proc dot_mv(aA: [] real, ax: [] real, ay: [] real) {
   const n = aA.shape[0];
   assert (aA.shape == (n,n));
   assert (ax.shape == (n,));
   assert (ay.shape == (n,));
   ref A = aA.reindex({1..n,1..n});
   ref x = ax.reindex({1..n});
   ref y = ay.reindex({1..n});
   forall i in 1..n do {
      var sum = 0.0;
      for j in 1..n do {
         sum += A[i,j]*x[j];
      }
      y[i] = sum;
   }
}
