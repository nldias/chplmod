// =============================================================================
// ==> unititme: a unit for cpu time
// =============================================================================
use Time only stopwatch;
use Random only fillRandom;
use smatrix only dot_mm;
const n = 250;
proc utime(): real {
   var
      a,
      b,
      c: [1..n,1..n] real;
   var runtime: stopwatch;
   for k in 1..100 do {
      runtime.start();
      fillRandom(a);
      fillRandom(b);
      dot_mm(a,b,c);
      runtime.stop();
   }
   return runtime.elapsed();
}
