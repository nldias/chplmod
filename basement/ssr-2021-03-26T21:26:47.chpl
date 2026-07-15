// ===================================================================
// ssr is a set of search and sort functions
//
// Nelson Luis Dias
// 19900000 (circa)
// 20060421 (today)
// 2021-03-20T14:44:19 now this is today: Chapel!
// ===================================================================

// -------------------------------------------------------------------
// --> indxsort: sorts an array of floats, using the heap algorithm,
// and sorting by index, stored in indx
//
// adapted from heapsort (hpsort) in Numerical Recipes
// -------------------------------------------------------------------
proc indxsort(
ref ax: [],
ref aindx: [] int
) {
// -------------------------------------------------------------------
// painful reindexing
// -------------------------------------------------------------------
   assert (ax.rank == 1);
   var n = ax.size;
   ref x = ax.reindex(1..n);
   assert (aindx.rank == 1);
   ref indx = aindx.reindex(1..n);
   var qq: ax.eltType;        // to be truly generic
// -------------------------------------------------------------------
// index is started from 1 to n
// -------------------------------------------------------------------
   for j in 1..n do {
      indx[j] = j;
   }
// -------------------------------------------------------------------
// beginning of algorithm proper
// -------------------------------------------------------------------
   var indxt: int;
   var l = n/2 + 1 ;
   var ir = n ;
   while true do {
      if (l > 1) then {
         l -= 1;
         indxt = indx[l];
         qq = x[indxt];
      }
      else {
         indxt=indx[ir];
         qq = x[indxt];
         indx[ir] = indx[1];
         ir -= 1;
         if ir == 1 then {
            indx[1] = indxt ;
            break ;
         }
      }
      var i = l;
      var j = l+1;
      while j <= ir do {
         if ((j < ir) && (x[indx[j]] < x[indx[j+1]])) then {
            j += 1;
         }
         if (qq < x[indx[j]]) then {
            indx[i] = indx[j];
            i = j;
            j = j+j;
         }
         else {
            j=ir+1;
         }
      }
      indx[i] = indxt;
   }
// -------------------------------------------------------------------
// we still need to fix the indices
// -------------------------------------------------------------------
   const del = ax.indices.first - 1;
   indx += del;
   return;
}

// -------------------------------------------------------------------
// --> countval: count how many times val occurs in x
// -------------------------------------------------------------------
proc countval(val: ?tv, x: [] ?tx): int {
// --------------------------------------------------------------------
// be careful with empty arrays
// --------------------------------------------------------------------
   assert (x.rank == 1);
   assert (tx == tv);
   var n = x.size;
   if n == 0 then {
      halt("ssr-->countval: empty array");
   }
   var count = 0;
// -------------------------------------------------------------------
// be careful with NANs
// -------------------------------------------------------------------
   if tv == real && isnan(val) then {
      for e in x do {
         if isnan(e) then count += 1;
      }
   }
   else {
      for e in x do {
         if e == val then count +=1 ;
      }
   }
   return count;
}

// -------------------------------------------------------------------
// butsum: sum all elements != val in array
// -------------------------------------------------------------------
proc butsum(val: ?tv, x: [] ?tx): tv {
   assert (x.rank == 1);
   assert (tx == tv);
   var n = x.size;
   if n == 0 then {
      halt("ssr-->countval: empty array");
   }
   var sum: tv = 0;
   if tv == real && isnan(val) then {
      for e in x do {
         if !isnan(e) then sum += e;
      }
   }
   else {
      for e in x do {
         if e != val then sum += e;
      }
   }
   return sum;
}
// -------------------------------------------------------------------
// --> whereval: where val occurs in x
//
// important: whereval *always* returns a 0-based array
// -------------------------------------------------------------------
use dgrow;
proc whereval(val, x: []): [] int {
// --------------------------------------------------------------------
// be careful with empty arrays
// --------------------------------------------------------------------
   assert (x.rank == 1);
   assert (x.eltType == val.type);
   var n = x.size;
   if n == 0 then {
      halt("--> whereval: empty array");
   }
   const xf = x.domain.first; // try to be agnostic
   const xl = x.domain.last;  // try to be agnostic
   const m = max(n/10,2);     // guess that 10% of x elements == val
   var dw = {0..#m};          // ... but at least 2
   var ww: [dw] int;          // where they are
   var ct = 0;                // count how many
   if val.type == real && isnan(val) then {  // be careful with NANs
      for i in xf..xl do {
         if isnan(x[i]) then {
            dgrow(ct,dw);
            ww[ct] = i;
            ct += 1;
         }
      }
   }
   else {                                    // just find them
      for i in xf..xl do {
         if x[i] == val then {
            dgrow(ct,dw);
            ww[ct] = i;
            ct += 1;
         }
      }
   }
   dw = {0..#ct};                   // adjust domain
   return ww;
}
// -------------------------------------------------------------------
// --> diff: calculates the 1st discrete difference of a 1D array
// bool--int "magic" is used!
//
// important: diff *always* returns a 0-based array
// -------------------------------------------------------------------
proc diff(ref ax: [] ?tx) {
   assert (ax.rank == 1);     // must be 1D
   var n = ax.size;           // count elements
   ref x = ax.reindex(1..n);  // reindex
// -------------------------------------------------------------------
// return an array with a domain that is compatible with ax's domain
// -------------------------------------------------------------------   
   type td;
   if tx == bool then {
      td = int;
   }
   else {
      td = tx;
   }
   var dd = {0..#(n-1)};      // always return a 0-based array
   var dx: [dd] td;           // the return array
   dx = x[2..n] - x[1..n-1];  // differentiated
   return dx;                 // end of the story
}

// -------------------------------------------------------------------
// --> linspace: my equivalent of a (simple!) numpy linspace. returns
// a 0-based 1D array with n linearly interpolated values including,
// and between, start and stop.
// -------------------------------------------------------------------
proc linspace(
   start: real,               // the first value
   stop: real,                // the last value
   n: int                     // how many do you want?
): [] real {
   var x: [0..#n] real;
   assert( n > 1);
   var dx = (stop-start)/(n-1);
   forall i in 0..#n do {     // this is a parallel algorithm!
      x[i] = start + i*dx;
   }
   return x;
}

// -------------------------------------------------------------------
// --> flip: flips a 1D array
// -------------------------------------------------------------------
proc flip(ref ax: [] ) {
   assert(ax.rank == 1);
   var n = ax.size;
   ref x = ax.reindex(0..#n);
   for i in 0..n/2 do {
      x[i] <=> x[n-1-i];
   }
}
   
