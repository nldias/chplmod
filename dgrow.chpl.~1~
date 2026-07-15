// ===================================================================
// ==> dgrow provides a module to "grow" a domain on demand
// ===================================================================
use IO only stderr;
// -------------------------------------------------------------------
// --> dgrow: grow the first dimension of a domain d by fr if the
// index i is beyond the last element of the first dimension
// -------------------------------------------------------------------
proc dgrow(
   i: int,          // an index of the first dimension
   ref d:domain,    // the domain to be grown
   fr: real=1.5     // the growth factor (1.2 grows by 20%, etc.)
) {
  assert(fr >= 0.0);               // just in case
  var dranges = d.dims();          // the ranges that constitute d      @\label{lin:dgrow-dranges}@
  var nfirst = dranges(0).first;   // the first index of the first dim  @\label{lin:dgrow-first}@
  var nlast =  dranges(0).last;    // the last index of the first dim   @\label{lin:dgrow-last}@
  if i < nfirst then {             // just in case
     stderr.writef(
         "dgrow-->dgrow: i=%i is less than the lowest index\n",i);
     halt();
  }
  else if i <= nlast then {        // no need to grow: return!      
     return;
  }
  else if i > nlast + 1 then {     // we only want the "next" index     @\label{lin:dgrow-next}@
     stderr.writef(
        "dgrow-->dgrow: i=%i must be at most equal to %i\n",i, nlast+1);
     halt();
  }
  var dsize: real = d.shape(0);    // current size of first dimension   @\label{lin:dgrow-size}@
  dsize *= fr;                     // grow it
  var nsize = dsize:int;           // back to int
  nlast = nfirst + nsize - 1;      // calculate the new last
  dranges(0) = nfirst..nlast;      // the new range of the first dim    @\label{lin:dgrow-dranges0}@
  d = dranges;                     // resize the whole d
}
