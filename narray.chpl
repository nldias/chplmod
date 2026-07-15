// =============================================================================
// ==> narray: n-dimensional arrays.
// Still in development as of 2025-11-12T17:47:38
// =============================================================================
record narray {
   param ran;                  // the rank
   var dom: domain(ran);           // the domain                                 @\label{lin:narray-dom}@
   var arr: [dom] real;            // the array                                  @\label{lin:narray-arr}@
   // var vfirst = dom.first;              // the first index after reind
   // var vlast = dom.last;                // the last index after reind
   // var vdelta = 0;                      // the array shift
   proc size: int {                // the size of a narray
      return dom.size;
   }
   // proc ref reindex(
   //    const in dv: range(int)
   //    ) {
   //    assert ( dv.size == dom.size );
   //    vfirst = dv.first;
   //    vlast = dv.last;
   //    vdelta = vfirst - dom.first;
   // }
   inline proc ref this(in k: int...?n) ref {     // access arr[k]               @\label{lin:this-ark}@
      return arr[k];
   }
   iter ref these() ref {               // iterate over the whole narray         @\label{lin:these-narray}@
      for i in dom do {
         yield arr[i];
      }
   }
   iter ref these(                      // iterate over several ranges           @\label{lin:these-range}@
      in rak: range(int)...?r
      ) ref {  
      for k in zip(rak) do {
         yield arr[k];
      }
   }
   proc init=(const in rhs: narray(?)) {
     this.ran = rhs.ran;
     this.dom = rhs.dom;
     this.arr = rhs.arr;
   }
   // operator :(a: narray(?), type t: narray(a.ran)) {
   //    var v: t = a;
   //    return v;
   // }
   operator =(ref lhs: narray(?), rhs: narray(lhs.ran)) {
     lhs.dom = rhs.dom;
     lhs.arr = rhs.arr;
   }   
   proc init=(const in a: real) {
      this.ran = 1;
      this.dom = {1..1};
      this.arr = a;
   }
   // --------------------------------------------------------------------------
   // a cast operator is needed for initialization and assignment from reals
   // --------------------------------------------------------------------------
   operator :(a: real, type t: narray(?)) {
      var v: t = a;
      return v;
   }
   // --------------------------------------------------------------------------
   // assign real to narray
   // --------------------------------------------------------------------------
   operator =(ref lhs: narray(?), const in a: real) {
      lhs.arr = a;
   }
}
/*
   /*
   // --------------------------------------------------------------------------
   // initialize with an array
   // --------------------------------------------------------------------------
   //    proc init=(const ref a: [] real) {
   //       this.ran = a.rank;
   //       this.dom = a.domain;
   //       this.arr = a;
   //    }
   // // --------------------------------------------------------------------------
   // // a cast operator is needed for initialization and assignment from arrays
   // // --------------------------------------------------------------------------
   // operator :(a: [] real, type t: narray) {
   //    var v: narray(a.rank,a.domain);
   //    v.arr = a;
   //    return v;
   // }
   */
// -------------------------------------------------------------------
// operator overloading: *all* of the subsequent operators return
// narray. if you want to assign an operation between an array b and a
// narray c to an array a, just say: a = b + c.arr
//
// the domain of the narray returned is always the same as the domain of
// the lhs, except when the lhs is a scalar
// -------------------------------------------------------------------
   operator +=(ref rhs: narray(this.ran)) {    // this + narray                                  
      this.arr += rhs.arr;
   }
   operator -=(ref rhs: narray) {    // this - narray
      this.arr -= rhs.arr ;
   }
   operator +(lhs: real, rhs: narray): narray(rhs.ran) {    // real + narray              @\label{lin:sum-real-narray}@
      var r = new narray(rhs.ran,rhs.dom);
      r.arr = lhs + rhs.arr;
      return r;
   }
   operator +(lhs: narray, rhs: real): narray {    // narray + real              @\label{lin:sum-narray-real}@
      var r = new narray(lhs.dom);
      r.arr = lhs.arr + rhs;
      return r;
   }
   // operator +(lhs: [] real, rhs: narray): narray
   //    where lhs.rank == 1 {                  // [] real + narray           @\label{lin:sum-array-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs + rhs.arr;
   //    return r;
   // }
   // operator +(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray + [] real           @\label{lin:sum-narray-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.dom);
   //    r.arr = lhs.arr + rhs;
   //    return r;
   // }
   operator +(lhs: narray, rhs: narray): narray {     // narray + narray               @\label{lin:sum-narray-narray}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new narray(lhs.dom);
      r.arr = lhs.arr + rhs.arr ;
      return r;
   }
   operator -(lhs: real, rhs: narray): narray {    // real - narray              @\label{lin:sub-real-narray}@
      var r = new narray(rhs.dom);
      r.arr = lhs - rhs.arr ;
      return r;
   }
   operator -(lhs: narray, rhs: real): narray {    // narray - real              @\label{lin:sub-narray-real}@
      var r = new narray(lhs.dom);
      r.arr = lhs.arr - rhs;
      return r;
   }
   // operator -(lhs: [] real, rhs: narray): narray
   //    where lhs.rank == 1 {                  // [] real - narray           @\label{lin:sub-array-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs - rhs.arr ;
   //    return r;
   // }
   // operator -(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray - [] real           @\label{lin:sub-narray-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.dom);
   //    r.arr = lhs.arr - rhs;
   //    return r;
   // }
   operator -(lhs: narray, rhs: narray): narray {     // narray - narray               @\label{lin:sub-narray-narray}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new narray(lhs.dom);
      r.arr = lhs.arr - rhs.arr;
      return r;
   }
   operator *(lhs: real, rhs: narray): narray {    // real * narray              @\label{lin:mul-real-narray}@
      var r = new narray(rhs.dom);
      r.arr = lhs*rhs.arr;
      return r;
   }
   operator *(lhs: narray, rhs: real): narray {    // narray * real              @\label{lin:mul-narray-real}@
      var r = new narray(lhs.dom);
      r.arr = narray.arr*rhs;
      return r;
   }
   // operator *(lhs: [] real, rhs: narray): narray
   //    where lhs.rank == 1 {                  // [] real * narray           @\label{lin:mul-array-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs*rhs.arr;
   //    return r;
   // }
   // operator *(lhs: narray, rhs: real): narray
   //    where rhs.rank == 1 {                  // narray * [] real           @\label{lin:mul-narray-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.dom);
   //    r.arr = narray.arr*rhs;
   //    return r;
   // }
   operator *(lhs: narray, rhs: narray): narray {     // narray * narray               @\label{lin:mul-narray-narray}@
      assert (lhs.size == rhs.size);
      var r = new narray(lhs.dom);
      r.arr = lhs.arr*rhs.arr;
      return r;
   }
   operator /(lhs: real, rhs: narray): narray {    // real / narray              @\label{lin:div-real-narray}@
      var r = new narray(rhs.dom);
      r.arr = lhs/rhs.arr;
      return r;
   }
   operator /(lhs: narray, rhs: real): narray {    // narray / real              @\label{lin:div-narray-real}@
      var r = new narray(lhs.dom);
      r.arr = lhs.arr/rhs;
      return r;
   }
   // operator /(lhs: [] real, rhs: narray): narray
   //    where lhs.rank == 1 {                  // [] real / narray           @\label{lin:div-array-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs/rhs.arr;
   //    return r;
   // }
   // operator /(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray / [] real           @\label{lin:div-narray-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.dom);
   //    r.arr = lhs.arr/rhs;
   //    return r;
   // }
   operator /(lhs: narray, rhs: narray): narray {     // narray / narray               @\label{lin:div-narray-narray}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new narray(lhs.dom);
      r.arr = lhs.arr/rhs.arr;
      return r;
   }
   operator **(lhs: narray, rhs: int): narray {    // narray**int                @\label{lin:exp-narray-int}@
      var r = new narray(lhs.dom);
      r.arr = lhs.arr**rhs;
      return r;
   }
   operator **(lhs: real, rhs: narray): narray {   // real**narray               @\label{lin:exp-real-narray}@
      var r = new narray(rhs.dom);
      r.arr = lhs**rhs.arr;
      return r;
   }
   operator **(lhs: narray, rhs: real): narray {   // narray**real               @\label{lin:exp-narray-real}@
      var r = new narray(lhs.dom);
      r.arr = lhs.arr**rhs;
      return r;
   }
   operator **(lhs: [] real, rhs: narray): narray
      where lhs.rank == 1 {                  // [] real**narray            @\label{lin:exp-array-narray}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new narray(rhs.dom);
      r.arr = lhs**rhs.arr;
      return r;
   }
   // operator **(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray** [] real           @\label{lin:exp-narray-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.dom);
   //    r.arr = lhs.arr**rhs;
   //    return r;
   // }
   operator **(lhs: narray, rhs: narray): narray {    // exponentiation          @\label{lin:exp-narray-narray}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new narray(lhs.dom);
      r.arr = lhs.arr**rhs.arr;
      return r;
   }            
   proc compare(i,j) {
      return arr[i] - arr[j] ;
   }
}

*/