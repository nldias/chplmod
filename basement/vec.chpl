// ===================================================================
// ==> vec: a variable-length 1d array that holds its own domain
// ===================================================================
record vec {
   var dom: domain(1);                  // the 1d domain                @\label{lin:vec-dom}@
   var arr: [dom] real;                 // the array                    @\label{lin:vec-arr}@
   var vfirst = dom.first;              // the first index after reind
   var vlast = dom.last;                // the last index after reind
   var vdelta = 0;                      // the array shift
   proc size: int {                     // the size of a vec
      return dom.size;
   }
   proc ref reindex(
      const in dv: range(int)
      ) {
      assert ( dv.size == dom.size );
      vfirst = dv.first;
      vlast = dv.last;
      vdelta = vfirst - dom.first;
   }
   // proc ref this() ref {                // access the whole array       @\label{lin:this-arr}@
   //    return arr;
   // }
   proc ref this(in k: int) ref {       // access arr[k]                @\label{lin:this-ark}@
      return arr[k-vdelta];
   }
   // iter ref these() ref {               // iterate over vec             @\label{lin:these-vec}@
   //    for i in dom do {
   //       yield arr[i-vdelta];
   //    }
   // }
   // iter ref these(                      // iterate over a range         @\label{lin:these-range}@
   //    in rak: range(int)
   //    ) ref {  
   //    for k in rak do {
   //       yield arr[k-vdelta];
   //    }
   // }
// -------------------------------------------------------------------
// operator overloading: *all* of the subsequent operators return
// vec. if you want to assign an operation between an array b and a
// vec c to an array a, just say: a = b + c.arr
//
// the domain of the vec returned is always the same as the domain of
// the lhs, except when the lhs is a scalar
// -------------------------------------------------------------------
   operator +(lhs: real, rhs: vec): vec {    // real + vec              @\label{lin:sum-real-vec}@
      var r = new vec(rhs.dom);
      r.arr = lhs + rhs.arr;
      return r;
   }
   operator +(lhs: vec, rhs: real): vec {    // vec + real              @\label{lin:sum-vec-real}@
      var r = new vec(lhs.dom);
      r.arr = lhs.arr + rhs;
      return r;
   }
   // operator +(lhs: [] real, rhs: vec): vec
   //    where lhs.rank == 1 {                  // [] real + vec           @\label{lin:sum-array-vec}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.domain);
   //    r.arr = lhs + rhs.arr;
   //    return r;
   // }
   // operator +(lhs: vec, rhs: [] real): vec
   //    where rhs.rank == 1 {                  // vec + [] real           @\label{lin:sum-vec-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.dom);
   //    r.arr = lhs.arr + rhs;
   //    return r;
   // }
   operator +(lhs: vec, rhs: vec): vec {     // vec + vec               @\label{lin:sum-vec-vec}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new vec(lhs.dom);
      r.arr = lhs.arr + rhs.arr ;
      return r;
   }
   operator -(lhs: real, rhs: vec): vec {    // real - vec              @\label{lin:sub-real-vec}@
      var r = new vec(rhs.dom);
      r.arr = lhs - rhs.arr ;
      return r;
   }
   operator -(lhs: vec, rhs: real): vec {    // vec - real              @\label{lin:sub-vec-real}@
      var r = new vec(lhs.dom);
      r.arr = lhs.arr - rhs;
      return r;
   }
   // operator -(lhs: [] real, rhs: vec): vec
   //    where lhs.rank == 1 {                  // [] real - vec           @\label{lin:sub-array-vec}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.domain);
   //    r.arr = lhs - rhs.arr ;
   //    return r;
   // }
   // operator -(lhs: vec, rhs: [] real): vec
   //    where rhs.rank == 1 {                  // vec - [] real           @\label{lin:sub-vec-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.dom);
   //    r.arr = lhs.arr - rhs;
   //    return r;
   // }
   operator -(lhs: vec, rhs: vec): vec {     // vec - vec               @\label{lin:sub-vec-vec}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new vec(lhs.dom);
      r.arr = lhs.arr - rhs.arr;
      return r;
   }
   operator *(lhs: real, rhs: vec): vec {    // real * vec              @\label{lin:mul-real-vec}@
      var r = new vec(rhs.dom);
      r.arr = lhs*rhs.arr;
      return r;
   }
   operator *(lhs: vec, rhs: real): vec {    // vec * real              @\label{lin:mul-vec-real}@
      var r = new vec(lhs.dom);
      r.arr = vec.arr*rhs;
      return r;
   }
   // operator *(lhs: [] real, rhs: vec): vec
   //    where lhs.rank == 1 {                  // [] real * vec           @\label{lin:mul-array-vec}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.domain);
   //    r.arr = lhs*rhs.arr;
   //    return r;
   // }
   // operator *(lhs: vec, rhs: real): vec
   //    where rhs.rank == 1 {                  // vec * [] real           @\label{lin:mul-vec-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.dom);
   //    r.arr = vec.arr*rhs;
   //    return r;
   // }
   operator *(lhs: vec, rhs: vec): vec {     // vec * vec               @\label{lin:mul-vec-vec}@
      assert (lhs.size == rhs.size);
      var r = new vec(lhs.dom);
      r.arr = lhs.arr*rhs.arr;
      return r;
   }
   operator /(lhs: real, rhs: vec): vec {    // real / vec              @\label{lin:div-real-vec}@
      var r = new vec(rhs.dom);
      r.arr = lhs/rhs.arr;
      return r;
   }
   operator /(lhs: vec, rhs: real): vec {    // vec / real              @\label{lin:div-vec-real}@
      var r = new vec(lhs.dom);
      r.arr = lhs.arr/rhs;
      return r;
   }
   // operator /(lhs: [] real, rhs: vec): vec
   //    where lhs.rank == 1 {                  // [] real / vec           @\label{lin:div-array-vec}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.domain);
   //    r.arr = lhs/rhs.arr;
   //    return r;
   // }
   // operator /(lhs: vec, rhs: [] real): vec
   //    where rhs.rank == 1 {                  // vec / [] real           @\label{lin:div-vec-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.dom);
   //    r.arr = lhs.arr/rhs;
   //    return r;
   // }
   operator /(lhs: vec, rhs: vec): vec {     // vec / vec               @\label{lin:div-vec-vec}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new vec(lhs.dom);
      r.arr = lhs.arr/rhs.arr;
      return r;
   }
   operator **(lhs: vec, rhs: int): vec {    // vec**int                @\label{lin:exp-vec-int}@
      var r = new vec(lhs.dom);
      r.arr = lhs.arr**rhs;
      return r;
   }
   operator **(lhs: real, rhs: vec): vec {   // real**vec               @\label{lin:exp-real-vec}@
      var r = new vec(rhs.dom);
      r.arr = lhs**rhs.arr;
      return r;
   }
   operator **(lhs: vec, rhs: real): vec {   // vec**real               @\label{lin:exp-vec-real}@
      var r = new vec(lhs.dom);
      r.arr = lhs.arr**rhs;
      return r;
   }
   operator **(lhs: [] real, rhs: vec): vec
      where lhs.rank == 1 {                  // [] real**vec            @\label{lin:exp-array-vec}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new vec(rhs.dom);
      r.arr = lhs**rhs.arr;
      return r;
   }
   // operator **(lhs: vec, rhs: [] real): vec
   //    where rhs.rank == 1 {                  // vec** [] real           @\label{lin:exp-vec-array}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new vec(lhs.dom);
   //    r.arr = lhs.arr**rhs;
   //    return r;
   // }
   operator **(lhs: vec, rhs: vec): vec {    // exponentiation          @\label{lin:exp-vec-vec}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new vec(lhs.dom);
      r.arr = lhs.arr**rhs.arr;
      return r;
   }            
}
proc tovec(x: [] real): vec where x.rank == 1 {
   var v = new vec(x.domain);
   v.arr = x;
   return v;
}