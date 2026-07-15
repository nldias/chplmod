// ===================================================================
// ==> ada: attached domain arrays vec (1D) and mat (2D)
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
   proc compare(i,j) {
      return arr[i] - arr[j] ;
   }
}

record mat {
   var dom: domain(2);                  // the 2d domain                @\label{lin:vec-dom}@
   var arr: [dom] real;                 // the array                    @\label{lin:vec-arr}@
   var v0first = dom.dim(0).first;              // the first index after reind
   var v1first = dom.dim(1).first;
   var v0delta = 0;                      // the array shift
   var v1delta = 0;
   proc size: int {                     // the size of a mat
      return dom.size;
   }
   proc shape: 2*int {
      return dom.shape;
   }
   proc ref reindex(
      const in dv0: range(int),
      const in dv1: range(int)
      ) {
      assert ( dom.shape == (dv0.size,dv1.size) );
      v0first = dv0.first;
      v0delta = v0first - dom.dim(0).first;
      v1first = dv1.first;
      v1delta = v1first - dom.dim(1).first;
   }
   proc ref this(in k: int, in l: int) ref {       // access arr[k]                @\label{lin:this-ark}@
      return arr[k-v0delta,l-v1delta];
   }
// -------------------------------------------------------------------
// operator overloading: *all* of the subsequent operators return
// mat. if you want to assign an operation between an array b and a
// mat c to an array a, just say: a = b + c.arr
//
// the domain of the mat returned is always the same as the domain of
// the lhs, except when the lhs is a scalar
// -------------------------------------------------------------------
   operator +(lhs: real, rhs: mat): mat {    // real + mat              @\label{lin:sum-real-mat}@
      var r = new mat(rhs.dom);
      r.arr = lhs + rhs.arr;
      return r;
   }
   operator +(lhs: mat, rhs: real): mat {    // mat + real              @\label{lin:sum-mat-real}@
      var r = new mat(lhs.dom);
      r.arr = lhs.arr + rhs;
      return r;
   }
   operator +(lhs: mat, rhs: mat): mat {     // mat + mat               @\label{lin:sum-mat-mat}@
      assert (lhs.shape== rhs.shape);
      var r = new mat(lhs.dom);
      r.arr = lhs.arr + rhs.arr ;
      return r;
   }
   operator -(lhs: real, rhs: mat): mat {    // real - mat              @\label{lin:sub-real-mat}@
      var r = new mat(rhs.dom);
      r.arr = lhs - rhs.arr ;
      return r;
   }
   operator -(lhs: mat, rhs: real): mat {    // mat - real              @\label{lin:sub-mat-real}@
      var r = new mat(lhs.dom);
      r.arr = lhs.arr - rhs;
      return r;
   }
   operator -(lhs: mat, rhs: mat): mat {     // mat - mat               @\label{lin:sub-mat-mat}@
      assert (lhs.shape == rhs.shape);
      var r = new mat(lhs.dom);
      r.arr = lhs.arr - rhs.arr;
      return r;
   }
}



proc tovec(x: [] real): vec where x.rank == 1 {
   var v = new vec(x.domain);
   v.arr = x;
   return v;
}

proc tomat(x: [] real): vec where x.rank == 2 {
   var m = new mat(x.domain);
   m.arr = x;
   return m;
}