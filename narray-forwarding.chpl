// =============================================================================
// ==> narray-forwarding: attached domain arrays reborn
// =============================================================================
record vec {
   var dom: domain(1);                  // the rectangular rank 1 domain                
   forwarding var arr: [dom] real;      // the array
   // --------------------------------------------------------------------------
   // my version of reindex involves a lot of copying that would not be
   // necessary if Chapel had pointers.
   // --------------------------------------------------------------------------
   proc ref reindex(const in d2: domain(1)) {
      assert (d2.size == this.dom.size);
      var aux: [d2] real = this.arr;
      this.dom = d2;
      this.arr = aux;
   }
   // --------------------------------------------------------------------------
   // Initialize with a scalar: potentially, this allows to initialize with an
   // empty domain if you say
   //
   // var x: vec = 1.0;
   //
   // This is not a bug, but a feature, because it will signal the user that
   // something is not right.
   // --------------------------------------------------------------------------
   proc init=(const in a: real) {
       this.arr = a;
   }
   // --------------------------------------------------------------------------
   // a cast operator is needed for initialization and assignment to vecs
   // --------------------------------------------------------------------------
   operator :(a: real, type t: vec) {
      var v: vec = a;
      return v;
   }
   // --------------------------------------------------------------------------
   // assign real to vec
   // --------------------------------------------------------------------------
   operator =(ref lhs: vec, const in rhs: real) {
      lhs.arr = rhs;
   }
   // --------------------------------------------------------------------------
   // initialize with an array
   // --------------------------------------------------------------------------
   proc init=(const ref a: [] real) where a.rank == 1 {
      this.dom = a.domain;
      this.arr = a;
   }
   // --------------------------------------------------------------------------
   // A cast operator is needed for initialization and assignment to
   // vecs. Weirdly, the second argument is acutally the return
   // --------------------------------------------------------------------------
   operator :(a: [] real, type t: vec) where a.rank == 1 {
      var v: vec = a;
      return v;
   }
   // --------------------------------------------------------------------------
   // assign array to vec
   // --------------------------------------------------------------------------
   operator =(ref lhs: vec, const ref rhs: [] real) where rhs.rank == 1 {
      assert (lhs.size == rhs.size);
      lhs.arr = rhs;
   }
   // --------------------------------------------------------------------------
   // initialize with a vec
   // --------------------------------------------------------------------------
   proc init=(const ref a: vec) {
      this.dom = a.dom;
      this.arr = a.arr;
   }
   // --------------------------------------------------------------------------
   // assign vec to vec
   // --------------------------------------------------------------------------
   operator =(ref lhs: vec, const ref rhs: vec) {
      lhs = rhs;
   }
   // -------------------------------------------------------------------
   // operator overloading: += and -=
   // -------------------------------------------------------------------
   operator +=( const in rhs: real) {
      this.arr += rhs;
   }
   operator -=( const in rhs: real) {
      this.arr -= rhs;
   }
   operator +=(ref rhs: [] real) where rhs.rank == 1 {    // this + vec                                  
      this.arr += rhs;
   }
   operator -=(rhs: [] real) where rhs.rank == 1 {    // this - vec
      this.arr -= rhs;
   }
   operator +=(ref rhs: vec) {    // this + vec                                  
      this.arr += rhs.arr;
   }
   operator -=(rhs: vec) {    // this - vec
      this.arr -= rhs.arr ;
   }
   // -----------------------------------------------------------------------------
   // various sums returning [] real
   // -----------------------------------------------------------------------------   
   operator +(lhs: vec, rhs: vec) : [] real where rhs.rank == 1 { 
      const n = lhs.size;
      assert (n == rhs.size);
      return lhs.arr + rhs.arr;
   }
   operator +(lhs: vec, rhs: real) : [] real { 
      return lhs.arr + rhs;
   }
   operator +(lhs: real, rhs: vec) : [] real {
      return lhs + rhs.arr;
   }
   // -----------------------------------------------------------------------------
   // various subtractions returning [] real
   // -----------------------------------------------------------------------------   
   operator -(lhs: vec, rhs: vec) : [] real { 
      const n = lhs.size;
      assert (n == rhs.size);
      return lhs.arr - rhs.arr;
   }
   operator -(lhs: vec, rhs: real) : [] real { 
      return lhs.arr - rhs;
   }
   operator -(lhs: real, rhs: vec) : [] real {
      return lhs - rhs.arr;
   }
}
// -----------------------------------------------------------------------------
// subtraction between vec and [] real returning vec
// -----------------------------------------------------------------------------   
   operator -(lhs: vec, rhs: [] real) : vec where rhs.rank == 1 {
      assert (lhs.size == rhs.size);
      var ret = new vec(lhs.dom);
      ret.arr = lhs.arr - rhs;
      return ret;
   }
record mat {
   var dom: domain(2);                  // the rectangular rank 1 domain                
   forwarding var arr: [dom] real;      // the array                             
   // --------------------------------------------------------------------------
   // Initialize with a scalar: potentially, this allows to initialize with an
   // empty domain if you say
   //
   // var x: vec = 1.0;
   //
   // This is not a bug, but a feature, because it will signal the user that
   // something is not right.
   // --------------------------------------------------------------------------
   proc init=(const in a: real) {
       this.arr = a;
   }
   // --------------------------------------------------------------------------
   // a cast operator is needed for initialization and assignment to mats
   // --------------------------------------------------------------------------
   operator :(a: real, type t: mat) {
      var m: mat = a;
      return m;
   }
   // --------------------------------------------------------------------------
   // A cast operator is needed for initialization and assignment to
   // mats. Weirdly, the second argument is actually the return
   // --------------------------------------------------------------------------
   operator :(a: [] real, type t: mat) where a.rank == 2 {
      var m: mat = a;
      return m;
   }
   // --------------------------------------------------------------------------
   // assign array to mat
   // --------------------------------------------------------------------------
   operator =(ref lhs: mat, const ref rhs: [] real) where rhs.rank == 2 {
      assert (lhs.dom.shape == rhs.shape);
      lhs.arr = rhs;
   }
   // --------------------------------------------------------------------------
   // initialize with a mat
   // --------------------------------------------------------------------------
   proc init=(const ref a: mat) {
      this.dom = a.dom;
      this.arr = a.arr;
   }
   // --------------------------------------------------------------------------
   // assign real to mat
   // --------------------------------------------------------------------------
   operator =(ref lhs: mat, const in rhs: real) {
      lhs.arr = rhs;
   }
   // --------------------------------------------------------------------------
   // initialize with an array
   // --------------------------------------------------------------------------
   proc init=(const ref a: [] real) where a.rank == 2 {
      this.dom = a.domain;
      this.arr = a;
   }
   // --------------------------------------------------------------------------
   // a cast operator is needed for initialization and assignment to mats
   // --------------------------------------------------------------------------
   operator :(a: [] real, type t: mat) where a.rank == 2 {
      var v: mat = a;
      return v;
   }
   // --------------------------------------------------------------------------
   // assign array to mat
   // --------------------------------------------------------------------------
   operator =(ref lhs: mat, const in rhs: [] real) where rhs.rank == 2 {
      lhs.dom = rhs.domain;
      lhs.arr = rhs;
   }
// -------------------------------------------------------------------
// operator overloading: += and -=
// -------------------------------------------------------------------
   operator +=( const in rhs: real) {
      this.arr += rhs;
   }
   operator -=( const in rhs: real) {
      this.arr -= rhs;
   }
   operator +=(ref rhs: [] real) where rhs.rank == 1 {    // this + mat                                  
      this.arr += rhs;
   }
   operator -=(rhs: [] real) where rhs.rank == 1 {    // this - mat
      this.arr -= rhs;
   }
   operator +=(ref rhs: mat) {    // this + mat                                  
      this.arr += rhs.arr;
   }
   operator -=(rhs: mat) {    // this - mat
      this.arr -= rhs.arr ;
   }
// -----------------------------------------------------------------------------
// various sums returning [] real
// -----------------------------------------------------------------------------   
   operator +(lhs: mat, rhs: mat) : [] real { 
      const n = lhs.size;
      assert (n == rhs.size);
      return lhs.arr + rhs.arr;
   }
   operator +(lhs: mat, rhs: real) : [] real { 
      return lhs.arr + rhs;
   }
   operator +(lhs: real, rhs: mat) : [] real {
      return lhs + rhs.arr;
   }
// -----------------------------------------------------------------------------
// various subtractions returning [] real
// -----------------------------------------------------------------------------   
   operator -(lhs: mat, rhs: mat) : [] real { 
      const n = lhs.size;
      assert (n == rhs.size);
      return lhs.arr - rhs.arr;
   }
   operator -(lhs: mat, rhs: real) : [] real { 
      return lhs.arr - rhs;
   }
   operator -(lhs: real, rhs: mat) : [] real {
      return lhs - rhs.arr;
   }
/*   
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
   //    where lhs.rank == 1 {                  // [] real - narray           @\label{lin:sub-narray-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs - rhs.arr ;
   //    return r;
   // }
   // operator -(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray - [] real           @\label{lin:sub-narray-narray}@
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
   //    where lhs.rank == 1 {                  // [] real * narray           @\label{lin:mul-narray-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs*rhs.arr;
   //    return r;
   // }
   // operator *(lhs: narray, rhs: real): narray
   //    where rhs.rank == 1 {                  // narray * [] real           @\label{lin:mul-narray-narray}@
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
   //    where lhs.rank == 1 {                  // [] real / narray           @\label{lin:div-narray-narray}@
   //    const n = lhs.size;
   //    assert (n == rhs.size);
   //    var r = new narray(lhs.domain);
   //    r.arr = lhs/rhs.arr;
   //    return r;
   // }
   // operator /(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray / [] real           @\label{lin:div-narray-narray}@
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
      where lhs.rank == 1 {                  // [] real**narray            @\label{lin:exp-narray-narray}@
      const n = lhs.size;
      assert (n == rhs.size);
      var r = new narray(rhs.dom);
      r.arr = lhs**rhs.arr;
      return r;
   }
   // operator **(lhs: narray, rhs: [] real): narray
   //    where rhs.rank == 1 {                  // narray** [] real           @\label{lin:exp-narray-narray}@
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
*/
}
