proc split(const in x : real(?w)) // return the integral and fractional parts of x
{
   // IEEE 754 functionality
   param b = bias(x);
   param p = precision(x);
   // get IEEE 754 erncodings of x itself and |x|
   const _x = x.transmute(uint(w));
   const _a = abs(x).transmute(uint(w));
   // extract the biased exponent of x
   const e = _a >> (p - 1);
   // extract the negative bit of x and 
   // then convert it to a real(w) number
   // this yields an IEEE 754 signed zero
   const z = (_x - _a).transmute(real(w));
    
   if e < b then // |x| < 1.0
   {
      return (z, x);
   }   
   else if x != x then // x is a NaN
   {
      return (x, x);
   }   
   else if e < b + (p - 1) then // x might have a fractional component
   {
      // create a mask corresponding to the significand of x
      param _S = (~(0:uint(w))) >> (w - p + 1);
      // get the bits associated with the fractional component
      const _t = _x & (_S >> (e - b));
        
      if _t != 0 then // there is a fractional component
      {
         // strip those bits from the encoding of x 
         // and convert the net to a real(w) number
         const _r = (_x - _t).transmute(real(w));
         return (_r, x - _r);
      }
   }
   return (x, z);
}
//  IEEE 754 precision - p (also = the digits in the significand)
//  (handles T = real(w) where w = 16, 32, 64 and 128)

proc precision(type T) param where isReal(T)
{
   param w = numBits(T);
   param p = w + 2 * (w >> 7) - 3 * (w >> 5) - 5;

   return p:uint(w);
}
proc precision(x : real(?w)) param do return precision(real(w));

//  IEEE 754 exponent bias - b (also = the largest unbiased exponent)

proc bias(type T) param where isReal(T)
{
   param w = numBits(T);
   param _1 = 1:uint(w);
   param p = precision(T):int(64);

   return (_1 << (w - 1 - p)) - _1;
}
proc bias(x : real(?w)) param do return bias(real(w));

// -----------------------------------------------------------------------------
// It is a bit of black magic, but this thing calls the corresponding
// C functions correctly.
// -----------------------------------------------------------------------------

extern proc hypot(x: real, y: real): real;
extern proc hypotf(x: real(32), y: real(32)): real(32);


