// -----------------------------------------------------------------------------
// nchpl: utility stuff
// -----------------------------------------------------------------------------
use IO.FormattedIO;
inline proc reverse(
   r: range(?)) {
   return r.lowBound .. r.highBound by -r.stride;
}

inline proc span(
   const in ifirst: int,
   const in ilast: int,
   const in iinc: int =  1
   ) : range(idxType=int, bounds=boundKind.both, strides=strideKind.any) {
   if ( iinc > 0 ) then {
      return ifirst .. ilast by iinc ;
   }
   else if ( iinc < 0 ) then {
      return ilast .. ifirst by iinc ;
   }
   else {
      halt("nchpl-->span: incremment (iinc) cannot be zero\n");
   }
}

// --------------------------------------------------------------------
// --> writelnn: writeln done right
// --------------------------------------------------------------------
proc writelnn(args ...?n) {
   for x in args do {
      write(x," ");
   }
   writeln();
}

// -----------------------------------------------------------------------------
// --> clrscr: a primitive clear screen
// -----------------------------------------------------------------------------
proc clrscr() {
   var cls: bytes = "";
   cls.appendByteValues(0o33);
   cls = cls + b"[H" + cls + b"[2J"; // C string "\033[H\033[2J"
   write(cls);
}




