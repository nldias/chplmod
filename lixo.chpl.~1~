// ===================================================================
// ==> qdran: a library for random number generation
// ===================================================================
private var next: uint(32) = 1;    // next exists between calls
private param maxx = max(uint(32));// to normalize betwenn 0 and 1
// -------------------------------------------------------------------
// --> ranqd: generates the next pseudorandom number from current
// -------------------------------------------------------------------
proc ranqd(): uint(32) {
   next = next*1664525 + 1013904223;
   return next;
}
// -------------------------------------------------------------------
// --> uranqd: normalizes to real between 0.0 and 1.0 (generates a
//     "uniform" variable ~U(0,1) )
// -------------------------------------------------------------------
proc uranqd(): real {
   next = next*1664525 + 1013904223;
   return (next:real)/maxx ;
}
// -------------------------------------------------------------------
// --> seeqd: seeds the sequence, changing global variable next
// -------------------------------------------------------------------
proc seedqd(in seed: uint(32)) {
   next = seed;
}
