// -----------------------------------------------------------------------------
// visicow: because it is easier to remember names than numbers :-) 
// -----------------------------------------------------------------------------
private var domstr: domain(string,parSafe=false);    // initially empty
private var domint: domain(1,int);                   // initially empty
private var namvar: [domint] string ;                // the names
private var varnum: [domstr] int ;                   // the indices
private var linval: [domint] string;                 // line by line values
private var cowFS = "";                              // the cow field separator
// -----------------------------------------------------------------------------
// --> cowfence: define the cow "fence" (separator)
// -----------------------------------------------------------------------------
proc cowfence(in FS: string) {
   cowFS = FS;
}
// -----------------------------------------------------------------------------
// --> cowhead: get the header, split it, and fill the variable names namvar and
// the indices varnum
// -----------------------------------------------------------------------------
proc cowhead(in line: string) {
   if line[0] != "#" then {
      halt("visicow-->cowhead: first character of header line must be '#'");
   }
   if cowFS == "" then {
      var namv2 = (line[1..]).split();
      namv2[0] = namv2[0].strip("# \t\r\n");
      var nf = namv2.size;
      domint = {0..nf-1};
      namvar = namv2;
   }
   else {
      var namv2 = (line[1..]).split(cowFS);
      namv2[0] = namv2[0].strip("# \t\r\n");
      var nf = namv2.size;
      domint = {0..nf-1};
      namvar = namv2;
   }
   for i in domint do {
      namvar[i] = namvar[i].strip();
      domstr += namvar[i];
      varnum[namvar[i]] = i;
   }
}
// -----------------------------------------------------------------------------
// --> cowstable: split a line into values and puts them in the stable
// -----------------------------------------------------------------------------
inline proc cowstable(in s: string) {
   var t = s.strip();
   if cowFS == "" then {
      linval = t.split();
   }
   else {
      linval = t.split(cowFS);
   }
}
// -----------------------------------------------------------------------------
// --> cowfield: return one of the fields converted to type t
// -----------------------------------------------------------------------------
inline proc cowfield(in s: string, type t = real): t {
   if !( t==string || t==real || t== int) then {
      compilerError("viscow-->cowfield type "+t:string+" should not be here\n");
   }
   var sval = linval[varnum[s]];
   if sval != "" then {
      return linval[varnum[s]]:t; // assume it is a valid string
   }
   // --------------------------------------------------------------------------
   // The next best thing to an integer nan.
   // --------------------------------------------------------------------------
   else if t == int then { 
      return -99999999;
   }
   // ---------------------------------------------------------------------------
   // must be real or imag
   // ---------------------------------------------------------------------------
   else if t == real then {
      return nan;
   }
   else {
      halt("viscow-->cowfield type "+t:string+" should not be here at runtime\n");
   }      
}
