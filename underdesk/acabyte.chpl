var aca = "açafrãolafrário";
var ib = aca.find("ão");
writeln("ib = ",ib);
writeln(ib.type: string);
writeln(aca.this(ib));
writeln(aca.this(5));
proc string.nfind(need: string): int {
   var n = this.size;
   var s = need.size;
   var i = 0;
   var found = (this[i..#s] == need);
   while i+s < n && !found do {
      writeln("i, s, need, substr: ",i," ",i+s," ",need," ",this[i..#s]);
      i += 1;
      found = (this[i..#s] == need);
   }
   if found then return i;
   return -1;
}
writeln(aca);
writeln("io  ",aca.ifind("io"));
writeln("aça ",aca.ifind("aça"));
writeln("rio ",aca.ifind("rio"));
writeln("não ",aca.ifind("não"));
proc string.ifind(need: string): int {
   var n = this.size;
   var s = need.size;
   var i = 0;
   do {
//      writeln("i, s, need, substr: ",i," ",i+s," ",need," ",this[i..#s]);
      if (this[i..#s] == need) then return i;
      i += 1;
   } while i <= n-s ;
   return -1;
}
