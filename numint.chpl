// =============================================================================
// ==> numint: rotinas de integração numérica
// =============================================================================
// -----------------------------------------------------------------------------
// --> trapezio(n,a,b,f): integra f entre a e b com n trapézios
// -----------------------------------------------------------------------------
proc trapezio(
   const in n: int,
   const in a: real,
   const in b: real,
   f
) : real {
   const deltax = (b-a)/n;
   const Se = f(a) + f(b);    // define Se
   var Si = 0.0;              // inicializa Si
   for k in 1..n do {         // calcula Si
      var xk = a + k*deltax; 
      Si += f(xk);
   }
   var I = Se + 2*Si;         // cálculo de I
   I *= deltax;                
   I /= 2;                     
   return I;
}   
// -----------------------------------------------------------------------------
// --> trapepsilon(epsilon,a,b,f): calcula a integral de f entre a e b
//   com erro absoluto epsilon, de forma eficiente
// -----------------------------------------------------------------------------
proc trapepsilon(
   const in epsilon: real,
   const in a: real,
   const in b: real,
   f
) : (real,real) {
   var eps = 2*epsilon;            // estabelece um erro inicial grande
   var n = 1;                      // n é o número de trapézios
//   writeln(a," ",f(a));
//   writeln(b," ",f(b));
//   exit(1);
   const Se = f(a) + f(b);         // Se não muda
   var deltax = (b-a)/n;           // primeiro deltax
   var dx2 = deltax/2;             // primeiro deltax/2
   var Siv = 0.0;                  // Si "velho"
   var Iv = Se*dx2;                // I "velho"
   var In: real;                   // I "novo"
   while eps > epsilon do {        // executa o loop pelo menos uma vez   
      var Sin = 0.0;               // Si "novo"
      n *= 2;                      // dobra o número de trapézios
      deltax /= 2;                 // divide deltax por dois
      dx2 = deltax/2;              // idem para dx2
      for i in 1..n by 2 do {      // apenas os ímpares...
         var xi = a + i*deltax;    // pula os ptos já calculados!
         Sin += f(xi);             // soma sobre os novos ptos internos
      }
      Sin = Sin + Siv;             // aproveita todos os ptos já calculados
      In = (Se + 2*Sin)*dx2;       // I "novo"
      eps = abs(In - Iv);          // calcula o erro absoluto
      Siv = Sin;                   // atualiza Siv
      Iv = In;                     // atualiza Iv
   }
   return (In,eps);
}
// -----------------------------------------------------------------------------
// --> untrap: unevenly spaced points, trapezoids
// -----------------------------------------------------------------------------
proc untrap(
   ref ax: [] real,
   ref ay: [] real
   ) : real {
   assert(ax.rank == 1);     // make sure arrays are 1D
   assert(ay.rank == 1);
   const n = ax.size - 1;    // check their size
   assert(ay.size == n+1);
   ref x = ax.reindex(0..n);     // re-index
   ref y = ay.reindex(0..n);
   var I = 0.0;
   for k in 1..n do {
      var h = x[k] - x[k-1];
      var b = (y[k] + y[k-1])/2.0;
      I += b*h;
   }
   return I;
}
