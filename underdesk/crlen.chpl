// -----------------------------------------------------------------------------
// crlen: CRLE without radiation, and only for lakes
// 
// References:
// 
// Morton, F.I. (1983a), Operational estimates of areal evapotranspiration and
// their significance to the science and practice of hydrology, Journal of
// Hydrology, 66:1--76.
// 
// Morton, F.I. (1983b), Operational estimates of lake evaporation, Journal of
// Hydrology, 66:77--100.
// 
// Morton, F.I., Ricard, F. e Fogarasi, S., Operational Estimates of Areal
// Evapotranspiration and Lake Evaporation -- Program WREVAP, NHRI Paper No 24,
// National Hydrology Research Institute, Inland Waters Directorate, Ottawa,
// Canada, 1985.
// 
// Nelson Luis Dias
//
// 2021-04-30T09:46:25 a new star is born
// -----------------------------------------------------------------------------
const sigma = 5.670374419e-8;  // Stefan-Boltzmann constant, W/m2/K4
const eps = 0.97;              // water emissivity
var 
   fi,              // latitude da estacao
   pp0:             // relacao patm local / nivel mar           
   real;
// -----------------------------------------------------------------------------
// --> cr_ini: calcula o fator pp0
// -----------------------------------------------------------------------------
proc crle_ini ( 
   const in H: real,           // altitude (m)
   const in Latitude: real     // latitude (rad)
) {
   fi = Latitude ;
   writlen("entrei em cr_ini\n");
   pp0 = exp( 5.256 * log( (288.0 - 0.0065*H ) / 288.0 ) ) ;
}
// -----------------------------------------------------------------------------
// --> crle_dew: cálculo da temperatura de ponto de orvalho iterativamente por
// Newton-Raphson
// -----------------------------------------------------------------------------
inline proc crle_dew(
   const in Ta: real,          // em K
   const in ea: real           // em Pa
): real {
   const DeltaT = 0.01;
   var dsd: real;   // slope at Td           
   var esd: real;   // vapor pressure at Td  
   var f: real;               // e*(T) - ea            
   var Td = Ta;        // dew-point temperature 
   while ( abs(DeltaT) >= 0.01 ) {
      (esd,dsd) = svpd(Td,"Tetens");
      f = esd - ea ;
      DeltaT = -f/dsd ;
      Td = Td + DeltaT ;
   }
   return Td ;
}
// -----------------------------------------------------------------------------
// --> cr_evapo: obtencao da evaporacao potencial e evaporacao em lago
//  ----------------------------------------------------------------------------
proc crle_evapo( 
   const in Ta: real,     // temperatura do ar                      
   const in Td: real,     // temperatura de pto de orvalho          
   const in alb: real,    // albedo                                 
   const in Rs: real,     // radiacao solar (W/m^2)                 
   const in Ra: real,     // radiacao atmosferica (W/m^2)           
   const in ReA: real,    // radiacao emitida aa temp. do ar (W/m^2)
   out Tp: real,          // temperatura de equilibrio (K)          
   out Rla: real,         // radiacao liq. a temp. do ar (W/m^2)    
   out Rlp: real,         // radiacao liq. a temp. Tp (W/m^2)       
   out LEp: real,         // evaporacao potencial (W/m^2)
   out LEw: real          // evaporacao lago (W/m^2)
) {
// -----------------------------------------------------------------------------
// variáveis
// -----------------------------------------------------------------------------
   const gamap0 = 66.5;
   const fz = 25.0e-02;
   const b0 = 1.12 ;
   const b1 = 13.0;
   const b2 = 1.12;
   var 
      aux,             // auxiliar, limita B               
      B,               // perda liq rad onda longa         
      dsa,             // inclinacao da pressa sat vapor Ta
      ea,              // pressao de vapor                 
      esa,             // pressao sat vapor temp. Ta       
      uksi,            // fator de estabilidade            
      lambida,         // coeficiente de transf. calor     
      ft,              // coeficiente de transf. vapor     
      gamap,           // cte psicrometrica pressao local  
      esp,             // pressao sat vapor temp. Tp       
      dsp,             // inclinacao da pressao sat vapor  
      deltaTp:         // diferenca entre 2 est. de Tp
      real;
// -----------------------------------------------------------------------------
// 12 - estima a radiacao liquida a temp. do ar, o fator de estabilidade, o
// coeficiente de transf. de vapor e o coeficiente de transf. de calor
// -----------------------------------------------------------------------------
   B = ReA - eps*Ra ;
   aux = 0.03 * eps * sigma * Ta**4; 
   if (B < aux) then {  
      B = aux ;
   }
   (Rla) = ( 1.0 - alb ) * Rs - B ;
   if (Rla < 0.0) then {  
      Rla = 0.0 ;
   }
// -----------------------------------------------------------------------------
// constantes abaixo de 0C                                                      
// -----------------------------------------------------------------------------
   gamap = gamap0 * pp0 ;
   var gelo = (Ta < 273.15) ;
   if (gelo) then {
      gamap /= 1.15 ;
      fz *= 1.15 ;
   }
// -----------------------------------------------------------------------------
// umidade com sat                                                              
// -----------------------------------------------------------------------------
   ea = svp(Td, "Tetens") ;
   (esa,dsa) = svpd(Ta,"Tetens");
   uksi = 0.28*(1.0 + ea/esa) +
      dsa * Rla / ( ( gamap / sqrt(pp0) ) * b0 * fz * ( esa - ea ) )  ;
   if (uksi > 1.0) {  
      uksi = 1.0 ;
   }
   ft = fz * uksi / sqrt(pp0) ;
   lambida = gamap + 4.0*eps*sigma*(Ta*Ta*Ta) / ft ;
// --------------------------------------------------------------------------------------
// 13 - escolhe valores iniciais para Tp, esp e dsp iguais a Tar, esa e dsa e obtém os 
// valores finais a partir da solucao iterativa de convergencia rapida das equacoes de 
// transferencia de vapor e de balanco de energia
// --------------------------------------------------------------------------------------
   Tp = Ta ;
   esp = esa ;
   dsp = dsa ;
   deltaTp = 0.01 ;
   while (abs(deltaTp) >= 0.01) do {
// --------------------------------------------------------------------------------------
// incremento na temperatura de equilibrio pelo metodo de Newton-Raphson
// --------------------------------------------------------------------------------------
      deltaTp = ( Rla/ft + ea - esp + lambida*(Ta - Tp) ) / ( dsp + lambida );
// --------------------------------------------------------------------------------------
// pressao de vapor e inclinacao a temp. de equilibrio note como, embora eu esteja 
// calculando esp, dsp aa temperatura de equilibrio, o argumento gelo, que define o 
// calculo da saturacao sobre agua ou gelo, ee mantido fixo e definido pelo valor
// de gelo == (Ta<273.15) !
// ---------------------------------------------------------------------------------------
      Tp += deltaTp ;
      (esp,dsp) = svpd(Tp,"Tetens");
   }
// ---------------------------------------------------------------------------------------
// calcula a evaporacao potencial Ep,
// a radiacao liquida a temp. de equilibrio Rlp
// e a evaporacao em lago Ew
// ---------------------------------------------------------------------------------------
   LEp = Rla - lambida * ft * ( Tp - Ta );
   Rlp = LEp + gamap * ft * ( Tp - Ta );
   LEw = b1 + ( b2 / (1.0 + gamap/dsp) ) * Rlp;
} 
