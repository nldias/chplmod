// --------------------------------------------------------------------------------------
// CRLEET: Modelo de calculo de evaporacao em lagos e evapotranspiração usando
// a relacao complementar
// 
// Observacao: os numeros que acompanham os comentarios seguem exatamente a documentacao
// do CRLE apresentada em Morton, 1983.  Os simbolos nem sempre sao os mesmos, mas a
// listagem do significado de cada variavel ao lado de sua declaracao torna a "traducao"
// imediata.
// 
// Referencias:
// 
// Morton, F.I. (1983a), Operational estimates of areal evapotranspiration and their
// significance to the science and practice of hydrology, Journal of Hydrology, 66:1--76.
// 
// Morton, F.I. (1983b), Operational estimates of lake evaporation, Journal of Hydrology,
// 66:77--100.
// 
// Morton, F.I., Ricard, F. e Fogarasi, S., Operational Estimates of Areal
// Evapotranspiration and Lake Evaporation -- Program WREVAP, NHRI Paper No 24, National
// Hydrology Research Institute, Inland Waters Directorate, Ottawa, Canada, 1985.
// 
// Nelson Luis Dias
//
// 19 ago 1997: calculo da pressao de saturacao de vapor e algumas constantes modificado
// para reproduzir melhor o programa WREVAP
//
// 04 nov 2004: inclusão do cálculo de evapotranspiração: a variável evopt define os 
// valores de algumas constantes
// --------------------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>   
#include "crleet.h"
// ---------------------------------------------------------------------------------------
// constantes
// ---------------------------------------------------------------------------------------
const double sigma  = 5.67e-8 ;          // cte stefan-boltzmann           
// ---------------------------------------------------------------------------------------
//   variaveis globais
// ---------------------------------------------------------------------------------------
static int
   evopt ;          // 0 == lake, anything else is evapotranspiration
static double
   eps,             // emissividade/absortividade     
   PA,		    // precipitação média anual
   pp0,             // relacao patm local / nivel mar           
   fi;              // latitude da estacao
// ---------------------------------------------------------------------------------------
// --> cr_Lheat: calor latente de evaporação em J/kg/K
//
// References:
// Linsley, Kohler e Paulhus (1975) "Hydrology for Engineers", McGraw-Hill.
// ---------------------------------------------------------------------------------------
inline double  cr_Lheat( 
double T            // temperatura da água em K
)  {
   return (3145780.0 - 2361.355 * T) ;
} 
// ---------------------------------------------------------------------------------------
// --> cr_sat: pressao de saturacao de vapor d'agua, sobre agua ou gelo, formula CRLE
// ---------------------------------------------------------------------------------------
inline double cr_sat(
int gelo,           // se gelo == 1, usa ctes do gelo
double T            // Temperatura, K
) {
// ---------------------------------------------------------------------------------------
// converte Ta de K para C
// ---------------------------------------------------------------------------------------
   T -= 273.15 ;
   if ( gelo ) {
      return 611 * exp( 21.88*T/(T+265.5)) ;
   }
   else {
      return 611 * exp(17.27*T/(T+237.3)) ;   
   }
}   

// ---------------------------------------------------------------------------------------
// --> cr_satd: pressao de saturacao de vapor e sua derivada
// ---------------------------------------------------------------------------------------
inline void cr_satd(
int gelo,           // se gelo == 1, usa ctes do gelo
double T,	    // Temperatura, K
double *e,	    // pressao de saturacao de vapor
double *d	    // e sua derivada
) {
// ---------------------------------------------------------------------------------------
// converte Ta de K para C
// ---------------------------------------------------------------------------------------
   T -= 273.15 ;
// ---------------------------------------------------------------------------------------
// gelo ou água?
// ---------------------------------------------------------------------------------------
   if ( gelo ) {
      (*e) = 611 * exp( 21.88*T/(T+265.5)) ;
      (*d) = 5809.14 * (*e) / ((T+265.5)*(T+265.5)) ;
   }
   else {
      (*e) = 611 * exp( 17.27*T/(T+237.3)) ;
      (*d) = 4098.17 * (*e) / ((T+237.3)*(T+237.3)) ;
   }
   return ;
}   
// ---------------------------------------------------------------------------------------
// --> cr_dew: cálculo da temperatura de ponto de orvalho iterativamente por Newton-Raphson
// ---------------------------------------------------------------------------------------
inline double cr_dew(
double Ta,          // em K
double ea           // em Pa
) {
double
   DeltaT = 0.01,   // temperature increments
   dsd,             // slope at Td           
   esd,             // vapor pressure at Td  
   f,               // e*(T) - ea            
   Td = Ta ;        // dew-point temperature 
   while ( fabs(DeltaT) >= 0.01 ) {
      cr_satd(0,Td,&esd,&dsd) ;
      f = esd - ea ;
      DeltaT = -f/dsd ;
      Td = Td + DeltaT ;
   }
   return Td ;
}
// ---------------------------------------------------------------------------------------
// --> cr_emrad: conversao de graus para radianos
// ---------------------------------------------------------------------------------------
inline double cr_emrad( double xgra) {
   return( xgra * M_PI / 180.0 ) ;
}
// ---------------------------------------------------------------------------------------
// --> cr_emgra: conversao de radianos para graus
// ---------------------------------------------------------------------------------------
inline double cr_emgra( double xrad) {
   return( xrad * 180.0 / M_PI) ;
}
// ---------------------------------------------------------------------------------------
// --> cr_emdec: converte uma leitura em graus, minutos e segundos para graus
//             decimais
//
//    entrada: ggg.mmss  ( em graus (gg), minutos(mm) e segundos(ss) )
//    saida  : ggg.xxxx  ( em graus decimais )
// ---------------------------------------------------------------------------------------
inline double cr_emdec(double xmsc) {
double
   xint,            // parte inteira de um numero
   xdec;            // graus decimais ...

   modf(xmsc,&xint) ;
   xdec = xint ;
   xmsc = ( xmsc - xint ) * 100.0 ;
   modf(xmsc,&xint) ;
   xdec = xdec + xint / 60.0 ;
   xmsc = ( xmsc - xint ) * 100.0 ;
   xdec = xdec + xmsc / 3600.0 ;
   return(xdec) ;

}
// ---------------------------------------------------------------------------------------
// --> cr_ini: calcula o fator pp0
//
// ---------------------------------------------------------------------------------------
void cr_ini ( 
int et,             // which model? 0 = EL and 1 = ET
double H,           // altitude (m)
double Latitude,    // latitude (rad)
double PAin         // precipitação média anual, mm
) {
   fi = Latitude ;
   fprintf(stderr,"entrei em cr_ini\n");
   pp0 = exp( 5.256 * log( (288.0 - 0.0065*H ) / 288.0 ) ) ;
   evopt = et ;
   if (evopt) {
      eps = 0.92 ;
      PA = PAin ;
   }
   else {
      eps = 0.97 ;
   }
}
// ---------------------------------------------------------------------------------------
// --> cr_rad: determina a radiacao solar a partir dos dados sobre mes, latitude e 
// numero de horas de  brilho de sol
// ---------------------------------------------------------------------------------------
void cr_rad ( 
int i,             // mes do ano
double S,           // insolacao media diaria (0<=S<=1)
double Td,          // temperatura de ponto de orvalho (K)
double Ta,          // temperatura do ar (K)
double *alb,        // albedo da superficie
double *Rs,         // radiacao solar (W/m^2)
double *Ra,         // radiacao atmosferica
double *ReA         // radiacao emitida a temp Ta
) {
// --------------------------------------------------------------------------------------
// constantes
// --------------------------------------------------------------------------------------
const double rs0 = 1354.0 ;     // constante solar, W/m2
double
   azd,                // albedo para céu limpo
//   azmax,	       // limite para azd
   ea,                 // pressao de vapor no ar             
   esa,                // pressao de saturacao de vapor      
   delta,              // declinacao do sol                  
   Z,                  // angulo zenital ao meio-dia         
   cosZ,               // cosseno angulo zenital ao meio-dia 
   H,                  // angulo horario                     
   cosH,               // cosseno de H                       
   cosz,               // cosseno de Z                       
   eta,                // raio-vetor sol-terra               
   ar2,                // quadrado inverso dist sol-terra UA 
   Rsea,               // rad solar extra-atmosferica        
   Rseo,               // rad solar com ceu claro            
   azz,                // valor zenital alb ceu claro s/neve 
   az,                 // valor zenital alb ceu claro        
   a0,                 // albedo medio ceu claro             
   c0,c1,c2,           // variaveis auxiliares               
   j,                  // coeficiente de turbidez            
   tau,                // transm.ceus.claros.rad.solar.dir.  
   taua,               // parcela tau absorcao               
   taub,               // parcela tau absorcao+difusao       
   W,                  // agua precipitavel                  
   ro;                 // aumento rad atmosf. pres. nuvens   

// --------------------------------------------------------------------------------------
// 2 - calcula:
// ea, pressao vapor dagua no ar = pressao saturacao aa temp. ponto de orvalho Td
// esa, pressao saturacao aa temp. ar Ta
// note que ea é sempre calculada usando a curva de pressao de saturacao de vapor sobre 
// a agua (0), enquanto que esa depende da temperatura do ar estar acima ou abaixo de zero
// ---------------------------------------------------------------------------------------
//   fprintf(stderr,"2") ; 
   ea = cr_sat(0,Td) ;
   esa = cr_sat((Ta<273.15),Ta) ;
//   printf("Ta, ea, esa = %6.2f %6.2f, %6.2f\n",Ta-273.15,ea,esa) ; */
// ---------------------------------------------------------------------------------------
// 3 - calcula varios angulos e funcoes para obter uma  estimativa da radiacao solar 
// extra-atmosferica Rsea; estima o valor zenital do albedo para céu limpo, azd
// ---------------------------------------------------------------------------------------
   fprintf(stderr,"3") ; 
   delta = cr_emrad( 23.2 * sin( cr_emrad(29.5 * i - 94.0) ) ) ;
   Z = fi - delta ;
   cosZ = cos( Z ) ;
   if (cosZ < 0.001) {
       cosZ = 0.001 ;
   }
   cosH  = 1.0 - cosZ / ( cos(fi) * cos(delta) ) ;
   if (cosH < -1.0) {
      cosH = -1.0 ;
   }
   H = acos(cosH) ;
   cosz  = cosZ + ( sin(H)/H - 1.0 ) * cos(fi) * cos(delta) ;
   fprintf(stderr,"cosz=%f\n",cosz) ; 
   eta   = 1.0 + sin( cr_emrad(29.5 * i - 106.0 ) ) / 60.0 ;
   ar2   = 1.0 / ( eta * eta ) ;
   Rsea  = rs0 * ar2 * ( H/M_PI ) * cosz ;
// --------------------------------------------------------------------------------------
// 4 - o albedo para um ceu claro e fixado no modelo CRLE em 0.05
// Estima az e a0
// A divisao por 100.0 converte de Pa para mb
// --------------------------------------------------------------------------------------
 //   fprintf(stderr,"4") ; */
   if (evopt) {
// ---------------------------------------------------------------------------------------
// albedo para céu limpo, no CRAE 1983
// ---------------------------------------------------------------------------------------
      azd = 0.26 - 0.00012*PA*pow(pp0,0.5)*(1.0 + fabs(cr_emgra(fi)/42) + pow(fabs(cr_emgra(fi)/42.0),2));
      if (azd < 0.11 ) {
	 azd = 0.11 ;
      }
      else if (azd > 0.17 ) {
	 azd = 0.17 ;
      }
      azz = azd ;
// ---------------------------------------------------------------------------------------
// a restrição C-12a de Morton parece funcionar muito MAL (em condições muito úmidas produz 
// albedos negativos), e por este motivo está sendo retirada do programa
// ---------------------------------------------------------------------------------------
//      azmax = 0.5*(0.91 - ea/esa);
//      if (azz > azmax) {
//	 azz = azmax ;
//      }
   } 
   else {
      azz = 0.05 ;
   }
   c0  = ( esa - ea ) / 100.0 ;
   if (c0 < 0.0) {
      c0 = 0.0 ;
   }
   else if (c0 > 1.0) {
      c0 = 1.0 ;
   }
   az  = azz + ( 1.0 - c0*c0 ) * ( 0.34 - azz ) ;
   a0  = az * 
      ( exp(1.08) - exp(Z*2.16/M_PI)*(cosZ*2.16 / M_PI + sin(Z)) ) / (1.473 * (1.0 - sin(Z) )) ;
/*
   printf("azz = %f\n",azz) ;
   printf("az  = %f\n",az) ;
   printf("a0  = %f\n",a0) ;
   getchar() ;             
*/
// ---------------------------------------------------------------------------------------
//   5 - estima a agua precipitavel W e o coeficiente de turbidez j
// ---------------------------------------------------------------------------------------
 //   fprintf(stderr,"5") ; */
   W = ( ea/100.0) / ( 0.49 + ( Ta - 273.15 )/129.0 ) ;
   c1 = 21.0 - ( Ta - 273.15 ) ;
   if (c1 < 0.0) {   
      c1 = 0.0 ;
   }
   else if (c1 > 5.0) {   
      c1 = 5.0 ;
   }
   j = ( 0.5 + 2.5 * ( cosz * cosz ) ) * exp( c1 * ( pp0 - 1.0 ) ) ;
// ---------------------------------------------------------------------------------------
//  6 - calcula a transmissividade de ceus claros aa radiacao solar direta
// ---------------------------------------------------------------------------------------
 //   fprintf(stderr,"6") ; */
/*
   fprintf(stderr,"cosz = %f\n",cosz) ;
   fprintf(stderr,"pp0 = %f\n",pp0) ;
   fprintf(stderr," j  = %f\n",j);
   fprintf(stderr,"W  = %f\n",W) ;
*/   
   tau =     - 0.089 * exp( 0.75 * log( pp0 / cosz ) ) ;
   tau = tau - 0.083 * exp( 0.90 * log( j   / cosz ) ) ;
   tau = tau - 0.029 * exp( 0.60 * log( W   / cosz ) ) ;
   tau = exp(tau) ;
// ---------------------------------------------------------------------------------------
// 7 - estima a parte de tau que e resultado de absorcao
// ---------------------------------------------------------------------------------------
 //   printf("7") ; */
   taua = - 0.0415 * exp( 0.90 * log( j / cosz ) ) ;
   taub = taua ;
   taua = taua - sqrt(0.0029) * exp( 0.30 * log( W / cosz ) ) ;
   taua = exp( taua ) ;
   taub = taub - 0.029 * exp ( 0.60 * log ( W / cosz ) ) ;
   taub = exp( taub ) ;
   if (taua < taub) {   
      taua = taub ;
   }
// ---------------------------------------------------------------------------------------
// 8 - calcula a radiacao solar com ceu claro Rseo e a radiacao solar Rs
// ---------------------------------------------------------------------------------------
 //   printf("8") ; */
   Rseo = Rsea * tau * ( 1.0 + ( 1.0 - tau/taua ) * ( 1.0 + a0*tau ) ) ;
   (*Rs) = S * Rseo + ( 0.08 + 0.30*S) * ( 1.0 - S ) * Rsea ;
// ---------------------------------------------------------------------------------------
// 9 - estima o albedo medio: note que estou convertendo a fórmula com Z em graus de 
//     Morton para usar Z em radianos; por isto há um fator explícito de 180/PI.
// ---------------------------------------------------------------------------------------
//   printf("9") 
   (*alb) = a0 * ( S + ( 1.0 - S) * ( 1.0 - (Z*180.0/M_PI) / 330.0 ) ) ;
   fprintf(stderr,"*** %8.4lf %8.4lf %8.4lf\n",a0,S,Z);
// ---------------------------------------------------------------------------------------
// 10 - estima o aumento proporcional na radiacao atmosferica devido a nuvens
// ---------------------------------------------------------------------------------------
 //   printf("10") 
   c2 = 10.0 * ( ea / esa - S - 0.42 ) ;
   if (c2 < 0.0) {   
      c2 = 0.0 ;
   }
   else if (c2 > 1.0) {  
      c2 = 1.0 ;
   }
   ro = 0.18 * ( (1.0 - c2)*(1.0 - S)*(1.0 - S) +
               c2 * sqrt( 1.0 - S ) ) * pp0 ;
// ---------------------------------------------------------------------------------------
// 11 - Calcula a perda líquida de radiação em comprimentos de onda longos considerando 
// uma temperatura de superfície igual à temperatura do ar.
//   
// Obs: a perda liquida pode ser obtida como:  B = Rea - eps*Ra
// ---------------------------------------------------------------------------------------
//   printf("11") 
// ---------------------------------------------------------------------------------------
// radiacao de corpo negro
// ---------------------------------------------------------------------------------------
   (*ReA) = sigma * exp( 4.0 * log(Ta) ) ;
// ---------------------------------------------------------------------------------------
// radiacao atmosferica 
// ---------------------------------------------------------------------------------------
   (*Ra)  = (*ReA) * ( 0.71 + 0.007 * ( ea / 100.0 ) * pp0 ) * ( 1.0 + ro ) ;
// ---------------------------------------------------------------------------------------
// radiacao emitida pela superficie
// ---------------------------------------------------------------------------------------
   (*ReA) = eps * (*ReA) ;
} 


// ---------------------------------------------------------------------------------------
// --> cr_rad_auto: cálculo do balanço radiativo a partir de dados médios mensais de: 
//     temperatura, pressão de vapor, radiação solar incidente
// ---------------------------------------------------------------------------------------
void cr_rad_auto ( 
int i,             // mes do ano
double Rs,          // insolacao media diaria (0<=S<=1)
double ea,          // pressão de vapor média
double Ta,          // temperatura do ar (K)
double *alb,        // albedo da superficie
double *Ra,         // radiacao atmosferica
double *ReA         // radiacao emitida a temp Ta
) {
// --------------------------------------------------------------------------------------
// constantes
// --------------------------------------------------------------------------------------
const double rs0 = 1354.0 ;     // constante solar, W/m2
double
   azd,                // albedo para céu limpo
//   azmax,	       // limite para azd
   esa,                // pressao de saturacao de vapor      
   delta,              // declinacao do sol                  
   Z,                  // angulo zenital ao meio-dia         
   cosZ,               // cosseno angulo zenital ao meio-dia 
   H,                  // angulo horario                     
   cosH,               // cosseno de H                       
   cosz,               // cosseno de Z                       
   eta,                // raio-vetor sol-terra               
   ar2,                // quadrado inverso dist sol-terra UA 
   Rsea,               // rad solar extra-atmosferica        
   Rseo,               // rad solar com ceu claro            
   azz,                // valor zenital alb ceu claro s/neve 
   az,                 // valor zenital alb ceu claro        
   a0,                 // albedo medio ceu claro             
   c0,c1,c2,           // variaveis auxiliares               
   j,                  // coeficiente de turbidez            
   S,		       // duração do brilho de sol
   tau,                // transm.ceus.claros.rad.solar.dir.  
   taua,               // parcela tau absorcao               
   taub,               // parcela tau absorcao+difusao       
   W,                  // agua precipitavel                  
   ro;                 // aumento rad atmosf. pres. nuvens   

// --------------------------------------------------------------------------------------
// 2 - calcula:
// ea, pressao vapor dagua no ar = pressao saturacao aa temp. ponto de orvalho Td
// esa, pressao saturacao aa temp. ar Ta
// note que ea é sempre calculada usando a curva de pressao de saturacao de vapor sobre 
// a agua (0), enquanto que esa depende da temperatura do ar estar acima ou abaixo de zero
// ---------------------------------------------------------------------------------------
//   fprintf(stderr,"2") ; 
// ---------------------------------------------------------------------------------------
// em cr_rad_auto, a pressão média mensal de vapor é ***fornecida*** ao programa, e 
// portanto não precisa ser calculada a partir da temperatura de ponto de orvalho
// ---------------------------------------------------------------------------------------
   esa = cr_sat((Ta<273.15),Ta) ;
//   printf("Ta, ea, esa = %6.2f %6.2f, %6.2f\n",Ta-273.15,ea,esa) ; */
// ---------------------------------------------------------------------------------------
// 3 - calcula varios angulos e funcoes para obter uma  estimativa da radiacao solar 
// extra-atmosferica Rsea; estima o valor zenital do albedo para céu limpo, azd
// ---------------------------------------------------------------------------------------
   delta = cr_emrad( 23.2 * sin( cr_emrad(29.5 * i - 94.0) ) ) ;
   Z = fi - delta ;
   cosZ = cos( Z ) ;
   if (cosZ < 0.001) {
       cosZ = 0.001 ;
   }
   cosH  = 1.0 - cosZ / ( cos(fi) * cos(delta) ) ;
   if (cosH < -1.0) {
      cosH = -1.0 ;
   }
   H = acos(cosH) ;
   cosz  = cosZ + ( sin(H)/H - 1.0 ) * cos(fi) * cos(delta) ;
   eta   = 1.0 + sin( cr_emrad(29.5 * i - 106.0 ) ) / 60.0 ;
   ar2   = 1.0 / ( eta * eta ) ;
   Rsea  = rs0 * ar2 * ( H/M_PI ) * cosz ;
// --------------------------------------------------------------------------------------
// 4 - o albedo para um ceu claro e fixado no modelo CRLE em 0.05
// Estima az e a0
// A divisao por 100.0 converte de Pa para mb
// --------------------------------------------------------------------------------------
 //   fprintf(stderr,"4") ; */
   if (evopt) {
// ---------------------------------------------------------------------------------------
// albedo para céu limpo, no CRAE 1983
// ---------------------------------------------------------------------------------------
      azd = 0.26 - 0.00012*PA*pow(pp0,0.5)*(1.0 + fabs(cr_emgra(fi)/42) + pow(fabs(cr_emgra(fi)/42.0),2));
      if (azd < 0.11 ) {
	 azd = 0.11 ;
      }
      else if (azd > 0.17 ) {
	 azd = 0.17 ;
      }
      azz = azd ;
// ---------------------------------------------------------------------------------------
// a restrição C-12a de Morton parece funcionar muito MAL (em condições muito úmidas produz 
// albedos negativos), e por este motivo está sendo retirada do programa
// ---------------------------------------------------------------------------------------
//      azmax = 0.5*(0.91 - ea/esa);
//      if (azz > azmax) {
//	 azz = azmax ;
//	 printf("***** azz = %8.4lf\n",azz);
//      }
   } 
   else {
      azz = 0.05 ;
   }
   c0  = ( esa - ea ) / 100.0 ;
   if (c0 < 0.0) {
      c0 = 0.0 ;
   }
   else if (c0 > 1.0) {
      c0 = 1.0 ;
   }
   az  = azz + ( 1.0 - c0*c0 ) * ( 0.34 - azz ) ;
   a0 = az *(exp(1.08) - exp(Z*2.16/M_PI)*(2.16*cosZ/M_PI + sin(Z)))/(1.473*(1.0-sin(Z)));
// ---------------------------------------------------------------------------------------
//   5 - estima a agua precipitavel W e o coeficiente de turbidez j
// ---------------------------------------------------------------------------------------
 //   fprintf(stderr,"5") ; */
   W = ( ea/100.0) / ( 0.49 + ( Ta - 273.15 )/129.0 ) ;
   c1 = 21.0 - ( Ta - 273.15 ) ;
   if (c1 < 0.0) {   
      c1 = 0.0 ;
   }
   else if (c1 > 5.0) {   
      c1 = 5.0 ;
   }
   j = ( 0.5 + 2.5 * ( cosz * cosz ) ) * exp( c1 * ( pp0 - 1.0 ) ) ;
// ---------------------------------------------------------------------------------------
//  6 - calcula a transmissividade de ceus claros aa radiacao solar direta
// ---------------------------------------------------------------------------------------
 //   fprintf(stderr,"6") ; */
/*
   fprintf(stderr,"cosz = %f\n",cosz) ;
   fprintf(stderr,"pp0 = %f\n",pp0) ;
   fprintf(stderr," j  = %f\n",j);
   fprintf(stderr,"W  = %f\n",W) ;
*/   
   tau =     - 0.089 * exp( 0.75 * log( pp0 / cosz ) ) ;
   tau = tau - 0.083 * exp( 0.90 * log( j   / cosz ) ) ;
   tau = tau - 0.029 * exp( 0.60 * log( W   / cosz ) ) ;
   tau = exp(tau) ;
// ---------------------------------------------------------------------------------------
// 7 - estima a parte de tau que e resultado de absorcao
// ---------------------------------------------------------------------------------------
//   printf("7") ; */
   taua = - 0.0415 * exp( 0.90 * log( j / cosz ) ) ;
   taub = taua ;
   taua = taua - sqrt(0.0029) * exp( 0.30 * log( W / cosz ) ) ;
   taua = exp( taua ) ;
   taub = taub - 0.029 * exp ( 0.60 * log ( W / cosz ) ) ;
   taub = exp( taub ) ;
   if (taua < taub) {   
      taua = taub ;
   }
// ---------------------------------------------------------------------------------------
// 8 - calcula a radiacao solar com ceu claro Rseo 
// em cr_rad_auto, a radiação solar é conhecida, e a incógnita é a duração de brilho de
// sol S; implementa C-41 e C-41a
// ---------------------------------------------------------------------------------------
//   printf("8") ; */
   Rseo = Rsea * tau * ( 1.0 + ( 1.0 - tau/taua ) * ( 1.0 + a0*tau ) ) ;
   S = 0.53*Rs/(Rseo - 0.47*Rs);
   if ( S > 1.0 ) S = 1.0 ;
   if (S < 0.0 ) S = 0.0 ;
// ---------------------------------------------------------------------------------------
// o trecho abaixo deixa de fazer sentido
//
//   (*Rs) = S * Rseo + ( 0.08 + 0.30*S) * ( 1.0 - S ) * Rsea ;
// ---------------------------------------------------------------------------------------
// 9 - estima o albedo medio
// ---------------------------------------------------------------------------------------
   (*alb) = a0 * ( S + ( 1.0 - S) * ( 1.0 - (Z*180.0/M_PI) / 330.0 ) ) ;
// ---------------------------------------------------------------------------------------
// 10 - estima o aumento proporcional na radiacao atmosferica devido a nuvens
// ---------------------------------------------------------------------------------------
   c2 = 10.0 * ( ea / esa - S - 0.42 ) ;
   if (c2 < 0.0) {   
      c2 = 0.0 ;
   }
   else if (c2 > 1.0) {  
      c2 = 1.0 ;
   }
   ro = 0.18 * ( (1.0 - c2)*(1.0 - S)*(1.0 - S) +
               c2 * sqrt( 1.0 - S ) ) * pp0 ;
// ---------------------------------------------------------------------------------------
// 11 - Calcula a perda líquida de radiação em comprimentos de onda longos considerando 
// uma temperatura de superfície igual à temperatura do ar.
//   
// Obs: a perda liquida pode ser obtida como:  B = Rea - eps*Ra
// ---------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------
// radiacao de corpo negro
// ---------------------------------------------------------------------------------------
   (*ReA) = sigma * exp( 4.0 * log(Ta) ) ;
// ---------------------------------------------------------------------------------------
// radiacao atmosferica 
// ---------------------------------------------------------------------------------------
   (*Ra)  = (*ReA) * ( 0.71 + 0.007 * ( ea / 100.0 ) * pp0 ) * ( 1.0 + ro ) ;
// ---------------------------------------------------------------------------------------
// radiacao emitida pela superficie
// ---------------------------------------------------------------------------------------
   (*ReA) = eps * (*ReA) ;
} 

// ---------------------------------------------------------------------------------------
// --> cr_evapo: obtencao da evaporacao potencial e evaporacao em lago
//  --------------------------------------------------------------------------------------
void cr_evapo( 
double Td,     // temperatura do ar                      
double Ta,     // temperatura de pto de orvalho          
double alb,    // albedo                                 
double Rs,     // radiacao solar (W/m^2)                 
double Ra,     // radiacao atmosferica (W/m^2)           
double ReA,    // radiacao emitida aa temp. do ar (W/m^2)
double *Tp,    // temperatura de equilibrio (K)          
double *Rla,   // radiacao liq. a temp. do ar (W/m^2)    
double *Rlp,   // radiacao liq. a temp. Tp (W/m^2)       
double *Ep,    // evaporacao potencial (kg/m^2/s)        
double *Ew     // evaporacao lago (kg/m^2/s)             
) {
// ---------------------------------------------------------------------------------------
// variáveis
// ---------------------------------------------------------------------------------------
int
   gelo ;           // para temperaturas abaixo de 0C    
const double
   gamap0 = 66.5;
double
   L,
   fz,
   b0,
   b1,
   b2,
   aux,             // auxiliar, limita B               
   B,               // perda liq rad onda longa         
   dsa,             // inclinacao da pressa sat vapor Ta
   ea,              // pressao de vapor                 
   esa,             // pressao sat vapor temp. Ta       
   uksi,            // fator de estabilidade            
   lambda,          // coeficiente de transf. calor     
   ft,              // coeficiente de transf. vapor     
   gamap,           // cte psicrometrica pressao local  
   esp,             // pressao sat vapor temp. Tp       
   dsp,             // inclinacao da pressao sat vapor  
   deltaTp;         // diferenca entre 2 est. de Tp     
// ---------------------------------------------------------------------------------------
// 12a - valores ligeiramente diferentes dependendo da opção CRAE ou CRLE
// ---------------------------------------------------------------------------------------
   if (evopt) {
      fz = 28.0e-02 ;
      b0 = 1.0 ;
      b1 = 14.0 ;
      b2 = 1.20 ;
   }
   else { 
      fz = 25.0e-02;
      b0 = 1.12 ;
      b1 = 13.0,
      b2 = 1.12;
   }
// ---------------------------------------------------------------------------------------
// 12 - estima a radiacao liquida a temp. do ar, o fator de estabilidade, o coeficiente de 
// transf. de vapor e o coeficiente de transf. de calor
// ---------------------------------------------------------------------------------------
   B = ReA - eps*Ra ;
   aux = 0.03 * eps * sigma * pow(Ta,4.0) ;
   if (B < aux) {  
      B = aux ;
   }
   (*Rla) = ( 1.0 - alb ) * Rs - B ;
   if ((*Rla) < 0.0) {  
      (*Rla) = 0.0 ;
   }
// ---------------------------------------------------------------------------------------
// constantes abaixo de 0C                                                        19970819
// ---------------------------------------------------------------------------------------
   gamap = gamap0 * pp0 ;
   gelo = (Ta < 273.15) ;
   if (gelo) {
      gamap /= 1.15 ;
      fz *= 1.15 ;
   }
// ---------------------------------------------------------------------------------------
// umidade com sat                                                                19970819
// ---------------------------------------------------------------------------------------
   ea = cr_sat(0,Td) ;
   cr_satd(gelo,Ta,&esa,&dsa) ;   
   uksi = 0.28*(1.0 + ea/esa) +
      dsa * (*Rla) / ( ( gamap / sqrt(pp0) ) * 
      b0 * fz * ( esa - ea ) )  ;
   if (uksi > 1.0) {  
      uksi = 1.0 ;
   }
   ft = fz * uksi / sqrt(pp0) ;
   lambda = gamap + 4.0*eps*sigma*(Ta*Ta*Ta) / ft ;
// --------------------------------------------------------------------------------------
// 13 - escolhe valores iniciais para Tp, esp e dsp iguais a Tar, esa e dsa e obtem os 
// valores finais a partir da solucao iterativa de convergencia rapida das equacoes de 
// transferencia de vapor e de balanco de energia
// --------------------------------------------------------------------------------------
   (*Tp) = Ta ;
   esp = esa ;
   dsp = dsa ;
   deltaTp = 0.01 ;
   while (fabs(deltaTp) >= 0.01) {
// --------------------------------------------------------------------------------------
// incremento na temperatura de equilibrio pelo metodo de Newton-Raphson
// --------------------------------------------------------------------------------------
      deltaTp = ( (*Rla)/ft + ea - esp + lambda*(Ta - (*Tp)) ) / ( dsp + lambda );
// --------------------------------------------------------------------------------------
// pressao de vapor e inclinacao a temp. de equilibrio note como, embora eu esteja 
// calculando esp, dsp aa temperatura de equilibrio, o argumento gelo, que define o 
// calculo da saturacao sobre agua ou gelo, ee mantido fixo e definido pelo valor
// de gelo == (Ta<273.15) !
// ---------------------------------------------------------------------------------------
      (*Tp) += deltaTp ;
      cr_satd(gelo,(*Tp),&esp,&dsp) ;
   }
// ---------------------------------------------------------------------------------------
// calcula a evaporacao potencial Ep,
// a radiacao liquida a temp. de equilibrio Rlp
// e a evaporacao em lago Ew
// ---------------------------------------------------------------------------------------
   L = cr_Lheat(*Tp);
   (*Ep) = ( (*Rla) - lambda * ft * ( (*Tp) - Ta ) ) / L ;
   (*Rlp) = L * (*Ep) + gamap * ft * ( (*Tp) - Ta ) ;
   (*Ew) = ( b1 + ( b2 / (1.0 + gamap/dsp) ) * (*Rlp) ) / L ; 
} 
