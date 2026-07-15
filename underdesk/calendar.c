/* -------------------------------------------------------------------------
   Calendar: Trata objetos de um calendario

   Nelson Luis da Costa Dias
   21-abr-1988  (em PASCAL)
   19-jun-1991  (em C)
   03-jan-2001  (em C)
   ------------------------------------------------------------------------- */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "definitions.h"
#include "tools.h"
#include "calendar.h"

   const int
      NNDias[2][13] = { {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31},
         {0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31} } ;

   const int NDiasAno[2] = {365,366} ;

   const char NomMes[2][13][4] = { 
   { "000","JAN","FEV","MAR","ABR","MAI","JUN",
           "JUL","AGO","SET","OUT","NOV","DEZ" },
   { "000","JAN","FEB","MAR","APR","MAY","JUN",
           "JUL","AUG","SEP","OCT","NOV","DEC" } } ;
   int
      calingua = PORTUGUES ;

/* ---------------------------------------------------------------------
   --> Idioma: modifica a variavel global calingula
   --------------------------------------------------------------------- */
   void Idioma(int ling)
   begin
      if ( ling == 1 ) then 
         calingua = INGLES ;
      end
      else begin
         calingua = PORTUGUES ;
      end ;
   end ;

/* -------------------------------------------------------------------------
   --> Avanca: faz o dia, o mes e o ano avancarem em 1 dia
   ------------------------------------------------------------------------- */
void Avanca(int *ano, int *mes, int *dia) {
int
   bis = Bissexto(*ano) ;
   (*dia)++ ;
   if ( (*dia) > NNDias[bis][(*mes)] ) {
      (*dia) = 1 ;
      (*mes)++ ;
      if ( (*mes) > 12 ) {
	 (*mes) = 1 ;
	 (*ano)++ ;
      }
   }
   return ;
}
/* -------------------------------------------------------------------------
   --> AvancaH: faz o dia, o mes, o ano e a hora avancarem em 1 hora
   ------------------------------------------------------------------------- */
void AvancaH(int *ano, int *mes, int *dia, int *hora) {
int
   bis = Bissexto(*ano) ;
   (*hora)++ ;
   if ( (*hora) > 24 ) {
      (*hora) = 1 ;
      (*dia)++ ;
      if ( (*dia) > NNDias[bis][(*mes)] ) {
	 (*dia) = 1 ;
	 (*mes)++ ;
	 if ( (*mes) > 12 ) {
	    (*mes) = 1 ;
	    (*ano)++ ;
	 }
      }
   }
   return ;
}

/* -------------------------------------------------------------------------
   --> Bissexto: Diz se um ano e' bissexto ou nao

          Entrada:
             Ano      -- Numero do ano de 0 em diante

          Saida:
             Bissexto -- Um valor logico ( TRUE ou FALSE )
   ------------------------------------------------------------------------- */
   int Bissexto( int Ano )

   begin
      return ( ( ( Ano % 4 == 0 ) && ( Ano % 100 != 0 ) ) || ( Ano % 400 == 0 ) ) ;
   end /* Bissexto */ ;

/* -------------------------------------------------------------------------
   --> DataMDi: Dado o fato de o ano ser ou nao bissexto e um dia do ano
                entre 1 e 366, retorna o mes e o dia do mes

          Entrada:
             Bis   -- Se o ano e' ou nao bissexto
             DiAno -- Dia do ano

          Saida:
             Mes   -- Mes do ano
             Dia   -- Dia do mes
   ------------------------------------------------------------------------- */
   void DataMDi( int Bis,
                 int DiAno,
                 int *Mes,
                 int *Dia )

   begin

      *Mes = 1 ;
      while ( DiAno > NNDias[Bis][(*Mes)] ) loop
         DiAno -= NNDias[Bis][(*Mes)] ;
         (*Mes)++ ;
      end ;
      (*Dia) = DiAno ;

   end /* DataMDi */  ;

/* -----------------------------------------------------------------------
   --> DataCor: Dado o fato de o ano ser ou nao bissexto e um par mes e
                dia, retorna um dia corrido entre 1 e 366

          Entrada:
             Bis    -- Se o ano e' ou nao bissexto
             Mes    -- Mes do ano
             Dia    -- Dia do mes

          Saida:
             DiAno  -- Dia do ano
   ------------------------------------------------------------------------ */
   void DataCor( int Bis,
                 int Mes,
                 int Dia,
                 int *DiAno )
   begin
   unsigned int
      i ;

      *DiAno = Dia ;
      i = 1 ;
      while  ( i < Mes ) loop
         (*DiAno) += NNDias[Bis][i]  ;
         i++ ;
      end ;

   end /* DataCor */ ;

/* ----------------------------------------------------------------------
   --> To_DiadoAno: converte uma string no formato "dd-mes-an" em uma
                    estrutura de dados DiadoAno
   ---------------------------------------------------------------------- */
   int To_DiadoAno(DiadoAno *data, const char *dstr)   
   begin
   unsigned int 
      i ;
   char
      nomes[4] ;


      if ( strlen(dstr) != 9 ) then
         printf("To_DiadoAno --> string c/ data deve ter tamanho 9") ;
         return 1 ;
      end ;
      sscanf(dstr,"%hd",&(data->Dia));
      strncpy(nomes,&dstr[3],3) ;
      nomes[3] = '\0' ;
      i = 0 ;
      do loop
         i++ ;
      end
      while ( (i<=12) && (strcmp(nomes,NomMes[calingua][i]) != 0) ) ;
      if ( i == 13 ) then
         return 1 ;
      end 
      data->Mes = i ;
      data->Ano = atoi(&dstr[7]) + 1900 ;
      return 0 ;
   end /* To_DiadoAno */

/* ----------------------------------------------------------------------
   --> To_dstr: converte uma estrutura de dados DiadoAno em uma string
   ---------------------------------------------------------------------- */
   int To_dstrD(char *dstr, DiadoAno *data) 
   begin
   if ( sprintf(dstr,"%2d-%s-%4d",(int) data->Dia,
                             NomMes[calingua][(int) data->Mes],
                             data->Ano ) == 11 ) then
      return 0 ;
   end
   return -1  ;
   end ;
/* ----------------------------------------------------------------------
   --> To_dstrI: converte uma estrutura de dados Instante em uma string
   ---------------------------------------------------------------------- */
   int To_dstrI(char *dstr, Instante *x) 
   begin
   if ( sprintf(dstr,"%02d-%3s-%04d-%02d:%02d:%05.2lf",
                     (int) x->Dia,
                     NomMes[calingua][(int) x->Mes],
                     x->Ano,
                     (int) x->Hora,
                     (int) x->Min,
                     x->Seg ) == 23 ) then
      return 0 ;
   end ;
   return -1 ;
   end 
/* ----------------------------------------------------------------------
   --> To_dstrIn: converte uma estrutura de dados Instante em uma string
   ---------------------------------------------------------------------- */

int To_dstrIn(char *dstr, Instante *x) {
   if ( sprintf(dstr,"%04d-%02d-%02dT%02d:%02d:%05.2lf",
		x->Ano,
		(int) x->Mes,
                (int) x->Dia,
                (int) x->Hora,
                (int) x->Min,
		x->Seg ) == 22 ) {
      return 0 ;
   }
   return -1 ;
}


// ---------------------------------------------------------------------------------------
// --> To_dstrInISO: converte uma estrutura de dados Instante em uma string ISO 8601
//     a string de formatação dos segundos deve ser da forma, por exemplo:
//     %02d  -- gera algo do tipo 2008-02-15T10:24:15
//     %8.5f -- gera algo do tipo 2008-02-15T10:24:15.00039
// ---------------------------------------------------------------------------------------

void To_dstrInISO(char *dstr, Instante *x, int m, int n) {
   if ( n == 0) {
      double ss = floor(x->Seg) ;
      sprintf(dstr,"%04d-%02d-%02dT%02d:%02d:%0*.*f",
	      x->Ano,
	      (int) x->Mes,
	      (int) x->Dia,
	      (int) x->Hora,
	      (int) x->Min,
	      m,n,
	      ss) ;
   }
   else {
      sprintf(dstr,"%04d-%02d-%02dT%02d:%02d:%0*.*f",
	      x->Ano,
	      (int) x->Mes,
	      (int) x->Dia,
	      (int) x->Hora,
	      (int) x->Min,
	      m,n,
	      x->Seg) ;
   }
}


/* ----------------------------------------------------------------------
   --> To_Indstr: converte uma string em uma estrutura de dados Instante
   ---------------------------------------------------------------------- */

int To_Indstr(Instante *x, char *dstr) {
   if ( sscanf(dstr,"%04d%02hd%02hd-%02hd:%02hd:%05lf",
		&(x->Ano),
		&(x->Mes),
                &(x->Dia),
                &(x->Hora),
                &(x->Min),
		&(x->Seg) ) == 20 ) {
      return 0 ;
   }
   return -1 ;
}


// ---------------------------------------------------------------------------------------
// --> To_IndstrISO: converte uma string ISO 8601 em uma estrutura de dados Instante
// ---------------------------------------------------------------------------------------

int To_IndstrISO(Instante *x, char *dstr) {
   return sscanf(dstr,"%04d-%02hd-%02hdT%02hd:%02hd:%lf",
		 &(x->Ano),
		 &(x->Mes),
		 &(x->Dia),
		 &(x->Hora),
		 &(x->Min),
		 &(x->Seg) ) ;
}


/* -----------------------------------------------------------------------
   --> LI_To_Instante: converte um inteiro into do tipo HHMMSS
       em uma estrutura de dados Instante.  Dia, Mes e Ano sao deixados
       em "branco"
   ----------------------------------------------------------------------- */
void LI_To_Instante(Instante *x, int y) 
begin
double
   ft = (double) y,
   gt,
   frac ;

   ft /= 100.0 ;
   gt = floor(ft) ;
   frac = ft - gt ;
   x->Seg = (float) (frac * 100.0) ;

   ft = gt/100.0 ;
   gt = floor(ft) ;
   frac = ft - gt ;
   x->Min = (char) fround (frac * 100.0) ;


   x->Hora = (char) (int) (gt) ;

end /* LI_To_Instante */ ;

/* -------------------------------------------------------------------------
   --> inc_mes: soma n meses
   ------------------------------------------------------------------------- */
   void inc_mes(Instante* x, int n)
   begin
   int mmm ;
   int aaa ;

      mmm = x->Mes + n ;
      aaa = x->Ano ;
      if ( mmm > 0 ) then
         while (mmm > 12 ) loop
            mmm -= 12 ;
            aaa++ ;
         end ;
      end
      else begin
         while ( mmm <= 0 ) loop
            mmm += 12 ;
            aaa-- ;
         end ;
      end ;
      x->Mes = mmm ;
      x->Ano = aaa ;

   end /* inc_mes */ ;

/* -------------------------------------------------------------------------
   --> inc_dia: soma n dias
   ------------------------------------------------------------------------- */
   void inc_dia(Instante* x,int n)
   begin
   int
      bis ;                  /* booleano para bissexto                       */
   int
      mmu,
      ddu ;
   int
      ddd,                   /* dia (corrido) do ano                         */
      aaa ;                  /* ano                                          */

/* -------------------------------------------------------------------------
   pega os valores presentes: ano atual, se bissexto, e a data corrida atual
   ------------------------------------------------------------------------- */
      aaa = x->Ano ;
      bis = Bissexto(x->Ano) ;
      DataCor(bis,x->Mes,x->Dia,&ddu) ;
/* -------------------------------------------------------------------------
   soma o total de dias ; se o total ultrapassar o tamanho DESTE ano:
   subtrai o numero de dias deste ano e segue para o ano seguinte, ate' que
   o algoritmo caia num ano do futuro.
   ------------------------------------------------------------------------- */
      ddd = ddu + n ;
      if ( ddd > 0 ) then
         while ( ddd > NDiasAno[bis] ) loop
            ddd -= NDiasAno[bis] ;
            aaa++ ;
            bis = Bissexto(aaa) ;
         end ;
      end
      else begin
         while ( ddd <= 0 ) loop
            aaa-- ;
            bis = Bissexto(aaa) ;
            ddd += NDiasAno[bis] ;
         end ;
      end ;
/* -------------------------------------------------------------------------
   A atualizacao da data agora e' fa'cil!
   ------------------------------------------------------------------------- */
      ddu = ddd ;
      x->Ano = aaa ;
      DataMDi(bis,ddu,&mmu,&ddu) ;
      x->Mes = mmu ;
      x->Dia = ddu ;
      return ;

   end /* inc_dia */ ;

/* -------------------------------------------------------------------------
   --> inc_hora: soma n horas
   ------------------------------------------------------------------------- */
   void inc_hora(Instante* x,int n)
   begin
   int
      ddd,
      hhh ;

      hhh = x->Hora + n ;
      ddd = 0 ;
      if ( hhh >= 0 ) then
         while ( hhh > 23 ) loop
            hhh -= 24 ;
            ddd++ ;
         end ;
      end
      else begin
         while ( hhh < 0 ) loop
            hhh += 24 ;
            ddd-- ;
         end ;
      end ;
      inc_dia(x,ddd) ;
      x->Hora = hhh ;
      return ;

   end /* inc_hora */ ;


/* -------------------------------------------------------------------------
   --> inc_min: soma n minutos
   ------------------------------------------------------------------------- */
   void inc_min(Instante* x,int n)
   begin
   int
      hhh,
      min;

      min  = x->Min + n ;
      hhh = 0 ;
      if ( min >= 0 ) then
         while ( min > 59 ) loop
            min -= 60 ;
            hhh++ ;
         end ;
      end
      else begin
         while ( min < 0 ) loop
            min += 60 ;
            hhh-- ;
         end ;
      end ;
      inc_hora(x,hhh) ;
      x->Min = min ;
      return ;

   end /* inc_min */ ;

/* -------------------------------------------------------------------------
   --> inc_seg: soma f = n.cccc... segundos
   2008-01-06T13:02:28 old form of inc_sec with 59.0; compiled, but not used anymore
   ------------------------------------------------------------------------- */
   void old_inc_seg(Instante* x,double f)
   begin
   int
      min;
   double
      seg ;

      seg  = x->Seg + f ;
      min = 0 ;
      if ( seg >= 0.0 ) then
         while ( seg > 59.0 ) loop
            seg -= 60.0 ;
            min++ ;
         end ;
      end
      else begin
         while ( seg < 0.0 ) loop
            seg += 60.0 ;
            min-- ;
         end ;
      end ;
      inc_min(x,min) ;
      x->Seg = seg ;
      return ;

   end /* inc_seg */ 


// ----------------------------------------------------------------------------
//  --> inc_seg: soma f = n.cccc... segundos
//
// Nelson Luís Dias
// 20020704
// 20020704
// 2008-01-06T13:03:12 inc_seg now is the same as ninc_seg
// ----------------------------------------------------------------------------
void inc_seg(Instante* x,double f) {
int
   min;
double
   seg ;
// ----------------------------------------------------------------------------
// probing an algo more than 10 years later
// ----------------------------------------------------------------------------
   seg  = x->Seg + f ;
   min = 0 ;
// ----------------------------------------------------------------------------
// seg > 0 means that I am summing seconds
// ----------------------------------------------------------------------------
   if ( seg >= 0.0 ) {
      while ( seg >= 60.0 ) {
	 seg -= 60.0 ;
	 min++ ;
      }
   }
// ----------------------------------------------------------------------------
// seg < 0 means that I am subtracting seconds
// ----------------------------------------------------------------------------
   else {
      while ( seg < 0.0 ) {
	 seg += 60.0 ;
	 min-- ;
      }
   }
   inc_min(x,min) ;
   x->Seg = seg ;
   return ;
}

// ----------------------------------------------------------------------------
//  --> ninc_seg: soma f = n.cccc... segundos
//
// Nelson Luís Dias
// 20020704
// 20020704
// ----------------------------------------------------------------------------
void ninc_seg(Instante* x,double f) {
int
   min;
double
   seg ;
// ----------------------------------------------------------------------------
// probing an algo more than 10 years later
// ----------------------------------------------------------------------------
   seg  = x->Seg + f ;
   min = 0 ;
// ----------------------------------------------------------------------------
// seg > 0 means that I am summing seconds
// ----------------------------------------------------------------------------
   if ( seg >= 0.0 ) {
      while ( seg >= 60.0 ) {
	 seg -= 60.0 ;
	 min++ ;
      }
   }
// ----------------------------------------------------------------------------
// seg < 0 means that I am subtracting seconds
// ----------------------------------------------------------------------------
   else {
      while ( seg < 0.0 ) {
	 seg += 60.0 ;
	 min-- ;
      }
   }
   inc_min(x,min) ;
   x->Seg = seg ;
   return ;
}






/* -------------------------------------------------------------------------
   --> nsegdia: numero de segundos no dia desde meia-noite
   ------------------------------------------------------------------------- */
   double nsegdia(Instante *x)
   begin
      return( x->Hora*3600.0 + x->Min*60.0 + x->Seg ) ;
   end /* nsegdia */ ;
/* ----------------------------------------------------------
   --> packtime: hora, min, seg.cc em hhmmsscc
   ---------------------------------------------------------- */
int packtime(double hora, double minu, double segs) {
   return (int)
      floor(hora*1.0e6 + minu*1.0e4 + segs*1.0e2 + 0.5) ;
} 
/* -----------------------------------------------------------
   --> unpacktime: hhmmsscc em hora, min, seg
   ---------------------------------------------------------- */
void unpacktime(
int itime, 
double *hora, 
double *minu, 
double *segs) {
   *hora = floor(itime/1.0e6) ;
   *minu = floor((itime-(*hora)*1.0e6)/1.0e4) ;
   *segs = (itime-(*hora)*1.0e6-(*minu)*1.0e4)/100.0 ;      
   return ;
}
/* ----------------------------------------------------------
   --> plustime
   ---------------------------------------------------------- */
int plustime(int x, int y) {
double
   hora,hrx,hry,
   minu,mnx,mny,
   segs,sgx,sgy ;
   unpacktime(x,&hrx,&mnx,&sgx) ;
   unpacktime(y,&hry,&mny,&sgy) ;
   segs = sgx + sgy ;
   minu = mnx + mny ;
   hora = hrx + hry ;
   if ( segs > 59.0 ) {
      minu += floor(segs / 60) ;
      segs = fmod(segs,60) ;
   }
   if ( minu > 59.0 ) {
      hora += floor(minu / 60) ;
      minu = fmod(minu,60) ;
   }
   if ( hora > 23 ) {
      hora = fmod(hora,24) ; 
   }
   return packtime(hora,minu,segs) ;
}
      
// ---------------------------------------------------------------------------------------
//  --> InstanteCmp: compara dois instantes, a la strcmp
//   
// Nelson Luís Dias
// 20030705 
// 20030705
// ---------------------------------------------------------------------------------------
int InstanteCmp(
Instante *a,        // um instante
Instante *b,	    // outro instante
int prec	    // prec-ésimos de segundo na comparação
) {
int 
   c;
// -----------------------------------------------------------------
// primeiro compara os anos
// -----------------------------------------------------------------
   c = a->Ano - b->Ano ;
   if (c) {
      return(c);
   } else {
// -----------------------------------------------------------------
// anos iguais: compara os meses
// -----------------------------------------------------------------
      c = a->Mes - b->Mes ;
      if (c) {
	 return(c);
      } else {
	 c = a->Dia - b->Dia ;
	 if (c) {
	    return(c);
	 } else {
	    c = a->Hora - b->Hora ;
	    if (c) {
	       return(c);
	    } else {
	       c = a->Min - b->Min ;
	       if (c) {
		  return(c) ;
	       } 
	       else {
	       int 
	          ma,
		  mb ;
// ------------------------------------------------------------------------------------------------
// takes the precision to prec: it is essential to round
// ------------------------------------------------------------------------------------------------
	          ma = floor(a->Seg*prec+0.5) ;
	          mb = floor(b->Seg*prec+0.5) ;
	          return (ma - mb) ;
	       }
	    }
	 }
      }
   }
}



/* -----------------------------------------------------------------
   --> InstanteCmp: compara dois instantes, a la strcmp
   
   Nelson Luís Dias
   20030705 
   20030705
   ----------------------------------------------------------------- */
int InstanteCmpBug(
Instante *a,        // um instante
Instante *b,	    // outro instante
int prec	    // prec-ésimos de segundo na comparação
) {
int 
   c;
// -----------------------------------------------------------------
// primeiro compara os anos
// -----------------------------------------------------------------
   c = a->Ano - b->Ano ;
   if (c) {
      return(c);
   } else {
      printf("Anos iguais\n");
// -----------------------------------------------------------------
// anos iguais: compara os meses
// -----------------------------------------------------------------
      c = a->Mes - b->Mes ;
      if (c) {
	 return(c);
      } else {
	 printf("meses iguais\n");
	 c = a->Dia - b->Dia ;
	 if (c) {
	    return(c);
	 } else {
	    printf("dias iguais\n");
	    c = a->Hora - b->Hora ;
	    if (c) {
	       return(c);
	    } else {
	       printf("horas iguais\n");
	       c = a->Min - b->Min ;
	       if (c) {
		  return(c) ;
	       } 
	       else {
		  printf("minutos iguais\n");
	       int 
	          ma,
		  mb ;
// ------------------------------------------------------------------------------------------------
// takes the precision to prec
// ------------------------------------------------------------------------------------------------
	       printf("%20.10f\n",a->Seg) ;
	       printf("%20.10f\n",b->Seg) ;
	       ma = floor( a->Seg*prec+0.5)  ; printf("ma = %d\n",ma) ;
	       mb = floor( b->Seg*prec+0.5) ; printf("mb = %d\n",mb) ;
		  return (ma - mb) ;
	       }
	    }
	 }
      }
   }
}

/* end Calendar */



// printf("%04d\n",time[new].Ano);	  
// printf("%02d\n",time[new].Mes);	  
// printf("%02d\n",time[new].Dia);	  
// printf("%02d\n",time[new].Hora);	  
// printf("%02d\n",time[new].Min);	  
// printf("%04.2f\n",time[new].Seg);  
