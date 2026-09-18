#!/home/nldias/miniconda3/bin/python3
# -*- coding: iso-8859-1 -*-
# ------------------------------------------------------------------------------
# --< atmgas.py: propriedades dos gases atmosféricos
#
# Nelson Luís Dias
# 2008-06-13T17:26:02 when it all began
# 2010-10-11T10:44:50 incluindo o cálculo da viscosidade do ar
#                     em função da temperatura
# 2019-08-07T14:34:21 traduzindo tudo para Português
# 2020-04-02T16:51:10 translating things (and renaming) to English
# 2020-04-06T16:02:47 including all gases in COESA, "U.S. Standard Atmosphere,
#                     1976". U.S. Government Printing Office, U.S. Government
#                     Printing Office, 1976, Table 3, plus H2O
# 2020-04-10T08:44:30 automating the calculation of gas constants
# 2020-04-29T09:32:45 English has won; translating
# 2020-04-29T09:46:40 the exact formula for e from y (and p)
# 2020-05-12T18:03:38 a few more translations: sphum and rho_air
# ------------------------------------------------------------------------------
from numba import jit
from math import exp
# ------------------------------------------------------------------------------
# a table of molar masses in kg/mol
# ------------------------------------------------------------------------------
Mmol = {
   'DRY':  28.966413e-3,
   'N2' :  28.013400e-3,
   'O2' :  31.99880e-3,
   'AR' :  39.94800e-3,   
   'H2O':  18.0153e-3,
   'CO2':  44.00950e-3,
   'NE' :  20.1797e-3,
   'HE' :   4.002602e-3,
   'KR' :  83.798e-3,
   'XE' : 131.293e-3,
   'CH4':  16.04250e-3,
   'H2' :   2.01588e-3,
   'O3' :  47.99820e-3,
   'N2O':  44.01280e-3
}
# --------------------------------------------------------------------
# the universal gas constant
# --------------------------------------------------------------------
Ru  = 8.31446261815    # J/K/mol constante universal dos gases
# --------------------------------------------------------------------
# the individual gas constants are handy
# --------------------------------------------------------------------
Rgas = {}
for gas in Mmol:
   Rgas[gas] = Ru/Mmol[gas]
#   print ('%8.4f' % Rgas[gas])
pass
# ------------------------------------------------------------------------------
# a table of useful abbreviations
# ------------------------------------------------------------------------------
Rd  = Rgas['DRY']      # gas constant for dry air
Rv  = Rgas['H2O']      # gas constant for water vapor
Rc  = Rgas['CO2']      # gas constant for CO2

# RO2 = 259.81         # constante de gás para Oxigênio
# RHe = 2078.5         # constante de gás para Hélio
# Cpd = 1005.0         # calor específico a pressão constante para Ar seco
# Cpv = 1850.0         # calor específico a pressão constante para Vapor d'água


# ------------------------------------------------------------------------------
# --> x2rho: converts molar fraction to density
#
# ------------------------------------------------------------------------------
def x2rho(x,p,T,gas) :
   '''
   x2rho: converts x in molar fraction (volumetric concentration) (parts per thousand, 
   million, billion, etc.) to density (kg/m^3).  

   Example: if = 12 ppm of CO2 at an atmospheric pressure of 101325 Pa, and a
   temperatura of 300 K, do:
   rho = x2rho(12e-6,101325,300,'CO2') 

   gas needs to be one of the strings in the Mmol dictionary
   '''
# ------------------------------------------------------------------------------
# obtains the gas constant
# ------------------------------------------------------------------------------
   R = Rgas[gas]
   return (x*p)/(R*T)

# ------------------------------------------------------------------------------
# these quantities are needed for viscosity and diffusivity calculations
# ------------------------------------------------------------------------------
mu0 = 18.18e-6  # dynamical viscosity of air at 20o C (Pa * s)
T0m = 293.15    # 20o C = 293.15 K
C0m = 120.00    # C = 120 K
# --------------------------------------------------------------------
# viscosidade dinâmica do ar
# --------------------------------------------------------------------
def viscair(T):
   ''' 
   returns the viscosity of air mu as a function of temperature T (in
   K). mu in Pa s. From: Montgomery, R. B. Viscosity and thermal
   conductivity of air and diffusivity of water vapor in air J of
   Meteorology, 1947, 4, 193-196, and Wikipedia (somewhere...)
   '''
   return mu0*((T0m + C0m)/(T + C0m))*((T/T0m)**1.5)
# ------------------------------------------------------------------------------
# viscosidade cinemática do ar
# ------------------------------------------------------------------------------
def difmom(p,T,q):
   '''
   retorna a viscosidade cinemática do ar em função da temperatura T 
   T em K, da pressão atmosférica p em Pa e da umidade específica q
   em kg/kg. nu_u em m^2/s
   '''
   rho = rho_ar(p,T,q)
   return viscair(T)/rho
def difcal(p,T,q):
   '''
   retorna a difusividade térmica do ar em função da pressão, temperatura e
   umidade específica. Ver: Montgomery, R. B. Viscosity and thermal conductivity
   of air and diffusivity of water vapor in air J of Meteorology, 1947, 4,
   193-196.
   '''
   return (difmom(p,T,q)/0.711) 
def difvap(p,T,q):
   '''
   retorna a difusividade molecular do vapor d'água no ar em função da pressão,
   temperatura e umidade específica. Ver: Montgomery, R. B. Viscosity and
   thermal conductivity of air and diffusivity of water vapor in air J of
   Meteorology, 1947, 4, 193-196.
   '''
   return (difmom(p,T,q)/0.596)




# ------------------------------------------------------------------------------
# --> dens2pp: converte densidade em concentração molar/volumétrica
# ------------------------------------------------------------------------------
def dens2pp(rhoi,p,T,gas) :
   '''
dens2pp(rhoi): converte rhoi (kg/m^3) em concentração volumétrica (partes por mil, 
milhão, bilhão, etc.).

Exemplo: se rhoi = 16 g/m^3 de H2O a uma pressão atmosférica de 101325 Pa, e a uma 
   temperatura de 300 K faça:
   x = dens2pp(16e-3,101325,300,'H2O') 
   '''
# ------------------------------------------------------------------------------
# calcula a constante do gás
# ------------------------------------------------------------------------------
   R = Ru/Mmol[gas]   
   return (rhoi*R*T)/p

# ------------------------------------------------------------------------------
# --> pressão atmosférica em função da altitude (em uma atmosfera padrão)
# ------------------------------------------------------------------------------
def pressalt(alt) :
   '''
pressalt: pressão atmosférica (em Pa) em função da altitude (em m) em uma 
atmosfera padrão
   '''
   return 101325.0 * ((288 - 0.0065*alt)/288)**5.256

# --------------------------------------------------------------------
# the exact eqn for e from T, y (and p)
# --------------------------------------------------------------------
@jit(nopython=True)
def eexact(T,y,p):
   es = svp(T)
   return (y * es)/(1 + (y-1)*(es/p))
# --------------------------------------------------------------------
# specific humidity as a function of water vapor pressure and
# atmospheric pressure
# --------------------------------------------------------------------
def sphum(e,p) :
   '''
sphum(e,p): specific humidity (in kg/kg) as a function of water vapor pressure
and atmospheric pressure (both in Pa)
   '''
   return e / ( 1.608*(p-e) + e ) 

# ------------------------------------------------------------------------------
# specific heat of humid air
# ------------------------------------------------------------------------------
def spheat(q):
   return 1005.0 + 845.0*(q) 
# ------------------------------------------------------------------------------
# constante de gás do ar úmido
# ------------------------------------------------------------------------------
def R_ar(q):
   return q* Rv + ( 1.0 - q )*Rd
# ------------------------------------------------------------------------------
# humid air density, kg/m^3
# ------------------------------------------------------------------------------
def rho_air(p,T,q):
   return p/(R_ar(q)*T)
   

def tempot(p0,p,T,q):
   R = rar(q)
   cp = calesp(q)
   return T * (p0/p)**(R/cp)

# ------------------------------------------------------------------------------
# --> calor latente de evaporação em função da temperatura termodinâmica
# ------------------------------------------------------------------------------
def latent(T) :
   '''
   latent: Latent heat of evaporation (in J kg^{-1}) as a function of 
   thermodynamic temperature (in K)

   From Dake 1972 "Evaporative Cooling of a body of water", Water Resour Res,
   (8)1087--1091.
   '''
   return (3142689.0 - 2356.01 * T)

# --------------------------------------------------------------------
# saturation vapor pressure
# --------------------------------------------------------------------
#@jit(nopython=True)
def svp(T,eqn='Richards'):
   '''
--> svp(T): saturation vapor pressure of water and its derivative
as a function of the absolute temperature T (K)

input: T 
output (es)

es == saturation vapor pressure (Pa)

References:

Brutsaert, W. (1982) "Evaporation Into the Atmosphere". D. Reidel 
Publishing Company, Dordrecht, Holland.

Richards, J. M. (1971) "Simple expression for the saturation vapour 
pressure of wter in the range -50o to 140o", Journal of Physics D: 
Applied Physics, v. 4, n. 4, L15--L18

Richards, J. M. (1971) "Simple expression for the saturation vapour 
pressure of wter in the range -50o to 140o", Journal of Physics D:
Applied Physics, v. 4, n. 6, 876

Murray, F. M. (1966) "On the computation of saturation vapor pressure", 
Journal of Applied Meteorolgy v. 6, 203---204.
   '''
   if eqn == 'Richards' :
      a1 = 13.3185
      a2 = -1.9760
      a3 = -0.6445
      a4 = -0.1229
      tr = 1.0 - 373.15 / T ;
      aux1 = ((( a4*tr + a3)*tr + a2)*tr + a1)*tr ;
      return (101325.0 * exp(aux1))
   elif eqn == 'Tetens' :
      if T >= 273.15 :
         return 610.78*exp(17.2693882*(T - 273.16)/(T - 35.86))
      else:
         return 610.78*exp(21.8745584*(T - 273.16)/(T - 7.66))
      pass
   else :
      exit("svp(T,eqn): eqn must be: omit == 'Richards', 'Richards', or 'Tetens'")
   pass

def psvd(T) :           
   '''
--> psvd(T): saturation vapor pressure of water and its derivative
as a function of the absolute temperature T (K)

input: T 
output (es,ds)

es == saturation vapor pressure (Pa)
ds == derivative of es with respect to T (Pa/K)

References:

Brutsaert, W. (1982) "Evaporation Into the Atmosphere". D. Reidel 
Publishing Company, Dordrecht, Holland.
   '''
   a1 = 13.3185
   a2 = -1.9760
   a3 = -0.6445
   a4 = -0.1229
   b0 = 13.3185    # b0 = a1
   b1 = -3.9520    # b1 = 2.0 * a2
   b2 = -1.9335    # b2 = 3.0 * a3
   b3 = -0.4916    # b3 = 4.0 * a4
   tr = 1.0 - 373.15 / T ;
   aux1 = ((( a4*tr + a3)*tr + a2)*tr + a1)*tr 
   aux2 = (( b3*tr + b2)*tr + b1)*tr + b0 
   es = 101325.0 * exp(aux1) 
   ds = 373.15 * (es) * aux2 / (T*T) 
   return(es,ds)


