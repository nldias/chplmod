#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-
def trapezio(n,a,b,f):
   '''
   trapezio(n,a,b,f): integra f entre a e b com n trapézios
   '''
   deltax = (b-a)/n
   Se = f(a) + f(b)           # define Se
   Si = 0.0                   # inicializa Si
   for k in range(1,n):       # calcula Si
      xk = a + k*deltax 
      Si += f(xk)
   I = Se + 2*Si              # cálculo de I
   I *= deltax                
   I /= 2                     
   return I
def trapepsilonlento(epsilon,a,b,f):
   '''
   trapepsilonlento(epsilon,a,b,f): calcula a integral de f entre a e
   b com erro absoluto epsilon, de forma ineficiente
   '''
   eps = 2*epsilon            # estabelece um erro inicial grande
   n = 1                      # um único trapézio
   Iv = trapezio(1,a,b,f)     # primeira estimativa, "velha"
   while eps > epsilon:       # loop
      n *= 2                  # dobra o número de trapézios
      In = trapezio(n,a,b,f)  # estimativa "nova", recalculada do zero
      eps = abs(In - Iv)      # calcula o erro absoluto
      Iv = In                 # atualiza a estimativa "velha"
   return (In,eps)
def trapepsilon(epsilon,a,b,f):
   '''
   trapepsilon(epsilon,a,b,f): calcula a integral de f entre a e b
   com erro absoluto epsilon, de forma eficiente
   '''
   eps = 2*epsilon                 # estabelece um erro inicial grande
   n = 1                           # n é o número de trapézios
   Se = f(a) + f(b)                # Se não muda
   deltax = (b-a)/n                # primeiro deltax
   dx2 = deltax/2                  # primeiro deltax/2
   Siv = 0.0                       # Si "velho"
   Iv = Se*dx2                     # I "velho"
   while eps > epsilon:            # executa o loop pelo menos uma vez   
      Sin = 0.0                    # Si "novo"
      n *= 2                       # dobra o número de trapézios
      deltax /= 2                  # divide deltax por dois
      dx2 = deltax/2               # idem para dx2
      for i in range(1,n,2):       # apenas os ímpares...
         xi = a + i*deltax         # pula os ptos já calculados!
         Sin += f(xi)              # soma sobre os novos ptos internos
      Sin = Sin + Siv              # aproveita todos os ptos já
                                   # calculados
      In = (Se + 2*Sin)*dx2        # I "novo"
      eps = abs(In - Iv)           # calcula o erro absoluto
      Siv = Sin                    # atualiza Siv
      Iv = In                      # atualiza Iv
   return (In,eps)
