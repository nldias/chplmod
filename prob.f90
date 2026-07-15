       REAL FUNCTION PPND16 (P, IFAULT)
C
C      ALGORITHMAS241 APPL. STATIST. (1988) VOL. 37, NO. 3
C
C      PRODUCES THE NORMALDEVIATE Z CORRESPONDING TO A GIVEN LOWER 
C      TAIL AREA OF P; Z IS ACCURATE TO ABOUT 1 PART IN 10**16. 
C
C      THE HASH SUMS BELOW ARE THE SUMS OF THE MANTISSAS OF THE 
C      COEFFICIENTS. THEY ARE INCLUDED FOR USE IN CHECKING 
C      TRANSCRIPTION. 
C
       REAL ZERO, ONE, HALF, SPLIT1, SPLIT2, CONST1, CONST2,
       AO, Al, A2, A3, A4, A5, A6, A7, Bl, B2, B3, B4, B5, B6, B7, 
       CO, Cl, C2, C3, C4, C5, C6, C7, Dl, D2, D3, D4, D5, D6, D7, 
       EO, El, E2, E3, E4, E5, E6, E7, Fl, F2, F3, F4, F5, F6, F7, 
       P, Q, R 
C 
       PARAMETER(ZERO = 0.OEO, ONE = 1.OEO, HALF = ONE/2.OEO, 
       SPLIT1 = 0.425E0, SPLIT2 = 5.OEO, * CONST1 = 0.180625E0, CONST2 = 1.6E0) 
C 
C COEFFICIENTS FOR P CLOSE TO 1/2
C                          
       PARAMETER(AO = 3.38713 28727 96366 6080E0, 
                 Al = 1.33141 66789 17843 7745E2,
                 A2 = 1.97159 09503 06551 4427E3, 
                 A3 = 1.37316 93765 50946 1125E4, 
                 A4 = 4.59219 53931 54987 1457E4, 
                 A5 = 6.72657 70927 00870 0853E4, 
                 A6 = 3.34305 75583 58812 8105E4, 
                 A7 = 2.50908 09287 30122 6727E3, 
                 Bl = 4.23133 30701 60091 1252E1, 
                 B2 = 6.87187 00749 20579 0830E2, 
                 B3 = 5.39419 60214 24751 1077E3, 
                 B4 = 2.12137 94301 58659 5867E4, 
                 B5 = 3.93078 95800 09271 0610E4, 
                 B6 = 2.87290 85735 72194 2674E4, 
                 B7 = 5.22649 52788 52854 5610E3) 
C     HASH SUM AB    55.88319 28806 14901 4439 
C 
C     COEFFICIENTS FOR P NEITHER CLOSE TO 1/2 NOR 0 OR 1
       PARAMETER(CO = 1.42343 71107 49683 57734E0,
                 Cl = 4.63033 78461 56545 29590E0, 
                 C2 = 5.76949 72214 60691 40550E0,
                 C3 = 3.64784 83247 63204 60504E0,
                 C4 = 1.27045 82524 52368 38258E0, 
                 C5 = 2.41780 72517 74506 11770E-1, 
                 C6 = 2.27238 44989 26918 45833E-2, 
                 C7 = 7.74545 01427 83414 07640E-4, 
                 Dl = 2.05319 16266 37758 82187E0, 
                 D2 = 1.67638 48301 83803 84940E0, 
                 D3 = 6.89767 33498 51000 04550E-1, 
                 D4 = 1.48103 97642 74800 74590E-1, 
                 D5 = 1.51986 66563 61645 71966E-2, 
                 D6 = 5.47593 80849 95344 94600E-4, 
                 D7 = 1.05075 00716 44416 84324E-9) 
C HASH SUM CD        49.33206 50330 16102 89036 
C 
C COEFFICIENTS FOR P NEAR 0 OR 1 
C
       PARAMETER(EQ = 6.65790 46435 01103 77720E0,
           El = 5.46378 49111 64114 36990E0, 
           E2 = 1.78482 65399 17291 33580E0, 
           E3 = 2.96560 57182 85048 91230E-1, 
           E4 = 2.65321 89526 57612 30930E-2, 
           E5 = 1.24266 09473 88078 43860E-3, 
           E6 = 2.71155 55687 43487 57815E-5, 
           E7 = 2.01033 43992 92288 13265E-7, 
           Fl = 5.99832 20655 58879 37690E-1, 
           F2 = 1.36929 88092 27358 05310E-1,
           F3 = 1.48753 61290 85061 48525E-2,
           F4 = 7.86869 13114 56132 59100E-4, 
           F5 = 1.84631 83175 10054 68180E-5, 
           F6 = 1.42151 17583 16445 88870E-7, 
           F7 = 2.04426 31033 89939 78564E-15) 
C HASH SUM EF  47.52583 31754 92896 71629 
C
       IFAULT = 0
       Q = P - HALF
       IF (ABS(Q) .LE. SPLIT1) THEN
          R = CONST1 - Q * Q
          PPND16 = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3) 
           * R + A2) * R + Al) * R + AO) / (((((((B7 * R + B6) * R + B5) 
           * R + B4) * R + B3) * R + B2) * R + B1) * R + ONE)
          RETURN
       ELSE IF (Q .LT. 0) THEN
          R = P
       ELSE
          R = ONE - P
       ENDIF
       IF (R .LE. ZERO) THEN
          IFAULT = 1
          PPND16 = ZERO
          RETURN
       ENDIF
       R = SQRT(-LOG(R))
       IF (R .LE. SPLIT2) THEN
          R = R - CONST2
          PPND16 = (((((((C7 * R + C6) * R + C5) * R + C4) * R
           + C3) * R + C2) * R + Cl) * R + C0) / (((((((D7 * R
           + D6) * R + D5) * R + D4) * R + D3) * R + D2) * R 
               + D1) * R + ONE)
       ELSE
          R = R - SPLIT2
          PPND16 = (((((((E7 * R + E6) * R + E5) * R + E4) * R 
           + E3) * R + E2) * R + El) * R + EO) / (((((((F7 * R 
           + F6) * R + F5) * R + F4) * R + F3) * R + F2) * R 
               + Fl) * R + ONE)
       ENDIF       
       IF (Q .LT. 0) PPND16 = -PPND16
       RETURN
       ENDIF
       END
