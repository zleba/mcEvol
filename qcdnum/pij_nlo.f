

C     ========================================
      DOUBLE PRECISION FUNCTION PP1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27

      IMPLICIT NONE
      DOUBLE PRECISION X, DMB_DILOG
      INTEGER NF

      DOUBLE PRECISION CX2, C1PX, C1MX, CLX, CLX2, CL1MX, CL1PX, CPFFX
      DOUBLE PRECISION CPFFMX, CS3X, CS2X, AAA, BBB, CCC, PQQ, PQQB

      DOUBLE PRECISION c11s6, c16s9, c20s3, c2s3, c4s3, c4s9, c5s3
      DOUBLE PRECISION pi, cpi2s3, cpia

      pi     = 3.14159265359
      c11s6  = 11./6.
      c16s9  = 16./9.
      c20s3  = 20./3.
      c2s3   = 2./3.
      c4s3   = 4./3.
      c4s9   = 4./9.
      c5s3   = 5./3.
      cpi2s3 = pi**2/3.
      cpia   = 67./18. - cpi2s3/2.


      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)
 
C--   P_F(x) in FP eq. (4.52)
      AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
     +      - 5.*C1MX
C--   0.5*P_G(x) in FP eq. (4.53)     
      BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
C--   P_N(x) in FP eq. (4.54)   
      CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
C--   P_QQ(x) in FP eq. (4.50) 
      PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
      PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
 
      PP1SFUNC = PQQ + PQQB
 
      RETURN
      END

C     ========================================
      DOUBLE PRECISION FUNCTION PM1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio NPB175(1980)27


      IMPLICIT NONE
      DOUBLE PRECISION X, DMB_DILOG
      INTEGER NF

      DOUBLE PRECISION CX2, C1PX, C1MX, CLX, CLX2, CL1MX, CL1PX, CPFFX
      DOUBLE PRECISION CPFFMX, CS3X, CS2X, AAA, BBB, CCC, PQQ, PQQB

      DOUBLE PRECISION c11s6, c16s9, c20s3, c2s3, c4s3, c4s9, c5s3
      DOUBLE PRECISION pi, cpi2s3, cpia

      pi     = 3.14159265359

      c11s6  = 11./6.
      c16s9  = 16./9.
      c20s3  = 20./3.
      c2s3   = 2./3.
      c4s3   = 4./3.
      c4s9   = 4./9.
      c5s3   = 5./3.
      cpi2s3 = pi**2/3.
      cpia   = 67./18. - cpi2s3/2.



      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

C--   P_F(x) in FP eq. (4.52) 
      AAA = - CPFFX*2.*CLX*CL1MX - (2.*X+3./C1MX)*CLX - .5*C1PX*CLX2
     +      - 5.*C1MX
C--   0.5*P_G(x) in FP eq. (4.53)     
      BBB =   CPFFX*(.5*CLX2+C11S6*CLX+CPIA) + C1PX*CLX + C20S3*C1MX
C--   P_N(x) in FP eq. (4.54)      
      CCC = - CPFFX*C2S3*(C5S3+CLX) - C4S3*C1MX
C--   P_QQ(x) in FP eq. (4.50) 
      PQQ  = C16S9*AAA + 4.*BBB + C2S3*NF*CCC
C--   0.5*P_A(x) in FP eq. (4.55) and P_QQBAR(x) eq. (4.51)      
      PQQB = - C4S9 * ( CPFFMX*CS2X + C1PX*CLX + 2.*C1MX )
 
      PM1SFUNC = PQQ - PQQB
 
      RETURN
      END
      


C     ========================================
      DOUBLE PRECISION FUNCTION GG1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

      IMPLICIT NONE
      DOUBLE PRECISION X, DMB_DILOG
      INTEGER NF

      DOUBLE PRECISION CX2, C1PX, C1MX, CLX, CLX2,CL1MX,CL1PX,CL1MX2
      DOUBLE PRECISION CS3X, CS2X, CPGG, CMPGG, AAA, BBB, CCC

      DOUBLE PRECISION pi, c20s3, c2s3, c4s3, cpi2s3

      pi     = 3.14159265359
      c20s3  = 20./3.
      c2s3   = 2./3.
      c4s3   = 4./3.
      cpi2s3 = pi**2/3.



      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      CPGG   = 1./C1MX + 1./X -2. + X - CX2
      CMPGG  = 1./C1PX - 1./X -2. - X - CX2

      AAA   = -16.+ 8.*X+ C20S3*CX2 + C4S3/X + (-6.-10.*X)*CLX +
     +        (-2.)*C1PX*CLX2
      BBB   = 2.* C1MX +  26./9.*(CX2-1./X) - C4S3*C1PX*CLX -
     +        20./9.*CPGG
      CCC   = 27./2.*C1MX + 67./9.*(CX2-1./X)+(-25./3.+11./3.*x-
     +        44./3.*CX2)*CLX+4.*C1PX*CLX2+(67./9.-4.*CLX*CL1MX +
     +        CLX2-CPI2S3)*CPGG + 2.*CMPGG*CS2X

      GG1SFUNC = C2S3*NF*AAA + 1.5*NF*BBB + 9.* CCC

      RETURN
      END
      


C     ========================================
      DOUBLE PRECISION FUNCTION GF1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|                   


      IMPLICIT NONE
      DOUBLE PRECISION X, DMB_DILOG
      INTEGER NF

      DOUBLE PRECISION CX2, C1PX, C1MX, CLX, CLX2, CL1MX, CL1PX, CL1MX2
      DOUBLE PRECISION CPGFX, CPGFMX, CS3X, CS2X, AAA, DDD


      DOUBLE PRECISION c136s3, c14s9, c182s9, c2s3, c38s3, c40s9, c44s3
      DOUBLE PRECISION pi, cpi2s3, cpie, cpif

      pi     = 3.14159265359
      c136s3 = 136./3.
      c14s9  = 14./9.
      c182s9 = 182./9.
      c2s3   = 2./3.
      c38s3  = 38./3.
      c40s9  = 40./9.
      c44s3  = 44./3.
      cpi2s3 = pi**2/3.
      cpie   = 5. - cpi2s3
      cpif   = cpi2s3 - 218./9.



      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPGFX  = CX2 + C1MX**2
      CPGFMX = CX2 + C1PX**2
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA =   4. - 9.*X + (4.*X-1.)*CLX + (2.*X-1.)*CLX2
     +      + 4.*CL1MX
     +      + (2.*CLX-2.*CLX*CL1MX+CLX2-2.*CL1MX+CL1MX2+CPIE)
     +      * 2. * CPGFX
      DDD =   C182S9 + C14S9*X + C40S9/X + (C136S3*X-C38S3)*CLX
     +      - 4.*CL1MX - (2.+8.*X)*CLX2 + 2.*CS2X*CPGFMX
     +      + (C44S3*CLX-CLX2-2.*CL1MX2+4.*CL1MX+CPIF) * CPGFX

      GF1SFUNC = C2S3*NF*AAA + 1.5*NF*DDD

      RETURN
      END
      


C     ========================================
      DOUBLE PRECISION FUNCTION FG1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

*     Warning: 'Furmanski-Petronzio' notation where
*                          
*              |F'| = |FF GF| |F| 
*              |G'| = |FG GG| |G|
*                  
*     In ../src/qcdwfun.f we will swap to math notation
*
*              |F'| = |FF FG| |F| 
*              |G'| = |GF GG| |G|

c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c      include 'qconst.inc'

      IMPLICIT NONE
      DOUBLE PRECISION X, DMB_DILOG
      INTEGER NF

      DOUBLE PRECISION CX2, C1PX, C1MX, CLX, CLX2, CL1MX, CL1PX
	DOUBLE PRECISION CL1MX2, CPFGX, CPFGMX, CS3X, CS2X  
	DOUBLE PRECISION AAA, BBB, CCC, c16s9, c4s3, cpi2s3,pi


      pi     = 3.14159265359
      c16s9  = 16./9.
      c4s3   = 4./3.
      cpi2s3 = pi**2/3.


      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CL1MX2 = CL1MX**2
      CPFGX  = (1.+C1MX**2) / X
      CPFGMX = - (1.+C1PX**2) / X
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA   = -5./2.- 7./2.*X+(2.+7./2.*X)*CLX+(-1.+0.5*X)*CLX2 
     +        -2.*X*CL1MX+ (-3.*CL1MX-CL1MX2)*CPFGX
      BBB   = 28./9.+65./18.*X+44./9.*CX2+(-12.-5.*X-8./3.*CX2)*CLX+
     +        (4.+X)*CLX2+2.*X*CL1MX+ (-2.*CLX*CL1MX+0.5*CLX2+
     +        11./3.*CL1MX+CL1MX2-0.5*CPI2S3+0.5)*CPFGX+CPFGMX*CS2X
      CCC   = -C4S3*X- (20./9.+C4S3*CL1MX)*CPFGX

      FG1SFUNC = C16S9*AAA+4.*BBB+2./3.*NF*CCC

      RETURN
      END
      



C     ========================================
      DOUBLE PRECISION FUNCTION FF1SFUNC(X,NF)
C     ========================================

C--   Furmanski and Petronzio PLB97(1980)437 eq. (11)

      IMPLICIT NONE
      DOUBLE PRECISION X, DMB_DILOG
      INTEGER NF
      DOUBLE PRECISION CX2, C1PX, C1MX, CLX, CLX2, CL1MX, CL1PX, CPFFX
      DOUBLE PRECISION CPFFMX, CS3X, CS2X, AAA, BBB, CCC

      DOUBLE PRECISION pi, c10s9, c112s9, c11s6, c14s3, c16s3, c16s9
      DOUBLE PRECISION     c2s3, c40s3, c40s9, cpi2s3, cpia  


      pi     = 3.14159265359
	c10s9  = 10./9.
      c112s9 = 112./9.
      c11s6  = 11./6.
      c14s3  = 14./3.
      c16s3  = 16./3.
      c16s9  = 16./9.
      c2s3   = 2./3.
      c40s3  = 40./3.
      c40s9  = 40./9.
      cpi2s3 = pi**2/3.
      cpia   = 67./18. - cpi2s3/2.


      CX2    = X**2
      C1PX   = 1.+X
      C1MX   = 1.-X
      CLX    = LOG(X)
      CLX2   = CLX**2
      CL1MX  = LOG(C1MX)
      CL1PX  = LOG(C1PX)
      CPFFX  = (1.+CX2) / C1MX
      CPFFMX = (1.+CX2) / C1PX
      CS3X   = -DMB_DILOG(-X)
      CS2X   = .5*(CLX2-CPI2S3) + 2.*(CS3X-CLX*CL1PX)

      AAA = - CPFFX*CLX*(1.5+2.*CL1MX) + 2.*CPFFMX*CS2X
     +      - 1. + X + (.5-1.5*X)*CLX - .5*C1PX*CLX2
      BBB =   CPFFX*(C11S6*CLX+.5*CLX2+CPIA) - CPFFMX*CS2X
     +      + C14S3*C1MX
      CCC = - CPFFX*(C10S9+C2S3*CLX) + C40S9/X - 2.*C1PX*CLX2
     +      - C16S3 + C40S3*X + (10.*X+C16S3*CX2+2.)*CLX
     +      - C112S9*CX2

      FF1SFUNC = C16S9*AAA + 4.*BBB + C2S3*NF*CCC

      RETURN
      END
      



C     ======================================
      double precision function dmb_dilog(x)
C     ======================================

C--   Cernlib routine DDILOG C332
C--
C--   Taken from /afs/cern.ch/asis/share/cern/2001/src/...
C--           .../mathlib/gen/c/dilog64.F

      implicit double precision (a-h,o-z)

      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF
      DMB_DILOG=H
      RETURN
      END


!      program AHOJ
!      write(*,*) 'ahoj'

!      endprogram
