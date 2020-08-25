c
      Subroutine F2NMC(T,x,qsq,F2,err_f2)
! This subroutine returns a value for F2 from the NMC parametrisation
! in Phy Lett. b295 159-168 Proton and Deuteron F2 Structure Functions
! in deep inelastic muon scattering   P. Amaudruz et. al.
!


      real a(7,2)/-0.1011,2.562,0.4121,-0.518,5.967,-10.197,4.685,
     >            -0.0996,2.489,0.4684,-1.924,8.159,-10.893,4.535/
      real b(4,2)/0.364,-2.764,0.0150,0.0186,
     >            0.252,-2.713,0.0254,0.0299/
      real c(4,2)/-1.179,8.24,-36.36,47.76,
     >            -1.221,7.50,-30.49,40.23/
      real Lam/0.25/
      real Qsqo/20.0/
      real log_term
      real A_x,B_x,C_x,F2,err_f2,x,qsq
      integer T                 ! 1 = Proton
                                ! 2 = Deuteron

      A_x=x**a(1,t)*(1-x)**a(2,t)*(a(3,t)+a(4,t)*(1-x)+a(5,t)*(1-x)**2
     >    +a(6,t)*(1-x)**3+a(7,t)*(1-x)**4)
      B_x=b(1,t)+B(2,t)*x+b(3,t)/(x+b(4,t))
      C_x=c(1,t)*x+c(2,t)*x**2+c(3,t)*x**3+c(4,t)*x**4
      Log_term = alog(qsq/lam**2)/alog(qsqo/lam**2)
      F2=A_x*(log_term)**B_x*(1+C_x/qsq)
      return
      end
c
c
c
      Subroutine WISER_ALL_SIG(E0,P,THETA_DEG,RAD_LEN,TYPE,SIGMA)

!------------------------------------------------------------------------------
! Calculate pi,K,p  cross section for electron beam on a proton target
! IntegrateQs over function WISER_FIT using integration routine QUADMO
! E0         is electron beam energy, OR max of Brem spectra
! P,E       is scattered particle  momentum,energy
! THETA_DEG  is kaon angle in degrees
! RAD_LEN (%)is the radiation length of target, including internal
!                (typically 5% for SLAC, 2% for JLab)
!               = .5 *(target radiation length in %) + 5. (+2 for JLab)
!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA ***
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! SIGMA      is output cross section in nanobars/GeV-str
!------------------------------------------------------------------------------

      IMPLICIT NONE       
      REAL E0,P,THETA_DEG,RAD_LEN,SIGMA
      INTEGER TYPE
      COMMON/WISER_ALL/ E,P_COM,COST,P_T,TYPE_COM,PARTICLE ,M_X,U_MAN
      REAL E,P_COM,COST,P_T,M_X,U_MAN
      INTEGER TYPE_COM,PARTICLE
!  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
      REAL A5(6)/-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53/
      REAL A6(6)/-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45/
      REAL MASS2(3)/.019488, .2437, .8804/
      REAL MASS(3)/.1396, .4973, .9383/ 
      REAL MP/.9383/,  MP2/.8804/, RADDEG/.0174533/
      REAL  M_L,SIG_E
      REAL*8 E_GAMMA_MIN,WISER_ALL_FIT,QUADMO,E08,EPSILON/.003/
      EXTERNAL WISER_ALL_FIT                        
      INTEGER N,CHARGE
                                            
      P_COM = P
      TYPE_COM = TYPE
      PARTICLE = (TYPE+1)/2       ! 1= pi, 2= K, 3 =P
      CHARGE = TYPE -2*PARTICLE +2  ! 1 for + charge, 2 for - charge
      E08 =E0
                
      E =SQRT(MASS2(PARTICLE) + P**2)

      COST = COS(RADDEG * THETA_DEG)
      P_T = P * SIN(RADDEG * THETA_DEG)
      IF(TYPE.LE.4) THEN  !mesons
       IF(CHARGE.EQ.1) THEN   ! K+ n final state
        M_X = MP
       ELSE   ! K- K+ P final state
        M_X = MP+ MASS(PARTICLE)
       ENDIF
      ELSE  ! baryons 
       IF(CHARGE.EQ.1) THEN   ! pi p  final state
        M_X = MASS(1)  ! pion mass
       ELSE   ! P P-bar  P final state
        M_X = 2.*MP
       ENDIF
      ENDIF
      E_GAMMA_MIN = (M_X**2 -MASS2(PARTICLE ) -MP2+2.*MP*E)/
     >  (2.*(MP -E +P*COST))
!      WRITE(10,'(''E_GAMMA_MIN='',F10.2,''  p_t='',F8.2)')
!     >     E_GAMMA_MIN,P_T
!      E_GAMMA_MIN = MP *(E + MASS(PARTILCE))/(MP -P*(1.-COST))
      
      IF(E_GAMMA_MIN.GT..1) THEN !Kinematically allowed?
       M_L = SQRT(P_T**2 + MASS2(PARTICLE))    

       IF(TYPE.NE.5) THEN  ! everything but proton
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) *EXP(A6(TYPE) *P_T**2/E)
       ELSE ! proton

        U_MAN = ABS(MP2 + MASS2(PARTICLE) -2.*MP*E)
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) 
       ENDIF
       SIGMA = P**2/E * 1000. * RAD_LEN/100. *SIG_E 
      ELSE   ! Kinematically forbidden
       SIGMA = 0.
      ENDIF

      RETURN
      END

      REAL*8 FUNCTION WISER_ALL_FIT(E_GAMMA)

!---------------------------------------------------------
! Calculates  pi, k, p  cross section for gamma + p -> k
!  It is already divided by E_GAMMA, the bremstrulung spectra
! David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
! Can be called from WISER_SIG using integration routine QUADMO
! E,P are KAON energy and momentum
! P_t is KAON transverse momentum
! P_CM is KAON center of mass momentum
! P_CM_L is KAON center of mass longitudinal momentum
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! E_GAMMA is photon energy.
!             Steve Rock 2/21/96
!---------------------------------------------------------
                           
      IMPLICIT NONE       
      COMMON/WISER_ALL/ E,P,COST,P_T,TYPE,PARTICLE,M_X,U_MAN

      REAL  E,P,COST,P_T,M_X,U_MAN
      INTEGER  TYPE  !  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
      INTEGER PARTICLE   ! 1= pi, 2= K, 3 =P
!  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
      REAL A1(6)/566.,  486.,   368., 18.2,  1.33E5,  1.63E3 / 
      REAL A2(6)/829.,  115.,   1.91, 307.,  5.69E4, -4.30E3 / 
      REAL A3(6)/1.79,  1.77,   1.91, 0.98,  1.41,    1.79 / 
      REAL A4(6)/2.10,  2.18,   1.15, 1.83,   .72,    2.24 /
      REAL A6/1.90/,A7/-.0117/ !proton only
      REAL MASS2(3)/.019488, .2437, .8804/
      REAL MASS(3)/.1396, .4973, .9383/ 
      REAL MP2/.8804/,MP/.9383/, RADDEG/.0174533/
      REAL X_R, S,B_CM, GAM_CM,  P_CM
      REAL P_CM_MAX, P_CM_L
      REAL*8 E_GAMMA
                                            

!Mandlestam variables                                                
      S = MP2 + 2.* E_GAMMA * MP    

!Go to Center of Mass to get X_R
      B_CM = E_GAMMA/(E_GAMMA+MP)
      GAM_CM = 1./SQRT(1.-B_CM**2)
      P_CM_L = -GAM_CM *B_CM *E + 
     >          GAM_CM * P * COST
      P_CM = SQRT(P_CM_L**2 + P_T**2)  


      P_CM_MAX =SQRT (S +(M_X**2-MASS2(PARTICLE))**2/S 
     >    -2.*(M_X**2 +MASS2(PARTICLE)) )/2.
      X_R =  P_CM/P_CM_MAX   
       IF(X_R.GT.1.) THEN  ! Out of kinematic range
        WISER_ALL_FIT = 0.
       ELSEIF(TYPE.NE.5) THEN  ! not the proton
        WISER_ALL_FIT = (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)/E_GAMMA  
       ELSE ! special formula for proton
        WISER_ALL_FIT = ( (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)          /
     >   (1.+U_MAN)**(A6+A7*S) )/E_GAMMA  
       ENDIF
      
      RETURN
      END
      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,LOWER,UPPER,EPSLON,NLVL)   1717.   
      REAL*8 FUNCT,LOWER,UPPER,EPSLON                                   1718.   
         INTEGER NLVL                                                   1719.   
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,RETRN(50),I                  1720.   
      REAL*8 VALINT(50,2), MX(50), RX(50), FMX(50), FRX(50),            1721.   
     1   FMRX(50), ESTRX(50), EPSX(50)                                  1722.   
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,L,       1723.   
     1   AREA, ABAREA,   M, COEF, ROMBRG,   EPS                         1724.   
         LEVEL = 0                                                      1725.   
         NLVL = 0                                                       1726.   
         ABAREA = 0.0                                                   1727.   
         L = LOWER                                                      1728.   
         R = UPPER                                                      1729.   
         FL = FUNCT(L)                                                  1730.   
         FM = FUNCT(0.5*(L+R))                                          1731.   
         FR = FUNCT(R)                                                  1732.   
         EST = 0.0                                                      1733.   
         EPS = EPSLON                                                   1734.   
  100 LEVEL = LEVEL+1                                                   1735.   
      M = 0.5*(L+R)                                                     1736.   
      COEF = R-L                                                        1737.   
      IF(COEF.NE.0) GO TO 150                                           1738.   
         ROMBRG = EST                                                   1739.   
         GO TO 300                                                      1740.   
  150 FML = FUNCT(0.5*(L+M))                                            1741.   
      FMR = FUNCT(0.5*(M+R))                                            1742.   
      ESTL = (FL+4.0*FML+FM)*COEF                                       1743.   
      ESTR = (FM+4.0*FMR+FR)*COEF                                       1744.   
      ESTINT = ESTL+ESTR                                                1745.   
      AREA=DABS(ESTL)+DABS(ESTR)                                        1746.   
      ABAREA=AREA+ABAREA-DABS(EST)                                      1747.   
      IF(LEVEL.NE.MAXLVL) GO TO 200                                     1748.   
         NLVL = NLVL+1                                                  1749.   
         ROMBRG = ESTINT                                                1750.   
         GO TO 300                                                      1751.   
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                         1752.   
     1         (LEVEL.LT.MINLVL))  GO TO 400                            1753.   
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                             1754.   
  300    LEVEL = LEVEL-1                                                1755.   
         I = RETRN(LEVEL)                                               1756.   
         VALINT(LEVEL, I) = ROMBRG                                      1757.   
         GO TO (500, 600), I                                            1758.   
  400    RETRN(LEVEL) = 1                                               1759.   
         MX(LEVEL) = M                                                  1760.   
         RX(LEVEL) = R                                                  1761.   
         FMX(LEVEL) = FM                                                1762.   
         FMRX(LEVEL) = FMR                                              1763.   
         FRX(LEVEL) = FR                                                1764.   
         ESTRX(LEVEL) = ESTR                                            1765.   
         EPSX(LEVEL) = EPS                                              1766.   
         EPS = EPS/1.4                                                  1767.   
         R = M                                                          1768.   
         FR = FM                                                        1769.   
         FM = FML                                                       1770.   
         EST = ESTL                                                     1771.   
         GO TO 100                                                      1772.   
  500    RETRN(LEVEL) = 2                                               1773.   
         L = MX(LEVEL)                                                  1774.   
         R = RX(LEVEL)                                                  1775.   
         FL = FMX(LEVEL)                                                1776.   
         FM = FMRX(LEVEL)                                               1777.   
         FR = FRX(LEVEL)                                                1778.   
         EST = ESTRX(LEVEL)                                             1779.   
         EPS = EPSX(LEVEL)                                              1780.   
         GO TO 100                                                      1781.   
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                          1782.   
      IF(LEVEL.GT.1) GO TO 300                                          1783.   
      QUADMO = ROMBRG /12.0D0                                           1784.   
      RETURN                                                            1785.   
      END                                                               1786.   



