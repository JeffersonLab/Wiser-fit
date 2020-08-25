      program ratios
! get p, pi, or K cross sections from Wiser fit
! valid for 5<E<19 GeV and 5<theta<50 degrees
! from Steve Rock
      implicit none
      integer ie,it
       REAL E,EP,TH,THG,RL,AM,R,SIN2,COS2,X,Q2,ANU,WSQ
       REAL F2p,eF2p,PIE,PIKM,PIKP,SIGPIM,SIGPIP,SIGKAP,SIGKAM
       REAL RATE_E,DEP,RATE_PI
       REAL W1,W2,DSIGX
       REAL*8 DX,DQ2,DF2p,DeF2p
c
c initial conditions
c
      E=10.6      ! beam electron in GeV
      TH=30.0      ! polar angle in degree
      RL=1.0      ! radiation length of the target 1% from Peter 2% from Wiser
      AM=0.9383
      R=0.2
      DEP=0.5
c
      do it=1,31
       TH=4+it ! E'
      do ie=1,20
       EP=1.0+(ie-1)*DEP ! E'
c     define kinematic variables
c
       SIN2=SIN(TH*3.1415928/360.)**2
       COS2=1.0-SIN2
       Q2=4.*E*EP*SIN2
       ANU=E-EP
        X=Q2/2./AM/ANU
       WSQ=AM**2+2.*AM*ANU-Q2
          IF(X.LT.0.9.AND.WSQ.GT.1.5) THEN
       DX=X
       DQ2=Q2
            CALL F2nmc_new(DX,DQ2,'H',DF2p,DeF2p)
            F2p=DF2p
            eF2p=DeF2p
c            CALL F2nmc(1,X,Q2,F2p,eF2p)
c
            W2=F2p/ANU
            W1=W2*(1.+ANU**2/Q2)/(R+1.)
            DSIGX=5.2/E/E/SIN2/SIN2*(COS2*W2+2.*SIN2*W1)
            RATE_E=DSIGX*0.001*2.E9*1.E-33*2.*6.023E23*DEP
            THG=TH
            CALL WISER_ALL_SIG (E,EP,THG,RL,2,SIGPIM)
            CALL WISER_ALL_SIG (E,EP,THG,RL,1,SIGPIP)
            CALL WISER_ALL_SIG (E,EP,THG,RL,4,SIGKAM)
            CALL WISER_ALL_SIG (E,EP,THG,RL,3,SIGKAP)
            PIE=SIGPIM/DSIGX
	    PIKM=SIGKAM/SIGPIM
	    PIKP=SIGKAP/SIGPIP
            RATE_PI=RATE_E*PIE
c            RADCOR=EP/18.*1.5
c            RATE_E=RATE_E/RADCOR
c            RATE=RATE_E+RATE_PI
       WRITE(6,'(A,1X,9F9.4,5E12.3)') 'EP=',E,TH,EP,Q2,X,WSQ
     *,PIE,PIKP,PIKM,DSIGX,SIGPIP,SIGPIM,SIGKAP,SIGKAM
          ELSE
c        WRITE(6,'(1X,4F12.3,A)') E,EP,X,Q2,'--OUT of Range'      
          ENDIF

      enddo
      enddo

      stop
      end

