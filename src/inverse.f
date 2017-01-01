C ==============================================================
C The following source code is obtained from 
C R package "limSolve" ver 1.5.5.1 written by 
C Karline Soetaert [aut, cre], Karel Van den Meersche [aut], Dick van Oevelen [aut], LAPACK authors [cph].             
C ==============================================================

C ***************************
C
C **************************
C Karline: removed the write statement; where they also passed an integer value, 
C this no long is the case. Each removed statement is preceded by:
C KARLINE: REMOVED WRITE		 


C*********************************************************************
C LEAST DISTANCE SUBROUTINE
C*********************************************************************

      SUBROUTINE ldp(G,H,NUnknowns,NConstraints,NW,X,XNorm,W,xIndex,           &
     &               Mode,verbose, IsError)


      INTEGER           :: NUnknowns,NConstraints,NW
      DOUBLE PRECISION  :: G(NConstraints,NUnknowns) 
      DOUBLE PRECISION  :: H(NConstraints)

      DOUBLE PRECISION  :: X(NUnknowns)
      DOUBLE PRECISION  :: XNorm

      LOGICAL           :: IsError, verbose

      DOUBLE PRECISION  :: W(NW)
      INTEGER           :: xIndex(NConstraints)
      INTEGER           :: Mode
      INTEGER :: xLDPSucces,xLDPNoUnknownsOrEquations,                         &
     &       xLDPToomanyIterations,xLDPIncompatibleConstraints,                &
     &       xLDPUnsolvable


      PARAMETER (xLDPSucces = 1,xLDPNoUnknownsOrEquations     = 2,             &
     &       xLDPToomanyIterations= 3,xLDPIncompatibleConstraints   = 4,       &
     &       xLDPUnsolvable  = -1)


      CALL xLDP(G,NConstraints,NConstraints,NUnknowns,H,X,Xnorm,W,             &
     &          xINdex,Mode)

      IsError=.TRUE.
      if (mode == xLDPSucces) IsError = .FALSE.

      IF (verbose) THEN 
      SELECT CASE (mode)
      CASE (xLDPNoUnknownsOrEquations)
       CALL xMESSAGE ("No unknowns or equations")
      CASE (xLDPToomanyIterations)
       CALL xMESSAGE ("Too many iterations")
      CASE (xLDPIncompatibleConstraints)
       CALL xMESSAGE ("Incompatible constraints ")
      CASE (xLDPUnsolvable       )
       CALL xMESSAGE ("LDP problem unsolvable")
      END SELECT
      ENDIF

      RETURN
      END SUBROUTINE LDP



C****************************************************************************
C SOLVING LINEAR LEAST SQUARES WITH LINEAR INEQUALITIES and EQUALITIES
C****************************************************************************
  
C----------------------------------------------------------------------------
C  Solves the least squares model with equality and inequality constraints
C
C  Mathematically: minimise SUM(A*x-B)^2 
C  where:          E*x=F
C                  G*x>=h   x> 0
C  x(NumUnknowns), A(NumApproximate,NumUnknowns), B(NumApproximate)
C               E(NumEquations,NumUnknowns)     , F(NumEquations)
C               G(NumInequalities,NumUnknowns)  , H(NumInequalities)
C----------------------------------------------------------------------------

      SUBROUTINE lsei (NUnknowns,NEquations,NConstraints,NApproximate,         &
     &          A,B,E,F,G,H,X,mIP,mdW,mWS,IP,W,WS,lpr,ProgOpt,                 &
     &          verbose,IsError)
 
      IMPLICIT NONE

C The arrays and dimensions for the mass balance inversion

      INTEGER          :: NUnknowns,NEquations,NConstraints,NApproximate
      INTEGER          :: mIP,mdW,mWS,lpr
      LOGICAL          :: IsError, verbose
      DOUBLE PRECISION  ::  A (NApproximate,NUnknowns),                        &
     &                      B (NApproximate)          ,                        &
     &                      E (NEquations,NUnknowns)  ,                        &
     &                      F (NEquations)            ,                        &
     &                      G (NConstraints,NUnknowns),                        &
     &                      H (NConstraints)          ,                        &
     &                      X (NUnknowns)                     
C work arrays
      DOUBLE PRECISION :: W(MDW,NUnknowns+1),WS(mWS)
      INTEGER          :: IP(mip)

      INTEGER          :: I,J,K,MOde,ME, MA, MG,N
      DOUBLE PRECISION :: RNORME, RNORML,ProgOpt(lpr)

C---------------------------------------------------------------------
C a linear Least Squares solver with linear Equality and Inequality constraints: (LSEI)
C
C Find MIN ||A*X=B||
C
C and where
C    E*x=F, G*X>H
C    x> 0
C---------------------------------------------------------------------

      N  = NUnknowns
      ME = NEquations
      MA = NApproximate
      MG = NConstraints

C W is Working array with E,A,G   ,F,B,H 

      RNormE = 0.D0
      RNORML = 0.D0
      MODE   = 0

       DO I = 1, ME
         DO J = 1, N
           W(I,J) = E(I,J)
         ENDDO
         W(I,N+1) = F(I)
       ENDDO        

      K = ME
       DO I = 1, MA
         DO J = 1, N
           W(I+K,J) = A(I,J)
         ENDDO
         W(I+K,N+1) = B(I)
       ENDDO        
      K = ME+MA
      DO I = 1,NConstraints
         DO J = 1, N
           W(I+K,J) = G(I,J)
         ENDDO
         W(I+K,N+1) = H(I)
       ENDDO        
      K = K + NConstraints

c      ProgOpt(1) = 1.D0


    
C CALLING SOLVER!

        CALL xdLSEI(W,                                                   &     
     &              MDW,                                                 &     
     &              ME,                                                  & 
     &              MA,                                                  &
     &              MG,                                                  &
     &              N,                                                   &
     &              ProgOpt,                                             &     
     &              X,                                                   &     
     &              RNORME,                                              &     
     &              RNORML,                                              &     
     &              MODE,                                                &     
     &              WS,                                                  &
     &              IP)


      IF (verbose) THEN
       SELECT CASE (Mode) 
       CASE(1)
           CALL XMESSAGE ("LSEI error: equalities contradictory")

       CASE(2)
           CALL XMESSAGE ("LSEI error: inequalities contradictory")

       CASE(3)
           CALL XMESSAGE                                                  &
     &    ("LSEI error: equalities + inequalities contradictory")

       CASE(4)
           CALL XMESSAGE("LSEI error: wrong input")       
       END SELECT
      ENDIF
      IsError = .FALSE.
      IF (mode > 0) Iserror = .TRUE.
      RETURN

      END SUBROUTINE LSEI


C                c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<c
C                c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<c
C                c            ERROR HANDLING          c
C                c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c
C                c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>c


      SUBROUTINE XMESSAGE (String)

      CHARACTER (LEN=*)   :: String

C Check whether it is safe to write
        CALL rwarn(String)

      END SUBROUTINE XMESSAGE



C######################################################################
C 
C Routines for solving quadratic and linear programming problem
C
C SUBROUTINES FROM LINPACK LIBRARY
C
C MAIN ROUTINES ARE xLDP
C xLSEI and  xdbocls
C 
C######################################################################



C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C    Quadratic programming solver    C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C


C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 15, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C  made compatible with fortran 95 by karline Soetaert
C  added x-prefix
C  captured writing to screen -> XMESSAGE

C*************************************************************************C
C LEAST DISTANCE SUBROUTINE
C*************************************************************************C


      SUBROUTINE xLDP (G,MDG,M,N,H,X,XNORM,W,xINDEX,MODE)     
C
C  Algorithm LDP: LEAST DISTANCE PROGRAMMING
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1974 MAR 1, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------

      IMPLICIT NONE
      INTEGER  :: xLDPSucces ,xLDPNoUnknownsOrEquations 
      INTEGER  :: xLDPToomanyIterations,xLDPIncompatibleConstraints
      INTEGER  :: xLDPUnsolvable 
      PARAMETER(xLDPSucces= 1,xLDPNoUnknownsOrEquations     = 2,                 &
     &                      xLDPToomanyIterations         = 3,                   &
     &                      xLDPIncompatibleConstraints   = 4,                   &
     &                      xLDPUnsolvable                = -1)
C Number of unknowns
      INTEGER  :: M, MDG,N         
C Succes or failure
      INTEGER  :: MODE             
      INTEGER  :: xINDEX(*)  
C Constraints G*X>H ; W=workarray; xnorm=residual norm
      DOUBLE PRECISION :: G(MDG,*), H(*), X(*)
      DOUBLE PRECISION :: W(*)                
      DOUBLE PRECISION :: XNORM               

      INTEGER          :: I, IW, IWDUAL, IY, IZ, J, JF, NP1
      DOUBLE PRECISION :: xDIFF, FAC, RNORM
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------

      MODE=xLDPSucces

      IF (N.LE.0) THEN
        MODE=xLDPNoUnknownsOrEquations
        RETURN
      ENDIF

      DO J=1,N   
        X(J)=ZERO     
      ENDDO
      XNORM=ZERO

      IF (M.LE.0) THEN
        MODE=xLDPNoUnknownsOrEquations           
        RETURN
      ENDIF
C   
C     THE DECLARED DIMENSION OF W() MUST BE AT LEAST (N+1)*(M+2)+2*M.   
C   
C      FIRST (N+1)*M LOCS OF W()   =  MATRIX E FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR F FOR PROBLEM NNLS.
C       NEXT     N+1 LOCS OF W()   =  VECTOR Z FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR Y FOR PROBLEM NNLS.
C       NEXT       M LOCS OF W()   =  VECTOR WDUAL FOR PROBLEM NNLS.    
C     COPY G**T INTO FIRST N ROWS AND M COLUMNS OF E.   
C     COPY H**T INTO ROW N+1 OF E.  
C   
      IW=0  
      DO J=1,M   
          DO I=1,N   
               IW=IW+1   
               W(IW)=G(J,I)  
          ENDDO
          IW=IW+1   
          W(IW)=H(J) 
      ENDDO   
      JF=IW+1   
C                                STORE N ZEROS FOLLOWED BY A ONE INTO F.
      DO I=1,N   
          IW=IW+1   
          W(IW)=ZERO    
      ENDDO
      W(IW+1)=ONE   
C   
      NP1=N+1   
      IZ=IW+2   
      IY=IZ+NP1 
      IWDUAL=IY+M   
C   
      CALL xNNLS (W,NP1,NP1,M,W(JF),W(IY),RNORM,W(IWDUAL),W(IZ),         &
     &  xINDEX,MODE)  
C                      USE THE FOLLOWING RETURN IF UNSUCCESSFUL IN NNLS.
      IF (MODE.NE.xLDPSucces) RETURN 

C      IF (RNORM) 130,130,50     KARLINE changed
C  50  
      IF (RNORM .LE. 0) THEN
C karline:remove ????
        MODE=xLDPUnsolvable
        RETURN
      ENDIF

      FAC=ONE   
      IW=IY-1   
      DO I=1,M   
         IW=IW+1   
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
         FAC=FAC-H(I)*W(IW)
      ENDDO
C   
C      IF (DIFF(ONE+FAC,ONE)) 130,130,70 
C   70 FAC=ONE/FAC   

      IF (xDIFF(ONE+FAC,ONE) .LE. 0) THEN
        MODE=xLDPIncompatibleConstraints
        RETURN
      ENDIF
      FAC=ONE/FAC   

      DO J=1,N   
         IW=IY-1   
         DO I=1,M   
            IW=IW+1   
C                               HERE WE ARE USING THE SOLUTION VECTOR Y.
            X(J)=X(J)+G(I,J)*W(IW)
         ENDDO 
         X(J)=X(J)*FAC 
      ENDDO

      DO J=1,N  
         XNORM=XNORM+X(J)**2   
      ENDDO

      XNORM=sqrt(XNORM) 
      RETURN
      END SUBROUTINE xLDP



C*************************************************************************C

C     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE)
C   
C  Algorithm NNLS: NONNEGATIVE LEAST SQUARES
C   
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 15, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
C     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
C     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM   
C   
C                      A * X = B  SUBJECT TO X .GE. 0   
C     ------------------------------------------------------------------
C                     Subroutine Arguments
C
C     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
C                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
C                     MATRIX, A.           ON EXIT A() CONTAINS 
C                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
C                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
C                     THIS SUBROUTINE.  
C     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
C             TAINS Q*B.
C     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
C             CONTAIN THE SOLUTION VECTOR.  
C     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
C             RESIDUAL VECTOR.  
C     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
C             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
C             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
C     ZZ()     AN M-ARRAY OF WORKING SPACE.     
C     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
C                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
C                 P AND Z AS FOLLOWS..  
C   
C                 INDEX(1)   THRU INDEX(NSETP) = SET P.     
C                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
C                 IZ1 = NSETP + 1 = NPP1
C                 IZ2 = N   
C     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
C             MEANINGS. 
C             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
C             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
C                   EITHER M .LE. 0 OR N .LE. 0.
C             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. 
C   
C     ------------------------------------------------------------------
      SUBROUTINE xnnls (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) 
C     ------------------------------------------------------------------
      integer I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX, J, JJ, JZ, L
      integer M, MDA, MODE,N, NPP1, NSETP, RTNKEY
C     integer INDEX(N)  
C     double precision A(MDA,N), B(M), W(N), X(N), ZZ(M)   
      integer INDEX(*)  
      double precision A(MDA,*), B(*), W(*), X(*), ZZ(*)
C Karline: changed dummy -> dummy(1) to avoid warining line 495         
      double precision ALPHA, ASAVE, CC, xDIFF, DUMMY(1), FACTOR, RNORM
      double precision SM, SS, T, TEMP, TWO, UNORM, UP, WMAX
      double precision ZERO, ZTEST
      parameter(FACTOR = 0.01d0)
      parameter(TWO = 2.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      MODE=1
      IF (M .le. 0 .or. N .le. 0) then
         MODE=2
         RETURN
      endif
      ITER=0
      ITMAX=3*N 
C   
C                    INITIALIZE THE ARRAYS INDEX() AND X(). 
C   
      DO I=1,N   
         X(I)=ZERO     
         INDEX(I)=I    
      ENDDO
C   
      IZ2=N 
      IZ1=1 
      NSETP=0   
      NPP1=1
C                             &*****  MAIN LOOP BEGINS HERE  ******     
   30 CONTINUE  
C                  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
C                        OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.    
C   
      IF (IZ1 .GT.IZ2.OR.NSETP.GE.M) GO TO 350   
C   
C         COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
C   
      DO IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         SM=ZERO   
         DO L=NPP1,M
           SM=SM+A(L,J)*B(L)   
         ENDDO  
         W(J)=SM   
      ENDDO
C                                   FIND LARGEST POSITIVE W(J). 
   60 continue
      WMAX=ZERO 
      DO IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         IF (W(J) .gt. WMAX) then
            WMAX=W(J)     
            IZMAX=IZ  
         endif
      ENDDO
C   
C             IF WMAX .LE. 0. GO TO TERMINATION.
C             THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
C   
      IF (WMAX .le. ZERO) go to 350
      IZ=IZMAX  
      J=INDEX(IZ)   
C   
C     THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.    
C     BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID  
C     NEAR LINEAR DEPENDENCE.   
C   
      ASAVE=A(NPP1,J)   
      CALL xH12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)    
      UNORM=ZERO
      IF (NSETP .ne. 0) then
          DO L=1,NSETP   
            UNORM=UNORM+A(L,J)**2     
          ENDDO
      endif
      UNORM=sqrt(UNORM) 
      IF (xDIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM) .gt. ZERO) then
C   
C        COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
C        AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).    
C   
        DO L=1,M  
          ZZ(L)=B(L)    
        ENDDO
        CALL xH12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)   
        ZTEST=ZZ(NPP1)/A(NPP1,J)  
C   
C                                     SEE IF ZTEST IS POSITIVE  
C   
         IF (ZTEST .gt. ZERO) go to 140
      endif
C   
C     REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.  
C     RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
C     COEFFS AGAIN.     
C   
      A(NPP1,J)=ASAVE   
      W(J)=ZERO 
      GO TO 60  
C   
C     THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
C     SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER  
C     TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN   
C     COL J,  SET W(J)=0.   
C   
  140 continue
      DO L=1,M  
        B(L)=ZZ(L)    
      ENDDO
C   
      INDEX(IZ)=INDEX(IZ1)  
      INDEX(IZ1)=J  
      IZ1=IZ1+1 
      NSETP=NPP1
      NPP1=NPP1+1   
C   
      IF (IZ1 .le. IZ2) then
         DO JZ=IZ1,IZ2 
            JJ=INDEX(JZ)  
            CALL xH12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
         ENDDO
      endif
C   
      IF (NSETP .ne. M) then
         DO L=NPP1,M   
           A(L,J)=ZERO   
         ENDDO
      endif
C   
      W(J)=ZERO 
C                                SOLVE THE TRIANGULAR SYSTEM.   
C                                STORE THE SOLUTION TEMPORARILY IN ZZ().
      RTNKEY = 1
      GO TO 400 
  200 CONTINUE  
C   
C                       &*****  SECONDARY LOOP BEGINS HERE ******   
C   
C                          ITERATION COUNTER.   
C 
  210 continue  
      ITER=ITER+1   
      IF (ITER .gt. ITMAX) then
         MODE=3
C         write (*,'(/a)') ' NNLS quitting on iteration count.'
      CALL XMESSAGE ('error in LDP - NNLS quitting on iteration count.')
         GO TO 350 
      endif
C   
C                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.    
C                                  IF NOT COMPUTE ALPHA.    
C   
      ALPHA=TWO 
      DO IP=1,NSETP 
         L=INDEX(IP)   
         IF (ZZ(IP) .le. ZERO) then
            T=-X(L)/(ZZ(IP)-X(L))     
            IF (ALPHA .gt. T) then
               ALPHA=T   
               JJ=IP 
            endif
         endif
      ENDDO 
C   
C          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL   
C          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.   
C   
      IF (ALPHA.EQ.TWO) GO TO 330   
C   
C          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO   
C          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.    
C   
      DO IP=1,NSETP 
         L=INDEX(IP)   
         X(L)=X(L)+ALPHA*(ZZ(IP)-X(L)) 
      ENDDO
C   
C        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I  
C        FROM SET P TO SET Z.   
C   
      I=INDEX(JJ)   
  260 continue
      X(I)=ZERO 
C   
      IF (JJ .ne. NSETP) then
         JJ=JJ+1   
         DO 280 J=JJ,NSETP 
            II=INDEX(J)   
            INDEX(J-1)=II 
            CALL xG1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))   
            A(J,II)=ZERO  
            DO L=1,N  
               IF (L.NE.II) then
C
C                 Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))  
C
                  TEMP = A(J-1,L)
                  A(J-1,L) = CC*TEMP + SS*A(J,L)
                  A(J,L)   =-SS*TEMP + CC*A(J,L)
               endif
            ENDDO   
C
C                 Apply procedure G2 (CC,SS,B(J-1),B(J))   
C
            TEMP = B(J-1)
            B(J-1) = CC*TEMP + SS*B(J)    
            B(J)   =-SS*TEMP + CC*B(J)    
  280    continue
      endif
C
      NPP1=NSETP
      NSETP=NSETP-1     
      IZ1=IZ1-1 
      INDEX(IZ1)=I  
C   
C        SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
C        BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
C        IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY   
C        THAT ARE NONPOSITIVE WILL BE SET TO ZERO   
C        AND MOVED FROM SET P TO SET Z. 
C   
      DO JJ=1,NSETP 
         I=INDEX(JJ)   
         IF (X(I) .le. ZERO) go to 260
      ENDDO 
C   
C         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
C   
      DO I=1,M  
        ZZ(I)=B(I)    
      ENDDO
      RTNKEY = 2
      GO TO 400 
  320 CONTINUE  
      GO TO 210 
C                      &*****  END OF SECONDARY LOOP  ******
C   
  330 continue
      DO IP=1,NSETP 
        I=INDEX(IP)   
        X(I)=ZZ(IP)   
      ENDDO
C        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.  
      GO TO 30  
C   
C                        &*****  END OF MAIN LOOP  ******   
C   
C                        COME TO HERE FOR TERMINATION.  
C                     COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.    
C 
  350 continue  
      SM=ZERO   
      IF (NPP1 .le. M) then
         DO I=NPP1,M   
           SM=SM+B(I)**2 
         ENDDO
      else
         DO J=1,N  
           W(J)=ZERO     
         ENDDO
      endif
      RNORM=sqrt(SM)    
      RETURN
C   
C     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE     
C     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().     
C   
  400 continue
      DO L=1,NSETP  
         IP=NSETP+1-L  
         IF (L .ne. 1) then
            DO II=1,IP
               ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)   
            ENDDO
         endif
         JJ=INDEX(IP)  
         ZZ(IP)=ZZ(IP)/A(IP,JJ)    
      ENDDO
C
      go to (200, 320), RTNKEY
      END SUBROUTINE xNNLS


C*************************************************************************C

      double precision FUNCTION xDIFF(X,Y)
C
C  Function used in tests that depend on machine precision.
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 7, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C
      double precision X, Y
      xDIFF=X-Y  
      RETURN
      END FUNCTION xDIFF  



C*************************************************************************C

      SUBROUTINE xG1 (A,B,CTERM,STERM,SIG)   
C
C     COMPUTE ORTHOGONAL ROTATION MATRIX..  
C
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 12, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C   
C     COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))   
C                        (-S,C)         (-S,C)(B)   (   0          )    
C     COMPUTE SIG = SQRT(A**2+B**2) 
C        SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT 
C        SIG MAY BE IN THE SAME LOCATION AS A OR B .
C     ------------------------------------------------------------------
      double precision A, B, CTERM, ONE, SIG, STERM, XR, YR, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
C     ------------------------------------------------------------------
      if (abs(A) .gt. abs(B)) then
         XR=B/A
         YR=sqrt(ONE+XR**2)
         CTERM=sign(ONE/YR,A)
         STERM=CTERM*XR
         SIG=abs(A)*YR     
         RETURN
      endif

      if (B .ne. ZERO) then
         XR=A/B
         YR=sqrt(ONE+XR**2)
         STERM=sign(ONE/YR,B)
         CTERM=STERM*XR
         SIG=abs(B)*YR     
         RETURN
      endif

      SIG=ZERO  
      CTERM=ZERO  
      STERM=ONE   
      RETURN
      END SUBROUTINE xG1

C*************************************************************************C

C     SUBROUTINE H12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C   
C  CONSTRUCTION AND/OR APPLICATION OF A SINGLE   
C  HOUSEHOLDER TRANSFORMATION..     Q = I + U*(U**T)/B   
C   
C  The original version of this code was developed by
C  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory
C  1973 JUN 12, and published in the book
C  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974.
C  Revised FEB 1995 to accompany reprinting of the book by SIAM.
C     ------------------------------------------------------------------
C                     Subroutine Arguments
C
C     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
C            Householder transformation, or Algorithm H2 to apply a
C            previously constructed transformation.
C     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
C     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
C            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
C            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
C     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
C            vector.  IUE is the storage increment between elements.  
C            On exit when MODE = 1, U() and UP contain quantities
C            defining the vector U of the Householder transformation.
C            on entry with MODE = 2, U() and UP should contain
C            quantities previously computed with MODE = 1.  These will
C            not be modified during the entry with MODE = 2.   
C     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
C            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
C            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
C            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
C     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().  
C     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().  
C     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
C            NO OPERATIONS WILL BE DONE ON C(). 
C     ------------------------------------------------------------------
      SUBROUTINE xH12 (MODE,LPIVOT,L1,M,U,IUE,UP,C,ICE,ICV,NCV)  
C     ------------------------------------------------------------------
      integer I, I2, I3, I4, ICE, ICV, INCR, IUE, J
      integer L1, LPIVOT, M, MODE, NCV
      double precision B, C(*), CL, CLINV, ONE, SM
C     double precision U(IUE,M)
      double precision U(IUE,*)
      double precision UP
      parameter(ONE = 1.0d0)
C     ------------------------------------------------------------------
      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN    
      CL=abs(U(1,LPIVOT))   
      IF (MODE.EQ.2) GO TO 60   
C                            &***** CONSTRUCT THE TRANSFORMATION. ******
          DO J=L1,M  
           CL=MAX(abs(U(1,J)),CL)  
          ENDDO
C KARLINE      IF (CL) 130,130,20
      IF (CL .LE. 0) GOTO 130
      CLINV=ONE/CL  
      SM=(U(1,LPIVOT)*CLINV)**2   
      DO J=L1,M  
        SM=SM+(U(1,J)*CLINV)**2 
      ENDDO
      CL=CL*SQRT(SM)   
C KARLINE     IF (U(1,LPIVOT)) 50,50,40     
      IF (U(1,LPIVOT) .LE. 0) GOTO 50
      CL=-CL
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL    
      GO TO 70  
C            &***** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
C   
CKARLINE   60 IF (CL) 130,130,70
   60 IF (CL .LE. 0) GOTO 130
   70 IF (NCV.LE.0) RETURN  
      B= UP*U(1,LPIVOT)
C                       B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
C   
CKARLINE      IF (B) 80,130,130 
      IF (B .GE. 0) GOTO 130 
      B=ONE/B   
      I2=1-ICV+ICE*(LPIVOT-1)   
      INCR=ICE*(L1-LPIVOT)  
          DO 120 J=1,NCV
          I2=I2+ICV     
          I3=I2+INCR    
          I4=I3 
          SM=C(I2)*UP
           DO I=L1,M  
             SM=SM+C(I3)*U(1,I)
             I3=I3+ICE 
           ENDDO
CKARLINE          IF (SM) 100,120,100   
          IF (SM .EQ. 0) GOTO 120
          SM=SM*B   
          C(I2)=C(I2)+SM*UP
            DO I=L1,M 
              C(I4)=C(I4)+SM*U(1,I)
              I4=I4+ICE 
            ENDDO  
  120     CONTINUE  
  130 RETURN
      END SUBROUTINE xH12





C*************************************************************************C

      DOUBLE PRECISION FUNCTION xDNRM2 (N, DX, INCX)
C***BEGIN PROLOGUE  xDNRM2
C***PURPOSE  Compute the Euclidean length (L2 norm) of a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A3B
C***TYPE      DOUBLE PRECISION (SNRM2-S, xDNRM2-D, SCNRM2-C)
C***KEYWORDS  BLAS, EUCLIDEAN LENGTH, EUCLIDEAN NORM, L2,
C             LINEAR ALGEBRA, UNITARY, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C    xDNRM2  double precision result (zero if N .LE. 0)
C
C     Euclidean norm of the N-vector stored in DX with storage
C     increment INCX.
C     If N .LE. 0, return with result = 0.
C     If N .GE. 1, then INCX must be .GE. 1
C
C     Four phase method using two built-in constants that are
C     hopefully applicable to all machines.
C         CUTLO = maximum of  SQRT(U/EPS)  over all known machines.
C         CUTHI = minimum of  SQRT(V)      over all known machines.
C     where
C         EPS = smallest no. such that EPS + 1. .GT. 1.
C         U   = smallest positive no.   (underflow limit)
C         V   = largest  no.            (overflow  limit)
C
C     Brief outline of algorithm.
C
C     Phase 1 scans zero components.
C     move to phase 2 when a component is nonzero and .LE. CUTLO
C     move to phase 3 when a component is .GT. CUTLO
C     move to phase 4 when a component is .GE. CUTHI/M
C     where M = N for X() real and M = 2*N for complex.
C
C     Values for CUTLO and CUTHI.
C     From the environmental parameters listed in the IMSL converter
C     document the limiting values are as follows:
C     CUTLO, S.P.   U/EPS = 2**(-102) for  Honeywell.  Close seconds are
C                   Univac and DEC at 2**(-103)
C                   Thus CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 for Univac, Honeywell, and DEC.
C                   Thus CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) for Honeywell and DEC.
C                   Thus CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   same as S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C     DATA CUTLO, CUTHI /4.441E-16,  1.304E19/
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  xDNRM2
      INTEGER NEXT, N,NN,INCX,I,J
      DOUBLE PRECISION DX(*), CUTLO, CUTHI, HITEST, SUM, XMAX, ZERO,      &
     &                 ONE
      SAVE CUTLO, CUTHI, ZERO, ONE
      DATA ZERO, ONE /0.0D0, 1.0D0/
C
      DATA CUTLO, CUTHI /8.232D-11,  1.304D19/
C***FIRST EXECUTABLE STATEMENT  xDNRM2
      IF (N .GT. 0) GO TO 10
         XDNRM2  = ZERO
         GO TO 300
C
C   10 ASSIGN 30 TO NEXT               CHANGED INTO:
   10 NEXT = 30      
      SUM = ZERO
      NN = N * INCX
C
C                                                 BEGIN MAIN LOOP
C
      I = 1
C  20    GO TO NEXT,(30, 50, 70, 110)  CHANGED INTO:
   20 SELECT CASE (NEXT)

        CASE (30)
          GOTO 30
        CASE (50)
          GOTO 50
        CASE (70)
          GOTO 70
        CASE (110)
          GOTO 110
                   
      END SELECT
      
      
        
   30 IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C     ASSIGN 50 TO NEXT               CHANGED INTO:
      NEXT = 50

      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF (DX(I) .EQ. ZERO) GO TO 200
      IF (ABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
C
C      ASSIGN 70 TO NEXT              CHANGED INTO:
      NEXT = 70
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
C      ASSIGN 110 TO NEXT              CHANGED INTO:
      NEXT = 110
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = ABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF (ABS(DX(I)) .GT. CUTLO) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF (ABS(DX(I)) .LE. XMAX) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = ABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI / N
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
C KARLINE

C      DO 95 J = I,NN,INCX
C      IF (ABS(DX(J)) .GE. HITEST) GO TO 100
C   95    SUM = SUM + DX(J)**2

      DO J = I,NN,INCX
         IF (ABS(DX(J)) .GE. HITEST) GO TO 100
         SUM = SUM + DX(J)**2
      ENDDO 

      XDNRM2 = SQRT(SUM)
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF (I .LE. NN) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      XDNRM2 = XMAX * SQRT(SUM)
  300 CONTINUE
      RETURN
      END

C*************************************************************************C

C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C     SOLVING SYSTEMS OF EQUATIONS   C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C



C********************************************************************

      subroutine  xdscal(n,da,dx,incx)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C
C     scales a vector by a constant.
C     uses unrolled loops for increment equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
C***BEGIN PROLOGUE  xDSCAL
C***PURPOSE  Multiply a vector by a constant.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A6
C***TYPE      DOUBLE PRECISION (SSCAL-S, xDSCAL-D, CSCAL-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scale factor
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C
C     --Output--
C       DX  double precision result (unchanged if N.LE.0)
C
C     Replace double precision DX by double precision DA*DX.
C     For I = 0 to N-1, replace DX(IX+I*INCX) with  DA * DX(IX+I*INCX),
C     where IX = 1 if INCX .GE. 0, else IX = 1+(1-N)*INCX.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900821  Modified to correct problem with a negative increment.
C           (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  xDSCAL

      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
C
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
C
C        code for increment not equal to 1
C
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
C
C        code for increment equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end                                        
C of xDscal




C********************************************************************
      subroutine xdaxpy(n,da,dx,incx,dy,incy)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C***BEGIN PROLOGUE  xDAXPY
C***PURPOSE  Compute a constant times a vector plus a vector.
C***LIBRARY   SLATEC (BLAS)
C***CATEGORY  D1A7
C***TYPE      DOUBLE PRECISION (SAXPY-S, xDAXPY-D, CAXPY-C)
C***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
C***AUTHOR  Lawson, C. L., (JPL)
C           Hanson, R. J., (SNLA)
C           Kincaid, D. R., (U. of Texas)
C           Krogh, F. T., (JPL)
C***DESCRIPTION
C
C                B L A S  Subprogram
C    Description of Parameters
C
C     --Input--
C        N  number of elements in input vector(s)
C       DA  double precision scalar multiplier
C       DX  double precision vector with N elements
C     INCX  storage spacing between elements of DX
C       DY  double precision vector with N elements
C     INCY  storage spacing between elements of DY
C
C     --Output--
C       DY  double precision result (unchanged if N .LE. 0)
C
C     Overwrite double precision DY with double precision DA*DX + DY.
C     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY),
C     where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
C     defined in a similar way using INCY.
C
C***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
C                 Krogh, Basic linear algebra subprograms for Fortran
C                 usage, Algorithm No. 539, Transactions on Mathematical
C                 Software 5, 3 (September 1979), pp. 308-323.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   791001  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  xDAXPY

C
C     constant times a vector plus a vector.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
C
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C        code for unequal increments or equal increments
C          not equal to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end                                        
C of xdaxpy

C********************************************************************

      double precision function xddot(n,dx,incx,dy,incy)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C
C     forms the dot product of two vectors.
C     uses unrolled loops for increments equal to one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
C
      xddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C        code for unequal increments or equal increments
C          not equal to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      xddot = dtemp
      return
C
C        code for both increments equal to 1
C
C
C        clean-up loop
C
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +                &
     &   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 xddot = dtemp
      return
      end                                        
C of xDDOT



      subroutine xdswap (n,dx,incx,dy,incy)
C
C     interchanges two vectors.
C     uses unrolled loops for increments equal one.
C     jack dongarra, linpack, 3/11/78.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
C
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
C
C       code for unequal increments or equal increments not equal
C         to 1
C
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
C
C       code for both increments equal to 1
C
C
C       clean-up loop
C
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end


C***********************************************************************
      SUBROUTINE XDCOPY(N,DX,INCX,DY,INCY)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(*),DY(*)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END


C********************************************************************
      integer function xidamax(n,dx,incx)

C-------------------------------------------------------------------*
C Subroutine obtained from LINPACK (ftp-site netlib.att.com)        &
C-------------------------------------------------------------------*
C
C     finds the index of element having max. absolute value.
C     jack dongarra, linpack, 3/11/78.
C     modified 3/93 to return if incx .le. 0.
C     modified 12/3/93, array(1) declarations changed to array(*)
C
      double precision dx(*),dmax
      integer i,incx,ix,n
C
      xidamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      xidamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
C
C        code for increment not equal to 1
C
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         xidamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
C
C        code for increment equal to 1
C
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         xidamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end                                        
C of xidamax


C***********************************************************************

      DOUBLE PRECISION FUNCTION D1MACH (IDUM)
      INTEGER  :: IDUM
C-----------------------------------------------------------------------
C This routine computes:
C the unit roundoff of the machine when IDUM = 4
C the largest value                when IDUM = 2
C Unit roundoffis defined as the smallest positive machine number
C u such that  1.0 + u .ne. 1.0
C Largest value is simply imposed...
C Subroutines/functions called by D1MACH.. None
C-----------------------------------------------------------------------
      DOUBLE PRECISION U, COMP
      DOUBLE PRECISION :: Prec(4)  
      LOGICAL          :: First(4) 
      SAVE Prec, FIRST
      DATA FIRST /.TRUE.,.TRUE.,.TRUE.,.TRUE./
      DATA Prec /1.D-8,1.D-8,1.D-8,1.D-8/


      IF (Idum > 4 .OR. Idum < 0) THEN
C         Write (*,*) "Error in function D1MACH"
C         Write (*,*) "NOT DEFINED FOR IDUM = ", Idum
       CALL XMESSAGE("Error in function D1MACH-NOT DEFINED FOR IDUM  ") 
      ENDIF

      IF (First(Idum)) THEN 

       First(Idum) = .FALSE.

       SELECT CASE (IDUM)

        CASE (2)
C Very large number
         D1MACH = 1.D300

        CASE (4)
C Unit roundoff
         U = 1.0D0
 10      U = U*0.5D0
         COMP = 1.0D0 + U
         IF (COMP .NE. 1.0D0) GO TO 10
         D1MACH = U*2.0D0
        CASE Default
C         Write (*,*) "Error in function D1MACH"
C         Write (*,*) "NOT DEFINED FOR IDUM = ", Idum
         CALL XMESSAGE("Error in function D1MACH-NOT DEFINED FOR IDUM ")
       END SELECT

       PREC (Idum) = D1MACH

      ELSE

       D1mach = Prec(IDUM)

      ENDIF

      RETURN

      END FUNCTION d1Mach





C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<C
C    LEAST SQUARES with constraints  C
C                DLSEI               C                 
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>C


C   LINPACK routine





c                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
c                !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
c                !    LEAST SQUARES with constraints  !
c                !                DLSEI               !                 
c                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
c                !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>!


c   LINPACK routine

      SUBROUTINE xDLSEI (W, MDW, ME, MA, MG, N, PRGOPT, X, RNORME,       &
     &   RNORML, MODE, WS, IP)
c***BEGIN PROLOGUE  xDLSEI
c***PURPOSE  Solve a linearly constrained least squares problem with
c            equality and inequality constraints, and optionally compute
c            a covariance matrix.
c***LIBRARY   SLATEC
c***CATEGORY  K1A2A, D9
c***TYPE      DOUBLE PRECISION (LSEI-S, xDLSEI-D)
c***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
c             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
c             QUADRATIC PROGRAMMING
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     Abstract
c
c     This subprogram solves a linearly constrained least squares
c     problem with both equality and inequality constraints, and, if the
c     user requests, obtains a covariance matrix of the solution
c     parameters.
c
c     Suppose there are given matrices E, A and G of respective
c     dimensions ME by N, MA by N and MG by N, and vectors F, B and H of
c     respective lengths ME, MA and MG.  This subroutine solves the
c     linearly constrained least squares problem
c
c                   EX = F, (E ME by N) (equations to be exactly
c                                       satisfied)
c                   AX = B, (A MA by N) (equations to be
c                                       approximately satisfied,
c                                       least squares sense)
c                   GX .GE. H,(G MG by N) (inequality constraints)
c
c     The inequalities GX .GE. H mean that every component of the
c     product GX must be .GE. the corresponding component of H.
c
c     In case the equality constraints cannot be satisfied, a
c     generalized inverse solution residual vector length is obtained
c     for F-EX.  This is the minimal length possible for F-EX.
c
c     Any values ME .GE. 0, MA .GE. 0, or MG .GE. 0 are permitted.  The
c     rank of the matrix E is estimated during the computation.  We call
c     this value KRANKE.  It is an output parameter in IP(1) defined
c     below.  Using a generalized inverse solution of EX=F, a reduced
c     least squares problem with inequality constraints is obtained.
c     The tolerances used in these tests for determining the rank
c     of E and the rank of the reduced least squares problem are
c     given in Sandia Tech. Rept. SAND-78-1290.  They can be
c     modified by the user if new values are provided in
c     the option list of the array PRGOPT(*).
c
c     The user must dimension all arrays appearing in the call list..
c     W(MDW,N+1),PRGOPT(*),X(N),WS(2*(ME+N)+K+(MG+2)*(N+7)),IP(MG+2*N+2)
c     where K=MAX(MA+MG,N).  This allows for a solution of a range of
c     problems in the given working space.  The dimension of WS(*)
c     given is a necessary overestimate.  Once a particular problem
c     has been run, the output parameter IP(3) gives the actual
c     dimension required for that problem.
c
c     The parameters for xDLSEI( ) are
c
c     Input.. All TYPE REAL variables are DOUBLE PRECISION
c
c     W(*,*),MDW,   The array W(*,*) is doubly subscripted with
c     ME,MA,MG,N    first dimensioning parameter equal to MDW.
c                   For this discussion let us call M = ME+MA+MG.  Then
c                   MDW must satisfy MDW .GE. M.  The condition
c                   MDW .LT. M is an error.
c
c                   The array W(*,*) contains the matrices and vectors
c
c                                  (E  F)
c                                  (A  B)
c                                  (G  H)
c
c                   in rows and columns 1,...,M and 1,...,N+1
c                   respectively.
c
c                   The integers ME, MA, and MG are the
c                   respective matrix row dimensions
c                   of E, A and G.  Each matrix has N columns.
c
c     PRGOPT(*)    This real-valued array is the option vector.
c                  If the user is satisfied with the nominal
c                  subprogram features set
c
c                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
c
c                  Otherwise PRGOPT(*) is a linked list consisting of
c                  groups of data of the following form
c
c                  LINK
c                  KEY
c                  DATA SET
c
c                  The parameters LINK and KEY are each one word.
c                  The DATA SET can be comprised of several words.
c                  The number of items depends on the value of KEY.
c                  The value of LINK points to the first
c                  entry of the next group of data within
c                  PRGOPT(*).  The exception is when there are
c                  no more options to change.  In that
c                  case, LINK=1 and the values KEY and DATA SET
c                  are not referenced.  The general layout of
c                  PRGOPT(*) is as follows.
c
c               ...PRGOPT(1) = LINK1 (link to first entry of next group)
c               .  PRGOPT(2) = KEY1 (key to the option change)
c               .  PRGOPT(3) = data value (data value for this change)
c               .       .
c               .       .
c               .       .
c               ...PRGOPT(LINK1)   = LINK2 (link to the first entry of
c               .                       next group)
c               .  PRGOPT(LINK1+1) = KEY2 (key to the option change)
c               .  PRGOPT(LINK1+2) = data value
c               ...     .
c               .       .
c               .       .
c               ...PRGOPT(LINK) = 1 (no more options to change)
c
c                  Values of LINK that are nonpositive are errors.
c                  A value of LINK .GT. NLINK=100000 is also an error.
c                  This helps prevent using invalid but positive
c                  values of LINK that will probably extend
c                  beyond the program limits of PRGOPT(*).
c                  Unrecognized values of KEY are ignored.  The
c                  order of the options is arbitrary and any number
c                  of options can be changed with the following
c                  restriction.  To prevent cycling in the
c                  processing of the option array, a count of the
c                  number of options changed is maintained.
c                  Whenever this count exceeds NOPT=1000, an error
c                  message is printed and the subprogram returns.
c
c                  Options..
c
c                  KEY=1
c                         Compute in W(*,*) the N by N
c                  covariance matrix of the solution variables
c                  as an output parameter.  Nominally the
c                  covariance matrix will not be computed.
c                  (This requires no user input.)
c                  The data set for this option is a single value.
c                  It must be nonzero when the covariance matrix
c                  is desired.  If it is zero, the covariance
c                  matrix is not computed.  When the covariance matrix
c                  is computed, the first dimensioning parameter
c                  of the array W(*,*) must satisfy MDW .GE. MAX(M,N).
c
c                  KEY=10
c                         Suppress scaling of the inverse of the
c                  normal matrix by the scale factor RNORM**2/
c                  MAX(1, no. of degrees of freedom).  This option
c                  only applies when the option for computing the
c                  covariance matrix (KEY=1) is used.  With KEY=1 and
c                  KEY=10 used as options the unscaled inverse of the
c                  normal matrix is returned in W(*,*).
c                  The data set for this option is a single value.
c                  When it is nonzero no scaling is done.  When it is
c                  zero scaling is done.  The nominal case is to do
c                  scaling so if option (KEY=1) is used alone, the
c                  matrix will be scaled on output.
c
c                  KEY=2
c                         Scale the nonzero columns of the
c                         entire data matrix.
c                  (E)
c                  (A)
c                  (G)
c
c                  to have length one.  The data set for this
c                  option is a single value.  It must be
c                  nonzero if unit length column scaling
c                  is desired.
c
c                  KEY=3
c                         Scale columns of the entire data matrix
c                  (E)
c                  (A)
c                  (G)
c
c                  with a user-provided diagonal matrix.
c                  The data set for this option consists
c                  of the N diagonal scaling factors, one for
c                  each matrix column.
c
c                  KEY=4
c                         Change the rank determination tolerance for
c                  the equality constraint equations from
c                  the nominal value of SQRT(DRELPR).  This quantity can
c                  be no smaller than DRELPR, the arithmetic-
c                  storage precision.  The quantity DRELPR is the
c                  largest positive number such that T=1.+DRELPR
c                  satisfies T .EQ. 1.  The quantity used
c                  here is internally restricted to be at
c                  least DRELPR.  The data set for this option
c                  is the new tolerance.
c
c                  KEY=5
c                         Change the rank determination tolerance for
c                  the reduced least squares equations from
c                  the nominal value of SQRT(DRELPR).  This quantity can
c                  be no smaller than DRELPR, the arithmetic-
c                  storage precision.  The quantity used
c                  here is internally restricted to be at
c                  least DRELPR.  The data set for this option
c                  is the new tolerance.
c
c                  For example, suppose we want to change
c                  the tolerance for the reduced least squares
c                  problem, compute the covariance matrix of
c                  the solution parameters, and provide
c                  column scaling for the data matrix.  For
c                  these options the dimension of PRGOPT(*)
c                  must be at least N+9.  The Fortran statements
c                  defining these options would be as follows:
c
c                  PRGOPT(1)=4 (link to entry 4 in PRGOPT(*))
c                  PRGOPT(2)=1 (covariance matrix key)
c                  PRGOPT(3)=1 (covariance matrix wanted)
c
c                  PRGOPT(4)=7 (link to entry 7 in PRGOPT(*))
c                  PRGOPT(5)=5 (least squares equas.  tolerance key)
c                  PRGOPT(6)=... (new value of the tolerance)
c
c                  PRGOPT(7)=N+9 (link to entry N+9 in PRGOPT(*))
c                  PRGOPT(8)=3 (user-provided column scaling key)
c
c                  CALL xDCOPY (N, D, 1, PRGOPT(9), 1)  (Copy the N
c                    scaling factors from the user array D(*)
c                    to PRGOPT(9)-PRGOPT(N+8))
c
c                  PRGOPT(N+9)=1 (no more options to change)
c
c                  The contents of PRGOPT(*) are not modified
c                  by the subprogram.
c                  The options for WNNLS( ) can also be included
c                  in this array.  The values of KEY recognized
c                  by WNNLS( ) are 6, 7 and 8.  Their functions
c                  are documented in the usage instructions for
c                  subroutine WNNLS( ).  Normally these options
c                  do not need to be modified when using xDLSEI( ).
c
c     IP(1),       The amounts of working storage actually
c     IP(2)        allocated for the working arrays WS(*) and
c                  IP(*), respectively.  These quantities are
c                  compared with the actual amounts of storage
c                  needed by xDLSEI( ).  Insufficient storage
c                  allocated for either WS(*) or IP(*) is an
c                  error.  This feature was included in xDLSEI( )
c                  because miscalculating the storage formulas
c                  for WS(*) and IP(*) might very well lead to
c                  subtle and hard-to-find execution errors.
c
c                  The length of WS(*) must be at least
c
c                  LW = 2*(ME+N)+K+(MG+2)*(N+7)
c
c                  where K = max(MA+MG,N)
c                  This test will not be made if IP(1).LE.0.
c
c                  The length of IP(*) must be at least
c
c                  LIP = MG+2*N+2
c                  This test will not be made if IP(2).LE.0.
c
c     Output.. All TYPE REAL variables are DOUBLE PRECISION
c
c     X(*),RNORME,  The array X(*) contains the solution parameters
c     RNORML        if the integer output flag MODE = 0 or 1.
c                   The definition of MODE is given directly below.
c                   When MODE = 0 or 1, RNORME and RNORML
c                   respectively contain the residual vector
c                   Euclidean lengths of F - EX and B - AX.  When
c                   MODE=1 the equality constraint equations EX=F
c                   are contradictory, so RNORME .NE. 0.  The residual
c                   vector F-EX has minimal Euclidean length.  For
c                   MODE .GE. 2, none of these parameters is defined.
c
c     MODE          Integer flag that indicates the subprogram
c                   status after completion.  If MODE .GE. 2, no
c                   solution has been computed.
c
c                   MODE =
c
c                   0  Both equality and inequality constraints
c                      are compatible and have been satisfied.
c
c                   1  Equality constraints are contradictory.
c                      A generalized inverse solution of EX=F was used
c                      to minimize the residual vector length F-EX.
c                      In this sense, the solution is still meaningful.
c
c                   2  Inequality constraints are contradictory.
c
c                   3  Both equality and inequality constraints
c                      are contradictory.
c
c                   The following interpretation of
c                   MODE=1,2 or 3 must be made.  The
c                   sets consisting of all solutions
c                   of the equality constraints EX=F
c                   and all vectors satisfying GX .GE. H
c                   have no points in common.  (In
c                   particular this does not say that
c                   each individual set has no points
c                   at all, although this could be the
c                   case.)
c
c                   4  Usage error occurred.  The value
c                      of MDW is .LT. ME+MA+MG, MDW is
c                      .LT. N and a covariance matrix is
c                      requested, or the option vector
c                      PRGOPT(*) is not properly defined,
c                      or the lengths of the working arrays
c                      WS(*) and IP(*), when specified in
c                      IP(1) and IP(2) respectively, are not
c                      long enough.
c
c     W(*,*)        The array W(*,*) contains the N by N symmetric
c                   covariance matrix of the solution parameters,
c                   provided this was requested on input with
c                   the option vector PRGOPT(*) and the output
c                   flag is returned with MODE = 0 or 1.
c
c     IP(*)         The integer working array has three entries
c                   that provide rank and working array length
c                   information after completion.
c
c                      IP(1) = rank of equality constraint
c                              matrix.  Define this quantity
c                              as KRANKE.
c
c                      IP(2) = rank of reduced least squares
c                              problem.
c
c                      IP(3) = the amount of storage in the
c                              working array WS(*) that was
c                              actually used by the subprogram.
c                              The formula given above for the length
c                              of WS(*) is a necessary overestimate.
c                              If exactly the same problem matrices
c                              are used in subsequent executions,
c                              the declared dimension of WS(*) can
c                              be reduced to this output value.
c     User Designated
c     Working Arrays..
c
c     WS(*),IP(*)              These are respectively type real
c                              and type integer working arrays.
c                              Their required minimal lengths are
c                              given above.
c
c***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
c                 linear least squares problems with equality and
c                 nonnegativity constraints, Report SAND77-0552, Sandia
c                 Laboratories, June 1978.
c               K. H. Haskell and R. J. Hanson, Selected algorithms for
c                 the linearly constrained least squares problem - a
c                 users guide, Report SAND78-1290, Sandia Laboratories,
c                 August 1979.
c               K. H. Haskell and R. J. Hanson, An algorithm for
c                 linear least squares problems with equality and
c                 nonnegativity constraints, Mathematical Programming
c                 21 (1981), pp. 98-118.
c               R. J. Hanson and K. H. Haskell, Two algorithms for the
c                 linearly constrained least squares problem, ACM
c                 Transactions on Mathematical Software, September 1982.
c***ROUTINES CALLED  D1MACH, xDASUM, xDAXPY, xDCOPY, xDDOT, xDH12, DLSI,
c                    xDNRM2, xDSCAL, xDSWAP, xXERMSG
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890618  Completely restructured and extensively revised (WRB & RWC)
c   890831  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
c   900510  Convert XERRWV calls to xXERMSG calls.  (RWC)
c   900604  DP version created from SP version.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  xDLSEI
      INTEGER IP(3), MA, MDW, ME, MG, MODE, N
      DOUBLE PRECISION PRGOPT(*), RNORME, RNORML, W(MDW,*), WS(*), X(*)
c
      EXTERNAL D1MACH, xDASUM, xDAXPY, xDCOPY,xDDOT,xDH12,DLSI,xDNRM2,               &
     &   xDSCAL, xDSWAP, xXERMSG
      DOUBLE PRECISION D1MACH, xDASUM, xDDOT, xDNRM2
c KARLINE: 
      DOUBLE PRECISION DRELPR, ENORM, FNORM, GAM, RB, RN, RNMAX, SIZE,             &
     &   SN, SNMAX, T, TAU, UJ, UP, VJ, XNORM, XNRME
      INTEGER I, IMAX, J, JP1, K, KEY, KRANKE, LAST, LCHK, LINK, M,                &
     &   MAPKE1, MDEQC, MEND, MEP1, N1, N2, NEXT, NLINK, NOPT, NP1,                &
     &   NTIMES
      LOGICAL COV, FIRST
      CHARACTER*8 XERN1, XERN2, XERN3, XERN4
      SAVE FIRST, DRELPR
c
      DATA FIRST /.TRUE./
c***FIRST EXECUTABLE STATEMENT  xDLSEI
c
c     Set the nominal tolerance used in the code for the equality
c     constraint equations.
c
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
      TAU = SQRT(DRELPR)
c
c     Check that enough storage was allocated in WS(*) and IP(*).
c
      MODE = 4
      IF (MIN(N,ME,MA,MG) .LT. 0) THEN
c         WRITE (XERN1, '(I8)') N
c         WRITE (XERN2, '(I8)') ME
c         WRITE (XERN3, '(I8)') MA
c         WRITE (XERN4, '(I8)') MG
C KARLINE: REMOVED WRITE		 
         CALL rwarn ('LSEI: THE VARIABLES N, ME,MA, MG MUST BE>0')
c         CALL xXERMSG ('SLATEC', 'LSEI', 'ALL OF THE VARIABLES N, ME,'//            &
c     &      ' MA, MG MUST BE .GE. 0 ENTERED ROUTINE WITH' //                        &
c     &      ' N  = ' // XERN1 //                                                    &
c     &      ' ME = ' // XERN2 //                                                    &
c     &      ' MA = ' // XERN3 //                                                    &
c     &      ' MG = ' // XERN4, 2, 1)
         RETURN
      ENDIF
c
      IF (IP(1).GT.0) THEN
         LCHK = 2*(ME+N) + MAX(MA+MG,N) + (MG+2)*(N+7)
         IF (IP(1).LT.LCHK) THEN
C KARLINE: REMOVED WRITE		 
         CALL rwarn ('LSEI: insufficient storage')
c            WRITE (XERN1, '(I8)') LCHK
c            CALL xXERMSG ('SLATEC', 'xDLSEI', 'INSUFFICIENT STORAGE ' //             &
c     &         'ALLOCATED FOR WS(*), NEED LW = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
c
      IF (IP(2).GT.0) THEN
         LCHK = MG + 2*N + 2
         IF (IP(2).LT.LCHK) THEN
C KARLINE: REMOVED WRITE		 
         CALL rwarn ('LSEI: insufficient storage')
c            WRITE (XERN1, '(I8)') LCHK
c            CALL xXERMSG ('SLATEC', 'xDLSEI', 'INSUFFICIENT STORAGE ' //             &
c     &         'ALLOCATED FOR IP(*), NEED LIP = ' // XERN1, 2, 1)
            RETURN
         ENDIF
      ENDIF
c
c     Compute number of possible right multiplying Householder
c     transformations.
c
      M = ME + MA + MG
      IF (N.LE.0 .OR. M.LE.0) THEN
         MODE = 0
         RNORME = 0
         RNORML = 0
         RETURN
      ENDIF
c
      IF (MDW.LT.M) THEN
        CALL xXERMSG ('SLATEC', 'xDLSEI', 'MDW.LT.ME+MA+MG IS AN ERROR',            &
     &      2, 1)
         RETURN
      ENDIF
c
      NP1 = N + 1
      KRANKE = MIN(ME,N)
      N1 = 2*KRANKE + 1
      N2 = N1 + N
c
c     Set nominal values.
c
c     The nominal column scaling used in the code is
c     the identity scaling.
c
      CALL XDCOPYSC (N, 1.D0, WS(N1), 1)
c
c     No covariance matrix is nominally computed.
c
      COV = .FALSE.
c
c     Process option vector.
c     Define bound for number of options to change.
c
      NOPT = 1000
      NTIMES = 0
c
c     Define bound for positive values of LINK.
c
      NLINK = 100000
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.EQ.0 .OR. LINK.GT.NLINK) THEN
         CALL xXERMSG('SLATEC','xDLSEI','THE OPTION VECTOR IS UNDEFINED'            &
     &   ,2,1)
         RETURN
      ENDIF
c
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
            CALL xXERMSG ('SLATEC','xDLSEI',                                         &
     &         'THE LINKS IN THE OPTION VECTOR ARE CYCLING.', 2, 1)
            RETURN
         ENDIF
c
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) THEN
            COV = PRGOPT(LAST+2) .NE. 0.D0
         ELSEIF (KEY.EQ.2 .AND. PRGOPT(LAST+2).NE.0.D0) THEN
            DO 110 J = 1,N
               T = xDNRM2(M,W(1,J),1)
               IF (T.NE.0.D0) T = 1.D0/T
               WS(J+N1-1) = T
  110       CONTINUE
         ELSEIF (KEY.EQ.3) THEN
            CALL xDCOPY (N, PRGOPT(LAST+2), 1, WS(N1), 1)
         ELSEIF (KEY.EQ.4) THEN
            TAU = MAX(DRELPR,PRGOPT(LAST+2))
         ENDIF
c
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
         CALL xXERMSG ('SLATEC', 'xDLSEI',                                           &
     &      'THE OPTION VECTOR IS UNDEFINED', 2, 1)
            RETURN
         ENDIF
c
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
c
      DO 120 J = 1,N
         CALL xDSCAL (M, WS(N1+J-1), W(1,J), 1)
  120 CONTINUE
c
      IF (COV .AND. MDW.LT.N) THEN
         CALL xXERMSG ('SLATEC', 'xDLSEI',                                           &
     &      'MDW .LT. N WHEN COV MATRIX NEEDED, IS AN ERROR', 2, 1)
         RETURN
      ENDIF
c
c     Problem definition and option vector OK.
c
      MODE = 0
c
c     Compute norm of equality constraint matrix and right side.
c
      ENORM = 0.D0
      DO 130 J = 1,N
         ENORM = MAX(ENORM,xDASUM(ME,W(1,J),1))
  130 CONTINUE
c
      FNORM = xDASUM(ME,W(1,NP1),1)
      SNMAX = 0.D0
      RNMAX = 0.D0
      DO 150 I = 1,KRANKE
c
c        Compute maximum ratio of vector lengths. Partition is at
c        column I.
c
         DO 140 K = I,ME
            SN = xDDOT(N-I+1,W(K,I),MDW,W(K,I),MDW)
            RN = xDDOT(I-1,W(K,1),MDW,W(K,1),MDW)
            IF (RN.EQ.0.D0 .AND. SN.GT.SNMAX) THEN
               SNMAX = SN
               IMAX = K
            ELSEIF (K.EQ.I .OR. SN*RNMAX.GT.RN*SNMAX) THEN
               SNMAX = SN
               RNMAX = RN
               IMAX = K
            ENDIF
  140    CONTINUE
c
c        Interchange rows if necessary.
c
         IF (I.NE.IMAX) CALL xDSWAP (NP1, W(I,1), MDW, W(IMAX,1), MDW)
         IF (SNMAX.GT.RNMAX*TAU**2) THEN
c
c        Eliminate elements I+1,...,N in row I.
c
            CALL xDH12 (1, I, I+1, N, W(I,1), MDW, WS(I), W(I+1,1), MDW,             &    
     &                1, M-I)
         ELSE
            KRANKE = I - 1
            GO TO 160
         ENDIF
  150 CONTINUE
c
c     Save diagonal terms of lower trapezoidal matrix.
c
  160 CALL xDCOPY (KRANKE, W, MDW+1, WS(KRANKE+1), 1)
c
c     Use Householder transformation from left to achieve
c     KRANKE by KRANKE upper triangular form.
c
      IF (KRANKE.LT.ME) THEN
         DO 170 K = KRANKE,1,-1
c
c           Apply transformation to matrix cols. 1,...,K-1.
c
         CALL xDH12 (1, K, KRANKE+1, ME, W(1,K), 1, UP, W, 1, MDW,K-1)
c
c           Apply to rt side vector.
c
        CALL xDH12 (2, K, KRANKE+1, ME, W(1,K), 1, UP, W(1,NP1),1,1, 1)
  170    CONTINUE
      ENDIF
c
c     Solve for variables 1,...,KRANKE in new coordinates.
c
      CALL xDCOPY (KRANKE, W(1, NP1), 1, X, 1)
      DO 180 I = 1,KRANKE
         X(I) = (X(I)-xDDOT(I-1,W(I,1),MDW,X,1))/W(I,I)
  180 CONTINUE
c
c     Compute residuals for reduced problem.
c
      MEP1 = ME + 1
      RNORML = 0.D0
      DO 190 I = MEP1,M
         W(I,NP1) = W(I,NP1) - xDDOT(KRANKE,W(I,1),MDW,X,1)
         SN = xDDOT(KRANKE,W(I,1),MDW,W(I,1),MDW)
         RN = xDDOT(N-KRANKE,W(I,KRANKE+1),MDW,W(I,KRANKE+1),MDW)
         IF (RN.LE.SN*TAU**2 .AND. KRANKE.LT.N)                                     &
     &      CALL XDCOPYSC (N-KRANKE, 0.D0, W(I,KRANKE+1), MDW)
  190 CONTINUE
c
c     Compute equality constraint equations residual length.
c
      RNORME = xDNRM2(ME-KRANKE,W(KRANKE+1,NP1),1)
c
c     Move reduced problem data upward if KRANKE.LT.ME.
c
      IF (KRANKE.LT.ME) THEN
         DO 200 J = 1,NP1
            CALL xDCOPY (M-ME, W(ME+1,J), 1, W(KRANKE+1,J), 1)
  200    CONTINUE
      ENDIF
c
c     Compute solution of reduced problem.
c
      CALL DLSI(W(KRANKE+1, KRANKE+1), MDW, MA, MG, N-KRANKE, PRGOPT,               &
     &         X(KRANKE+1), RNORML, MODE, WS(N2), IP(2))
c
c     Test for consistency of equality constraints.
c
      IF (ME.GT.0) THEN
         MDEQC = 0
         XNRME = xDASUM(KRANKE,W(1,NP1),1)
         IF (RNORME.GT.TAU*(ENORM*XNRME+FNORM)) MDEQC = 1
         MODE = MODE + MDEQC
c
c        Check if solution to equality constraints satisfies inequality
c        constraints when there are no degrees of freedom left.
c
         IF (KRANKE.EQ.N .AND. MG.GT.0) THEN
            XNORM = xDASUM(N,X,1)
            MAPKE1 = MA + KRANKE + 1
            MEND = MA + KRANKE + MG
            DO 210 I = MAPKE1,MEND
               SIZE = xDASUM(N,W(I,1),MDW)*XNORM + ABS(W(I,NP1))
               IF (W(I,NP1).GT.TAU*SIZE) THEN
                  MODE = MODE + 2
                  GO TO 290
               ENDIF
  210       CONTINUE
         ENDIF
      ENDIF
c
c     Replace diagonal terms of lower trapezoidal matrix.
c
      IF (KRANKE.GT.0) THEN
         CALL xDCOPY (KRANKE, WS(KRANKE+1), 1, W, MDW+1)
c
c        Reapply transformation to put solution in original coordinates.
c
         DO 220 I = KRANKE,1,-1
            CALL xDH12 (2, I, I+1, N, W(I,1), MDW, WS(I), X, 1, 1, 1)
  220    CONTINUE
c
c        Compute covariance matrix of equality constrained problem.
c
         IF (COV) THEN
            DO 270 J = MIN(KRANKE,N-1),1,-1
               RB = WS(J)*W(J,J)
               IF (RB.NE.0.D0) RB = 1.D0/RB
               JP1 = J + 1
               DO 230 I = JP1,N
                  W(I,J) = RB*xDDOT(N-J,W(I,JP1),MDW,W(J,JP1),MDW)
  230          CONTINUE
c
               GAM = 0.5D0*RB*xDDOT(N-J,W(JP1,J),1,W(J,JP1),MDW)
               CALL xDAXPY (N-J, GAM, W(J,JP1), MDW, W(JP1,J), 1)
               DO 250 I = JP1,N
                  DO 240 K = I,N
                     W(I,K) = W(I,K) + W(J,I)*W(K,J) + W(I,J)*W(J,K)
                     W(K,I) = W(I,K)
  240             CONTINUE
  250          CONTINUE
               UJ = WS(J)
               VJ = GAM*UJ
               W(J,J) = UJ*VJ + UJ*VJ
               DO 260 I = JP1,N
                  W(J,I) = UJ*W(I,J) + VJ*W(J,I)
  260          CONTINUE
               CALL xDCOPY (N-J, W(J, JP1), MDW, W(JP1,J), 1)
  270       CONTINUE
         ENDIF
      ENDIF
c
c     Apply the scaling to the covariance matrix.
c
      IF (COV) THEN
         DO 280 I = 1,N
            CALL xDSCAL (N, WS(I+N1-1), W(I,1), MDW)
            CALL xDSCAL (N, WS(I+N1-1), W(1,I), 1)
  280    CONTINUE
      ENDIF
c
c     Rescale solution vector.
c
  290 IF (MODE.LE.1) THEN
         DO 300 J = 1,N
            X(J) = X(J)*WS(N1+J-1)
  300    CONTINUE
      ENDIF
c
      IP(1) = KRANKE
      IP(3) = IP(3) + 2*KRANKE + N
      RETURN
      END





      SUBROUTINE DLSI (W, MDW, MA, MG, N, PRGOPT, X, RNORM, MODE, WS,               &
     &   IP)
c***BEGIN PROLOGUE  DLSI
c***SUBSIDIARY
c***PURPOSE  Subsidiary to xDLSEI
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (LSI-S, DLSI-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c***DESCRIPTION
c
c     This is a companion subprogram to xDLSEI.  The documentation for
c     xDLSEI has complete usage instructions.
c
c     Solve..
c              AX = B,  A  MA by N  (least squares equations)
c     subject to..
c
c              GX.GE.H, G  MG by N  (inequality constraints)
c
c     Input..
c
c      W(*,*) contains  (A B) in rows 1,...,MA+MG, cols 1,...,N+1.
c                       (G H)
c
c     MDW,MA,MG,N
c              contain (resp) var. dimension of W(*,*),
c              and matrix dimensions.
c
c     PRGOPT(*),
c              Program option vector.
c
c     OUTPUT..
c
c      X(*),RNORM
c
c              Solution vector(unless MODE=2), length of AX-B.
c
c      MODE
c              =0   Inequality constraints are compatible.
c              =2   Inequality constraints contradictory.
c
c      WS(*),
c              Working storage of dimension K+N+(MG+2)*(N+7),
c              where K=MAX(MA+MG,N).
c      IP(MG+2*N+1)
c              Integer working storage
c
c***ROUTINES CALLED  D1MACH, xDASUM, xDAXPY, xDCOPY, xDDOT, xDH12, xDHFTI,
c                    DLPDP, xDSCAL, xDSWAP
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890618  Completely restructured and extensively revised (WRB & RWC)
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900328  Added TYPE section.  (WRB)
c   900604  DP version created from SP version.  (RWC)
c   920422  Changed CALL to xDHFTI to include variable MA.  (WRB)
c***END PROLOGUE  DLSI  - karline: added RNORMV(1) to avoid warning when calling xdhfti     
      INTEGER IP(*), MA, MDW, MG, MODE, N
      DOUBLE PRECISION PRGOPT(*), RNORM, W(MDW,*),WS(*),X(*),RNORMV(1)
c
      EXTERNAL D1MACH, xDASUM,xDAXPY,xDCOPY,xDDOT,xDH12,xDHFTI,DLPDP,               &
     &   xDSCAL, xDSWAP
      DOUBLE PRECISION D1MACH, xDASUM, xDDOT
c
      DOUBLE PRECISION ANORM, DRELPR, FAC, GAM, RB, TAU, TOL, XNORM
      INTEGER I, J, K, KEY, KRANK, KRM1, KRP1, L, LAST, LINK, M, MAP1,              &
     &   MDLPDP, MINMAN, N1, N2, N3, NEXT, NP1, MDB 
      LOGICAL COV, FIRST, SCLCOV
c
      SAVE DRELPR, FIRST
      DATA FIRST /.TRUE./
c
c***FIRST EXECUTABLE STATEMENT  DLSI
c
c     Set the nominal tolerance used in the code.
c
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
      TOL = SQRT(DRELPR)
c
      MODE = 0
      RNORM = 0.D0
      M = MA + MG
      NP1 = N + 1
      KRANK = 0
      IF (N.LE.0 .OR. M.LE.0) GO TO 370
c
c     To process option vector.
c
      COV = .FALSE.
      SCLCOV = .TRUE.
      LAST = 1
      LINK = PRGOPT(1)
c
  100 IF (LINK.GT.1) THEN
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.1) COV = PRGOPT(LAST+2) .NE. 0.D0
         IF (KEY.EQ.10) SCLCOV = PRGOPT(LAST+2) .EQ. 0.D0
         IF (KEY.EQ.5) TOL = MAX(DRELPR,PRGOPT(LAST+2))
         NEXT = PRGOPT(LINK)
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
c
c     Compute matrix norm of least squares equations.
c
      ANORM = 0.D0
      DO 110 J = 1,N
         ANORM = MAX(ANORM,xDASUM(MA,W(1,J),1))
  110 CONTINUE
c
c     Set tolerance for xDHFTI( ) rank test.
c
      TAU = TOL*ANORM
c
c     Compute Householder orthogonal decomposition of matrix.
c
      CALL XDCOPYSC (N, 0.D0, WS, 1)
      CALL xDCOPY (MA, W(1, NP1), 1, WS, 1)
      K = MAX(M,N)
      MINMAN = MIN(MA,N)
      N1 = K + 1
      N2 = N1 + N
c      CALL xDHFTI (W, MDW, MA, N, WS, MA, 1, TAU, KRANK, RNORM, WS(N2),  WS(N1), IP)
c KARLINE:ADDED both next sentences ...
        RNORMV(1) = RNORM       
        MDB = MAX(MA,N)
      CALL xDHFTI (W, MDW, MA, N, WS, MDB, 1,TAU,KRANK, RNORMV, WS(N2),             &
     &           WS(N1), IP)   
      RNORM = RNORMV(1)   ! and this one added as well
c and changed that...
      FAC = 1.D0
      GAM = MA - KRANK
      IF (KRANK.LT.MA .AND. SCLCOV) FAC = RNORM**2/GAM
c
c     Reduce to DLPDP and solve.
c
      MAP1 = MA + 1
c
c     Compute inequality rt-hand side for DLPDP.
c
      IF (MA.LT.M) THEN
         IF (MINMAN.GT.0) THEN
            DO 120 I = MAP1,M
               W(I,NP1) = W(I,NP1) - xDDOT(N,W(I,1),MDW,WS,1)
  120       CONTINUE
c
c           Apply permutations to col. of inequality constraint matrix.
c
            DO 130 I = 1,MINMAN
               CALL xDSWAP (MG, W(MAP1,I), 1, W(MAP1,IP(I)), 1)
  130       CONTINUE
c
c           Apply Householder transformations to constraint matrix.
c
            IF (KRANK.GT.0 .AND. KRANK.LT.N) THEN
               DO 140 I = KRANK,1,-1
                  CALL xDH12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),             &
     &                      W(MAP1,1), MDW, 1, MG)
  140          CONTINUE
            ENDIF
c
c           Compute permuted inequality constraint matrix times r-inv.
c
            DO 160 I = MAP1,M
               DO 150 J = 1,KRANK
                W(I,J) = (W(I,J)-xDDOT(J-1,W(1,J),1,W(I,1),MDW))/W(J,J)
  150          CONTINUE
  160       CONTINUE
         ENDIF
c
c        Solve the reduced problem with DLPDP algorithm,
c        the least projected distance problem.
c
         CALL DLPDP(W(MAP1,1), MDW, MG, KRANK, N-KRANK, PRGOPT, X,                  &
     &             XNORM, MDLPDP, WS(N2), IP(N+1))
c
c        Compute solution in original coordinates.
c
         IF (MDLPDP.EQ.1) THEN
            DO 170 I = KRANK,1,-1
               X(I) = (X(I)-xDDOT(KRANK-I,W(I,I+1),MDW,X(I+1),1))/W(I,I)
  170       CONTINUE
c
c           Apply Householder transformation to solution vector.
c
            IF (KRANK.LT.N) THEN
               DO 180 I = 1,KRANK
                  CALL xDH12 (2, I, KRANK+1, N, W(I,1), MDW, WS(N1+I-1),             &
     &                      X, 1, 1, 1)
  180          CONTINUE
            ENDIF
c
c           Repermute variables to their input order.
c
            IF (MINMAN.GT.0) THEN
               DO 190 I = MINMAN,1,-1
                  CALL xDSWAP (1, X(I), 1, X(IP(I)), 1)
  190          CONTINUE
c
c              Variables are now in original coordinates.
c              Add solution of unconstrained problem.
c
               DO 200 I = 1,N
                  X(I) = X(I) + WS(I)
  200          CONTINUE
c
c              Compute the residual vector norm.
c
               RNORM = SQRT(RNORM**2+XNORM**2)
            ENDIF
         ELSE
            MODE = 2
         ENDIF
      ELSE
         CALL xDCOPY (N, WS, 1, X, 1)
      ENDIF
c
c     Compute covariance matrix based on the orthogonal decomposition
c     from xDHFTI( ).
c
      IF (.NOT.COV .OR. KRANK.LE.0) GO TO 370
      KRM1 = KRANK - 1
      KRP1 = KRANK + 1
c
c     Copy diagonal terms to working array.
c
      CALL xDCOPY (KRANK, W, MDW+1, WS(N2), 1)
c
c     Reciprocate diagonal terms.
c
      DO 210 J = 1,KRANK
         W(J,J) = 1.D0/W(J,J)
  210 CONTINUE
c
c     Invert the upper triangular QR factor on itself.
c
      IF (KRANK.GT.1) THEN
         DO 230 I = 1,KRM1
            DO 220 J = I+1,KRANK
               W(I,J) = -xDDOT(J-I,W(I,I),MDW,W(I,J),1)*W(J,J)
  220       CONTINUE
  230    CONTINUE
      ENDIF
c
c     Compute the inverted factor times its transpose.
c
      DO 250 I = 1,KRANK
         DO 240 J = I,KRANK
            W(I,J) = xDDOT(KRANK+1-J,W(I,J),MDW,W(J,J),MDW)
  240    CONTINUE
  250 CONTINUE
c
c     Zero out lower trapezoidal part.
c     Copy upper triangular to lower triangular part.
c
      IF (KRANK.LT.N) THEN
         DO 260 J = 1,KRANK
            CALL xDCOPY (J, W(1,J), 1, W(J,1), MDW)
  260    CONTINUE
c
         DO 270 I = KRP1,N
            CALL XDCOPYSC (I, 0.D0, W(I,1), MDW)
  270    CONTINUE
c
c        Apply right side transformations to lower triangle.
c
         N3 = N2 + KRP1
         DO 330 I = 1,KRANK
            L = N1 + I
            K = N2 + I
            RB = WS(L-1)*WS(K-1)
c
c           If RB.GE.0.D0, transformation can be regarded as zero.
c
            IF (RB.LT.0.D0) THEN
               RB = 1.D0/RB
c
c              Store unscaled rank one Householder update in work array.
c
               CALL XDCOPYSC (N, 0.D0, WS(N3), 1)
               L = N1 + I
               K = N3 + I
               WS(K-1) = WS(L-1)
c
               DO 280 J = KRP1,N
                  WS(N3+J-1) = W(I,J)
  280          CONTINUE
c
               DO 290 J = 1,N
                  WS(J) = RB*(xDDOT(J-I,W(J,I),MDW,WS(N3+I-1),1)+                    &
     &                    xDDOT(N-J+1,W(J,J),1,WS(N3+J-1),1))
  290          CONTINUE
c
               L = N3 + I
               GAM = 0.5D0*RB*xDDOT(N-I+1,WS(L-1),1,WS(I),1)
               CALL xDAXPY (N-I+1, GAM, WS(L-1), 1, WS(I), 1)
               DO 320 J = I,N
                  DO 300 L = 1,I-1
                     W(J,L) = W(J,L) + WS(N3+J-1)*WS(L)
  300             CONTINUE
c
                  DO 310 L = I,J
                     W(J,L) = W(J,L) + WS(J)*WS(N3+L-1)+WS(L)*WS(N3+J-1)
  310             CONTINUE
  320          CONTINUE
            ENDIF
  330    CONTINUE
c
c        Copy lower triangle to upper triangle to symmetrize the
c        covariance matrix.
c
         DO 340 I = 1,N
            CALL xDCOPY (I, W(I,1), MDW, W(1,I), 1)
  340    CONTINUE
      ENDIF
c
c     Repermute rows and columns.
c
      DO 350 I = MINMAN,1,-1
         K = IP(I)
         IF (I.NE.K) THEN
            CALL xDSWAP (1, W(I,I), 1, W(K,K), 1)
            CALL xDSWAP (I-1, W(1,I), 1, W(1,K), 1)
            CALL xDSWAP (K-I-1, W(I,I+1), MDW, W(I+1,K), 1)
            CALL xDSWAP (N-K, W(I, K+1), MDW, W(K, K+1), MDW)
         ENDIF
  350 CONTINUE
c
c     Put in normalized residual sum of squares scale factor
c     and symmetrize the resulting covariance matrix.
c
      DO 360 J = 1,N
         CALL xDSCAL (J, FAC, W(1,J), 1)
         CALL xDCOPY (J, W(1,J), 1, W(J,1), MDW)
  360 CONTINUE
c
  370 IP(1) = KRANK
      IP(2) = N + MAX(M,N) + (MG+2)*(N+7)
      RETURN
      END





      SUBROUTINE DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE,                  &
     &   RNORM, IDOPE, DOPE, DONE)
c***BEGIN PROLOGUE  DWNLIT
c***SUBSIDIARY
c***PURPOSE  Subsidiary to DWNNLS
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (WNLIT-S, DWNLIT-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     This is a companion subprogram to DWNNLS( ).
c     The documentation for DWNNLS( ) has complete usage instructions.
c
c     Note  The M by (N+1) matrix W( , ) contains the rt. hand side
c           B as the (N+1)st col.
c
c     Triangularize L1 by L1 subsystem, where L1=MIN(M,L), with
c     col interchanges.
c
c***SEE ALSO  DWNNLS
c***ROUTINES CALLED  xDCOPY, xDH12, xDROTM, xDROTMG, xDSCAL, xDSWAP, DWNLT1,
c                    DWNLT2, DWNLT3, xIDAMAX
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890618  Completely restructured and revised.  (WRB & RWC)
c   890620  Revised to make WNLT1, WNLT2, and WNLT3 subroutines.  (RWC)
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900328  Added TYPE section.  (WRB)
c   900604  DP version created from SP version. .  (RWC)
c***END PROLOGUE  DWNLIT
      INTEGER IDOPE(*), IPIVOT(*), ITYPE(*), L, M, MDW, N
      DOUBLE PRECISION DOPE(*), H(*), RNORM, SCALE(*), W(MDW,*)
      LOGICAL DONE
c
      EXTERNAL xDCOPY, xDH12, xDROTM, xDROTMG, xDSCAL, xDSWAP, DWNLT1,              &
     &   DWNLT2, DWNLT3, xIDAMAX
      INTEGER xIDAMAX
      LOGICAL DWNLT2
c
      DOUBLE PRECISION ALSQ, AMAX, EANORM, FACTOR, HBAR, RN, SPARAM(5),             &
     &   T, TAU
      INTEGER I, I1, IMAX, IR, J, J1, JJ, JP, KRANK, L1, LB, LEND, ME,              &
     &   MEND, NIV, NSOLN
      LOGICAL INDEP, RECALC
c
c***FIRST EXECUTABLE STATEMENT  DWNLIT
      ME    = IDOPE(1)
      NSOLN = IDOPE(2)
      L1    = IDOPE(3)
c
      ALSQ   = DOPE(1)
      EANORM = DOPE(2)
      TAU    = DOPE(3)
c
      LB     = MIN(M-1,L)
      RECALC = .TRUE.
      RNORM  = 0.D0
      KRANK  = 0
c
c     We set FACTOR=1.0 so that the heavy weight ALAMDA will be
c     included in the test for column independence.
c
      FACTOR = 1.D0
      LEND = L
      DO 180 I=1,LB
c
c        Set IR to point to the I-th row.
c
         IR = I
         MEND = M
         CALL DWNLT1 (I, LEND, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,            &
     &                W)
c
c        Update column SS and find pivot column.
c
         CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
c
c        Perform column interchange.
c        Test independence of incoming column.
c
  130    IF (DWNLT2(ME, MEND, IR, FACTOR, TAU, SCALE, W(1,I))) THEN
c
c           Eliminate I-th column below diagonal using modified Givens
c           transformations applied to (A B).
c
c           When operating near the ME line, use the largest element
c           above it as the pivot.
c
            DO 160 J=M,I+1,-1
               JP = J-1
               IF (J.EQ.ME+1) THEN
                  IMAX = ME
                  AMAX = SCALE(ME)*W(ME,I)**2
                  DO 150 JP=J-1,I,-1
                     T = SCALE(JP)*W(JP,I)**2
                     IF (T.GT.AMAX) THEN
                        IMAX = JP
                        AMAX = T
                     ENDIF
  150             CONTINUE
                  JP = IMAX
               ENDIF
c
               IF (W(J,I).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(JP), SCALE(J), W(JP,I), W(J,I),                &
     &                         SPARAM)
                  W(J,I) = 0.D0
                  CALL xDROTM (N+1-I, W(JP,I+1), MDW, W(J,I+1), MDW,                 &
     &                        SPARAM)
               ENDIF
  160       CONTINUE
         ELSE IF (LEND.GT.I) THEN
c
c           Column I is dependent.  Swap with column LEND.
c           Perform column interchange,
c           and find column in remaining set with largest SS.
c
            CALL DWNLT3 (I, LEND, M, MDW, IPIVOT, H, W)
            LEND = LEND - 1
            IMAX = xIDAMAX(LEND-I+1, H(I), 1) + I - 1
            HBAR = H(IMAX)
            GO TO 130
         ELSE
            KRANK = I - 1
            GO TO 190
         ENDIF
  180 CONTINUE
      KRANK = L1
c
  190 IF (KRANK.LT.ME) THEN
         FACTOR = ALSQ
         DO 200 I=KRANK+1,ME
            CALL XDCOPYSC (L, 0.D0, W(I,1), MDW)
  200    CONTINUE
c
c        Determine the rank of the remaining equality constraint
c        equations by eliminating within the block of constrained
c        variables.  Remove any redundant constraints.
c
         RECALC = .TRUE.
         LB = MIN(L+ME-KRANK, N)
         DO 270 I=L+1,LB
            IR = KRANK + I - L
            LEND = N
            MEND = ME
            CALL DWNLT1 (I, LEND, ME, IR, MDW, RECALC, IMAX, HBAR, H,               &
     &                   SCALE, W)
c
c           Update col ss and find pivot col
c
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
c
c           Perform column interchange
c           Eliminate elements in the I-th col.
c
            DO 240 J=ME,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL xDROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),               &
     &                        SPARAM)
                  W(J,I) = 0.D0
                  CALL xDROTM (N+1-I, W(J-1,I+1), MDW,W(J,I+1), MDW,                 &
     &                        SPARAM)
               ENDIF
  240       CONTINUE
c
c           I=column being eliminated.
c           Test independence of incoming column.
c           Remove any redundant or dependent equality constraints.
c
            IF (.NOT.DWNLT2(ME, MEND, IR, FACTOR,TAU,SCALE,W(1,I))) THEN
               JJ = IR
               DO 260 IR=JJ,ME
                  CALL XDCOPYSC (N, 0.D0, W(IR,1), MDW)
                  RNORM = RNORM + (SCALE(IR)*W(IR,N+1)/ALSQ)*W(IR,N+1)
                  W(IR,N+1) = 0.D0
                  SCALE(IR) = 1.D0
c
c                 Reclassify the zeroed row as a least squares equation.
c
                  ITYPE(IR) = 1
  260          CONTINUE
c
c              Reduce ME to reflect any discovered dependent equality
c              constraints.
c
               ME = JJ - 1
               GO TO 280
            ENDIF
  270    CONTINUE
      ENDIF
c
c     Try to determine the variables KRANK+1 through L1 from the
c     least squares equations.  Continue the triangularization with
c     pivot element W(ME+1,I).
c
  280 IF (KRANK.LT.L1) THEN
         RECALC = .TRUE.
c
c        Set FACTOR=ALSQ to remove effect of heavy weight from
c        test for column independence.
c
         FACTOR = ALSQ
         DO 350 I=KRANK+1,L1
c
c           Set IR to point to the ME+1-st row.
c
            IR = ME+1
            LEND = L
            MEND = M
            CALL DWNLT1 (I, L, M, IR, MDW, RECALC, IMAX, HBAR, H, SCALE,            &
     &                   W)
c
c           Update column SS and find pivot column.
c
            CALL DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
c
c           Perform column interchange.
c           Eliminate I-th column below the IR-th element.
c
            DO 320 J=M,IR+1,-1
               IF (W(J,I).NE.0.D0) THEN
                 CALL xDROTMG (SCALE(J-1), SCALE(J), W(J-1,I), W(J,I),               &
     &                        SPARAM)
                  W(J,I) = 0.D0
                  CALL xDROTM (N+1-I, W(J-1,I+1),  MDW, W(J,I+1), MDW,               &
     &                        SPARAM)
               ENDIF
  320       CONTINUE
c
c           Test if new pivot element is near zero.
c           If so, the column is dependent.
c           Then check row norm test to be classified as independent.
c
            T = SCALE(IR)*W(IR,I)**2
            INDEP = T .GT. (TAU*EANORM)**2
            IF (INDEP) THEN
               RN = 0.D0
               DO 340 I1=IR,M
                  DO 330 J1=I+1,N
                     RN = MAX(RN, SCALE(I1)*W(I1,J1)**2)
  330             CONTINUE
  340          CONTINUE
               INDEP = T .GT. RN*TAU**2
            ENDIF
c
c           If independent, swap the IR-th and KRANK+1-th rows to
c           maintain the triangular form.  Update the rank indicator
c           KRANK and the equality constraint pointer ME.
c
            IF (.NOT.INDEP) GO TO 360
            CALL xDSWAP(N+1, W(KRANK+1,1), MDW, W(IR,1), MDW)
            CALL xDSWAP(1, SCALE(KRANK+1), 1, SCALE(IR), 1)
c
c           Reclassify the least square equation as an equality
c           constraint and rescale it.
c
            ITYPE(IR) = 0
            T = SQRT(SCALE(KRANK+1))
            CALL xDSCAL(N+1, T, W(KRANK+1,1), MDW)
            SCALE(KRANK+1) = ALSQ
            ME = ME+1
            KRANK = KRANK+1
  350    CONTINUE
      ENDIF
c
c     If pseudorank is less than L, apply Householder transformation.
c     from right.
c
  360 IF (KRANK.LT.L) THEN
         DO 370 J=KRANK,1,-1
            CALL xDH12 (1, J, KRANK+1, L, W(J,1), MDW, H(J), W, MDW, 1,              &
     &                J-1)
  370    CONTINUE
      ENDIF
c
      NIV = KRANK + NSOLN - L
      IF (L.EQ.N) DONE = .TRUE.
c
c     End of initial triangularization.
c
      IDOPE(1) = ME
      IDOPE(2) = KRANK
      IDOPE(3) = NIV
      RETURN
      END





      SUBROUTINE DWNLSM (W, MDW, MME, MA, N, L, PRGOPT, X, RNORM, MODE,             &
     &   IPIVOT, ITYPE, WD, H, SCALE, Z, TEMP, D)
c***BEGIN PROLOGUE  DWNLSM
c***SUBSIDIARY
c***PURPOSE  Subsidiary to DWNNLS
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (WNLSM-S, DWNLSM-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     This is a companion subprogram to DWNNLS.
c     The documentation for DWNNLS has complete usage instructions.
c
c     In addition to the parameters discussed in the prologue to
c     subroutine DWNNLS, the following work arrays are used in
c     subroutine DWNLSM  (they are passed through the calling
c     sequence from DWNNLS for purposes of variable dimensioning).
c     Their contents will in general be of no interest to the user.
c
c     Variables of type REAL are DOUBLE PRECISION.
c
c         IPIVOT(*)
c            An array of length N.  Upon completion it contains the
c         pivoting information for the cols of W(*,*).
c
c         ITYPE(*)
c            An array of length M which is used to keep track
c         of the classification of the equations.  ITYPE(I)=0
c         denotes equation I as an equality constraint.
c         ITYPE(I)=1 denotes equation I as a least squares
c         equation.
c
c         WD(*)
c            An array of length N.  Upon completion it contains the
c         dual solution vector.
c
c         H(*)
c            An array of length N.  Upon completion it contains the
c         pivot scalars of the Householder transformations performed
c         in the case KRANK.LT.L.
c
c         SCALE(*)
c            An array of length M which is used by the subroutine
c         to store the diagonal matrix of weights.
c         These are used to apply the modified Givens
c         transformations.
c
c         Z(*),TEMP(*)
c            Working arrays of length N.
c
c         D(*)
c            An array of length N that contains the
c         column scaling for the matrix (E).
c                                       (A)
c
c***SEE ALSO  DWNNLS
c***ROUTINES CALLED  D1MACH, xDASUM, xDAXPY, xDCOPY, xDH12, xDNRM2, xDROTM,
c                    xDROTMG, xDSCAL, xDSWAP, DWNLIT, xIDAMAX, xXERMSG
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890618  Completely restructured and revised.  (WRB & RWC)
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
c   900328  Added TYPE section.  (WRB)
c   900510  Fixed an error message.  (RWC)
c   900604  DP version created from SP version.  (RWC)
c   900911  Restriction on value of ALAMDA included.  (WRB)
c***END PROLOGUE  DWNLSM
      INTEGER IPIVOT(*), ITYPE(*), L, MA, MDW, MME, MODE, N
      DOUBLE PRECISION D(*), H(*), PRGOPT(*), RNORM, SCALE(*), TEMP(*),             &
     &   W(MDW,*), WD(*), X(*), Z(*)
c
      EXTERNAL D1MACH,xDASUM,xDAXPY,xDCOPY,xDH12,xDNRM2,xDROTM,xDROTMG,             &
     &   xDSCAL, xDSWAP, DWNLIT, xIDAMAX, xXERMSG
      DOUBLE PRECISION D1MACH, xDASUM, xDNRM2
      INTEGER xIDAMAX
c
      DOUBLE PRECISION ALAMDA, ALPHA, ALSQ, AMAX, BLOWUP, BNORM,                    &
     &   DOPE(3), DRELPR, EANORM, FAC, SM, SPARAM(5), T, TAU, WMAX, Z2,             &
     &   ZZ
      INTEGER I, IDOPE(3), IMAX, ISOL, ITEMP, ITER, ITMAX, IWMAX, J,                &
     &   JCON, JP, KEY, KRANK, L1, LAST, LINK, M, ME, NEXT, NIV, NLINK,             &  
     &   NOPT, NSOLN, NTIMES
      LOGICAL DONE, FEASBL, FIRST, HITCON, POS
c
      SAVE DRELPR, FIRST
      DATA FIRST /.TRUE./
c***FIRST EXECUTABLE STATEMENT  DWNLSM
c
c     Initialize variables.
c     DRELPR is the precision for the particular machine
c     being used.  This logic avoids resetting it every entry.
c
      IF (FIRST) DRELPR = D1MACH(4)
      FIRST = .FALSE.
c
c     Set the nominal tolerance used in the code.
c
      TAU = SQRT(DRELPR)
c
      M = MA + MME
      ME = MME
      MODE = 2
c
c     To process option vector
c
      FAC = 1.D-4
c
c     Set the nominal blow up factor used in the code.
c
      BLOWUP = TAU
c
c     The nominal column scaling used in the code is
c     the identity scaling.
c
      CALL XDCOPYSC (N, 1.D0, D, 1)
c
c     Define bound for number of options to change.
c
      NOPT = 1000
c
c     Define bound for positive value of LINK.
c
      NLINK = 100000
      NTIMES = 0
      LAST = 1
      LINK = PRGOPT(1)
      IF (LINK.LE.0 .OR. LINK.GT.NLINK) THEN
         CALL xXERMSG ('SLATEC', 'DWNLSM',                                          &
     &      'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
         RETURN
      ENDIF
c
  100 IF (LINK.GT.1) THEN
         NTIMES = NTIMES + 1
         IF (NTIMES.GT.NOPT) THEN
         CALL xXERMSG ('SLATEC', 'DWNLSM',                                          &
     &      'IN DWNNLS, THE LINKS IN THE OPTION VECTOR ARE CYCLING.',               &
     &      3, 1)
            RETURN
         ENDIF
c
         KEY = PRGOPT(LAST+1)
         IF (KEY.EQ.6 .AND. PRGOPT(LAST+2).NE.0.D0) THEN
            DO 110 J = 1,N
               T = xDNRM2(M,W(1,J),1)
               IF (T.NE.0.D0) T = 1.D0/T
               D(J) = T
  110       CONTINUE
         ENDIF
c
         IF (KEY.EQ.7) CALL xDCOPY (N, PRGOPT(LAST+2), 1, D, 1)
         IF (KEY.EQ.8) TAU = MAX(DRELPR,PRGOPT(LAST+2))
         IF (KEY.EQ.9) BLOWUP = MAX(DRELPR,PRGOPT(LAST+2))
c
         NEXT = PRGOPT(LINK)
         IF (NEXT.LE.0 .OR. NEXT.GT.NLINK) THEN
            CALL xXERMSG ('SLATEC', 'DWNLSM',                                       &
     &         'IN DWNNLS, THE OPTION VECTOR IS UNDEFINED', 3, 1)
            RETURN
         ENDIF
c
         LAST = LINK
         LINK = NEXT
         GO TO 100
      ENDIF
c
      DO 120 J = 1,N
         CALL xDSCAL (M, D(J), W(1,J), 1)
  120 CONTINUE
c
c     Process option vector
c
      DONE = .FALSE.
      ITER = 0
      ITMAX = 3*(N-L)
      MODE = 0
      NSOLN = L
      L1 = MIN(M,L)
c
c     Compute scale factor to apply to equality constraint equations.
c
      DO 130 J = 1,N
         WD(J) = xDASUM(M,W(1,J),1)
  130 CONTINUE
c
      IMAX = xIDAMAX(N,WD,1)
      EANORM = WD(IMAX)
      BNORM = xDASUM(M,W(1,N+1),1)
      ALAMDA = EANORM/(DRELPR*FAC)
c
c     On machines, such as the VAXes using D floating, with a very
c     limited exponent range for double precision values, the previously
c     computed value of ALAMDA may cause an overflow condition.
c     Therefore, this code further limits the value of ALAMDA.
c
      ALAMDA = MIN(ALAMDA,SQRT(D1MACH(2)))
c
c     Define scaling diagonal matrix for modified Givens usage and
c     classify equation types.
c
      ALSQ = ALAMDA**2
      DO 140 I = 1,M
c
c        When equation I is heavily weighted ITYPE(I)=0,
c        else ITYPE(I)=1.
c
         IF (I.LE.ME) THEN
            T = ALSQ
            ITEMP = 0
         ELSE
            T = 1.D0
            ITEMP = 1
         ENDIF
         SCALE(I) = T
         ITYPE(I) = ITEMP
  140 CONTINUE
c
c     Set the solution vector X(*) to zero and the column interchange
c     matrix to the identity.
c
      CALL XDCOPYSC (N, 0.D0, X, 1)
      DO 150 I = 1,N
         IPIVOT(I) = I
  150 CONTINUE
c
c     Perform initial triangularization in the submatrix
c     corresponding to the unconstrained variables.
c     Set first L components of dual vector to zero because
c     these correspond to the unconstrained variables.
c
      CALL XDCOPYSC (L, 0.D0, WD, 1)
c
c     The arrays IDOPE(*) and DOPE(*) are used to pass
c     information to DWNLIT().  This was done to avoid
c     a long calling sequence or the use of COMMON.
c
      IDOPE(1) = ME
      IDOPE(2) = NSOLN
      IDOPE(3) = L1
c
      DOPE(1) = ALSQ
      DOPE(2) = EANORM
      DOPE(3) = TAU
      CALL DWNLIT (W, MDW, M, N, L, IPIVOT, ITYPE, H, SCALE, RNORM,                 &
     &            IDOPE, DOPE, DONE)
      ME    = IDOPE(1)
      KRANK = IDOPE(2)
      NIV   = IDOPE(3)
c
c     Perform WNNLS algorithm using the following steps.
c
c     Until(DONE)
c        compute search direction and feasible point
c        when (HITCON) add constraints
c        else perform multiplier test and drop a constraint
c        fin
c     Compute-Final-Solution
c
c     To compute search direction and feasible point,
c     solve the triangular system of currently non-active
c     variables and store the solution in Z(*).
c
c     To solve system
c     Copy right hand side into TEMP vector to use overwriting method.
c
  160 IF (DONE) GO TO 330
      ISOL = L + 1
      IF (NSOLN.GE.ISOL) THEN
         CALL xDCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 170 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
c
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL xDAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  170    CONTINUE
      ENDIF
c
c     Increment iteration counter and check against maximum number
c     of iterations.
c
      ITER = ITER + 1
      IF (ITER.GT.ITMAX) THEN
         MODE = 1
         DONE = .TRUE.
      ENDIF
c
c     Check to see if any constraints have become active.
c     If so, calculate an interpolation factor so that all
c     active constraints are removed from the basis.
c
      ALPHA = 2.D0
      HITCON = .FALSE.
      DO 180 J = L+1,NSOLN
         ZZ = Z(J)
         IF (ZZ.LE.0.D0) THEN
            T = X(J)/(X(J)-ZZ)
            IF (T.LT.ALPHA) THEN
               ALPHA = T
               JCON = J
            ENDIF
            HITCON = .TRUE.
         ENDIF
  180 CONTINUE
c
c     Compute search direction and feasible point
c
      IF (HITCON) THEN
c
c        To add constraints, use computed ALPHA to interpolate between
c        last feasible solution X(*) and current unconstrained (and
c        infeasible) solution Z(*).
c
         DO 190 J = L+1,NSOLN
            X(J) = X(J) + ALPHA*(Z(J)-X(J))
  190    CONTINUE
         FEASBL = .FALSE.
c
c        Remove column JCON and shift columns JCON+1 through N to the
c        left.  Swap column JCON into the N th position.  This achieves
c        upper Hessenberg form for the nonactive constraints and
c        leaves an upper Hessenberg matrix to retriangularize.
c
  200    DO 210 I = 1,M
            T = W(I,JCON)
            CALL xDCOPY (N-JCON, W(I, JCON+1), MDW, W(I, JCON), MDW)
            W(I,N) = T
  210    CONTINUE
c
c        Update permuted index vector to reflect this shift and swap.
c
         ITEMP = IPIVOT(JCON)
         DO 220 I = JCON,N - 1
            IPIVOT(I) = IPIVOT(I+1)
  220    CONTINUE
         IPIVOT(N) = ITEMP
c
c        Similarly permute X(*) vector.
c
         CALL xDCOPY (N-JCON, X(JCON+1), 1, X(JCON), 1)
         X(N) = 0.D0
         NSOLN = NSOLN - 1
         NIV = NIV - 1
c
c        Retriangularize upper Hessenberg matrix after adding
c        constraints.
c
         I = KRANK + JCON - L
         DO 230 J = JCON,NSOLN
            IF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.0) THEN
c
c              Zero IP1 to I in column J
c
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),              &
     &                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,                &
     &                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.1) THEN
c
c              Zero IP1 to I in column J
c
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),              &
     &                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,                &
     &                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.1 .AND. ITYPE(I+1).EQ.0) THEN
               CALL xDSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
               CALL xDSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
               ITEMP = ITYPE(I+1)
               ITYPE(I+1) = ITYPE(I)
               ITYPE(I) = ITEMP
c
c              Swapped row was formerly a pivot element, so it will
c              be large enough to perform elimination.
c              Zero IP1 to I in column J.
c
               IF (W(I+1,J).NE.0.D0) THEN
                  CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J), W(I+1,J),              &
     &                         SPARAM)
                  W(I+1,J) = 0.D0
                  CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,                &
     &                        SPARAM)
               ENDIF
            ELSEIF (ITYPE(I).EQ.0 .AND. ITYPE(I+1).EQ.1) THEN
               IF (SCALE(I)*W(I,J)**2/ALSQ.GT.(TAU*EANORM)**2) THEN
c
c                 Zero IP1 to I in column J
c
                  IF (W(I+1,J).NE.0.D0) THEN
                     CALL xDROTMG (SCALE(I), SCALE(I+1), W(I,J),                     &
     &                            W(I+1,J), SPARAM)
                     W(I+1,J) = 0.D0
                     CALL xDROTM (N+1-J, W(I,J+1), MDW, W(I+1,J+1), MDW,             &
     &                           SPARAM)
                  ENDIF
               ELSE
                  CALL xDSWAP (N+1, W(I,1), MDW, W(I+1,1), MDW)
                  CALL xDSWAP (1, SCALE(I), 1, SCALE(I+1), 1)
                  ITEMP = ITYPE(I+1)
                  ITYPE(I+1) = ITYPE(I)
                  ITYPE(I) = ITEMP
                  W(I+1,J) = 0.D0
               ENDIF
            ENDIF
            I = I + 1
  230    CONTINUE
c
c        See if the remaining coefficients in the solution set are
c        feasible.  They should be because of the way ALPHA was
c        determined.  If any are infeasible, it is due to roundoff
c        error.  Any that are non-positive will be set to zero and
c        removed from the solution set.
c
         DO 240 JCON = L+1,NSOLN
            IF (X(JCON).LE.0.D0) GO TO 250
  240    CONTINUE
         FEASBL = .TRUE.
  250    IF (.NOT.FEASBL) GO TO 200
      ELSE
c
c        To perform multiplier test and drop a constraint.
c
         CALL xDCOPY (NSOLN, Z, 1, X, 1)
         IF (NSOLN.LT.N) CALL XDCOPYSC (N-NSOLN, 0.D0, X(NSOLN+1), 1)
c
c        Reclassify least squares equations as equalities as necessary.
c
         I = NIV + 1
  260    IF (I.LE.ME) THEN
            IF (ITYPE(I).EQ.0) THEN
               I = I + 1
            ELSE
               CALL xDSWAP (N+1, W(I,1), MDW, W(ME,1), MDW)
               CALL xDSWAP (1, SCALE(I), 1, SCALE(ME), 1)
               ITEMP = ITYPE(I)
               ITYPE(I) = ITYPE(ME)
               ITYPE(ME) = ITEMP
               ME = ME - 1
            ENDIF
            GO TO 260
         ENDIF
c
c        Form inner product vector WD(*) of dual coefficients.
c
         DO 280 J = NSOLN+1,N
            SM = 0.D0
            DO 270 I = NSOLN+1,M
               SM = SM + SCALE(I)*W(I,J)*W(I,N+1)
  270       CONTINUE
            WD(J) = SM
  280    CONTINUE
c
c        Find J such that WD(J)=WMAX is maximum.  This determines
c        that the incoming column J will reduce the residual vector
c        and be positive.
c
  290    WMAX = 0.D0
         IWMAX = NSOLN + 1
         DO 300 J = NSOLN+1,N
            IF (WD(J).GT.WMAX) THEN
               WMAX = WD(J)
               IWMAX = J
            ENDIF
  300    CONTINUE
         IF (WMAX.LE.0.D0) GO TO 330
c
c        Set dual coefficients to zero for incoming column.
c
         WD(IWMAX) = 0.D0
c
c        WMAX .GT. 0.D0, so okay to move column IWMAX to solution set.
c        Perform transformation to retriangularize, and test for near
c        linear dependence.
c
c        Swap column IWMAX into NSOLN-th position to maintain upper
c        Hessenberg form of adjacent columns, and add new column to
c        triangular decomposition.
c
         NSOLN = NSOLN + 1
         NIV = NIV + 1
         IF (NSOLN.NE.IWMAX) THEN
            CALL xDSWAP (M, W(1,NSOLN), 1, W(1,IWMAX), 1)
            WD(IWMAX) = WD(NSOLN)
            WD(NSOLN) = 0.D0
            ITEMP = IPIVOT(NSOLN)
            IPIVOT(NSOLN) = IPIVOT(IWMAX)
            IPIVOT(IWMAX) = ITEMP
         ENDIF
c
c        Reduce column NSOLN so that the matrix of nonactive constraints
c        variables is triangular.
c
         DO 320 J = M,NIV+1,-1
            JP = J - 1
c
c           When operating near the ME line, test to see if the pivot
c           element is near zero.  If so, use the largest element above
c           it as the pivot.  This is to maintain the sharp interface
c           between weighted and non-weighted rows in all cases.
c
            IF (J.EQ.ME+1) THEN
               IMAX = ME
               AMAX = SCALE(ME)*W(ME,NSOLN)**2
               DO 310 JP = J - 1,NIV,-1
                  T = SCALE(JP)*W(JP,NSOLN)**2
                  IF (T.GT.AMAX) THEN
                     IMAX = JP
                     AMAX = T
                  ENDIF
  310          CONTINUE
               JP = IMAX
            ENDIF
c
            IF (W(J,NSOLN).NE.0.D0) THEN
               CALL xDROTMG (SCALE(JP), SCALE(J), W(JP,NSOLN),                       &
     &                      W(J,NSOLN), SPARAM)
               W(J,NSOLN) = 0.D0
               CALL xDROTM (N+1-NSOLN, W(JP,NSOLN+1), MDW, W(J,NSOLN+1),             &
     &                     MDW, SPARAM)
            ENDIF
  320    CONTINUE
c
c        Solve for Z(NSOLN)=proposed new value for X(NSOLN).  Test if
c        this is nonpositive or too large.  If this was true or if the
c        pivot term was zero, reject the column as dependent.
c
         IF (W(NIV,NSOLN).NE.0.D0) THEN
            ISOL = NIV
            Z2 = W(ISOL,N+1)/W(ISOL,NSOLN)
            Z(NSOLN) = Z2
            POS = Z2 .GT. 0.D0
            IF (Z2*EANORM.GE.BNORM .AND. POS) THEN
               POS = .NOT. (BLOWUP*Z2*EANORM.GE.BNORM)
            ENDIF
c
c           Try to add row ME+1 as an additional equality constraint.
c           Check size of proposed new solution component.
c           Reject it if it is too large.
c
         ELSEIF (NIV.LE.ME .AND. W(ME+1,NSOLN).NE.0.D0) THEN
            ISOL = ME + 1
            IF (POS) THEN
c
c              Swap rows ME+1 and NIV, and scale factors for these rows.
c
               CALL xDSWAP (N+1, W(ME+1,1), MDW, W(NIV,1), MDW)
               CALL xDSWAP (1, SCALE(ME+1), 1, SCALE(NIV), 1)
               ITEMP = ITYPE(ME+1)
               ITYPE(ME+1) = ITYPE(NIV)
               ITYPE(NIV) = ITEMP
               ME = ME + 1
            ENDIF
         ELSE
            POS = .FALSE.
         ENDIF
c
         IF (.NOT.POS) THEN
            NSOLN = NSOLN - 1
            NIV = NIV - 1
         ENDIF
         IF (.NOT.(POS.OR.DONE)) GO TO 290
      ENDIF
      GO TO 160
c
c     Else perform multiplier test and drop a constraint.  To compute
c     final solution.  Solve system, store results in X(*).
c
c     Copy right hand side into TEMP vector to use overwriting method.
c
  330 ISOL = 1
      IF (NSOLN.GE.ISOL) THEN
         CALL xDCOPY (NIV, W(1,N+1), 1, TEMP, 1)
         DO 340 J = NSOLN,ISOL,-1
            IF (J.GT.KRANK) THEN
               I = NIV - NSOLN + J
            ELSE
               I = J
            ENDIF
c
            IF (J.GT.KRANK .AND. J.LE.L) THEN
               Z(J) = 0.D0
            ELSE
               Z(J) = TEMP(I)/W(I,J)
               CALL xDAXPY (I-1, -Z(J), W(1,J), 1, TEMP, 1)
            ENDIF
  340    CONTINUE
      ENDIF
c
c     Solve system.
c
      CALL xDCOPY (NSOLN, Z, 1, X, 1)
c
c     Apply Householder transformations to X(*) if KRANK.LT.L
c
      IF (KRANK.LT.L) THEN
         DO 350 I = 1,KRANK
            CALL xDH12 (2, I, KRANK+1, L, W(I,1), MDW, H(I), X, 1, 1, 1)
  350    CONTINUE
      ENDIF
c
c     Fill in trailing zeroes for constrained variables not in solution.
c
      IF (NSOLN.LT.N) CALL XDCOPYSC (N-NSOLN, 0.D0, X(NSOLN+1), 1)
c
c     Permute solution vector to natural order.
c
      DO 380 I = 1,N
         J = I
  360    IF (IPIVOT(J).EQ.I) GO TO 370
         J = J + 1
         GO TO 360
c
  370    IPIVOT(J) = IPIVOT(I)
         IPIVOT(I) = J
         CALL xDSWAP (1, X(J), 1, X(I), 1)
  380 CONTINUE
c
c     Rescale the solution using the column scaling.
c
      DO 390 J = 1,N
         X(J) = X(J)*D(J)
  390 CONTINUE
c
      DO 400 I = NSOLN+1,M
         T = W(I,N+1)
         IF (I.LE.ME) T = T/ALAMDA
         T = (SCALE(I)*T)*T
         RNORM = RNORM + T
  400 CONTINUE
c
      RNORM = SQRT(RNORM)
      RETURN
      END





      SUBROUTINE DWNLT1 (I, LEND, MEND, IR, MDW, RECALC, IMAX, HBAR, H,             &
     &   SCALE, W)
c***BEGIN PROLOGUE  DWNLT1
c***SUBSIDIARY
c***PURPOSE  Subsidiary to WNLIT
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (WNLT1-S, DWNLT1-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     To update the column Sum Of Squares and find the pivot column.
c     The column Sum of Squares Vector will be updated at each step.
c     When numerically necessary, these values will be recomputed.
c
c***SEE ALSO  DWNLIT
c***ROUTINES CALLED  xIDAMAX
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
c   900604  DP version created from SP version.  (RWC)
c***END PROLOGUE  DWNLT1
      INTEGER I, IMAX, IR, LEND, MDW, MEND
      DOUBLE PRECISION H(*), HBAR, SCALE(*), W(MDW,*)
      LOGICAL RECALC
c
      EXTERNAL xIDAMAX
      INTEGER xIDAMAX
c
      INTEGER J, K
c
c***FIRST EXECUTABLE STATEMENT  DWNLT1
      IF (IR.NE.1 .AND. (.NOT.RECALC)) THEN
c
c        Update column SS=sum of squares.
c
         DO 10 J=I,LEND
            H(J) = H(J) - SCALE(IR-1)*W(IR-1,J)**2
   10    CONTINUE
c
c        Test for numerical accuracy.
c
         IMAX = xIDAMAX(LEND-I+1, H(I), 1) + I - 1
         RECALC = (HBAR+1.E-3*H(IMAX)) .EQ. HBAR
      ENDIF
c
c     If required, recalculate column SS, using rows IR through MEND.
c
      IF (RECALC) THEN
         DO 30 J=I,LEND
            H(J) = 0.D0
            DO 20 K=IR,MEND
               H(J) = H(J) + SCALE(K)*W(K,J)**2
   20       CONTINUE
   30    CONTINUE
c
c        Find column with largest SS.
c
         IMAX = xIDAMAX(LEND-I+1, H(I), 1) + I - 1
         HBAR = H(IMAX)
      ENDIF
      RETURN
      END





      LOGICAL FUNCTION DWNLT2 (ME, MEND, IR, FACTOR, TAU, SCALE, WIC)
c***BEGIN PROLOGUE  DWNLT2
c***SUBSIDIARY
c***PURPOSE  Subsidiary to WNLIT
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (WNLT2-S, DWNLT2-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     To test independence of incoming column.
c
c     Test the column IC to determine if it is linearly independent
c     of the columns already in the basis.  In the initial tri. step,
c     we usually want the heavy weight ALAMDA to be included in the
c     test for independence.  In this case, the value of FACTOR will
c     have been set to 1.E0 before this procedure is invoked.
c     In the potentially rank deficient problem, the value of FACTOR
c     will have been set to ALSQ=ALAMDA**2 to remove the effect of the
c     heavy weight from the test for independence.
c
c     Write new column as partitioned vector
c           (A1)  number of components in solution so far = NIV
c           (A2)  M-NIV components
c     And compute  SN = inverse weighted length of A1
c                  RN = inverse weighted length of A2
c     Call the column independent when RN .GT. TAU*SN
c
c***SEE ALSO  DWNLIT
c***ROUTINES CALLED  (NONE)
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
c   900604  DP version created from SP version.  (RWC)
c***END PROLOGUE  DWNLT2
      DOUBLE PRECISION FACTOR, SCALE(*), TAU, WIC(*)
      INTEGER IR, ME, MEND
c
      DOUBLE PRECISION RN, SN, T
      INTEGER J
c
c***FIRST EXECUTABLE STATEMENT  DWNLT2
      SN = 0.E0
      RN = 0.E0
      DO 10 J=1,MEND
         T = SCALE(J)
         IF (J.LE.ME) T = T/FACTOR
         T = T*WIC(J)**2
c
         IF (J.LT.IR) THEN
            SN = SN + T
         ELSE
            RN = RN + T
         ENDIF
   10 CONTINUE
      DWNLT2 = RN .GT. SN*TAU**2
      RETURN
      END





      SUBROUTINE DWNLT3 (I, IMAX, M, MDW, IPIVOT, H, W)
c***BEGIN PROLOGUE  DWNLT3
c***SUBSIDIARY
c***PURPOSE  Subsidiary to WNLIT
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (WNLT3-S, DWNLT3-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     Perform column interchange.
c     Exchange elements of permuted index vector and perform column
c     interchanges.
c
c***SEE ALSO  DWNLIT
c***ROUTINES CALLED  xDSWAP
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890620  Code extracted from WNLIT and made a subroutine.  (RWC))
c   900604  DP version created from SP version.  (RWC)
c***END PROLOGUE  DWNLT3
      INTEGER I, IMAX, IPIVOT(*), M, MDW
      DOUBLE PRECISION H(*), W(MDW,*)
c
      EXTERNAL xDSWAP
c
      DOUBLE PRECISION T
      INTEGER ITEMP
c
c***FIRST EXECUTABLE STATEMENT  DWNLT3
      IF (IMAX.NE.I) THEN
         ITEMP        = IPIVOT(I)
         IPIVOT(I)    = IPIVOT(IMAX)
         IPIVOT(IMAX) = ITEMP
c
         CALL xDSWAP(M, W(1,IMAX), 1, W(1,I), 1)
c
         T       = H(IMAX)
         H(IMAX) = H(I)
         H(I)    = T
      ENDIF
      RETURN
      END





      SUBROUTINE DWNNLS (W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE,              &
     &   IWORK, WORK)
c***BEGIN PROLOGUE  DWNNLS
c***PURPOSE  Solve a linearly constrained least squares problem with
c            equality constraints and nonnegativity constraints on
c            selected variables.
c***LIBRARY   SLATEC
c***CATEGORY  K1A2A
c***TYPE      DOUBLE PRECISION (WNNLS-S, DWNNLS-D)
c***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
c             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
c             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c     Abstract
c
c     This subprogram solves a linearly constrained least squares
c     problem.  Suppose there are given matrices E and A of
c     respective dimensions ME by N and MA by N, and vectors F
c     and B of respective lengths ME and MA.  This subroutine
c     solves the problem
c
c               EX = F, (equations to be exactly satisfied)
c
c               AX = B, (equations to be approximately satisfied,
c                        in the least squares sense)
c
c               subject to components L+1,...,N nonnegative
c
c     Any values ME.GE.0, MA.GE.0 and 0.LE. L .LE.N are permitted.
c
c     The problem is reposed as problem DWNNLS
c
c               (WT*E)X = (WT*F)
c               (   A)    (   B), (least squares)
c               subject to components L+1,...,N nonnegative.
c
c     The subprogram chooses the heavy weight (or penalty parameter) WT.
c
c     The parameters for DWNNLS are
c
c     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
c
c     W(*,*),MDW,  The array W(*,*) is double subscripted with first
c     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
c                  discussion let us call M = ME + MA.  Then MDW
c                  must satisfy MDW.GE.M.  The condition MDW.LT.M
c                  is an error.
c
c                  The array W(*,*) contains the matrices and vectors
c
c                       (E  F)
c                       (A  B)
c
c                  in rows and columns 1,...,M and 1,...,N+1
c                  respectively.  Columns 1,...,L correspond to
c                  unconstrained variables X(1),...,X(L).  The
c                  remaining variables are constrained to be
c                  nonnegative. The condition L.LT.0 or L.GT.N is
c                  an error.
c
c     PRGOPT(*)    This double precision array is the option vector.
c                  If the user is satisfied with the nominal
c                  subprogram features set
c
c                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
c
c                  Otherwise PRGOPT(*) is a linked list consisting of
c                  groups of data of the following form
c
c                  LINK
c                  KEY
c                  DATA SET
c
c                  The parameters LINK and KEY are each one word.
c                  The DATA SET can be comprised of several words.
c                  The number of items depends on the value of KEY.
c                  The value of LINK points to the first
c                  entry of the next group of data within
c                  PRGOPT(*).  The exception is when there are
c                  no more options to change.  In that
c                  case LINK=1 and the values KEY and DATA SET
c                  are not referenced. The general layout of
c                  PRGOPT(*) is as follows.
c
c               ...PRGOPT(1)=LINK1 (link to first entry of next group)
c               .  PRGOPT(2)=KEY1 (key to the option change)
c               .  PRGOPT(3)=DATA VALUE (data value for this change)
c               .       .
c               .       .
c               .       .
c               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
c               .                       next group)
c               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
c               .  PRGOPT(LINK1+2)=DATA VALUE
c               ...     .
c               .       .
c               .       .
c               ...PRGOPT(LINK)=1 (no more options to change)
c
c                  Values of LINK that are nonpositive are errors.
c                  A value of LINK.GT.NLINK=100000 is also an error.
c                  This helps prevent using invalid but positive
c                  values of LINK that will probably extend
c                  beyond the program limits of PRGOPT(*).
c                  Unrecognized values of KEY are ignored.  The
c                  order of the options is arbitrary and any number
c                  of options can be changed with the following
c                  restriction.  To prevent cycling in the
c                  processing of the option array a count of the
c                  number of options changed is maintained.
c                  Whenever this count exceeds NOPT=1000 an error
c                  message is printed and the subprogram returns.
c
c                  OPTIONS..
c
c                  KEY=6
c                         Scale the nonzero columns of the
c                  entire data matrix
c                  (E)
c                  (A)
c                  to have length one. The DATA SET for
c                  this option is a single value.  It must
c                  be nonzero if unit length column scaling is
c                  desired.
c
c                  KEY=7
c                         Scale columns of the entire data matrix
c                  (E)
c                  (A)
c                  with a user-provided diagonal matrix.
c                  The DATA SET for this option consists
c                  of the N diagonal scaling factors, one for
c                  each matrix column.
c
c                  KEY=8
c                         Change the rank determination tolerance from
c                  the nominal value of SQRT(SRELPR).  This quantity
c                  can be no smaller than SRELPR, The arithmetic-
c                  storage precision.  The quantity used
c                  here is internally restricted to be at
c                  least SRELPR.  The DATA SET for this option
c                  is the new tolerance.
c
c                  KEY=9
c                         Change the blow-up parameter from the
c                  nominal value of SQRT(SRELPR).  The reciprocal of
c                  this parameter is used in rejecting solution
c                  components as too large when a variable is
c                  first brought into the active set.  Too large
c                  means that the proposed component times the
c                  reciprocal of the parameter is not less than
c                  the ratio of the norms of the right-side
c                  vector and the data matrix.
c                  This parameter can be no smaller than SRELPR,
c                  the arithmetic-storage precision.
c
c                  For example, suppose we want to provide
c                  a diagonal matrix to scale the problem
c                  matrix and change the tolerance used for
c                  determining linear dependence of dropped col
c                  vectors.  For these options the dimensions of
c                  PRGOPT(*) must be at least N+6.  The FORTRAN
c                  statements defining these options would
c                  be as follows.
c
c                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
c                  PRGOPT(2)=7 (user-provided scaling key)
c
c                  CALL xDCOPY(N,D,1,PRGOPT(3),1) (copy the N
c                  scaling factors from a user array called D(*)
c                  into PRGOPT(3)-PRGOPT(N+2))
c
c                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
c                  PRGOPT(N+4)=8 (linear dependence tolerance key)
c                  PRGOPT(N+5)=... (new value of the tolerance)
c
c                  PRGOPT(N+6)=1 (no more options to change)
c
c
c     IWORK(1),    The amounts of working storage actually allocated
c     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
c                  respectively.  These quantities are compared with
c                  the actual amounts of storage needed for DWNNLS( ).
c                  Insufficient storage allocated for either WORK(*)
c                  or IWORK(*) is considered an error.  This feature
c                  was included in DWNNLS( ) because miscalculating
c                  the storage formulas for WORK(*) and IWORK(*)
c                  might very well lead to subtle and hard-to-find
c                  execution errors.
c
c                  The length of WORK(*) must be at least
c
c                  LW = ME+MA+5*N
c                  This test will not be made if IWORK(1).LE.0.
c
c                  The length of IWORK(*) must be at least
c
c                  LIW = ME+MA+N
c                  This test will not be made if IWORK(2).LE.0.
c
c     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
c
c     X(*)         An array dimensioned at least N, which will
c                  contain the N components of the solution vector
c                  on output.
c
c     RNORM        The residual norm of the solution.  The value of
c                  RNORM contains the residual vector length of the
c                  equality constraints and least squares equations.
c
c     MODE         The value of MODE indicates the success or failure
c                  of the subprogram.
c
c                  MODE = 0  Subprogram completed successfully.
c
c                       = 1  Max. number of iterations (equal to
c                            3*(N-L)) exceeded. Nearly all problems
c                            should complete in fewer than this
c                            number of iterations. An approximate
c                            solution and its corresponding residual
c                            vector length are in X(*) and RNORM.
c
c                       = 2  Usage error occurred.  The offending
c                            condition is noted with the error
c                            processing subprogram, xXERMSG( ).
c
c     User-designated
c     Working arrays..
c
c     WORK(*)      A double precision working array of length at least
c                  M + 5*N.
c
c     IWORK(*)     An integer-valued working array of length at least
c                  M+N.
c
c***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
c                 linear least squares problems with equality and
c                 nonnegativity constraints, Report SAND77-0552, Sandia
c                 Laboratories, June 1978.
c               K. H. Haskell and R. J. Hanson, Selected algorithms for
c                 the linearly constrained least squares problem - a
c                 users guide, Report SAND78-1290, Sandia Laboratories,
c                 August 1979.
c               K. H. Haskell and R. J. Hanson, An algorithm for
c                 linear least squares problems with equality and
c                 nonnegativity constraints, Mathematical Programming
c                 21 (1981), pp. 98-118.
c               R. J. Hanson and K. H. Haskell, Two algorithms for the
c                 linearly constrained least squares problem, ACM
c                 Transactions on Mathematical Software, September 1982.
c               c. L. Lawson and R. J. Hanson, Solving Least Squares
c                 Problems, Prentice-Hall, Inc., 1974.
c***ROUTINES CALLED  DWNLSM, xXERMSG
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890618  Completely restructured and revised.  (WRB & RWC)
c   891006  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
c   900510  Convert XERRWV calls to xXERMSG calls, change Prologue
c           comments to agree with WNNLS.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  DWNNLS
      INTEGER IWORK(*), L, L1, L2, L3, L4, L5, LIW, LW, MA, MDW, ME,                &
     &     MODE, N
      DOUBLE PRECISION  PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
      CHARACTER*8 XERN1
c***FIRST EXECUTABLE STATEMENT  DWNNLS
      MODE = 0
      IF (MA+ME.LE.0 .OR. N.LE.0) RETURN
c
      IF (IWORK(1).GT.0) THEN
         LW = ME + MA + 5*N
         IF (IWORK(1).LT.LW) THEN
C KARLINE: REMOVED WRITE		 
         CALL rwarn ('LSEI: insufficient storage')
		 
C            WRITE (XERN1, '(I8)') LW
C            CALL xXERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //            &
C     &         'ALLOCATED FOR WORK(*), NEED LW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
c
      IF (IWORK(2).GT.0) THEN
         LIW = ME + MA + N
         IF (IWORK(2).LT.LIW) THEN
C KARLINE: REMOVED WRITE		 
         CALL rwarn ('LSEI: insufficient storage')
		 
C            WRITE (XERN1, '(I8)') LIW
C            CALL xXERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' //            &
C     &         'ALLOCATED FOR IWORK(*), NEED LIW = ' // XERN1, 2, 1)
            MODE = 2
            RETURN
         ENDIF
      ENDIF
c
      IF (MDW.LT.ME+MA) THEN
         CALL xXERMSG ('SLATEC', 'DWNNLS',                                          &
     &      'THE VALUE MDW.LT.ME+MA IS AN ERROR', 2, 1)
         MODE = 2
         RETURN
      ENDIF
c
      IF (L.LT.0 .OR. L.GT.N) THEN
         CALL xXERMSG ('SLATEC', 'DWNNLS',                                          &
     &      'L.GE.0 .AND. L.LE.N IS REQUIRED', 2, 1)
         MODE = 2
         RETURN
      ENDIF
c
c     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
c     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
c     REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).
c
      L1 = N + 1
      L2 = L1 + N
      L3 = L2 + ME + MA
      L4 = L3 + N
      L5 = L4 + N
c
      CALL DWNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK,              &
     &            IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3),                 &
     &            WORK(L4), WORK(L5))
      RETURN
      END



      SUBROUTINE xDROTM (N,DX,INCX,DY,INCY,DPARAM)
c
c     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
c
c     (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
c     (DY**T)
c
c     DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
c     LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
c     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
c
c     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
c
c       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
c     H=(          )    (          )    (          )    (          )
c       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
c     SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
c
      DOUBLE PRECISION DFLAG,DH12,DH22,DX,TWO,Z,DH11,DH21,                          &
     & DPARAM,DY,W,ZERO
      INTEGER N,INCX,INCY,NSTEPS,I,KX,KY
      DIMENSION DX(*),DY(*),DPARAM(5)
      DATA ZERO,TWO/0.D0,2.D0/
c
      DFLAG=DPARAM(1)
      IF(N .LE. 0 .OR.(DFLAG+TWO.EQ.ZERO)) GO TO 140
          IF(.NOT.(INCX.EQ.INCY.AND. INCX .GT.0)) GO TO 70
c
               NSTEPS=N*INCX
               IF(DFLAG) 50,10,30
   10          CONTINUE
               DH12=DPARAM(4)
               DH21=DPARAM(3)
                    DO 20 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W+Z*DH12
                    DY(I)=W*DH21+Z
   20               CONTINUE
               GO TO 140
   30          CONTINUE
               DH11=DPARAM(2)
               DH22=DPARAM(5)
                    DO 40 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z
                    DY(I)=-W+DH22*Z
   40               CONTINUE
               GO TO 140
   50          CONTINUE
               DH11=DPARAM(2)
               DH12=DPARAM(4)
               DH21=DPARAM(3)
               DH22=DPARAM(5)
                    DO 60 I=1,NSTEPS,INCX
                    W=DX(I)
                    Z=DY(I)
                    DX(I)=W*DH11+Z*DH12
                    DY(I)=W*DH21+Z*DH22
   60               CONTINUE
               GO TO 140
   70     CONTINUE
          KX=1
          KY=1
          IF(INCX .LT. 0) KX=1+(1-N)*INCX
          IF(INCY .LT. 0) KY=1+(1-N)*INCY
c
          IF(DFLAG)120,80,100
   80     CONTINUE
          DH12=DPARAM(4)
          DH21=DPARAM(3)
               DO 90 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W+Z*DH12
               DY(KY)=W*DH21+Z
               KX=KX+INCX
               KY=KY+INCY
   90          CONTINUE
          GO TO 140
  100     CONTINUE
          DH11=DPARAM(2)
          DH22=DPARAM(5)
               DO 110 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z
               DY(KY)=-W+DH22*Z
               KX=KX+INCX
               KY=KY+INCY
  110          CONTINUE
          GO TO 140
  120     CONTINUE
          DH11=DPARAM(2)
          DH12=DPARAM(4)
          DH21=DPARAM(3)
          DH22=DPARAM(5)
               DO 130 I=1,N
               W=DX(KX)
               Z=DY(KY)
               DX(KX)=W*DH11+Z*DH12
               DY(KY)=W*DH21+Z*DH22
               KX=KX+INCX
               KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
          RETURN
          END





      SUBROUTINE xDROTMG (DD1,DD2,DX1,DY1,DPARAM)
c
c     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
c     THE SECOND COMPONENT OF THE 2-VECTOR  (DSQRT(DD1)*DX1,DSQRT(DD2)*
c     DY2)**T.
c     WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
c
c     DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
c
c       (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
c     H=(          )    (          )    (          )    (          )
c       (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
c     LOCATIONS 2-4 OF DPARAM CONTAIN DH11, DH21, DH12, AND DH22
c     RESPECTIVELY. (VALUES OF 1.D0, -1.D0, OR 0.D0 IMPLIED BY THE
c     VALUE OF DPARAM(1) ARE NOT STORED IN DPARAM.)
c
c     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
c     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
c     OF DD1 AND DD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
c
      DOUBLE PRECISION GAM,ONE,RGAMSQ,DD2,DH11,DH21,DPARAM,DP2,                     &
     & DQ2,DU,DY1,ZERO,GAMSQ,DD1,DFLAG,DH12,DH22,DP1,DQ1,                           &
     & DTEMP,DX1,TWO
      INTEGER IGO
      DIMENSION DPARAM(5)
c
      DATA ZERO,ONE,TWO /0.D0,1.D0,2.D0/
      DATA GAM,GAMSQ,RGAMSQ/4096.D0,16777216.D0,5.9604645D-8/
      IF(.NOT. DD1 .LT. ZERO) GO TO 10
c       GO ZERO-H-D-AND-DX1..
          GO TO 60
   10 CONTINUE
c     CASE-DD1-NONNEGATIVE
      DP2=DD2*DY1
      IF(.NOT. DP2 .EQ. ZERO) GO TO 20
          DFLAG=-TWO
          GO TO 260
c     REGULAR-CASE..
   20 CONTINUE
      DP1=DD1*DX1
      DQ2=DP2*DY1
      DQ1=DP1*DX1
c
      IF(.NOT. DABS(DQ1) .GT. DABS(DQ2)) GO TO 40
          DH21=-DY1/DX1
          DH12=DP2/DP1
c
          DU=ONE-DH12*DH21
c
          IF(.NOT. DU .LE. ZERO) GO TO 30
c         GO ZERO-H-D-AND-DX1..
               GO TO 60
   30     CONTINUE
               DFLAG=ZERO
               DD1=DD1/DU
               DD2=DD2/DU
               DX1=DX1*DU
c         GO SCALE-CHECK..
               GO TO 100
   40 CONTINUE
          IF(.NOT. DQ2 .LT. ZERO) GO TO 50
c         GO ZERO-H-D-AND-DX1..
               GO TO 60
   50     CONTINUE
               DFLAG=ONE
               DH11=DP1/DP2
               DH22=DX1/DY1
               DU=ONE+DH11*DH22
               DTEMP=DD2/DU
               DD2=DD1/DU
               DD1=DTEMP
               DX1=DY1*DU
c         GO SCALE-CHECK
               GO TO 100
c     PROCEDURE..ZERO-H-D-AND-DX1..
   60 CONTINUE
          DFLAG=-ONE
          DH11=ZERO
          DH12=ZERO
          DH21=ZERO
          DH22=ZERO
c
          DD1=ZERO
          DD2=ZERO
          DX1=ZERO
c         RETURN..
          GO TO 220
c     PROCEDURE..FIX-H..
   70 CONTINUE
      IF(.NOT. DFLAG .GE. ZERO) GO TO 90
c
          IF(.NOT. DFLAG .EQ. ZERO) GO TO 80
          DH11=ONE
          DH22=ONE
          DFLAG=-ONE
          GO TO 90
   80     CONTINUE
          DH21=-ONE
          DH12=ONE
          DFLAG=-ONE
   90 CONTINUE
C      GO TO IGO,(120,150,180,210)  CHANGED INTO
      SELECT CASE (IGO)

        CASE (120)
          GOTO 120
        CASE (150)
          GOTO 150
        CASE (180)
          GOTO 180
        CASE (210)
          GOTO 210
                   
      END SELECT
      

C karline: IGO DOES NOT HAVE A VALUE ?
c     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
          IF(.NOT. DD1 .LE. RGAMSQ) GO TO 130
               IF(DD1 .EQ. ZERO) GO TO 160
C               ASSIGN 120 TO IGO  CHANGED INTO:
               IGO = 120
c              FIX-H..
               GO TO 70
  120          CONTINUE
               DD1=DD1*GAM**2
               DX1=DX1/GAM
               DH11=DH11/GAM
               DH12=DH12/GAM
          GO TO 110
  130 CONTINUE
  140     CONTINUE
          IF(.NOT. DD1 .GE. GAMSQ) GO TO 160
C               ASSIGN 150 TO IGO   CAHNGED INTO
               IGO = 150
c              FIX-H..
               GO TO 70
  150          CONTINUE
               DD1=DD1/GAM**2
               DX1=DX1*GAM
               DH11=DH11*GAM
               DH12=DH12*GAM
          GO TO 140
  160 CONTINUE
  170     CONTINUE
          IF(.NOT. DABS(DD2) .LE. RGAMSQ) GO TO 190
               IF(DD2 .EQ. ZERO) GO TO 220
C               ASSIGN 180 TO IGO  CHANGED INTO
               IGO = 180
c              FIX-H..
               GO TO 70
  180          CONTINUE
               DD2=DD2*GAM**2
               DH21=DH21/GAM
               DH22=DH22/GAM
          GO TO 170
  190 CONTINUE
  200     CONTINUE
          IF(.NOT. DABS(DD2) .GE. GAMSQ) GO TO 220
C               ASSIGN 210 TO IGO  CHANGED INTO
               IGO = 210

c              FIX-H..
               GO TO 70
  210          CONTINUE
               DD2=DD2/GAM**2
               DH21=DH21*GAM
               DH22=DH22*GAM
          GO TO 200
  220 CONTINUE
          IF(DFLAG)250,230,240
  230     CONTINUE
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               GO TO 260
  240     CONTINUE
               DPARAM(2)=DH11
               DPARAM(5)=DH22
               GO TO 260
  250     CONTINUE
               DPARAM(2)=DH11
               DPARAM(3)=DH21
               DPARAM(4)=DH12
               DPARAM(5)=DH22
  260 CONTINUE
          DPARAM(1)=DFLAG
          RETURN
      END



      SUBROUTINE DLPDP (A, MDA, M, N1, N2, PRGOPT, X, WNORM, MODE, WS,              &
     &   IS)
c***BEGIN PROLOGUE  DLPDP
c***SUBSIDIARY
c***PURPOSE  Subsidiary to xDLSEI
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (LPDP-S, DLPDP-D)
c***AUTHOR  Hanson, R. J., (SNLA)
c           Haskell, K. H., (SNLA)
c***DESCRIPTION
c
c  **** Double Precision version of LPDP ****
c     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
c     where N=N1+N2.  This is a slight overestimate for WS(*).
c
c     Determine an N1-vector W, and
c               an N2-vector Z
c     which minimizes the Euclidean length of W
c     subject to G*W+H*Z .GE. Y.
c     This is the least projected distance problem, LPDP.
c     The matrices G and H are of respective
c     dimensions M by N1 and M by N2.
c
c     Called by subprogram DLSI( ).
c
c     The matrix
c                (G H Y)
c
c     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
c
c     The solution (W) is returned in X(*).
c                  (Z)
c
c     The value of MODE indicates the status of
c     the computation after returning to the user.
c
c          MODE=1  The solution was successfully obtained.
c
c          MODE=2  The inequalities are inconsistent.
c
c***SEE ALSO  xDLSEI
c***ROUTINES CALLED  xDCOPY, xDDOT, xDNRM2, xDSCAL, DWNNLS
c***REVISION HISTORY  (YYMMDD)
c   790701  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900328  Added TYPE section.  (WRB)
c   910408  Updated the AUTHOR section.  (WRB)
c***END PROLOGUE  DLPDP
c
      INTEGER I, IS(*), IW, IX, J, L, M, MDA, MODE, MODEW, N, N1, N2,               &
     &     NP1
      DOUBLE PRECISION A(MDA,*), xDDOT, xDNRM2, FAC, ONE,                             &
     &     PRGOPT(*), RNORM, SC, WNORM, WS(*), X(*), YNORM, ZERO
      SAVE ZERO, ONE, FAC
      DATA ZERO,ONE /0.0D0,1.0D0/, FAC /0.1D0/
c***FIRST EXECUTABLE STATEMENT  DLPDP
      N = N1 + N2
      MODE = 1
      IF (M .GT. 0) GO TO 20
         IF (N .LE. 0) GO TO 10
            X(1) = ZERO
            CALL xDCOPY(N,X,0,X,1)
   10    CONTINUE
         WNORM = ZERO
      GO TO 200
   20 CONTINUE
c        BEGIN BLOCK PERMITTING ...EXITS TO 190
            NP1 = N + 1
c
c           SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
            DO 40 I = 1, M
               SC = xDNRM2(N,A(I,1),MDA)
               IF (SC .EQ. ZERO) GO TO 30
                  SC = ONE/SC
                  CALL xDSCAL(NP1,SC,A(I,1),MDA)
   30          CONTINUE
   40       CONTINUE
c
c           SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
            YNORM = xDNRM2(M,A(1,NP1),1)
            IF (YNORM .EQ. ZERO) GO TO 50
               SC = ONE/YNORM
               CALL xDSCAL(M,SC,A(1,NP1),1)
   50       CONTINUE
c
c           SCALE COLS OF MATRIX H.
            J = N1 + 1
   60       IF (J .GT. N) GO TO 70
               SC = xDNRM2(M,A(1,J),1)
               IF (SC .NE. ZERO) SC = ONE/SC
               CALL xDSCAL(M,SC,A(1,J),1)
               X(J) = SC
               J = J + 1
            GO TO 60
   70       CONTINUE
            IF (N1 .LE. 0) GO TO 130
c
c              COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
               IW = 0
               DO 80 I = 1, M
c
c                 MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
                  CALL xDCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
c
c                 MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
                  CALL xDCOPY(N1,A(I,1),MDA,WS(IW+1),1)
                  IW = IW + N1
c
c                 MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
   80          CONTINUE
               WS(IW+1) = ZERO
               CALL xDCOPY(N,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N
               WS(IW+1) = ONE
               IW = IW + 1
c
c              SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U.GE.0.  THE
c              MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
c              F = TRANSPOSE OF (0,...,0,1).
               IX = IW + 1
               IW = IW + M
c
c              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
c              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,NP1,N2,NP1-N2,M,0,PRGOPT,WS(IX),RNORM,                &
     &                     MODEW,IS,WS(IW+1))
c
c              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
               SC = ONE - xDDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*ABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)                 &
     &            GO TO 110
                  SC = ONE/SC
                  DO 90 J = 1, N1
                     X(J) = SC*xDDOT(M,A(1,J),1,WS(IX),1)
   90             CONTINUE
c
c                 COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS
c                 VECTOR.
                  DO 100 I = 1, M
                     A(I,NP1) = A(I,NP1) - xDDOT(N1,A(I,1),MDA,X,1)
  100             CONTINUE
               GO TO 120
  110          CONTINUE
                  MODE = 2
c        .........EXIT
                  GO TO 190
  120          CONTINUE
  130       CONTINUE
            IF (N2 .LE. 0) GO TO 180
c
c              COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
               IW = 0
               DO 140 I = 1, M
                  CALL xDCOPY(N2,A(I,N1+1),MDA,WS(IW+1),1)
                  IW = IW + N2
                  WS(IW+1) = A(I,NP1)
                  IW = IW + 1
  140          CONTINUE
               WS(IW+1) = ZERO
               CALL xDCOPY(N2,WS(IW+1),0,WS(IW+1),1)
               IW = IW + N2
               WS(IW+1) = ONE
               IW = IW + 1
               IX = IW + 1
               IW = IW + M
c
c              SOLVE RV=S SUBJECT TO V.GE.0.  THE MATRIX R =(TRANSPOSE
c              OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
c              OF (0,...,0,1)).
c
c              DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF
c              DWNNLS( ).
               IS(1) = 0
               IS(2) = 0
               CALL DWNNLS(WS,N2+1,0,N2+1,M,0,PRGOPT,WS(IX),RNORM,MODEW,            &
     &                     IS,WS(IW+1))
c
c              COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
               SC = ONE - xDDOT(M,A(1,NP1),1,WS(IX),1)
               IF (ONE + FAC*ABS(SC) .EQ. ONE .OR. RNORM .LE. ZERO)                 &
     &            GO TO 160
                  SC = ONE/SC
                  DO 150 J = 1, N2
                     L = N1 + J
                     X(L) = SC*xDDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150             CONTINUE
               GO TO 170
  160          CONTINUE
                  MODE = 2
c        .........EXIT
                  GO TO 190
  170          CONTINUE
  180       CONTINUE
c
c           ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
            CALL xDSCAL(N,YNORM,X,1)
            WNORM = xDNRM2(N1,X,1)
  190    CONTINUE
  200 CONTINUE
      RETURN
      END




      SUBROUTINE xDH12 (MODE, LPIVOT, L1, M, U, IUE, UP, C, ICE, ICV,                &
     &   NCV)
c***BEGIN PROLOGUE  DH12
c***SUBSIDIARY
c***PURPOSE  Subsidiary to xDHFTI, xDLSEI and DWNNLS
c***LIBRARY   SLATEC
c***TYPE      DOUBLE PRECISION (H12-S, DH12-D)
c***AUTHOR  (UNKNOWN)
c***DESCRIPTION
c
c      &** DOUBLE PRECISION VERSION OF H12 ******
c
c     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
c     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
c
c     Construction and/or application of a single
c     Householder transformation..     Q = I + U*(U**T)/B
c
c     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
c     LPIVOT is the index of the pivot element.
c     L1,M   If L1 .LE. M   the transformation will be constructed to
c            zero elements indexed from L1 through M.   If L1 GT. M
c            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
c     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
c                   IUE is the storage increment between elements.
c                                       On exit from H1 U() and UP
c                   contain quantities defining the vector U of the
c                   Householder transformation.   On entry to H2 U()
c                   and UP should contain quantities previously computed
c                   by H1.  These will not be modified by H2.
c     C()    On entry to H1 or H2 C() contains a matrix which will be
c            regarded as a set of vectors to which the Householder
c            transformation is to be applied.  On exit C() contains the
c            set of transformed vectors.
c     ICE    Storage increment between elements of vectors in C().
c     ICV    Storage increment between vectors in C().
c     NCV    Number of vectors in C() to be transformed. If NCV .LE. 0
c            no operations will be done on C().
c
c***SEE ALSO  xDHFTI, xDLSEI, DWNNLS
c***ROUTINES CALLED  xDAXPY, xDDOT, xDSWAP
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   890831  Modified array declarations.  (WRB)
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900328  Added TYPE section.  (WRB)
c   900911  Added xDDOT to DOUBLE PRECISION statement.  (WRB)
c***END PROLOGUE  DH12
      INTEGER I, I2, I3, I4, ICE, ICV, INCR, IUE, J, KL1, KL2, KLP,                 &
     &     L1, L1M1, LPIVOT, M, MML1P2, MODE, NCV
      DOUBLE PRECISION B, C, CL, CLINV, ONE, UL1M1, SM, U, UP, xDDOT
      DIMENSION U(IUE,*), C(*)
c     BEGIN BLOCK PERMITTING ...EXITS TO 140
c***FIRST EXECUTABLE STATEMENT  DH12
         ONE = 1.0D0
c
c     ...EXIT
         IF (0 .GE. LPIVOT .OR. LPIVOT .GE. L1 .OR. L1 .GT. M) GO TO 140
         CL = ABS(U(1,LPIVOT))
         IF (MODE .EQ. 2) GO TO 40
c           &***** CONSTRUCT THE TRANSFORMATION. ******
            DO 10 J = L1, M
               CL = MAX(ABS(U(1,J)),CL)
   10       CONTINUE
            IF (CL .GT. 0.0D0) GO TO 20
c     .........EXIT
               GO TO 140
   20       CONTINUE
            CLINV = ONE/CL
            SM = (U(1,LPIVOT)*CLINV)**2
            DO 30 J = L1, M
               SM = SM + (U(1,J)*CLINV)**2
   30       CONTINUE
            CL = CL*SQRT(SM)
            IF (U(1,LPIVOT) .GT. 0.0D0) CL = -CL
            UP = U(1,LPIVOT) - CL
            U(1,LPIVOT) = CL
         GO TO 50
   40    CONTINUE
c        &***** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
c
         IF (CL .GT. 0.0D0) GO TO 50
c     ......EXIT
            GO TO 140
   50    CONTINUE
c     ...EXIT
         IF (NCV .LE. 0) GO TO 140
         B = UP*U(1,LPIVOT)
c        B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
c
         IF (B .LT. 0.0D0) GO TO 60
c     ......EXIT
            GO TO 140
   60    CONTINUE
         B = ONE/B
         MML1P2 = M - L1 + 2
         IF (MML1P2 .LE. 20) GO TO 80
            L1M1 = L1 - 1
            KL1 = 1 + (L1M1 - 1)*ICE
            KL2 = KL1
            KLP = 1 + (LPIVOT - 1)*ICE
            UL1M1 = U(1,L1M1)
            U(1,L1M1) = UP
            IF (LPIVOT .NE. L1M1) CALL xDSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
            DO 70 J = 1, NCV
               SM = xDDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
               SM = SM*B
               CALL xDAXPY(MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
               KL1 = KL1 + ICV
   70       CONTINUE
            U(1,L1M1) = UL1M1
c     ......EXIT
            IF (LPIVOT .EQ. L1M1) GO TO 140
            KL1 = KL2
            CALL xDSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
         GO TO 130
   80    CONTINUE
            I2 = 1 - ICV + ICE*(LPIVOT - 1)
            INCR = ICE*(L1 - LPIVOT)
            DO 120 J = 1, NCV
               I2 = I2 + ICV
               I3 = I2 + INCR
               I4 = I3
               SM = C(I2)*UP
               DO 90 I = L1, M
                  SM = SM + C(I3)*U(1,I)
                  I3 = I3 + ICE
   90          CONTINUE
               IF (SM .EQ. 0.0D0) GO TO 110
                  SM = SM*B
                  C(I2) = C(I2) + SM*UP
                  DO 100 I = L1, M
                     C(I4) = C(I4) + SM*U(1,I)
                     I4 = I4 + ICE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
      RETURN
      END




      SUBROUTINE xDHFTI (A, MDA, M, N, B, MDB, NB, TAU, KRANK, RNORM, H,             &
     &   G, IP)
c***BEGIN PROLOGUE  xDHFTI
c***PURPOSE  Solve a least squares problem for banded matrices using
c            sequential accumulation of rows of the data matrix.
c            Exactly one right-hand side vector is permitted.
c***LIBRARY   SLATEC
c***CATEGORY  D9
c***TYPE      DOUBLE PRECISION (HFTI-S, xDHFTI-D)
c***KEYWORDS  CURVE FITTING, LEAST SQUARES
c***AUTHOR  Lawson, C. L., (JPL)
c           Hanson, R. J., (SNLA)
c***DESCRIPTION
c
c     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
c
c     This subroutine solves a linear least squares problem or a set of
c     linear least squares problems having the same matrix but different
c     right-side vectors.  The problem data consists of an M by N matrix
c     A, an M by NB matrix B, and an absolute tolerance parameter TAU
c     whose usage is described below.  The NB column vectors of B
c     represent right-side vectors for NB distinct linear least squares
c     problems.
c
c     This set of problems can also be written as the matrix least
c     squares problem
c
c                       AX = B,
c
c     where X is the N by NB solution matrix.
c
c     Note that if B is the M by M identity matrix, then X will be the
c     pseudo-inverse of A.
c
c     This subroutine first transforms the augmented matrix (A B) to a
c     matrix (R C) using premultiplying Householder transformations with
c     column interchanges.  All subdiagonal elements in the matrix R are
c     zero and its diagonal elements satisfy
c
c                       ABS(R(I,I)).GE.ABS(R(I+1,I+1)),
c
c                       I = 1,...,L-1, where
c
c                       L = MIN(M,N).
c
c     The subroutine will compute an integer, KRANK, equal to the number
c     of diagonal terms of R that exceed TAU in magnitude. Then a
c     solution of minimum Euclidean length is computed using the first
c     KRANK rows of (R C).
c
c     To be specific we suggest that the user consider an easily
c     computable matrix norm, such as, the maximum of all column sums of
c     magnitudes.
c
c     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
c     norm of B), it is suggested that TAU be set approximately equal to
c     EPS*(norm of A).
c
c     The user must dimension all arrays appearing in the call list..
c     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
c     permits the solution of a range of problems in the same array
c     space.
c
c     The entire set of parameters for xDHFTI are
c
c     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
c
c     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
c                       matrix A of the least squares problem AX = B.
c                       The first dimensioning parameter of the array
c                       A(*,*) is MDA, which must satisfy MDA.GE.M
c                       Either M.GE.N or M.LT.N is permitted.  There
c                       is no restriction on the rank of A.  The
c                       condition MDA.LT.M is considered an error.
c
c     B(*),MDB,NB       If NB = 0 the subroutine will perform the
c                       orthogonal decomposition but will make no
c                       references to the array B(*).  If NB.GT.0
c                       the array B(*) must initially contain the M by
c                       NB matrix B of the least squares problem AX =
c                       B.  If NB.GE.2 the array B(*) must be doubly
c                       subscripted with first dimensioning parameter
c                       MDB.GE.MAX(M,N).  If NB = 1 the array B(*) may
c                       be either doubly or singly subscripted.  In
c                       the latter case the value of MDB is arbitrary
c                       but it should be set to some valid integer
c                       value such as MDB = M.
c
c                       The condition of NB.GT.1.AND.MDB.LT. MAX(M,N)
c                       is considered an error.
c
c     TAU               Absolute tolerance parameter provided by user
c                       for pseudorank determination.
c
c     H(*),G(*),IP(*)   Arrays of working space used by xDHFTI.
c
c     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
c
c     A(*,*)            The contents of the array A(*,*) will be
c                       modified by the subroutine. These contents
c                       are not generally required by the user.
c
c     B(*)              On return the array B(*) will contain the N by
c                       NB solution matrix X.
c
c     KRANK             Set by the subroutine to indicate the
c                       pseudorank of A.
c
c     RNORM(*)          On return, RNORM(J) will contain the Euclidean
c                       norm of the residual vector for the problem
c                       defined by the J-th column vector of the array
c                       B(*,*) for J = 1,...,NB.
c
c     H(*),G(*)         On return these arrays respectively contain
c                       elements of the pre- and post-multiplying
c                       Householder transformations used to compute
c                       the minimum Euclidean length solution.
c
c     IP(*)             Array in which the subroutine records indices
c                       describing the permutation of column vectors.
c                       The contents of arrays H(*),G(*) and IP(*)
c                       are not generally required by the user.
c
c***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
c                 Problems, Prentice-Hall, Inc., 1974, Chapter 14.
c***ROUTINES CALLED  D1MACH, xDH12, xXERMSG
c***REVISION HISTORY  (YYMMDD)
c   790101  DATE WRITTEN
c   890531  Changed all specific intrinsics to generic.  (WRB)
c   891006  Cosmetic changes to prologue.  (WRB)
c   891006  REVISION DATE from Version 3.2
c   891214  Prologue converted to Version 4.0 format.  (BAB)
c   900315  CALLs to XERROR changed to CALLs to xXERMSG.  (THJ)
c   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC)
c   920501  Reformatted the REFERENCES section.  (WRB)
c***END PROLOGUE  xDHFTI
      INTEGER I, II, IOPT, IP(*), IP1, J, JB, JJ, K, KP1, KRANK, L,                 &
     &     LDIAG, LMAX, M, MDA, MDB, N, NB, NERR
      DOUBLE PRECISION A, B, D1MACH, DZERO, FACTOR,                                 &
     &     G, H, HMAX, RELEPS, RNORM, SM, SM1, SZERO, TAU, TMP
      DIMENSION A(MDA,*),B(MDB,*),H(*),G(*),RNORM(*)
      SAVE RELEPS
      DATA RELEPS /0.D0/
c     BEGIN BLOCK PERMITTING ...EXITS TO 360
c***FIRST EXECUTABLE STATEMENT  xDHFTI
         IF (RELEPS.EQ.0.D0) RELEPS = D1MACH(4)
         SZERO = 0.0D0
         DZERO = 0.0D0
         FACTOR = 0.001D0
c
         K = 0
         LDIAG = MIN(M,N)
         IF (LDIAG .LE. 0) GO TO 350
c           BEGIN BLOCK PERMITTING ...EXITS TO 130
c              BEGIN BLOCK PERMITTING ...EXITS TO 120
                  IF (MDA .GE. M) GO TO 10
                     NERR = 1
                     IOPT = 2
                     CALL xXERMSG ('SLATEC', 'xDHFTI',                               &
     &                  'MDA.LT.M, PROBABLE ERROR.',                                &
     &                  NERR, IOPT)
c     ...............EXIT
                     GO TO 360
   10             CONTINUE
c
                  IF (NB .LE. 1 .OR. MAX(M,N) .LE. MDB) GO TO 20
                     NERR = 2
                     IOPT = 2
                     CALL xXERMSG ('SLATEC', 'xDHFTI',                               &
     &                  'MDB.LT.MAX(M,N).AND.NB.GT.1. PROBABLE ERROR.',             &
     &                  NERR, IOPT)
c     ...............EXIT
                     GO TO 360
   20             CONTINUE
c
                  DO 100 J = 1, LDIAG
c                    BEGIN BLOCK PERMITTING ...EXITS TO 70
                        IF (J .EQ. 1) GO TO 40
c
c                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
c                          ..
                           LMAX = J
                           DO 30 L = J, N
                              H(L) = H(L) - A(J-1,L)**2
                              IF (H(L) .GT. H(LMAX)) LMAX = L
   30                      CONTINUE
c                    ......EXIT
                           IF (FACTOR*H(LMAX) .GT. HMAX*RELEPS) GO TO 70
   40                   CONTINUE
c
c                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
c                       ..
                        LMAX = J
                        DO 60 L = J, N
                           H(L) = 0.0D0
                           DO 50 I = J, M
                              H(L) = H(L) + A(I,L)**2
   50                      CONTINUE
                           IF (H(L) .GT. H(LMAX)) LMAX = L
   60                   CONTINUE
                        HMAX = H(LMAX)
   70                CONTINUE
c                    ..
c                     LMAX HAS BEEN DETERMINED
c
c                     DO COLUMN INTERCHANGES IF NEEDED.
c                    ..
                     IP(J) = LMAX
                     IF (IP(J) .EQ. J) GO TO 90
                        DO 80 I = 1, M
                           TMP = A(I,J)
                           A(I,J) = A(I,LMAX)
                           A(I,LMAX) = TMP
   80                   CONTINUE
                        H(LMAX) = H(J)
   90                CONTINUE
c
c                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
c                     AND B.
c                    ..
                     CALL xDH12(1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,              &
     &                         N-J)
                     CALL xDH12(2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
  100             CONTINUE
c
c                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE,
c                  TAU.
c                 ..
                  DO 110 J = 1, LDIAG
c              ......EXIT
                     IF (ABS(A(J,J)) .LE. TAU) GO TO 120
  110             CONTINUE
                  K = LDIAG
c           ......EXIT
                  GO TO 130
  120          CONTINUE
               K = J - 1
  130       CONTINUE
            KP1 = K + 1
c
c           COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
c
            IF (NB .LT. 1) GO TO 170
            DO 160 JB = 1, NB
               TMP = SZERO
               IF (M .LT. KP1) GO TO 150
               DO 140 I = KP1, M
                  TMP = TMP + B(I,JB)**2
  140          CONTINUE
  150          CONTINUE
               RNORM(JB) = SQRT(TMP)
  160       CONTINUE
  170       CONTINUE
c           SPECIAL FOR PSEUDORANK = 0
            IF (K .GT. 0) GO TO 210
               IF (NB .LT. 1) GO TO 200
               DO 190 JB = 1, NB
                  DO 180 I = 1, N
                     B(I,JB) = SZERO
  180             CONTINUE
  190          CONTINUE
  200          CONTINUE
            GO TO 340
  210       CONTINUE
c
c               IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
c               DECOMPOSITION OF FIRST K ROWS.
c              ..
               IF (K .EQ. N) GO TO 230
                  DO 220 II = 1, K
                     I = KP1 - II
                     CALL xDH12(1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  220             CONTINUE
  230          CONTINUE
c
c
               IF (NB .LT. 1) GO TO 330
               DO 320 JB = 1, NB
c
c                  SOLVE THE K BY K TRIANGULAR SYSTEM.
c                 ..
                  DO 260 L = 1, K
                     SM = DZERO
                     I = KP1 - L
                     IP1 = I + 1
                     IF (K .LT. IP1) GO TO 250
                     DO 240 J = IP1, K
                        SM = SM + A(I,J)*B(J,JB)
  240                CONTINUE
  250                CONTINUE
                     SM1 = SM
                     B(I,JB) = (B(I,JB) - SM1)/A(I,I)
  260             CONTINUE
c
c                  COMPLETE COMPUTATION OF SOLUTION VECTOR.
c                 ..
                  IF (K .EQ. N) GO TO 290
                     DO 270 J = KP1, N
                        B(J,JB) = SZERO
  270                CONTINUE
                     DO 280 I = 1, K
                        CALL xDH12(2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,              &
     &                            MDB,1)
  280                CONTINUE
  290             CONTINUE
c
c                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
c                   COLUMN INTERCHANGES.
c                 ..
                  DO 310 JJ = 1, LDIAG
                     J = LDIAG + 1 - JJ
                     IF (IP(J) .EQ. J) GO TO 300
                        L = IP(J)
                        TMP = B(L,JB)
                        B(L,JB) = B(J,JB)
                        B(J,JB) = TMP
  300                CONTINUE
  310             CONTINUE
  320          CONTINUE
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
c        ..
c         THE SOLUTION VECTORS, X, ARE NOW
c         IN THE FIRST  N  ROWS OF THE ARRAY B(,).
c
         KRANK = K
  360 CONTINUE
      RETURN
      END



C***********************************************************************
C Karline: added as adaptation from xdcopy, to avoid warning of rank mismatch
C***********************************************************************
      SUBROUTINE XDCOPYSC(N,DX,DY,INCY)
C
C     COPIES A scalar, X, TO A VECTOR, Y. 
C
      DOUBLE PRECISION DX, DY(*)
      INTEGER I,INCY,IY,N
C
      IF(N.LE.0)RETURN

      IY = 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1

      DO 10 I = 1,N
        DY(IY) = DX
        IY = IY + INCY
   10 CONTINUE
      RETURN
      END

