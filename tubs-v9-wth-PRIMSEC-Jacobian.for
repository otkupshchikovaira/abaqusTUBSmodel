C==================================================================
C  UMAT with real PRCREEP / SECREEP / TECREEP / VPFLOW
C  - IMPLICIT NONE everywhere
C  - supports NTENS = 3,4,6 (uses only first NTENS components)
c  - all material constants are PARAMETER inside subroutines
C  - safety tiny to avoid div/0
C - check connection of STRESS and STATEV for primary creep strains (STATEV(K+6) used by original)
C==================================================================

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
      IMPLICIT NONE

C --- standard UMAT arguments (explicit types)
      CHARACTER*80 CMNAME
      INTEGER NDI, NSHR, NTENS, NSTATV, NPROPS
      INTEGER NOEL, NPT, LAYER, KSPT, KSTEP, KINC
      REAL*8 STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS,NTENS)
      REAL*8 SSE, SPD, SCD, RPL, DDSDDT(*), DRPLDE(*), DRPLDT(*)
      REAL*8 STRAN(NTENS), DSTRAN(NTENS), TIME(2), DTIME, TEMP, DTEMP
      REAL*8 PREDEF(*), DPRED(*), PROPS(*), COORDS(*)
      REAL*8 DROT(3,3), PNEWDT, CELENT, DFGRD0(3,3), DFGRD1(3,3)

C --- locals
      INTEGER I,J,K,LSTR
      REAL*8 DDS(6,6), DELSTRAN(6), DSTRESS(6),TERTSTRAN(6)
      REAL*8 PRIMSTRAN(6), PRIMCREEPRATE(6), SECCREEPRATE(6)
      REAL*8 TERTCREEPRATE(6), VISCOPLASTICRATE(6), DQP(6), DQF(6)
      REAL*8 DDQ(6),DQ(6),DFP(6),DSFP(6),DSFS(6),M1P(6,6),M1S(6,6)
      REAL*8 M2P(6,6),M2S(6,6), AP(6,6), AS(6,6),KSUM(6,6)
      REAL*8 EFFSIGMA, EFFPRIMEPS, PRESS, TEPSVOL,DENOM1,F_active
      REAL*8 E, NU, EBULK3, EG2, EG, ELAM,DENOM,FP,FS,FS1,DOUTER
      REAL*8 tiny
      REAL*8 S(6), PS(3)
      REAL*8 alphaP, alphaS, alphaEff
      PARAMETER (E = 30.0D9, NU = 0.3D0)
      PARAMETER (tiny = 1.0D-20)
      external dgetrf, dgetrs
C --- material constants (internal parametrization)


      NDI  = 3
      NSHR = 3
      LSTR = 1   ! 1 = stresses

C     ???????? ??????? ??????????:
C     PS(1), PS(2), PS(3) � ?????? ??????? ???????? ??????????
C --- elastic moduli
      EBULK3 = E / (1.0D0 - 2.0D0*NU)
      EG2 = E / (1.0D0 + NU)
      EG  = EG2 / 2.0D0
      ELAM = (EBULK3 - EG2) / 3.0D0

C --- zero local arrays (6-size storage)
      DO I=1,6
         DO J=1,6
            DDS(I,J) = 0.0D0
            DDQ(I,J) = 0.0D0
            M1P(I,J) = 0.0D0
            M1S(I,J) = 0.0D0
            M2P(I,J) = 0.0D0
            M2S(I,J) = 0.0D0
            AP(I,J) = 0.0D0
            AS(I,J) = 0.0D0
            KSUM(I,J) = 0.0D0
         END DO
         DELSTRAN(I) = 0.0D0
         DSTRESS(I)  = 0.0D0
         PRIMSTRAN(I) = 0.0D0
         TERTSTRAN(I)=0.0D0
         PRIMCREEPRATE(I) = 0.0D0
         SECCREEPRATE(I) = 0.0D0
         TERTCREEPRATE(I) = 0.0D0
         VISCOPLASTICRATE(I) = 0.0D0
         DFP(I) = 0.0D0
         DFS(I) = 0.0D0
         DSFP(I) = 0.0D0
         DSFS(I) = 0.0D0
         DQ(I) = 0.0D0
         DQP(I) = 0.0D0
         DQF(I) = 0.0D0
         STRESS1(I)=0.0D0
      END DO

C --- form elastic stiffness (only relevant components set)
      CALL FORM_ELASTIC_DDS_USER(DDS, NTENS, EG2, ELAM, EG)

C --- compute effective stress and pressure (safe for NTENS < 6)
      CALL COMPUTE_EFFECTIVE_USER(STRESS, NTENS, EFFSIGMA, PRESS)
C --- get primary strain from STATEV if available (STATEV(K+6) used by original)
      DO K = 1, NTENS
         IF ( (K+6) .LE. NSTATV ) THEN
            PRIMSTRAN(K) = STATEV(K+6)
            TERTSTRAN(K) = STATEV(K+6+NTENS)
            STRESS1(K)=STATEV(6+2*NTENS+3+K)
         ELSE
            PRIMSTRAN(K) = 0.0D0
            TERTSTRAN(K) = 0.0D0
         END IF
      END DO
      CALL COMPUTE_EFFECTIVE_EPS_USER(PRIMSTRAN, NTENS, EFFPRIMEPS)

C --- volumetric strain
C      IF (NTENS .GE. 3) THEN
         TEPSVOL = (TERTSTRAN(1) + TERTSTRAN(2) + TERTSTRAN(3)) / 3.0D0
C      ELSE
C         EPSVOL = 0.0D0
C      END IF

      CALL SPRINC(STRESS, PS, 1, NDI, NSHR)
      PS(1)=PS(1)
      PS(2)=PS(2)
      PS(3)=PS(3)
C --- call creep / viscoplastic subroutines (they fill 6-component rate arrays)
      CALL PRCREEP_REAL(EFFSIGMA, EFFPRIMEPS, DTIME, PRIMCREEPRATE, 
     1 DQP, STRESS, NTENS,DENOM,FP)
      DENOM1=DENOM
      CALL SECREEP_REAL(EFFSIGMA, DTIME, SECCREEPRATE, STRESS, NTENS,FS)
      FS1=FS
      CALL TECREEP_REAL(TEPSVOL, PRESS, EFFSIGMA, TERTCREEPRATE, 
     1 STRESS, NTENS)
      CALL VPFLOW_REAL(TEPSVOL, PRESS, EFFSIGMA, VISCOPLASTICRATE,
     1 STRESS, NTENS, DTIME,F_active,PS)
      F_active=F_active
C=====================================================================
C     MULTIPLICATIOM for DDSDDE
C=====================================================================
      CALL DDSDDEDERIVATIVES(DDQ, DQ, STRESS1, DENOM, NTENS)
      CALL PRIMARYFUNCTIONS(DOUTER,STRESS1,NTENS,EFFSIGMA,
     1 EFFPRIMEPS,DENOM,FP,DFP)
      CALL SECONDARYFUNCTIONS(DOUTER,STRESS1,NTENS,EFFSIGMA,
     1 DENOM,FS,DFS)    
C      CE*DFP/CE*DFS
      CALL MAT6x6_VEC6(DDS, DFP, DSFP)
      CALL MAT6x6_VEC6(DDS, DFS, DSFS)
      CALL VEC6_VEC6(DSFP, DQ, M1P)
      CALL VEC6_VEC6(DSFS, DQ, M1S)
      CALL DGEMM('N','N',6,6,6,FP, DDS,6, DDQ,6, 0.0D0, M2P,6)
      CALL DGEMM('N','N',6,6,6,FS, DDS,6, DDQ,6, 0.0D0, M2S,6)
      DO I = 1, 6
       DO J = 1, 6
          AP(I,J) = 1/NP*(M1P(I,J) + M1P(I,J))
          AS(I,J) = 1/NS*(M1S(I,J) + M1S(I,J))
          KSUM(I,J) = DTIME*(AP(I,J)+AS(I,J))
       END DO
       KSUM(I,I)=M(I,I)+1.0D0
      END DO
      CALL dgetrf(6, 6, KSUM, 6, ipiv, info)

      IF (info .ne. 0) THEN
        WRITE(*,*) 'WARNING: KSUM is singular or badly conditioned'
        ! ????? ????? ????????? dt ??? ???????? regularization
      END IF
      RHS = DDS
      CALL dgetrs('N', 6, 6, KSUM, 6, ipiv, RHS, 6, info)
      IF (info .ne. 0) THEN
        WRITE(*,*) 'ERROR solving linear system'
      END IF
      DDSDDE = RHS
      CALL MAT6x6_VEC6(DDSDDE, DSTRAN, DSTRESS1)      
C=====================================================================
C     UPDATE STRAINS AND STRESSES
C=====================================================================

C --- build elastic increment: DELSTRAN = DSTRAN - creep_rates * DTIME
      DO K = 1, NTENS
C         SECCREEPRATE(K)=0.0
         VISCOPLASTICRATE(K)=0.0
         TERTCREEPRATE(K)=0.0
C         PRIMCREEPRATE(K)=0.0
         DELSTRAN(K) = DSTRAN(K) - (PRIMCREEPRATE(K)+SECCREEPRATE(K)+
     1    TERTCREEPRATE(K)+VISCOPLASTICRATE(K))*DTIME
      END DO

C --- compute stress increment DSTRESS = DDS * DELSTRAN
      CALL KMLT1_USER(DDS, DELSTRAN, DSTRESS, NTENS)

CC --- update stress
      DO K = 1, NTENS
         STRESS(K) = STRESS(K) + DSTRESS(K)
         STRESS1(K) = STRESS1(K) + DSTRESS1(K)
      END DO
CC --- return elastic DDS as DDSDDE (consistent tangent approximate)
      DO I = 1, NTENS
         DO J = 1, NTENS
            DDSDDE(I,J) = DDS(I,J)
         END DO
      END DO
C --- update STATEV primary creep strains (store PRIMSTRAN incrementally as in original)
      DO K = 1, NTENS
         IF ( (K+6) .LE. NSTATV ) THEN
            STATEV(K+6) = PRIMSTRAN(K) + PRIMCREEPRATE(K)*DTIME
            STATEV(K+6+NTENS)=TERTSTRAN(K)+TERTCREEPRATE(K)*DTIME
            STATEV(6+2*NTENS+3+K)=STRESS1(K)
         END IF
      END DO
      DO K=1,3
          STATEV(K)=PS(K)
      END DO
      STATEV(6+2*NTENS+1)=DSTRESS(1)/DSTRESS1(1)
      STATEV(6+2*NTENS+2)=DSTRESS(2)/DSTRESS1(2)
      STATEV(6+2*NTENS+3)=DSTRESS(3)/DSTRESS1(3)
      RETURN
      END

C==================================================================
C  FORM_ELASTIC_DDS_USER: build elastic matrix for NTENS = 3,4,6
C==================================================================
      SUBROUTINE FORM_ELASTIC_DDS_USER(DDS, NTENS, EG2, ELAM, EG)
      IMPLICIT NONE
      INTEGER NTENS, I, J
      REAL*8 DDS(6,6), EG2, ELAM, EG

      DO I = 1, 6
         DO J = 1,6
            DDS(I,J) = 0.0D0
         END DO
      END DO

      IF (NTENS .EQ. 3) THEN
         DDS(1,1) = EG2 + ELAM
         DDS(1,2) = ELAM
         DDS(1,3) = 0.0D0
         DDS(2,1) = ELAM
         DDS(2,2) = EG2 + ELAM
         DDS(2,3) = 0.0D0
         DDS(3,1) = 0.0D0
         DDS(3,2) = 0.0D0
         DDS(3,3) = EG
      ELSE IF (NTENS .EQ. 4) THEN
         DO I = 1, 3
            DO J = 1, 3
               DDS(I,J) = ELAM
            END DO
            DDS(I,I) = EG2 + ELAM
         END DO
         DDS(4,4) = EG
      ELSE
         DO I = 1, 3
            DO J = 1, 3
               DDS(I,J) = ELAM
            END DO
            DDS(I,I) = EG2 + ELAM
         END DO
         DDS(4,4) = EG
         DDS(5,5) = EG
         DDS(6,6) = EG
      END IF

      RETURN
      END

C==================================================================
C  COMPUTE_EFFECTIVE_USER: effective stress (von Mises-like norm)
C==================================================================
      SUBROUTINE COMPUTE_EFFECTIVE_USER(STRESS, NTENS, EFFSIGMA, PRESS)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 STRESS(*), EFFSIGMA, PRESS, arg, tiny
      PARAMETER (tiny = 1.0D-20)

      arg = 0.0D0
      IF (NTENS .GE. 3) THEN
         arg = 0.5D0*((STRESS(1)-STRESS(2))**2 + (STRESS(2)-
     1   STRESS(3))**2 + (STRESS(1)-STRESS(3))**2)
      END IF
      IF (NTENS .GE. 4) arg = arg +
     1    3.0D0*( ( (NTENS .GE. 4) * STRESS(4)**2 ) )
      IF (NTENS .GE. 5) arg = arg + 3.0D0*( ( (NTENS .GE. 5) *
     1 STRESS(5)**2 ) )
      IF (NTENS .GE. 6) arg = arg + 3.0D0*( ( (NTENS .GE. 6) *
     1 STRESS(6)**2 ) )
      EFFSIGMA = SQRT(MAX(arg, tiny))
      PRESS = (STRESS(1) + STRESS(2) + STRESS(3)) / 3.0D0
      RETURN
      END

C==================================================================
C  COMPUTE_EFFECTIVE_EPS_USER: effective strain measure for primary creep
C==================================================================
      SUBROUTINE COMPUTE_EFFECTIVE_EPS_USER(PRIMSTRAN, NTENS,
     1 EFFPRIMEPS)
      IMPLICIT NONE
      INTEGER NTENS
      REAL*8 PRIMSTRAN(*), EFFPRIMEPS, arg, tiny
      PARAMETER (tiny = 1.0D-20)

      arg = 0.0D0
      IF (NTENS .GE. 1) arg = arg + (PRIMSTRAN(1))**2
      IF (NTENS .GE. 2) arg = arg + (PRIMSTRAN(2))**2
      IF (NTENS .GE. 3) arg = arg + (PRIMSTRAN(3))**2
      IF (NTENS .GE. 4) arg = arg + 0.5D0*(PRIMSTRAN(4))**2
      IF (NTENS .GE. 5) arg = arg + 0.5D0*(PRIMSTRAN(5))**2
      IF (NTENS .GE. 6) arg = arg + 0.5D0*(PRIMSTRAN(6))**2
      EFFPRIMEPS = SQRT(MAX(2.0D0/3.0D0 * arg, tiny))
      RETURN
      END

C==================================================================
C  KMLT1_USER: multiply matrix (6x6) by vector (6) -> vector (6)
C==================================================================
      SUBROUTINE KMLT1_USER(DM1, DM2, DM, NTENS)
      IMPLICIT NONE
      INTEGER I, K, NTENS
      REAL*8 DM1(6,6), DM2(6), DM(6), x
      DO I = 1, NTENS
         x = 0.0D0
         DO K = 1, NTENS
            x = x + DM1(I,K) * DM2(K)
         END DO
         DM(I) = x
      END DO
      RETURN
      END

C==================================================================
C  PRCREEP_REAL: primary creep rate (real implementation)
C==================================================================
      SUBROUTINE PRCREEP_REAL(EFFSIGMA, EFFPRIMEPS, DTIME,
     1 PRIMCREEPRATE, DQP, STRESS, NTENS,DENOM,FP)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 EFFSIGMA, EFFPRIMEPS, DTIME
      REAL*8 PRIMCREEPRATE(6), DQP(6), STRESS(*)
      REAL*8 Ep, MEXP, NP, FP, arg, DENOM, s1,s2,s3,s4,s5,s6, tiny

      PARAMETER (Ep = 110.0D6, MEXP = 2.0D0, NP = 550.0D6, 
     1 tiny = 1.0D-20)

      FP = Ep * ( (EFFSIGMA / Ep)**MEXP - EFFPRIMEPS )
      IF (FP .LT. 0.0D0) FP = 0.0D0

      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF

      arg = 0.5D0 * ( (s1 - s2)**2 + (s2 - s3)**2 + (s1 - s3)**2) +
     1 3.0D0*(s4**2 + s5**2 + s6**2)
      DENOM = SQRT(MAX(arg, tiny))

      DQP(1) = 0.5D0*(2.0D0*s1 - s2 - s3) / DENOM
      DQP(2) = 0.5D0*(-s1 + 2.0D0*s2 - s3) / DENOM
      DQP(3) = 0.5D0*(-s2 + 2.0D0*s3 - s1) / DENOM
      DQP(4) = 3.0D0 * s4 / DENOM
      DQP(5) = 3.0D0 * s5 / DENOM
      DQP(6) = 3.0D0 * s6 / DENOM

      DO K = 1, 6
         PRIMCREEPRATE(K) = (FP / NP) * DQP(K)
      END DO

      RETURN
      END

C==================================================================
C  SECREEP_REAL: secondary creep rate
C==================================================================
      SUBROUTINE SECREEP_REAL(EFFSIGMA, DTIME, SECCREEPRATE,
     1 STRESS, NTENS,FS)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 EFFSIGMA, DTIME
      REAL*8 SECCREEPRATE(6), DQP(6), STRESS(*)
      REAL*8 p0, NEXP, NS, arg, DENOM, s1,s2,s3,s4,s5,s6, tiny,FS

      PARAMETER (p0 = 1.0D6, NEXP = 1.0D0, NS = 5.0D11, tiny = 1.0D-20)
      FS = p0*(EFFSIGMA/p0)**NEXP

      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF

      arg = 0.5D0 * ( (s1 - s2)**2 + (s2 - s3)**2 + (s1 - s3)**2) +
     1 3.0D0*(s4**2 + s5**2 + s6**2)
      DENOM = SQRT(MAX(arg, tiny))

      DQP(1) = 0.5D0*(2.0D0*s1 - s2 - s3) / DENOM
      DQP(2) = 0.5D0*(-s1 + 2.0D0*s2 - s3) / DENOM
      DQP(3) = 0.5D0*(-s2 + 2.0D0*s3 - s1) / DENOM
      DQP(4) = 3.0D0 * s4 / DENOM
      DQP(5) = 3.0D0 * s5 / DENOM
      DQP(6) = 3.0D0 * s6 / DENOM

      DO K = 1, 6
         SECCREEPRATE(K) = ( FS / NS ) * DQP(K)
      END DO

      RETURN
      END

C==================================================================
C  TECREEP_REAL: tertiary creep rate
C==================================================================
      SUBROUTINE TECREEP_REAL(TEPSVOL, PRESS, EFFSIGMA, TERTCREEPRATE,
     1 STRESS, NTENS)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 TEPSVOL, PRESS, EFFSIGMA
      REAL*8 TERTCREEPRATE(6), DQT(6), STRESS(*)
      REAL*8 SIGF0, NTpar, PHIf, PSI, SIGMAF, FT, arg, DENOM, s1,s2,s3,
     1 s4,s5,s6, tiny,M

      PARAMETER (SIGF0 = 15.0D6, NTpar = 2500.0D6, PHIf = 0.35D0, 
     1 PSI = 0.35D0, tiny = 1.0D-20,M=400.0D6)

      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF

      SIGMAF = SIGF0 + M*TEPSVOL  ! placeholder for M*EPSVOL if desired
      FT = (1.0D0/3.0D0*(3.0D0 - SIN(PHIf))*EFFSIGMA/(1.0D0 - 
     1 SIN(PHIf)) - 2.0D0*SIN(PHIf)*PRESS/(1.0D0 - SIN(PHIf)) - SIGMAF)
      IF (FT .LT. 0.0D0) FT = 0.0D0

      arg = 0.5D0 * ( (s1 - s2)**2 + (s2 - s3)**2 + (s1 - s3)**2) +
     1 3.0D0*(s4**2 + s5**2 + s6**2 )
      DENOM = 2.0D0 * SQRT(MAX(arg, tiny))

      DQT(1) = (1.0D0/12.0D0)*(3.0D0 - SIN(PSI))*(2.0D0*s1 - s2 - 
     1 s3)/DENOM - (1.0D0/3.0D0)*SIN(PSI)
      DQT(2) = (1.0D0/12.0D0)*(3.0D0 - SIN(PSI))*(-s1 + 2.0D0*s2 -
     1 s3)/DENOM - (1.0D0/3.0D0)*SIN(PSI)
      DQT(3) = (1.0D0/12.0D0)*(3.0D0 - SIN(PSI))*(-s2 + 2.0D0*s3 - 
     1 s1)/DENOM - (1.0D0/3.0D0)*SIN(PSI)
      DQT(4) = 0.5D0*(3.0D0 - SIN(PSI))*s4/DENOM
      DQT(5) = 0.5D0*(3.0D0 - SIN(PSI))*s5/DENOM
      DQT(6) = 0.5D0*(3.0D0 - SIN(PSI))*s6/DENOM

      DO K = 1, 6
         TERTCREEPRATE(K) = (FT / NTpar) * DQT(K)
      END DO

      RETURN
      END

C==================================================================
C  VPFLOW_REAL: viscoplastic flow rate
C==================================================================
      SUBROUTINE VPFLOW_REAL(TEPSVOL, PRESS, EFFSIGMA, VISCOPLASTICRATE,
     1 STRESS, NTENS, DTIME,F_active,PS)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 TEPSVOL, PRESS, EFFSIGMA, DTIME
      REAL*8 VISCOPLASTICRATE(6), STRESS(*),PS(3)
      REAL*8 F_active, F_tmp, SIGP_MAX, SIGTENSILE, PHI, PHIresidual
      REAL*8 SIGCresidual, NFpar, arg, DENOM, s1,s2,s3,s4,s5,s6, tiny
      REAL*8 DQF(6), SIGMAC,SIGC0,Ncoef

      PARAMETER (SIGTENSILE = 2.0D6, PHI = 0.35D0, 
     1 PHIresidual = 0.1745D0,SIGCresidual = 5.0D6, 
     2 NFpar = 2500.0D6, Ncoef = 400.0D6, SIGC0 = 25.0D6,
     3 tiny = 1.0D-20)

      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF

C----- compute compressive strength SIGMAC (eq. like in user's code)
      SIGMAC = Ncoef * TEPSVOL + SIGC0

C----- maximum principal (approx) � use normal components
      SIGP_MAX = PS(1)
      IF (PS(2) .GT. SIGP_MAX) SIGP_MAX = PS(2)
      IF (PS(3) .GT. SIGP_MAX) SIGP_MAX = PS(3)

      IF (SIGP_MAX .GT. SIGTENSILE) THEN
         F_active = SIGP_MAX
      ELSE
         F_tmp = ( (3.0D0 - SIN(PHI)) / 3.0D0 ) * EFFSIGMA / (1.0D0 -
     1   SIN(PHI)) - 2.0D0 * SIN(PHI) * PRESS / (1.0D0 - SIN(PHI)) - 
     2    SIGMAC
         IF (F_tmp .LE. 0.0D0) THEN
            F_active = 0.0D0
         ELSE
            F_active = ( (3.0D0 - SIN(PHIresidual)) / 3.0D0 ) *
     1       EFFSIGMA / (1.0D0 - SIN(PHIresidual)) - 
     2                  2.0D0 * SIN(PHIresidual) * PRESS / (1.0D0 - 
     3                  SIN(PHIresidual)) - SIGCresidual
            IF (F_active .LT. 0.0D0) F_active = 0.0D0
         END IF
      END IF

C----- flow direction (deviatoric)
      arg = 0.5D0 * ( (s1 - s2)**2 + (s2 - s3)**2 + (s1 - s3)**2) +
     1 3.0D0*(s4**2 + s5**2 + s6**2)
      DENOM = 2.0D0 * SQRT(MAX(arg, tiny))

      DQF(1) = 0.5D0*(2.0D0*s1 - s2 - s3) / DENOM
      DQF(2) = 0.5D0*(-s1 + 2.0D0*s2 - s3) / DENOM
      DQF(3) = 0.5D0*(-s2 + 2.0D0*s3 - s1) / DENOM
      DQF(4) = 3.0D0 * s4 / DENOM
      DQF(5) = 3.0D0 * s5 / DENOM
      DQF(6) = 3.0D0 * s6 / DENOM

      DO K = 1, 6
         VISCOPLASTICRATE(K) = (F_active / NFpar) * DQF(K)
      END DO

      RETURN
      END
C==================================================================
C  MULTIPLICATION MATRIX(6,6) TO VECTOR(6)
C==================================================================
      SUBROUTINE MAT6x6_VEC6(A, x, y)
        IMPLICIT NONE
        DOUBLE PRECISION A(6,6), x(6), y(6)
        INTEGER i, j
        DO i = 1, 6
        y(i) = 0.0D0
        DO j = 1, 6
            y(i) = y(i) + A(i,j) * x(j)
        END DO
        END DO
      RETURN
      END
C==================================================================
C  outer product VECTOR(6) TO VECTOR(6)
C==================================================================      
      SUBROUTINE VEC6_VEC6(a, b, C)
      IMPLICIT NONE
      DOUBLE PRECISION a(6), b(6), C(6,6)
      INTEGER i, j
      DO i = 1, 6
       DO j = 1, 6
          C(i,j) = a(i) * b(j)
       END DO
      END DO
      RETURN
      END
C==================================================================
C  End of UMAT with creep/viscoplastic implementations
C     Jacobian matrix under
C==================================================================
C==================================================================
C  DERIVATIVES FOR JACOBIAN
C==================================================================
      SUBROUTINE DDSDDEDERIVATIVES(DDQP, DQP, STRESS, NTENS)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 DOUTER,DENOM3
      REAL*8 DDQP(6,6), DQP(6), STRESS(*)
      REAL*8 arg, DENOM, s1,s2,s3,s4,s5,s6, tiny

      PARAMETER (tiny = 1.0D-20)

      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF

      arg = 0.5D0 * ( (s1 - s2)**2 + (s2 - s3)**2 + (s1 - s3)**2) +
     1 3.0D0*(s4**2 + s5**2 + s6**2)
      DENOM = SQRT(MAX(arg, tiny))
      DENOM3 = (SQRT(MAX(arg, tiny)))**3

      DQ(1) = 0.5D0*(2.0D0*s1 - s2 - s3) / DENOM
      DQ(2) = 0.5D0*(-s1 + 2.0D0*s2 - s3) / DENOM
      DQ(3) = 0.5D0*(-s2 + 2.0D0*s3 - s1) / DENOM
      DQ(4) = 3.0D0 * s4 / DENOM
      DQ(5) = 3.0D0 * s5 / DENOM
      DQ(6) = 3.0D0 * s6 / DENOM  
      
      DDQ(1,1) = - (2D0*s1 - s2 - s3)**2 / (4D0*DENOM3)+  
     1  1D0 / DENOM
      DDQ(1,2) = - (2D0*s1 - s2 - s3)*( -s1 + 2D0*s2 - s3 ) / 
     1 (4D0*DENOM3)- 1D0/(2D0*DENOM)
      DDQ(1,3) = - (2D0*s1 - s2 - s3)*( -s2 + 2D0*s3 - s1 ) /
     1 (4D0*DENOM3)- 1D0/(2D0*DENOM)
      DDQ(1,4) = - 3D0*(2D0*s1 - s2 - s3)*s4 / (2D0*DENOM3)
      DDQ(1,5) = - 3D0*(2D0*s1 - s2 - s3)*s5 / (2D0*DENOM3)
      DDQ(1,6) = - 3D0*(2D0*s1 - s2 - s3)*s6 / (2D0*DENOM3)
      
      DDQ(2,1) = - (2D0*s1 - s2 - s3)*( -s1 + 2D0*s2 - s3 ) /
     1 (4D0*DENOM3) - 1D0/(2D0*DENOM)
      DDQ(2,2) = - ( -s1 + 2D0*s2 - s3 )**2 / (4D0*DENOM3)+
     1 1D0 / DENOM
      DDQ(2,3) = - ( -s1 + 2D0*s2 - s3 )*( -s2 + 2D0*s3 - s1 ) / 
     1 (4D0*DENOM3)- 1D0/(2D0*DENOM)
      DDQ(2,4) = - 3D0*( -s1 + 2D0*s2 - s3 )*s4 / (2D0*DENOM3)
      DDQ(2,5) = - 3D0*( -s1 + 2D0*s2 - s3 )*s5 / (2D0*DENOM3)
      DDQ(2,6) = - 3D0*( -s1 + 2D0*s2 - s3 )*s6 / (2D0*DENOM3)
      
      DDQ(3,1) = - (2D0*s1 - s2 - s3)*( -s2 + 2D0*s3 - s1 ) / 
     1 (4D0*DENOM3) - 1D0/(2D0*DENOM)
      DDQ(3,2) = - ( -s1 + 2D0*s2 - s3 )*( -s2 + 2D0*s3 - s1 ) /
     1 (4D0*DENOM3) - 1D0/(2D0*DENOM)
      DDQ(3,3) = - ( -s2 + 2D0*s3 - s1 )**2 / (4D0*DENOM3)+
     1 1D0 / DENOM
      DDQ(3,4) = - 3D0*( -s2 + 2D0*s3 - s1 )*s4 / (2D0*DENOM3)
      DDQ(3,5) = - 3D0*( -s2 + 2D0*s3 - s1 )*s5 / (2D0*DENOM3)
      DDQ(3,6) = - 3D0*( -s2 + 2D0*s3 - s1 )*s6 / (2D0*DENOM3)

      DDQ(4,1) = -3D0*(2D0*s1 - s2 - s3)*s4 / (2D0*DENOM3)
      DDQ(4,2) = -3D0*( -s1 + 2D0*s2 - s3 )*s4 / (2D0*DENOM3)
      DDQ(4,3) = -3D0*( -s2 + 2D0*s3 - s1 )*s4 / (2D0*DENOM3)
      DDQ(4,4) = - 9D0*s4*s4 / (DENOM3)+ 3D0/DENOM
      DDQ(4,5) = - 9D0*s4*s5 / (DENOM3)
      DDQ(4,6) = - 9D0*s4*s6 / (DENOM3)
      
      DDQ(5,1) = -3D0*(2D0*s1 - s2 - s3)*s5 / (2D0*DENOM3)
      DDQ(5,2) = -3D0*( -s1 + 2D0*s2 - s3 )*s5 / (2D0*DENOM3)
      DDQ(5,3) = -3D0*( -s2 + 2D0*s3 - s1 )*s5 / (2D0*DENOM3)
      DDQ(5,4) = - 9D0*s4*s5 / (DENOM3)
      DDQ(5,5) = - 9D0*s5*s5 / (DENOM3)+ 3D0/DENOM
      DDQ(5,6) = - 9D0*s5*s6 / (DENOM3)
      
      DDQ(6,1) = -3D0*(2D0*s1 - s2 - s3)*s6 / (2D0*DENOM3)
      DDQ(6,2) = -3D0*( -s1 + 2D0*s2 - s3 )*s6 / (2D0*DENOM3)
      DDQ(6,3) = -3D0*( -s2 + 2D0*s3 - s1 )*s6 / (2D0*DENOM3)
      DDQ(6,4) = - 9D0*s4*s6 / (DENOM3)
      DDQ(6,5) = - 9D0*s5*s6 / (DENOM3)
      DDQ(6,6) = - 9D0*s6*s6 / (DENOM3)+ 3D0/DENOM
      RETURN
      END
    
C==================================================================
C  PRCREEP FUNCTIONS FOR JACOBIAN
C==================================================================  
      SUBROUTINE PRIMARYFUNCTIONS(DOUTER,STRESS,NTENS,EFFSIGMA,
     1 EFFPRIMEPS,DENOM,FP,DFP)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 EFFSIGMA, EFFPRIMEPS, DTIME, tiny
      REAL*8 STRESS(*),DFP(6)
      REAL*8 Ep, MEXP, NP, FP, arg, DENOM, s1,s2,s3,s4,s5,s6

      PARAMETER (tiny = 1.0D-20)

      FP = Ep * ( (EFFSIGMA / Ep)**MEXP - EFFPRIMEPS )
      IF (FP .LT. 0.0D0) FP = 0.0D0
      
      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF      
      DOUTER = Ep*MEXP*((EFFSIGMA)**(MEXP-1))/ ((Ep)**MEXP)
      DFP(1) = DOUTER* 0.5D0*(2.0D0*s1 - s2 - s3) / DENOM
      DFP(2) = DOUTER* 0.5D0*(-s1 + 2.0D0*s2 - s3) / DENOM
      DFP(3) = DOUTER* 0.5D0*(-s2 + 2.0D0*s3 - s1) / DENOM
      DFP(4) = DOUTER* 3.0D0 * s4 / DENOM
      DFP(5) = DOUTER* 3.0D0 * s5 / DENOM
      DFP(6) = DOUTER* 3.0D0 * s6 / DENOM  
      RETURN
      END
C==================================================================
C  SECCREEP FUNCTIONS FOR JACOBIAN
C==================================================================  
      SUBROUTINE SECONDARYFUNCTIONS(DOUTER,STRESS,NTENS,EFFSIGMA,
     1 DENOM,FS,DFS)
      IMPLICIT NONE
      INTEGER NTENS, K
      REAL*8 EFFSIGMA, EFFPRIMEPS, DTIME, tiny
      REAL*8 STRESS(*),DFS(6)
      REAL*8 p0, NEXP, NS, arg, DENOM, s1,s2,s3,s4,s5,s6, tiny,FS

      PARAMETER (p0 = 1.0D6, NEXP = 1.0D0, NS = 5.0D11, tiny = 1.0D-20)
      
      FS = p0*(EFFSIGMA/p0)**NEXP
      
      s1 = STRESS(1)
      s2 = STRESS(2)
      s3 = STRESS(3)
      IF (NTENS .GE. 4) THEN
         s4 = STRESS(4)
      ELSE
         s4 = 0.0D0
      END IF
      IF (NTENS .GE. 5) THEN
         s5 = STRESS(5)
      ELSE
         s5 = 0.0D0
      END IF
      IF (NTENS .GE. 6) THEN
         s6 = STRESS(6)
      ELSE
         s6 = 0.0D0
      END IF
      DOUTER = p0*NEXP*((EFFSIGMA)**(NEXP-1))/ ((P0)**NEXP)
      DFS(1) = DOUTER* 0.5D0*(2.0D0*s1 - s2 - s3) / DENOM
      DFS(2) = DOUTER* 0.5D0*(-s1 + 2.0D0*s2 - s3) / DENOM
      DFS(3) = DOUTER* 0.5D0*(-s2 + 2.0D0*s3 - s1) / DENOM
      DFS(4) = DOUTER* 3.0D0 * s4 / DENOM
      DFS(5) = DOUTER* 3.0D0 * s5 / DENOM
      DFS(6) = DOUTER* 3.0D0 * s6 / DENOM  
      RETURN
      END