!**********************************************************************
      MODULE UNWAM

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      IMPLICIT NONE

      LOGICAL :: APPLY_DXP_CORR = .TRUE.
!!!! JB:
!!!! USE_EXACT_FORMULA_SPHERICAL_AREA = .TRUE. DOES NOT SEEM TO WORK !!!
!!! I need to invertigate further
      LOGICAL :: USE_EXACT_FORMULA_SPHERICAL_AREA = .FALSE.
      LOGICAL :: LNANINFCHK = .TRUE.

!     THE FOLLOWING FLAG ARE PART OF THE INPUT NAMELIST
!     SEE *MPUSERIN*
      LOGICAL :: LIMPLICIT
      LOGICAL :: REFRACTION_IMPL
      LOGICAL :: FREQ_SHIFT_IMPL
      LOGICAL :: SOURCE_IMPL
      LOGICAL :: LNONL
      LOGICAL :: BLOCK_GAUSS_SEIDEL
      LOGICAL :: LLIMT
      LOGICAL :: L_SOLVER_NORM
      LOGICAL :: LCHKCONV

      INTEGER(KIND=JWIM), ALLOCATABLE :: JA_IE(:,:,:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: IA(:), JA(:)
      INTEGER(KIND=JWIM), ALLOCATABLE :: POSI(:,:), I_DIAG(:)
      INTEGER(KIND=JWIM) :: NNZ

      REAL(KIND=JWRU), PARAMETER :: YPMAX=87.5_JWRU

      REAL(KIND=JWRU), ALLOCATABLE :: ASPAR_JAC(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: B_JAC(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: CAD_THE(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: CAS_SIG(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: U_JACOBI(:,:,:)
      REAL(KIND=JWRU), ALLOCATABLE :: DS_INCR(:)
      REAL(KIND=JWRU) :: DT4D, DT4F, DDIR
      !
      ! Solver thresholds
      !
      INTEGER(KIND=JWRU) :: maxiter = 100
      REAL(KIND=JWRU) :: solverthr
      REAL(KIND=JWRU) :: pmin

      CONTAINS 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PROPAG_UNWAM(FL1, FL3)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        ------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!       Driver routine for the fluctuation splitting schemes 
!       Parallelization over freq. and directions unsing OpenMP untill dom. decomp. is not ready 
!       References: Roland, 2008, Roland et al. 2006, Zanke et al. 2006, Csik et al. 2002 
!                   Hubbard & Roe, 2000, Struis et al. 1998 

!     Externals.  EXPLICIT_N_SCHEME, EXPLICIT_LFPSI_SCHEME, EXPLICIT_PSI_SCHEME
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
!      USE yowpd
      USE YOWTEST  , ONLY : IU06
      USE yownodepool, ONLY : NP, iplg, ipgl
      USE yowdatapool, only: myrank
      USE yowpd, only: MNP=>npa
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP

      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(INOUT)  :: FL1(NINF-1:NSUP,NANG,NFRE)
      REAL(KIND=JWRB), INTENT(OUT)    :: FL3(NINF-1:NSUP,NANG,NFRE)
 
      INTEGER(KIND=JWIM) :: IS, ID, IP

#ifdef DEBUG
      CALL COHERENCY_ERROR_3D(FL1, "testing the function FL1")
      CALL COHERENCY_ERROR_3D(FL3, "testing the function FL3")
#endif
      IF (LIMPLICIT) THEN
        CALL IMPLICIT_N_SCHEME_BLOCK(FL1, FL3)
        Print *, 'Leaving IMPLICIT_N_SCHEME_BLOCK'
      ELSE
        IF (LVECTOR) THEN 
          CALL EXPLICIT_N_SCHEME_VECTOR(FL1,FL3)
        ELSE
!$OMP     PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(ID,IS)
          DO ID = 1, NANG
            DO IS = 1, MSC
              CALL EXPLICIT_N_SCHEME(IS,ID,FL1(:,ID,IS),FL3(:,ID,IS))
             !CALL EXPLICIT_PSI_SCHEME(IS,ID,FL1(:,ID,IS),FL3(:,ID,IS))
             !CALL EXPLICIT_LF_SCHEME(IS,ID,FL1(:,ID,IS),FL3(:,ID,IS))
            END DO
          END DO
!$OMP     END PARALLEL DO
        ENDIF
      ENDIF

      END SUBROUTINE PROPAG_UNWAM
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FLUCTCFL(IS, ID, DTMAX)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  EXPLICIT_N_SCHEME, EXPLICIT_LFPSI_SCHEME, EXPLICIT_PSI_SCHEME
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
      USE yowpd, only: MNE=>ne, INE, MNP=>npa

!     --------------------------------------------------------------
      IMPLICIT NONE

      REAL(KIND=JWRU) :: K(3,MNE)
      REAL(KIND=JWRB), INTENT(OUT) :: DTMAX
      INTEGER(KIND=JWIM) :: IS, ID
      INTEGER(KIND=JWIM) :: I, J, I1, I2, I3
      INTEGER(KIND=JWIM) :: IP, IE, POS
      REAL(KIND=JWRU) :: KSUM, KMAX, LAMBDA(2)
      REAL(KIND=JWRU) :: DTMAX_EXP, DTMAX_GLOBAL_EXP
      REAL(KIND=JWRU) :: REST, C(2,MNP)
      REAL(KIND=JWRB) :: DIFRU, USOC

      DTMAX_GLOBAL_EXP = 10.E14_JWRU
      CALL CADVXY(IS,ID,C)
      DO IE = 1, MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
        LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
        K(1,IE)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
        K(2,IE)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
        K(3,IE)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
      END DO
      J = 0
      DO IP = 1, MNP
        KSUM = 0.0_JWRU
        KMAX = 0.0_JWRU
        DO I = 1, CCON(IP)
          J = J + 1
          IE    = IE_CELL(J)
          POS   = POS_CELL(J)
          KSUM  = KSUM + MAX(K(POS,IE),0.0_JWRU)
          IF ( ABS(K(POS,IE)) > KMAX ) KMAX = ABS(K(POS,IE))
        END DO
        IF (KSUM > 0.0_JWRU) THEN
          DTMAX_EXP = SI(IP)/KSUM
        ELSE
          DTMAX_EXP = 10.E14_JWRU
        END IF
        IF (DTMAX_GLOBAL_EXP>DTMAX_EXP) DTMAX_GLOBAL_EXP = DTMAX_EXP
      END DO
      DTMAX = REAL(DTMAX_GLOBAL_EXP)
      REST  = ABS(MOD(DT4A/DTMAX_GLOBAL_EXP,1.0_JWRU))
      IF (REST > THR .AND. REST < 0.5_JWRU) THEN
        ITER_EXP(IS,ID) = ABS(NINT(DT4A/DTMAX_GLOBAL_EXP)) + 1
      ELSE
        ITER_EXP(IS,ID) = ABS(NINT(DT4A/DTMAX_GLOBAL_EXP))
      END IF
      END SUBROUTINE FLUCTCFL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME  ( IS, ID, UOLD, UNEW )

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Narrow Stencil (N) scheme using contour integration method (Csik et al.) for flux conservation 
!                 1st order time and space, monot1.d0, conservative, quasi-positive
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------
         USE MPL_MPIF
         USE YOWUNPOOL
         use yowdatapool, only: myrank
         USE yowpd, only: MNE=>ne, MNP=>npa, INE, exchange, comm
         USE YOWMPP   , ONLY : NINF, NSUP, IRANK
         USE YOWSTAT  , ONLY : IDELPRO
         USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

!     -------------------------------------------------------------------------

         IMPLICIT NONE

         INTEGER(KIND=JWIM), INTENT(IN)    :: IS,ID
         REAL(KIND=JWRB), INTENT(IN)   :: UOLD(NINF-1:NSUP)
         REAL(KIND=JWRB), INTENT(OUT)  :: UNEW(NINF-1:NSUP)
!
! local integer
!
         INTEGER(KIND=JWIM) :: IP, IE, IT
         INTEGER(KIND=JWIM) :: I1, I2, I3
         INTEGER(KIND=JWIM) :: NI(3),K
!
! local double
!
         REAL(KIND=JWRU) :: FT
         REAL(KIND=JWRU) :: UTILDE
         REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP, DTMAX_EXP, DTMAX_GLOBAL_EXP_LOC
         REAL(KIND=JWRU) :: REST
         REAL(KIND=JWRU) :: LAMBDA(2), DT4AI
         REAL(KIND=JWRU) :: FL11,FL12,FL21,FL22,FL31,FL32
         REAL(KIND=JWRU) :: KTMP(3)
         REAL(KIND=JWRU) :: KKSUM(MNP), ST(MNP), N(MNE), U3(3), ST3(3)
         REAL(KIND=JWRU) :: C(2,MNP), DTSI(MNP), CFLXY
         REAL(KIND=JWRU) :: FLALL(3,MNE), U(MNP)
         REAL(KIND=JWRU) :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(KIND=JWRU) :: KELEM(3,MNE)
!
! local parameter
!
         REAL(KIND=JWRU) :: TMP
         REAL(KIND=JWRB) :: ZHOOK_HANDLE

         INTEGER(KIND=JWIM) :: ierr

!     -------------------------------------------------------------------------

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('EXPLICIT_N_SCHEME',0,ZHOOK_HANDLE)
#endif

!
! set time step
!
         DT4A = REAL(IDELPRO,JWRU)
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
         CALL CADVXY(IS,ID,C)

!???JB is the call needed ???
!!!!!         call EXCHANGE(C(1,:))
!???JB is the call needed ???
!!!!         call EXCHANGE(C(2,:))

!         if(IS==1 .and. ID==1) then
!           write(DBG%FHNDL,*) "C"
!           do ip=1, np_global
!             if(ipgl(ip) /= 0) then
!               write(DBG%FHNDL,*) ip, C(:,ipgl(ip))
!             endif
!             if(ghostgl(ip) /= 0) then
!               write(DBG%FHNDL,*) ip, C(:,np+ghostgl(ip))
!             endif
!           end do
!         endif

!         DO IP = 1, MNP
!           C(:,IP) = C(:,IP) * IOBPD(ID,IP)
!         END DO
!
!        Calculate K-Values and contour based quantities ...
!

         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)

            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
!            WRITE(DBG%FHNDL,*) 'LAMBDA', LAMBDA
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
!            WRITE(DBG%FHNDL,'(A10,6D15.4)') 'TEST IEN', IEN(:,IE)
            KTMP  = KELEM(:,IE)
            TMP   = SUM(MIN(0.0_JWRU,KTMP))
            N(IE) = - 1.0_JWRU/MIN(-THR,TMP)
            KELEM(:,IE) = MAX(0.0_JWRU,KTMP)
!            WRITE(DBG%FHNDL,'(A20,3I10,3F15.10)') 'TESTING KELEM' ,IS, ID, IE, KELEM(:,IE)
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = 2.0_JWRU*FL11+FL12
            FL112 = 2.0_JWRU*FL12+FL11
            FL211 = 2.0_JWRU*FL21+FL22
            FL212 = 2.0_JWRU*FL22+FL21
            FL311 = 2.0_JWRU*FL31+FL32
            FL312 = 2.0_JWRU*FL32+FL31
            FLALL(1,IE) = (FL311 + FL212) * ONESIXTH + KELEM(1,IE)
            FLALL(2,IE) = (FL111 + FL312) * ONESIXTH + KELEM(2,IE)
            FLALL(3,IE) = (FL211 + FL112) * ONESIXTH + KELEM(3,IE)
         END DO


         IF (LCALC) THEN ! If the current field or water level --> new CFL 
           KKSUM = 0.0_JWRU
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
           END DO
           !write(DBG%FHNDL, *) "KKSUM", KKSUM

! THOMAS find global dt max/min
           DTMAX_GLOBAL_EXP = LARGE
           DTMAX_GLOBAL_EXP_LOC = LARGE
           DO IP = 1, MNP
             DTMAX_EXP = SI(IP)/MAX(SMALL, KKSUM(IP))
             IF (LCFL) THEN
               CFLCXY(1,IP) = MAX(CFLCXY(1,IP), C(1,IP))
               CFLCXY(2,IP) = MAX(CFLCXY(2,IP), C(2,IP))
               CFLCXY(3,IP) = MAX(CFLCXY(3,IP), DT4A/DTMAX_EXP)
             END IF
             DTMAX_GLOBAL_EXP_LOC = MIN(DTMAX_GLOBAL_EXP_LOC,DTMAX_EXP)

!             if(DTMAX_GLOBAL_EXP < 0.1) then
!               WRITE(DBG%FHNDL,*) 'LOCAL',IP,SI(IP),KKSUM(IP),DTMAX_EXP
!               WRITE(DBG%FHNDL,*) DTMAX_GLOBAL_EXP
!             endif
           END DO

           CALL MPI_ALLREDUCE(DTMAX_GLOBAL_EXP_LOC, DTMAX_GLOBAL_EXP,    &
     &                        1, MPI_REAL8, MPI_MIN, COMM, IERR)


           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,1.0_JWRU))
           IF (REST .LT. THR) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           ELSE IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           END IF
!           WRITE(DBG%FHNDL,'(2I5,3F15.8,I10)') IS, ID, DT4A, CFLXY, 
!     &                            DTMAX_GLOBAL_EXP, ITER_EXP(IS,ID)
         END IF

         DT4AI    = DT4A/REAL(ITER_EXP(IS,ID),JWRU)
         DTSI(:)  = DT4AI/SI(:)

         U(1:MNP) = REAL(UOLD(1:MNP),JWRU)
         CALL EXCHANGE(U)
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
         DO IT = 1, ITER_EXP(IS,ID)
           ST = 0.0_JWRU ! Init. ... only used over the residual nodes see IP loop
           DO IE = 1, MNE
             NI     = INE(:,IE)
             U3     = U(NI) 
             UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3*IOBPD(ID,NI)))
             !UTILDE = N(IE)*(DOT_PRODUCT(FLALL(:,IE),U3))
             ST(NI) = ST(NI)*IOBPD(ID,NI)+KELEM(:,IE)*(U3 - UTILDE) 
             !ST(NI) = ST(NI)+KELEM(:,IE)*(U3 - UTILDE)
           END DO
           U = MAX(0.0_JWRU,U-DTSI*ST*REAL(IOBWB,JWRU))!*REAL(IOBPD(ID,:),JWRU)
           CALL EXCHANGE(U)
         END DO  ! ----> End Iteration

         UNEW(NINF-1)=0.0_JWRB
         UNEW(1:MNP) = REAL(U(1:MNP),JWRB)


#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('EXPLICIT_N_SCHEME',1,ZHOOK_HANDLE)
#endif

      END SUBROUTINE EXPLICIT_N_SCHEME

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXCHANGE_FOR_FL1_FL3_SL (AC)

      USE YOWUNPOOL, only : MSC, MDC
      USE yowpd, only: MNP=>npa, INE, exchange

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP

      IMPLICIT NONE

      REAL(KIND=JWRB), intent(inout) :: AC(NINF-1:NSUP,NANG,NFRE)
      integer(KIND=JWIM) :: IS, ID, IP
      REAL(KIND=JWRU) :: ACexch(MNP)
      REAL(KIND=JWRB)    :: ACtest(MNP)

      DO is = 1 , msc
        DO id = 1 , mdc
          DO ip = 1 , mnp
            ACexch(ip) = REAL(AC(ip,id,is),JWRU)
          ENDDO
          CALL exchange(ACexch)
          DO ip = 1 , mnp
            AC(ip,id,is) = REAL(ACexch(ip),JWRB)
            ACtest(ip) = REAL(ACexch(ip),JWRB)
          ENDDO
        ENDDO
      ENDDO
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'After IS/ID loop'
      CALL COHERENCY_ERROR_3D(AC, "testing AC")
      FLUSH(740+MyRankGlobal)
#endif
      END SUBROUTINE EXCHANGE_FOR_FL1_FL3_SL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_PSI_SCHEME  ( IS, ID, UOLD, UNEW )

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  PSI (positive streamline invariant) CRD (contour approach) scheme accoridng to Struijs et al. Csik et al. and Roland, 2008
!                 2nd order in space on regular meshes in steady state, better than 1st order in space on irregular meshes. Best scheme according 
!                 to me ... 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL
      USE yowpd, only: MNE=>ne, MNP=>npa, INE, exchange
      USE YOWMPP   , ONLY : NINF, NSUP
      USE YOWSTAT  , ONLY : IDELPRO
      use yowexchangemodule

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN)    :: IS,ID
      REAL(KIND=JWRB), INTENT(IN)   :: UOLD(NINF-1:NSUP)
      REAL(KIND=JWRB), INTENT(OUT)  :: UNEW(NINF-1:NSUP)
!
! local integer
!
      INTEGER(KIND=JWIM) :: IP, IE, IT
      INTEGER(KIND=JWIM) :: I1, I2, I3
      INTEGER(KIND=JWIM) :: NI(3),K
!
! local double
!
      REAL(KIND=JWRU) :: FT
      REAL(KIND=JWRU) :: UTILDE
      REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP, DTMAX_EXP
      REAL(KIND=JWRU) :: REST
      REAL(KIND=JWRU) :: LAMBDA(2), DT4AI
      REAL(KIND=JWRU) :: FL11,FL12,FL21,FL22,FL31,FL32
      REAL(KIND=JWRU) :: KTMP(3)
      REAL(KIND=JWRU) :: KKSUM(MNP), ST(MNP), N(MNE), U3(3), ST3(3)
      REAL(KIND=JWRU) :: C(2,MNP), DTSI(MNP), CFLXY
      REAL(KIND=JWRU) :: FLALL(3,MNE), U(MNP)
      REAL(KIND=JWRU) :: FL111, FL112, FL211, FL212, FL311, FL312
      REAL(KIND=JWRU) :: KELEM(3,MNE)
      REAL(KIND=JWRU) :: THETA_L(3)
      REAL(KIND=JWRU) :: BET1(3), BETAHAT(3)
!
! local parameter
!
      REAL(KIND=JWRU) :: TMP
!
! set time step
!
      DT4A = REAL(IDELPRO,JWRU)
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
      CALL CADVXY(IS,ID,C)
!
!        Calculate K-Values and contour based quantities ...
!
      DO IE = 1, MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
        LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
        KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
        KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
        KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
        KTMP  = KELEM(:,IE)
        TMP   = SUM(MIN(0.0_JWRU,KTMP))
        N(IE) = - 1.0_JWRU/MIN(-THR,TMP)
        KELEM(:,IE) = MAX(0.0_JWRU,KTMP)
        FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
        FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
        FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
        FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
        FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
        FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
        FL111 = 2.0_JWRU*FL11+FL12
        FL112 = 2.0_JWRU*FL12+FL11
        FL211 = 2.0_JWRU*FL21+FL22
        FL212 = 2.0_JWRU*FL22+FL21
        FL311 = 2.0_JWRU*FL31+FL32
        FL312 = 2.0_JWRU*FL32+FL31
        FLALL(1,IE) = FL311 + FL212
        FLALL(2,IE) = FL111 + FL312
        FLALL(3,IE) = FL211 + FL112
      END DO

      IF (LCALC) THEN ! If the current field or water level --> new CFL 
        KKSUM = 0.0_JWRU
        DO IE = 1, MNE
          NI = INE(:,IE)
          KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
        END DO
        DTMAX_GLOBAL_EXP = LARGE
        DO IP = 1, MNP
          DTMAX_EXP = SI(IP)/MAX(SMALL,KKSUM(IP))
          DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
        END DO
        CFLXY = DT4A/DTMAX_GLOBAL_EXP
        REST  = ABS(MOD(CFLXY,1.0_JWRU))
        IF (REST .LT. THR) THEN
          ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
        ELSE IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
          ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
        ELSE
          ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
        END IF
      END IF

      DT4AI    = DT4A/REAL(ITER_EXP(IS,ID),JWRU)
      DTSI(:)  = DT4AI/SI(:)

      U(1:MNP) = REAL(UOLD(1:MNP),JWRU)
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
      DO IT = 1, ITER_EXP(IS,ID)
        ST = 0.0_JWRU
        DO IE = 1, MNE
          NI   = INE(:,IE)
          FT   = -ONESIXTH*DOT_PRODUCT(U(NI)*IOBPD(ID,NI),FLALL(:,IE))
          UTILDE = N(IE)*(DOT_PRODUCT(KELEM(:,IE),U(NI)*IOBPD(ID,NI))-FT)
          THETA_L(:) = KELEM(:,IE) * (U(NI) - UTILDE)
          IF (ABS(FT) .GT. 0.0_JWRU) THEN
            BET1(:) = THETA_L(:)/FT
            IF (ANY( BET1 .LT. 0.0_JWRU) ) THEN
              BETAHAT(1)= BET1(1) + 0.5_JWRU * BET1(2)
              BETAHAT(2)= BET1(2) + 0.5_JWRU * BET1(3)
              BETAHAT(3)= BET1(3) + 0.5_JWRU * BET1(1)
              BET1(1)= MAX(0.0_JWRU,MIN(BETAHAT(1),1.0_JWRU-BETAHAT(2),1.0_JWRU))
              BET1(2)= MAX(0.0_JWRU,MIN(BETAHAT(2),1.0_JWRU-BETAHAT(3),1.0_JWRU))
              BET1(3)= MAX(0.0_JWRU,MIN(BETAHAT(3),1.0_JWRU-BETAHAT(1),1.0_JWRU))
              THETA_L(:) = FT * BET1
            END IF
          ELSE
            THETA_L(:) = 0.0_JWRU
          END IF
          ST(NI) = ST(NI)*IOBPD(ID,NI) + THETA_L ! the 2nd term are the theta values of each node ...
        END DO
        U = MAX(0.0_JWRU,U-DTSI*ST*REAL(IOBWB,JWRU))!*REAL(IOBPD(ID,:),JWRU)
      END DO  ! ----> End Iteration
      CALL EXCHANGE(U)
      UNEW(NINF-1)=0.0_JWRB
      UNEW(1:MNP) = REAL(U(1:MNP),JWRB)
      END SUBROUTINE EXPLICIT_PSI_SCHEME

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_LFPSI_SCHEME(IS,ID,UOLD,UNEW)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Law-Friedrich Flux Corrected Transport blended PSI schemes based on Hubbard & Roe (2000), Roland (2008)
!                 True 2nd order space/time scheme on narrow stencil. A bit expensive, needs GSE correction if mesh fine and directioan resolution 
!                 coarse 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

         USE YOWUNPOOL
         USE yowpd, only: MNE=>ne, MNP=>npa, INE, exchange
         USE YOWMPP   , ONLY : NINF, NSUP
         USE YOWSTAT  , ONLY : IDELPRO

         IMPLICIT NONE

         INTEGER(KIND=JWIM), INTENT(IN)    :: IS,ID
         REAL(KIND=JWRB), INTENT(IN)   :: UOLD(NINF-1:NSUP)
         REAL(KIND=JWRB), INTENT(OUT)  :: UNEW(NINF-1:NSUP)
!
! local integer
!
         INTEGER(KIND=JWIM) :: IP, IE, IT
         INTEGER(KIND=JWIM) :: I1, I2, I3, K
         INTEGER(KIND=JWIM) :: NI(3)
!
! local double
!
         REAL(KIND=JWRU) :: FT
         REAL(KIND=JWRU) :: UTILDE

         REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP, DTMAX_EXP

         REAL(KIND=JWRU) :: REST
         REAL(KIND=JWRU) :: TMP(3), TMP1

         REAL(KIND=JWRU) :: LAMBDA(2), DT4AI
         REAL(KIND=JWRU) :: BET1(3), BETAHAT(3), BL

         REAL(KIND=JWRU) :: FL11,FL12,FL21,FL22,FL31,FL32

         REAL(KIND=JWRU) :: THETA_L(3,MNE),THETA_H(3),THETA_ACE(3,MNE),UTMP(3)
         REAL(KIND=JWRU) :: WII(2,MNP), UL(MNP), USTARI(2,MNP), KTMP(3)

         REAL(KIND=JWRU) :: KKSUM(MNP), ST(MNP)
         REAL(KIND=JWRU) :: PM(MNP), PP(MNP), UIM(MNP), UIP(MNP)

         REAL(KIND=JWRU) :: C(2,MNP), U(MNP), DTSI(MNP), CFLXY, N(MNE)
         REAL(KIND=JWRU) :: FL111, FL112, FL211, FL212, FL311, FL312
         REAL(KIND=JWRU) :: KELEM(3,MNE), FLALL(3,MNE)
!
! local parameter
!
         BL = 0.0_JWRU
!
! set time step
!
         DT4A = REAL(IDELPRO,JWRU)
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
         CALL CADVXY(IS,ID,C)
!
!        Calculate K-Values and contour based quantities ...
!
         DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
!            WRITE(DBG%FHNDL,*) 'LAMBDA', LAMBDA
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
!            WRITE(DBG%FHNDL,'(A10,6D15.4)') 'TEST IEN', IEN(:,IE)
            KTMP  = KELEM(:,IE)
            TMP1   = SUM(MIN(0.0_JWRU,KTMP))
            N(IE) = - 1.0_JWRU/MIN(-THR,TMP1)
            KELEM(:,IE) = MAX(0.0_JWRU,KTMP)
!            WRITE(DBG%FHNDL,'(A20,3I10,3F15.10)') 'TESTING KELEM', 
!     &                                    IS, ID, IE, KELEM(:,IE)
            FL11  = C(1,I2) * IEN(1,IE) + C(2,I2) * IEN(2,IE)
            FL12  = C(1,I3) * IEN(1,IE) + C(2,I3) * IEN(2,IE)
            FL21  = C(1,I3) * IEN(3,IE) + C(2,I3) * IEN(4,IE)
            FL22  = C(1,I1) * IEN(3,IE) + C(2,I1) * IEN(4,IE)
            FL31  = C(1,I1) * IEN(5,IE) + C(2,I1) * IEN(6,IE)
            FL32  = C(1,I2) * IEN(5,IE) + C(2,I2) * IEN(6,IE)
            FL111 = 2.0_JWRU*FL11+FL12
            FL112 = 2.0_JWRU*FL12+FL11
            FL211 = 2.0_JWRU*FL21+FL22
            FL212 = 2.0_JWRU*FL22+FL21
            FL311 = 2.0_JWRU*FL31+FL32
            FL312 = 2.0_JWRU*FL32+FL31
            FLALL(1,IE) = FL311 + FL212
            FLALL(2,IE) = FL111 + FL312
            FLALL(3,IE) = FL211 + FL112
         END DO

         IF (LCALC) THEN ! If the current field or water level --> new CFL 
           KKSUM = 0.0_JWRU
           DO IE = 1, MNE
             NI = INE(:,IE)
             KKSUM(NI) = KKSUM(NI) + KELEM(:,IE)
           END DO
           DTMAX_GLOBAL_EXP = LARGE
           DO IP = 1, MNP
             DTMAX_EXP = SI(IP)/MAX(SMALL,KKSUM(IP))
             DTMAX_GLOBAL_EXP = MIN ( DTMAX_GLOBAL_EXP, DTMAX_EXP)
             !WRITE(DBG%FHNDL,*)'LOCAL DT',IP,SI(IP),KKSUM(IP),DTMAX_EXP
           END DO
           CFLXY = DT4A/DTMAX_GLOBAL_EXP
           REST  = ABS(MOD(CFLXY,1.0_JWRU))
           IF (REST .LT. THR) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           ELSE IF (REST .GT. THR .AND. REST .LT. 0.5_JWRU) THEN
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY)) + 1
           ELSE
             ITER_EXP(IS,ID) = ABS(NINT(CFLXY))
           END IF
!           WRITE(DBG%FHNDL,'(2I5,3F15.8,I10)') IS, ID, DT4A, CFLXY, 
!     &                            DTMAX_GLOBAL_EXP, ITER_EXP(IS,ID)
         END IF

         DT4AI    = REAL(DT4A,JWRU)/REAL(ITER_EXP(IS,ID),JWRU)
         DTSI(:)  = DT4AI/SI(:)

         U(1:MNP)  = REAL(UOLD(1:MNP),JWRU)
!
!  Loop over all sub time steps, all quantities in this loop depend on the solution U itself !!!
!
         DO IT = 1, ITER_EXP(IS,ID)
!
! Element loop
!
           ST = 0.0_JWRU
           PM = 0.0_JWRU
           PP = 0.0_JWRU
           DO IE = 1, MNE
             NI   = INE(:,IE)
        FT   = -ONESIXTH*DOT_PRODUCT(U(NI)*IOBPD(ID,NI),FLALL(:,IE))
        UTILDE = N(IE)*(DOT_PRODUCT(KELEM(:,IE),U(NI)*IOBPD(ID,NI))-FT)
             THETA_L(:,IE) = KELEM(:,IE) * (U(NI) - UTILDE)
             IF (ABS(FT) .GT. 0.0_JWRU) THEN
               BET1(:) = THETA_L(:,IE)/FT
               IF (ANY( BET1 .LT. 0.0_JWRU) ) THEN
                 BETAHAT(1)= BET1(1) + 0.5_JWRU * BET1(2)
                 BETAHAT(2)= BET1(2) + 0.5_JWRU * BET1(3)
                 BETAHAT(3)= BET1(3) + 0.5_JWRU * BET1(1)
                 BET1(1)= MAX(0.0_JWRU,MIN(BETAHAT(1),1.0_JWRU-BETAHAT(2),1.0_JWRU))
                 BET1(2)= MAX(0.0_JWRU,MIN(BETAHAT(2),1.0_JWRU-BETAHAT(3),1.0_JWRU))
                 BET1(3)= MAX(0.0_JWRU,MIN(BETAHAT(3),1.0_JWRU-BETAHAT(1),1.0_JWRU))
                 THETA_L(:,IE) = FT * BET1
               END IF
             ELSE
               THETA_L(:,IE) = 0.0_JWRU
             END IF
              ST(NI)  = ST(NI) + THETA_L(:,IE)
              !THETA_H = (ONETHIRD+DT4AI/(2.*TRIA(IE)) * KELEM(:,IE) )*FT ! LAX
              THETA_H = (ONETHIRD+TWOTHIRD*KELEM(:,IE)/SUM(MAX(0.0_JWRU,KELEM(:,IE))))*FT  ! CENTRAL
              THETA_ACE(:,IE) = THETA_H-THETA_L(:,IE)
              PP(NI) =  PP(NI) + MAX( 0.0_JWRU, -THETA_ACE(:,IE)) * DTSI(NI)
              PM(NI) =  PM(NI) + MIN( 0.0_JWRU, -THETA_ACE(:,IE)) * DTSI(NI)
            END DO
!
            !U = MAX(0.0_JWRU,U-DTSI*ST)*REAL(IOBPD(ID,:),JWRU)
            UL = MAX(0.0_JWRU,U-DTSI*ST*REAL(IOBWB,JWRU))

            USTARI(1,:) = MAX(UL,U) 
            USTARI(2,:) = MIN(UL,U) 

            UIP = 0.
            UIM = 0.
            DO IE = 1, MNE
              NI = INE(:,IE)
              UIP(NI) = MAX (UIP(NI), MAXVAL( USTARI(1,NI) ))
              UIM(NI) = MIN (UIM(NI), MINVAL( USTARI(2,NI) ))
            END DO

            WII(1,:) = MIN(1.0_JWRU,(UIP-UL)/MAX( REAL(SMALL,JWRU),PP))
            WII(2,:) = MIN(1.0_JWRU,(UIM-UL)/MIN(-REAL(SMALL,JWRU),PM))

            ST = 0.0_JWRU
            DO IE = 1, MNE
               I1 = INE(1,IE)
               I2 = INE(2,IE)
               I3 = INE(3,IE)
               IF (THETA_ACE(1,IE) .LT. 0.0_JWRU) THEN
                 TMP(1) = WII(1,I1)
               ELSE
                 TMP(1) = WII(2,I1)
               END IF
               IF (THETA_ACE(2,IE) .LT. 0.0_JWRU) THEN
                 TMP(2) = WII(1,I2)
               ELSE
                 TMP(2) = WII(2,I2)
               END IF
               IF (THETA_ACE(3,IE) .LT. 0.0_JWRU) THEN
                 TMP(3) = WII(1,I3)
               ELSE
                 TMP(3) = WII(2,I3)
               END IF
               TMP1 = MINVAL(TMP)
               ST(I1) = ST(I1) + THETA_ACE(1,IE) * TMP1! * (1.d0 - BL) + BL * THETA_L(1,IE)
               ST(I2) = ST(I2) + THETA_ACE(2,IE) * TMP1! * (1.d0 - BL) + BL * THETA_L(2,IE)
               ST(I3) = ST(I3) + THETA_ACE(3,IE) * TMP1! * (1.d0 - BL) + BL * THETA_L(3,IE)
            END DO
            U = MAX(0.0_JWRU,UL-DTSI*ST*REAL(IOBWB,JWRU))!*REAL(IOBPD(ID,:),JWRU)
         END DO  ! ----> End Iteration
         CALL EXCHANGE(U)
         UNEW(NINF-1)=0.0_JWRB
         UNEW(1:MNP) = REAL(U(1:MNP),JWRB)

!         IF (LADVTEST) THEN
!           WRITE(4001)  RTIME
!           WRITE(4001) (SNGL(C(1,IP)), SNGL(C(2,IP)), AC2(IP,1,1), IP = 1, MNP)
!           CALL CHECKCONS(U,SUMAC2)
!           IF (MINVAL(SNGL(U)) .LT. MINTEST) MINTEST = MINVAL(SNGL(U))
!           WRITE (*,*) 'VOLUMES AT T0, T1 and T2',SUMACt0, SUMAC1, SUMAC2, MINTEST
!           WRITE (*,*) 'VOLUME ERROR: TOTAL and ACTUAL', 100.0-((SUMACt0-SUMAC2)/SUMACt0)*&
!      &                                       100.0, 100.0-((SUMAC1-SUMAC2)/SUMAC1)*100.0
!         END IF

      END SUBROUTINE EXPLICIT_LFPSI_SCHEME

!**********************************************************************
!*
!**********************************************************************
      SUBROUTINE CADVXY(IS,ID,C)
      USE YOWUNPOOL, ONLY : LADVTEST, LCUR, LSPHE, CG, CURTXY, DEGRAD, REARTH
      USE YOWPD,     ONLY : MNP=>npa, XP=>x, YP=>y
      USE YOWFRED  , ONLY : COSTH, SINTH
      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN)  :: IS, ID
      REAL(KIND=JWRU), INTENT(OUT)  :: C(2,MNP)

      INTEGER(KIND=JWIM) :: IP

      REAL(KIND=JWRB) :: ZHOOK_HANDLE
      REAL(KIND=JWRB) :: DGRTH, DGRTHM1

! ----------------------------------------------------------------------
#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CADVXY',0,ZHOOK_HANDLE)
#endif

      DGRTH = DEGRAD*REARTH
      DGRTHM1 = 1.0_JWRU/(DGRTH)
      IF (LADVTEST) THEN
        C(1,:) =   YP
        C(2,:) = - XP
      ELSE
        DO IP = 1, MNP
          IF (LCUR) THEN
            C(1,IP) = CG(IP,IS)*SINTH(ID)+CURTXY(IP,1)
            C(2,IP) = CG(IP,IS)*COSTH(ID)+CURTXY(IP,2)
          ELSE
            C(1,IP) = CG(IP,IS)*SINTH(ID)
            C(2,IP) = CG(IP,IS)*COSTH(ID)
          END IF
          IF (LSPHE) THEN
            C(1,IP) = C(1,IP)*1.0_JWRU/(DGRTH*COS(MIN(ABS(YP(IP)),YPMAX)*DEGRAD))
            C(2,IP) = C(2,IP)*DGRTHM1
          END IF
        END DO
      END IF

#ifdef ECMWF
      IF (LHOOK) CALL DR_HOOK('CADVXY',1,ZHOOK_HANDLE)
#endif
      END SUBROUTINE CADVXY

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UNWAM_ARRAYS

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Intialize UnWam arrays 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      USE YOWUNPOOL 
      USE yowpd, only: MNE=>ne, MNP=>npa
      USE YOWSHAL  , ONLY : DEPTH
      IMPLICIT NONE
      IF(ALLOCATED(DEPTH)) DEALLOCATE(DEPTH)
      ALLOCATE( DEPTH(MNP,1))
      ALLOCATE( CCON(MNP) ); CCON = 0
      ALLOCATE( SI(MNP) ); SI = 0.0_JWRU
      ALLOCATE( ITER_EXP(MSC,MDC) ); ITER_EXP = 0
      ALLOCATE( TRIA(MNE) ); TRIA = 0.0_JWRU
      IF (LCFL) THEN
        ALLOCATE(CFLCXY(3,MNP)); CFLCXY = 0.0_JWRU
      END IF
      ALLOCATE( IEN(6,MNE) ); IEN = 0
      ALLOCATE( IOBP(MNP)  ); IOBP = 0
      ALLOCATE( IOBPD(MDC,MNP)); IOBPD = 0
      ALLOCATE( IOBWB(MNP) ); IOBWB = 1 ! for boundary nodes we set this to 0 since in this case we omit advection at these nodes ... 
      ALLOCATE( CG(MNP,MSC) ); CG = 0.0_JWRU
      ALLOCATE( CURTXY(2,MNP)); CURTXY = 0.0_JWRU
      END SUBROUTINE INIT_UNWAM_ARRAYS

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UNWAM

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Initalize Unwam ... 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF*

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------
      USE MPL_MPIF
      USE YOWUNPOOL, ONLY : BND, DBG, GRID, FILEDEF
      USE YOWSHAL  , ONLY : DEPTH    ,DEPTHA    ,TOOSHALLOW
      USE YOWPARAM , ONLY : NIBLO
      USE yowpd, only: MNP=>npa, z, comm, np_global, initPD, setDimSize
      USE YOWPARAM , ONLY : NANG, NFRE
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IERR, IP

      CALL setDimSize(NANG, NFRE)
      CALL SET_UNWAM_HANDLES ! set file handles
      CALL initPD("system.dat")
      CALL mpi_barrier(comm, ierr)
      NIBLO = np_global
      CALL INIT_UNWAM_ARRAYS
!JB do not allow negative depth !!!!
      DO IP = 1, MNP 
        DEPTH(IP,1) = MAX(Z(IP),REAL(TOOSHALLOW,JWRU))
      ENDDO

      CALL INIT_FLUCT        ! init fluctuation splitting stuff 
      CALL SET_IOBP          ! boundary point marker 
      CALL SET_IOBPD         ! boundary directional marker 
      CLOSE(BND%FHNDL)
      END SUBROUTINE INIT_UNWAM

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SPHERICAL_COORDINATE_DISTANCE(LON1, LON2, LAT1, LAT2, DIST)
!     Purpose.: computes the distance on a sphere of radius=1 between (LON1,LAT1) and (LON2,LAT2)
!     --------

!        Explicit arguments :  
!        --------------------   

!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------

!     Externals.
!     ----------

!     Reference.
!     ----------
!     http://en.wikipedia.org/wiki/Great-circle_distance

!     Author.
!     -------

!     Modifications.
!     --------------
!     --------------------------------------------------------------
      USE YOWPCONS , ONLY : RAD

      IMPLICIT NONE

      REAL(KIND=JWRU), INTENT(IN) :: LON1, LON2, LAT1, LAT2
      REAL(KIND=JWRU), INTENT(OUT) :: DIST
      REAL(KIND=JWRU) :: SLAT, SLON, C1, C2

      SLAT = SIN(0.5_JWRU*(LAT1-LAT2)*RAD)**2
      SLON = SIN(0.5_JWRU*(LON1-LON2)*RAD)**2
      C1=COS(LAT1*RAD)
      C2=COS(LAT2*RAD)
      DIST = SQRT(MAX(SLAT+C1*C2*SLON,0.0_JWRU))

      IF (DIST .ge. 1.0_JWRU) THEN
        DIST=0.0_JWRU
      ELSE
        DIST = 2.0_JWRU*ASIN(DIST) 
      END IF
      END SUBROUTINE SPHERICAL_COORDINATE_DISTANCE


      SUBROUTINE SPHERICAL_COORDINATE_AREA(LON1, LON2, LON3, LAT1, LAT2, LAT3, AREA)
      IMPLICIT NONE
      REAL(KIND=JWRU), INTENT(IN) :: LON1, LON2, LON3, LAT1, LAT2, LAT3
      REAL(KIND=JWRU), INTENT(OUT) :: AREA
      REAL(KIND=JWRU) :: DistA, DistB, DistC, DistS
      REAL(KIND=JWRU) :: eTan1, eTan2, eTan3, eTan4
      REAL(KIND=JWRU) :: eProd, sqrtProd
      CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON2, LAT1, LAT2, DistA)
      CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON3, LAT1, LAT3, DistB)
      CALL SPHERICAL_COORDINATE_DISTANCE(LON2, LON3, LAT2, LAT3, DistC)
      DistS=0.5_JWRU*(DistA + DistB + DistC)
      eTan1=tan(0.5_JWRU*DistS)
      eTan2=tan(0.5_JWRU*(DistS - DistA))
      eTan3=tan(0.5_JWRU*(DistS - DistB))
      eTan4=tan(0.5_JWRU*(DistS - DistC))
      eProd=eTan1*eTan2*eTan3*eTan4
      sqrtProd=SQRT(eProd)
      AREA=4.0_JWRU*ATAN(sqrtProd)
      END SUBROUTINE SPHERICAL_COORDINATE_AREA

      SUBROUTINE CORRECT_SINGLE_DXP(DXP)
      IMPLICIT NONE
      REAL(KIND=JWRU), INTENT(INOUT) :: DXP
      IF (DXP .LE. -180.0_JWRU) THEN
        DXP=DXP + 360.0_JWRU
      END IF
      IF (DXP .GE. 180.0_JWRU) THEN
        DXP=DXP - 360.0_JWRU
      END IF
      END SUBROUTINE CORRECT_SINGLE_DXP


      SUBROUTINE INIT_FLUCT
      USE YOWUNPOOL
      USE yowpd, only: MNE=>ne, MNP=>npa, INE, XP=>x, YP=>y, z, exchange
      USE YOWFRED  , ONLY : FR
      USE YOWPARAM , ONLY : NANG, NFRE
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: I, J, K
      INTEGER(KIND=JWIM) :: IP, IE, IS, ISTAT
      INTEGER(KIND=JWIM) :: POS, POS_J, POS_K
      INTEGER(KIND=JWIM) :: IP_I, IP_J, IP_K
      INTEGER(KIND=JWIM) :: I1, I2, I3, NI(3)
      INTEGER(KIND=JWIM) :: CHILF(MNP), COUNT_MAX
      INTEGER(KIND=JWIM) :: ITMP(MNP)
      INTEGER(KIND=JWIM) :: TMPINE
      INTEGER(KIND=JWIM) :: INEXT, IP_NEXT, IE_ADJ, nbMatch, ICON
      INTEGER(KIND=JWIM) :: IADJ, POS_PREV, IP_ADJ_PREV
      INTEGER(KIND=JWIM) :: IE2, POS_NEXT, IP_ADJ_NEXT
      INTEGER(KIND=JWIM) :: MAX_DEG
      INTEGER(KIND=JWIM) :: POS_TRICK(3,2)
      INTEGER(KIND=JWIM), ALLOCATABLE :: PTABLE(:,:)

      REAL(KIND=JWRU) :: P1(2), P2(2), P3(2)
      REAL(KIND=JWRU) :: R1(2), R2(2), R3(2)
      REAL(KIND=JWRU) :: N1(2), N2(2), N3(2)
      REAL(KIND=JWRU) :: TRIA03
      REAL(KIND=JWRU) :: DXP1, DXP2, DXP3 
      REAL(KIND=JWRU) :: DYP1, DYP2, DYP3
      REAL(KIND=JWRU) :: DBLTMP, WN, WVC, WVK, WVCG
      REAL(KIND=JWRU) :: TL1, TL2, TL3
      REAL(KIND=JWRU) :: TLMIN,TMPTLMIN
      REAL(KIND=JWRU) :: TLMAX,TMPTLMAX
      REAL(KIND=JWRU) :: AVETA
      REAL(KIND=JWRU) :: PROV1, PROV2, PROV3
      REAL(KIND=JWRU) :: AREA, AREA_RAD

      LOGICAL   :: LWRONG = .FALSE.
      LOGICAL   :: CHECK_COMBIN_ORIENT

      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      TLMIN = LARGE
      TLMAX = SMALL 
      AVETA = 0.0_JWRU
!
!     CALCULATE NORMAL VECTORS
! 
      DO IE = 1, MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        DXP1 = XP(I2) - XP(I1) ! dx dy for each element and node ...
        DYP1 = YP(I2) - YP(I1)
        DXP2 = XP(I3) - XP(I2)
        DYP2 = YP(I3) - YP(I2)
        DXP3 = XP(I1) - XP(I3)
        DYP3 = YP(I1) - YP(I3)
        IF (APPLY_DXP_CORR) THEN
          CALL CORRECT_SINGLE_DXP(DXP1)
          CALL CORRECT_SINGLE_DXP(DXP2)
          CALL CORRECT_SINGLE_DXP(DXP3)
        END IF
        IEN(1,IE) = - DYP2
        IEN(2,IE) =   DXP2
        IEN(3,IE) = - DYP3
        IEN(4,IE) =   DXP3
        IEN(5,IE) = - DYP1
        IEN(6,IE) =   DXP1
        DBLTMP = (DXP3*DYP1 - DYP3*DXP1)*0.5_JWRU
        IF (LSPHE .and. USE_EXACT_FORMULA_SPHERICAL_AREA) THEN
          CALL SPHERICAL_COORDINATE_AREA(XP(I1), XP(I2), XP(I3),    &
     &           YP(I1), YP(I2), YP(I3), AREA_RAD)
          AREA=AREA_RAD*RADDEG*RADDEG
        ELSE
          AREA=DBLTMP
        END IF
        TRIA(IE) = AREA
        IF (TRIA(IE) .LT. SMALL) THEN
          TMPINE = INE(2,IE)
          INE(2,IE) = INE(3,IE)
          INE(3,IE) = TMPINE
          I2 = INE(2,IE)
          I3 = INE(3,IE)
          TRIA(IE) = -1.0*TRIA(IE)
          PROV1=IEN(6,IE) ! DXP1
          PROV2=IEN(2,IE) ! DXP2
          PROV3=IEN(4,IE) ! DXP3
          IEN(6,IE)=-PROV3
          IEN(2,IE)=-PROV2
          IEN(4,IE)=-PROV1
          PROV1= - IEN(5,IE)
          PROV2= - IEN(1,IE)
          PROV3= - IEN(3,IE)
          IEN(1,IE) = PROV2
          IEN(3,IE) = PROV1
          IEN(5,IE) = PROV3
          LWRONG = .TRUE.     ! check element orientation and swap if needed ...
        END IF
        TL1 = SQRT(IEN(5,IE)**2 + IEN(6,IE)**2)  ! Edge length's 
        TL2 = SQRT(IEN(3,IE)**2 + IEN(4,IE)**2)
        TL3 = SQRT(IEN(1,IE)**2 + IEN(2,IE)**2)
        TMPTLMIN = MIN(TL1, TL2, TL3)                ! max. min. edge lengts of hte whole mesh 
        TMPTLMAX = MAX(TL1, TL2, TL3)
        IF (TLMIN > TMPTLMIN) THEN
          TLMIN = TMPTLMIN
        END IF
        IF (TLMAX < TMPTLMAX) THEN
          TLMAX = TMPTLMAX
        END IF
        AVETA = AVETA+TRIA(IE)
      END DO
      AVETA = AVETA / MNE ! Average edge area ...
!
! Calculate the max. number of connected elements count_max
!
      SI(:)   = 0.0_JWRU ! Median Dual Patch Area of each Node
      CCON(:) = 0        ! Number of connected Elements
      DO IE = 1 , MNE
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        CCON(I1) = CCON(I1) + 1
        CCON(I2) = CCON(I2) + 1
        CCON(I3) = CCON(I3) + 1
        TRIA03 = ONETHIRD * TRIA(IE)
        SI(I1) = SI(I1) + TRIA03
        SI(I2) = SI(I2) + TRIA03
        SI(I3) = SI(I3) + TRIA03
      ENDDO
      CALL exchange(SI)
      DO IP =1, MNP
        IF(SI(IP) .LT. SMALL) THEN
          write(DBG%FHNDL,*) 'SI(IP) NEG. OR LESS THAN THR. IP=', IP
          STOP 'SI NEG. OR LESS THAN THR'
        ENDIF
      END DO
      MAXMNECON  = MAXVAL(CCON)
      ALLOCATE(CELLVERTEX(MNP,MAXMNECON,2), stat=istat)
      IF (istat/=0) stop 'wwm_fluctsplit, allocate error 4'
      CELLVERTEX(:,:,:) = 0 ! Stores for each node the Elementnumbers of the connected Elements
      CHILF             = 0 ! and the Position of the position of the Node in the Element Index
      DO IE = 1, MNE
        DO J=1,3
          I = INE(J,IE)
          CHILF(I) = CHILF(I)+1
          CELLVERTEX(I,CHILF(I),1) = IE
          CELLVERTEX(I,CHILF(I),2) = J
        END DO
      ENDDO
!
!        Emulates loop structure and counts max. entries in the different pointers that have to be designed
!
      J = 0
      DO IP = 1, MNP
        DO I = 1, CCON(IP)
          J = J + 1
        END DO
      END DO
      COUNT_MAX = J ! Max. Number of entries in the pointers used in the calculations
      IF (COUNT_MAX.ne.3*MNE) THEN
        Print *, 'COUNT_MAX=', COUNT_MAX
        Print *, 'MNE=', MNE
        STOP 'Do Not Sleep Before solving the problem'
      ENDIF

      ALLOCATE (IE_CELL(COUNT_MAX), POS_CELL(COUNT_MAX), stat=istat)
      IF (istat/=0) STOP 'MEMORY ERROR AT unwam.F 1210'
      ALLOCATE (IE_CELL2(MNP,MAXMNECON), POS_CELL2(MNP,MAXMNECON), stat=istat)
      IF (istat/=0) STOP 'MEMORY ERROR AT unwam.F 1212'

      IE_CELL  = 0
      POS_CELL = 0
      IE_CELL2  = 0
      POS_CELL2 = 0

      J = 0
      DO IP = 1, MNP
        DO I = 1, CCON(IP)
          J = J + 1
          IE_CELL(J)      = CELLVERTEX(IP,I,1)
          POS_CELL(J)     = CELLVERTEX(IP,I,2)
          IE_CELL2(IP,I)  = CELLVERTEX(IP,I,1)
          POS_CELL2(IP,I) = CELLVERTEX(IP,I,2)
        END DO
      END DO
      DEALLOCATE(CELLVERTEX)

      CHECK_COMBIN_ORIENT=.TRUE.
      IF (CHECK_COMBIN_ORIENT) THEN
        DO IE=1,MNE
          DO I=1,3
            INEXT=POS_TRICK(I,1)
            IP=INE(I, IE)
            IP_NEXT=INE(I, IE)
            nbMatch=0
            IE_ADJ=-1
            DO ICON=1,CCON(IP)
              IE2=IE_CELL2(IP,ICON)
              IF (IE .ne. IE2) THEN
                POS=POS_CELL2(IP, ICON)
                POS_NEXT=POS_TRICK(POS,1)
                IP_ADJ_NEXT=INE(POS_NEXT,IE2)
                IF (IP_ADJ_NEXT .eq. IP_NEXT) THEN
                  Print *, 'Combinatorial orientability problem'
                  Print *, 'IE=', IE, ' IE2=', IE2
                  Print *, 'IP=', IP, ' IP_NEXT=', IP_NEXT
                  STOP
                END IF
                POS_PREV=POS_TRICK(POS,2)
                IP_ADJ_PREV=INE(POS_PREV,IE2)
                IF (IP_ADJ_PREV .eq. IP_NEXT) THEN
                  nbMatch=nbMatch+1
                  IE_ADJ=IE2
                END IF
              END IF
            END DO
            IF (nbMatch .gt. 1) THEN
              Print *, 'nbMatch is too large.'
              Print *, 'Should be 0 for boundary edge'
              Print *, 'Should be 1 for interior edges'
              Print *, 'nbMatch=', nbMatch
              STOP
            END IF
          END DO
        END DO
      END IF


      IF (LIMPLICIT) THEN
        ALLOCATE(PTABLE(COUNT_MAX,7), JA_IE(3,3,MNE), stat=istat)
        IF (istat/=0) stop 'wwm_fluctsplit, allocate error 6'
        PTABLE(:,:) = 0. ! Table storing some other values needed to design the sparse matrix pointers.
        J = 0
        DO IP = 1, MNP
          DO I = 1, CCON(IP)
            J = J + 1
            IE    = IE_CELL(J)
            POS   = POS_CELL(J)
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            IF (POS == 1) THEN
              POS_J = 2
              POS_K = 3
            ELSE IF (POS == 2) THEN
              POS_J = 3
              POS_K = 1
            ELSE
              POS_J = 1
              POS_K = 2
            END IF
            IP_I = IP
            IP_J = INE(POS_J,IE)
            IP_K = INE(POS_K,IE)
            PTABLE(J,1) = IP_I ! Node numbers of the connected elements
            PTABLE(J,2) = IP_J
            PTABLE(J,3) = IP_K
            PTABLE(J,4) = POS  ! Position of the nodes in the element index
            PTABLE(J,5) = POS_J
            PTABLE(J,6) = POS_K
            PTABLE(J,7) = IE   ! Element numbers same as IE_CELL
          END DO
        END DO
!
! Count number of nonzero entries in the matrix ...
! Basically, each connected element may have two off-diagonal
! contribution and one diagonal related to the connected vertex itself ...
!
        J = 0
        NNZ = 0
        DO IP = 1, MNP
          ITMP(:) = 0
          DO I = 1, CCON(IP)
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            ITMP(IP)   = 1
            ITMP(IP_J) = 1
            ITMP(IP_K) = 1
          END DO
          NNZ = NNZ + SUM(ITMP)
        END DO
!
! Allocate sparse matrix pointers using the Compressed Sparse Row Format CSR ... this is now done only of MNP nodes
! The next step is to do it for the whole Matrix MNP * MSC * MDC
! see ...:x
!
! JA Pointer according to the convention in my thesis see p. 123
! IA Pointer according to the convention in my thesis see p. 123
        ALLOCATE (JA(NNZ), IA(MNP+1), POSI(3,COUNT_MAX), I_DIAG(MNP), ASPAR_JAC(NANG,NFRE,NNZ), stat=istat)
        IF (istat/=0) stop 'wwm_fluctsplit, allocate error 6'
        JA = 0
        IA = 0
        POSI = 0
! Points to the position of the matrix entry in the mass matrix
! according to the CSR matrix format see p. 124
        J = 0
        K = 0
        IA(1) = 1
        MAX_DEG=0
        DO IP = 1, MNP ! Run through all rows
          ITMP=0
          DO I = 1, CCON(IP) ! Check how many entries there are ...
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            ITMP(IP)   = 1
            ITMP(IP_J) = 1
            ITMP(IP_K) = 1
          END DO
          IADJ=0
          DO I = 1, MNP ! Run through all columns
            IF (ITMP(I) .GT. 0) THEN
              K = K + 1
              IF (I .ne. IP) THEN
                IADJ=IADJ + 1
              END IF
              JA(K) = I
            END IF
          END DO
          IF (IADJ .gt. MAX_DEG) THEN
            MAX_DEG=IADJ
          END IF
          IA(IP + 1) = K + 1
        END DO


        J = 0
        DO IP = 1, MNP
          DO I = 1, CCON(IP)
            J = J + 1
            IP_J  = PTABLE(J,2)
            IP_K  = PTABLE(J,3)
            DO K = IA(IP), IA(IP+1) - 1
              IF (IP   == JA(K)) POSI(1,J)  = K
              IF (IP   == JA(K)) I_DIAG(IP) = K
              IF (IP_J == JA(K)) POSI(2,J)  = K
              IF (IP_K == JA(K)) POSI(3,J)  = K
            END DO
          END DO
        END DO

        J=0
        DO IP=1,MNP
          DO I = 1, CCON(IP)
            J=J+1
            IE    =  IE_CELL(J)
            POS   =  POS_CELL(J)
            I1    =  POSI(1,J)
            I2    =  POSI(2,J)
            I3    =  POSI(3,J)
            JA_IE(POS,1,IE) = I1
            JA_IE(POS,2,IE) = I2
            JA_IE(POS,3,IE) = I3
          END DO
        END DO
        !
        ! Arrays for Jacobi-Gauss-Seidel solver
        !
        DEALLOCATE(PTABLE)
      ENDIF
!
!     Group velocities
!
      DO IS = 1, MSC
        DO IP = 1, MNP 
          CALL WAVEKCG8(z(IP), FR(IS)*PI2, WN, WVC, WVK, WVCG)
          CG(IP,IS) = WVCG
        END DO 
      END DO
      END SUBROUTINE INIT_FLUCT

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBP

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Estimate boundary nodes of an abirtrary unstructured mesh ...
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE yowpd, only: MNP=>npa, MNE=>ne, ipgl, exchange,np_global, INE

        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: I, IWILD(MNP)
        INTEGER(KIND=JWIM) :: I1, I2, I3, IE, IP, ID, IFSTAT, IP_global
        INTEGER(KIND=JWIM), POINTER :: STATUS(:)
        INTEGER(KIND=JWIM), POINTER :: COLLECTED(:)
        INTEGER(KIND=JWIM), POINTER :: NEXTVERT(:)
        INTEGER(KIND=JWIM), POINTER :: PREVVERT(:)
        INTEGER(KIND=JWIM) :: ISFINISHED, INEXT, IPREV
        INTEGER(KIND=JWIM) :: IPNEXT, IPPREV, ZNEXT, ITMP
        INTEGER(KIND=JWIM) :: IOBPcopy(MNP), eDiff, nbDiff

        REAL(KIND=JWRB) :: BNDTMP
        REAL(KIND=JWRU) :: x1, y1, x2, y2
        REAL(KIND=JWRU) :: EVX, EVY
        REAL(KIND=JWRU) :: eDet1, eDet2
        REAL(KIND=JWRU) :: ATMP, BTMP
        REAL(KIND=JWRU) :: rtemp(MNP)

        LOGICAL :: LFLIVE

!
! open and read boundary nodes file ...
!
        IOBP    = 0
        IOBPD   = 0

        OPEN(BND%FHNDL, FILE = BND%FNAME, STATUS = 'OLD')
        CALL RHEADER_NODE(BND%FHNDL,ITMP,ITMP)
        IP = 1
        DO IP_global = 1, np_global
          IP = ipgl(IP_global)
          IF(IP /= 0) THEN
            READ(BND%FHNDL,*) ITMP, BNDTMP, BNDTMP, BNDTMP
            IF (BNDTMP .GT. 0.) IOBP(IP) = INT(BNDTMP)
            IF (.NOT.  LBCWA) THEN
              IF (IOBP(IP) .EQ. 2) IOBP(IP) = 1
            END IF
          END IF
        END DO

        IWBMNP = 0
        DO IP = 1, MNP
          IF (IOBP(IP) == 2) IWBMNP = IWBMNP + 1 ! Local number of boundary nodes ...
        END DO
!
! map boundary nodes ... needed later for the decomposition ...
!
!        ALLOCATE( IWBNDLC(IWBMNP) ) 
!        IWBMNP = 0
!        DO IP = 1, MNP
!          IF (IOBP(IP) == 2) THEN
!            IWBMNP = IWBMNP + 1
!            IWBNDLC(IWBMNP) = IP ! Stores local wave boundary index 
!          END IF
!        END DO
!
! find islands and domain boundary ....
!
        ALLOCATE(STATUS(MNP))
        ALLOCATE(COLLECTED(MNP))
        ALLOCATE(PREVVERT(MNP))
        ALLOCATE(NEXTVERT(MNP))

        STATUS(:) = 0

        DO IE=1,MNE
          DO I=1,3
            IF (I.EQ.1) THEN
              IPREV=3
            ELSE
              IPREV=I-1
            END IF
            IF (I.EQ.3) THEN
              INEXT=1
            ELSE
              INEXT=I+1
            END IF
            IP=INE(I,IE)
            IPNEXT=INE(INEXT,IE)
            IPPREV=INE(IPREV,IE)
            IF (STATUS(IP).EQ.0) THEN
              STATUS(IP)=1
              PREVVERT(IP)=IPPREV
              NEXTVERT(IP)=IPNEXT
            END IF
          END DO
        END DO

        STATUS(:)=0
        DO
          COLLECTED(:)=0
          DO IE=1,MNE
            DO I=1,3
              IF (I.EQ.1) THEN
                IPREV=3
              ELSE
                IPREV=I-1
              END IF
              IF (I.EQ.3) THEN
                INEXT=1
              ELSE
                INEXT=I+1
              END IF
              IP=INE(I,IE)
              IPNEXT=INE(INEXT,IE)
              IPPREV=INE(IPREV,IE)
              IF (STATUS(IP).eq.0) THEN
                ZNEXT=NEXTVERT(IP)
                IF (ZNEXT.eq.IPPREV) THEN
                  COLLECTED(IP)=1
                  NEXTVERT(IP)=IPNEXT
                  IF (NEXTVERT(IP).eq.PREVVERT(IP)) THEN
                    STATUS(IP)=1
                  END IF
                END IF
              END IF
            END DO
          END DO
          ISFINISHED=1
          DO IP=1,MNP
!    Correction of MDS, begin
            IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
              STATUS(IP)=-1
            END IF
!    Correction of MDS, end
            IF (STATUS(IP).eq.0) THEN
              ISFINISHED=0
            END IF
          END DO
          IF (ISFINISHED.eq.1) THEN
            EXIT
          END IF
        END DO
        DO IP=1,MNP
          IF (STATUS(IP).eq.-1 .AND. IOBP(IP) .EQ. 0) THEN
            IOBP(IP)=1
          END IF
        END DO

        DEALLOCATE(STATUS)
        DEALLOCATE(COLLECTED)
        DEALLOCATE(NEXTVERT)
        DEALLOCATE(PREVVERT)
!
! allocate wave boundary arrays ... 
!
!        IF (LINHOM) THEN
!          OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
!          ALLOCATE( WBAC(MSC,MDC,IWBMNP) ); WBAC = 0.
!          IF (LBINTER) THEN ! For time interpolation 
!            ALLOCATE( WBACOLD(MSC,MDC,IWBMNP) ); WBACOLD = 0.
!            ALLOCATE( WBACNEW(MSC,MDC,IWBMNP) ); WBACNEW = 0.
!            ALLOCATE( DSPEC(MSC,MDC,IWBMNP) ); DSPEC   = 0.
!          END IF
!        ELSE
!          OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')
!          ALLOCATE( WBAC(MSC,MDC,1) ); WBAC = 0.
!          IF (LBINTER) THEN
!            ALLOCATE( WBACOLD(MSC,MDC,1) ); WBACOLD = 0.
!            ALLOCATE( WBACNEW(MSC,MDC,1) ); WBACNEW = 0.
!            ALLOCATE( DSPEC(MSC,MDC,1) ); DSPEC   = 0.
!          END IF
!        ENDIF ! LINHOM

!        WRITE(DBG%FHNDL,'("+TRACE...",A,I10)') 
!     &                  'Number of Active Wave Boundary Nodes', IWBMNP
!
! write test output ... 
!
!        WRITE(DBG%FHNDL,*) 'BOUNDARY MAPPING'
        DO IP = 1, MNP
!           WRITE(DBG%FHNDL,*) IP, IOBP(IP)
        END DO

        IOBPcopy=IOBP
        call exchange(IOBPcopy)

        rtemp = REAL(IOBP(:),JWRU)
        call exchange(rtemp)
        IOBP(:) = INT(rtemp)

!          call exchange(IOBP)

        nbDiff=0
        DO IP=1,MNP
          IF (IOBP(IP) .ne. IOBPcopy(IP)) THEN
            nbDiff=nbDiff+1
!            WRITE(DBG%FHNDL,*) 'IP/IOBP/IOBPcopy=', IP, 
!     &           IOBP(IP), IOBPcopy(IP)
          END IF
        END DO
!        WRITE(DBG%FHNDL,*) 'nbDiff=', nbDiff
      END SUBROUTINE SET_IOBP

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_IOBPD

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Estimate spectral direction that are pointing into the domain from the boundary ...
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL 
        USE yowpd, only: MNP=>npa, INE, MNE=>ne, exchange
        USE YOWFRED, ONLY : DFIM, COSTH, SINTH

        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: I1, I2, I3, IE, IP, ID, NI(3)
        INTEGER(KIND=JWIM) :: I, IWILD(MNP)
        INTEGER(KIND=JWIM), POINTER :: STATUS(:)
        INTEGER(KIND=JWIM), POINTER :: COLLECTED(:)
        INTEGER(KIND=JWIM), POINTER :: NEXTVERT(:)
        INTEGER(KIND=JWIM), POINTER :: PREVVERT(:)
        INTEGER(KIND=JWIM) :: ISFINISHED, INEXT, IPREV
        INTEGER(KIND=JWIM) :: IPNEXT, IPPREV, ZNEXT

        REAL(KIND=JWRU) :: DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
        REAL(KIND=JWRU) :: x1, y1, x2, y2
        REAL(KIND=JWRU) :: EVX, EVY
        REAL(KIND=JWRU) :: eDet1, eDet2
        REAL(KIND=JWRU) :: rtemp(MNP)
!
! SET IOBPD ...
!
        DO IE=1,MNE
          I1   =   INE(1,IE)
          I2   =   INE(2,IE)
          I3   =   INE(3,IE)
          DXP1 =   IEN(6,IE)
          DYP1 = - IEN(5,IE)
          DXP2 =   IEN(2,IE)
          DYP2 = - IEN(1,IE)
          DXP3 =   IEN(4,IE)
          DYP3 = - IEN(3,IE)
!2do ... modifly wave direction by currents ...
          DO ID=1,MDC
            EVX=SINTH(ID)
            EVY=COSTH(ID)
            DO I=1,3
               IF (I.eq.1) THEN
                  x1=   DXP1
                  y1=   DYP1
                  x2= - DXP3
                  y2= - DYP3
                  IP=   I1
               END IF
               IF (I.eq.2) THEN
                  x1 =   DXP2
                  y1 =   DYP2
                  x2 = - DXP1
                  y2 = - DYP1
                  IP =   I2
               END IF
               IF (I.eq.3) THEN
                  x1 =   DXP3
                  y1 =   DYP3
                  x2 = - DXP2
                  y2 = - DYP2
                  IP =   I3
               END IF
               eDet1 = SMALL-x1*EVY+y1*EVX
               eDet2 = SMALL+x2*EVY-y2*EVX
               IF ((eDet1.gt.0.0_JWRU).and.(eDet2.gt.0.0_JWRU)) THEN
                  IOBPD(ID,IP)=1
               END IF
            END DO
          END DO
        END DO

        DO IP = 1, MNP
          IF ( LBCWA ) THEN
            IF ( IOBP(IP) == 2 ) THEN
              IOBWB(IP) = 0
              IOBPD(:,IP) = 1
            ENDIF
          END IF
          IF ( IOBP(IP) == 3 ) THEN ! If Neumann boundary condition is given set IOBP to 4
            IOBPD(:,IP) = 1 ! Update Neumann nodes ...
          END IF
        END DO

        DO ID=1, MDC
          rtemp = IOBPD(ID,:)
          call exchange(rtemp)
          IOBPD(ID,:) = rtemp
        END DO

        RETURN
      END SUBROUTINE SET_IOBPD

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_SYSTEM

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Check system size and allocate 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL, ONLY : GRID, FILEDEF!, MNP, MNE
        USE yowpd, only: MNE=>ne, INE, MNP=>npa
        USE YOWPARAM , ONLY : NIBLO

        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: NBOUND, NDOMAIN
        INTEGER(KIND=JWIM) :: I, ITMP, NGESAMT

        OPEN(GRID%FHNDL, FILE = GRID%FNAME, STATUS='old')
        CALL RHEADER_NODE(GRID%FHNDL,NBOUND,NDOMAIN)
        MNP = NBOUND + NDOMAIN 
        NIBLO = MNP
        DO I = 1, MNP
          READ(GRID%FHNDL, *) 
        END DO 
        CALL RHEADER_ELEMENT(GRID%FHNDL,MNE)
        CLOSE(GRID%FHNDL)
  
      END SUBROUTINE CHECK_SYSTEM

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UNWAM_OUT(IHANDLE)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Write everything out to disk ...

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE YOWSHAL  , ONLY : DEPTH
        use yowpd, only: MNP=>npa, MNE=>ne, XP=>x, YP=>y, INE
        IMPLICIT NONE

        INTEGER(KIND=JWIM), INTENT(IN) :: IHANDLE

!        WRITE(DBG%FHNDL,*) MNP, MNE, MDC, MSC

        WRITE(IHANDLE) MNP, MNE, MSC, MDC
        WRITE(IHANDLE) XP
        WRITE(IHANDLE) YP 
        WRITE(IHANDLE) DEPTH
        WRITE(IHANDLE) CCON
        WRITE(IHANDLE) SI
        WRITE(IHANDLE) TRIA
        WRITE(IHANDLE) INE
        WRITE(IHANDLE) IEN
        WRITE(IHANDLE) CG
        WRITE(IHANDLE) IOBP
        WRITE(IHANDLE) IOBPD
        WRITE(IHANDLE) IOBWB

      END SUBROUTINE UNWAM_OUT

!**********************************************************************
!!*                                                                   *
!**********************************************************************
      SUBROUTINE UNWAM_IN(IHANDLE)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Read everything from disk ...

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE YOWSHAL  , ONLY : DEPTH
!  you can not fill pdlib array direct.
!         use yowpd, only: MNP=>npa, MNE=>ne, XP=>x, YP=>y, INE
        IMPLICIT NONE

        INTEGER(KIND=JWIM), INTENT(IN) :: IHANDLE
! 
!         READ(IHANDLE) MNP, MNE, MSC, MDC

        CALL INIT_UNWAM_ARRAYS
!         READ(IHANDLE) XP
!         READ(IHANDLE) YP 
        READ(IHANDLE) DEPTH
        READ(IHANDLE) CCON
        READ(IHANDLE) SI
        READ(IHANDLE) TRIA
!         READ(IHANDLE) INE
        READ(IHANDLE) IEN
        READ(IHANDLE) CG
        READ(IHANDLE) IOBP
        READ(IHANDLE) IOBPD
        READ(IHANDLE) IOBWB

      END SUBROUTINE UNWAM_IN

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_UNWAM_HANDLES

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Set file handles using ECMWF tool ....

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

        USE YOWUNPOOL
        USE YOWTEST  , ONLY : IU06
        USE YOWMPP   , ONLY : IRANK, NPROC
        IMPLICIT NONE

        INTEGER(KIND=JWIM) :: I_GET_UNIT

!        GRID%FNAME = 'system.dat'
!        GRID%FHNDL = I_GET_UNIT(IU06, GRID%FNAME, 'r', 'f', 0) 

        BND%FNAME = 'sysbnd.dat'
        BND%FHNDL = I_GET_UNIT(IU06, BND%FNAME, 'w', 'f', 0)

        DBG%FNAME = 'unwamdbg.%p.dat'
        CALL EXPAND_STRING(IRANK,NPROC,0,0,DBG%FNAME,1)
        DBG%FHNDL = I_GET_UNIT(IU06, DBG%FNAME, 'w', 'f', 0)

        IF(LLUNBINOUT) THEN
          XFN_HS%FNAME = 'erghs.bin'
          XFN_HS%FHNDL = I_GET_UNIT(IU06, XFN_HS%FNAME, 'w', 'u', 0)

          XFN_TM%FNAME = 'ergtm.bin'
          XFN_TM%FHNDL = I_GET_UNIT(IU06, XFN_TM%FNAME, 'w', 'u', 0)

          XFN_TEST%FNAME = 'ergtest.bin'
          XFN_TEST%FHNDL = I_GET_UNIT(IU06, XFN_TEST%FNAME, 'w', 'u', 0)
        ENDIF

      END SUBROUTINE SET_UNWAM_HANDLES

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT(X,Y,Z,XP,YP,Wi,Zi,LSAME)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Interpolate any point within 1.d0 element based on linear shape functions ...
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE

      REAL(KIND=JWRU), INTENT(IN)  :: X(3), Y(3), Z(3)
      REAL(KIND=JWRU), INTENT(IN)  :: XP, YP
      REAL(KIND=JWRU), INTENT(OUT) :: Zi
      REAL(KIND=JWRU), INTENT(OUT) :: WI(3)

      REAL(KIND=JWRU)               :: y1,y2,y3,x1,x2,x3,z1,z2,z3
      REAL(KIND=JWRU), SAVE         :: A,B,C,D
      REAL(KIND=JWRU), PARAMETER    :: THR = TINY(1.0_JWRU)

      LOGICAL, INTENT(IN)  :: LSAME

      IF (.NOT. LSAME) THEN
        x1 = X(1); x2 = X(2); x3 = X(3)
        y1 = Y(1); y2 = Y(2); y3 = Y(3)
        z1 = Z(1); z2 = Z(2); z3 = Z(3)
        A = y1*(z2 - z3)  +  y2*(z3 - z1) +  y3*(z1 - z2)
        B = z1*(x2 - x3)  +  z2*(x3 - x1) +  z3*(x1 - x2)
        C = x1*(y2 - y3)  +  x2*(y3 - y1) +  x3*(y1 - y2)
        D = -A*x1 - B*y1 - C*z1
        IF (ABS(C) .GT. THR ) THEN
          WI(1) = -A/C
          WI(2) = -B/C
          WI(3) = -D/C
        ELSE
          WI    = 0.0_JWRU
        END IF 
      END IF
      Zi = REAL(WI(1) * REAL(XP,JWRU) + WI(2) * REAL(YP,JWRU) + WI(3))

      END SUBROUTINE INTELEMENT

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVEKCG8(DEP, SIGIN, WN, WVC, WVK, WVCG)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     Estimate the max. integration time step and amount of iterations ...

!     Externals.  Estimate group vel. 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------
         USE YOWUNPOOL, ONLY : SMALL
         USE YOWPCONS , ONLY : G

         IMPLICIT NONE

         REAL(KIND=JWRU),   INTENT(IN)  :: DEP
         REAL(KIND=JWRU), INTENT(IN)  :: SIGIN
         REAL(KIND=JWRU), INTENT(OUT) :: WVC, WVK, WVCG, WN
         REAL(KIND=JWRU) :: SGDLS , AUX1, AUX2
         REAL(KIND=JWRU) :: WKDEP, DMIN
! 
!AR: I put the mindepth for kcg to 0.1d0
!
         DMIN = 0.1_JWRU

         IF (SIGIN .LT. SMALL) THEN
            WN = 0.0_JWRU
            WVK=10.0_JWRU
            WVCG=0.0_JWRU
            RETURN
         END IF

         IF (DEP > DMIN) THEN
            SGDLS = SIGIN*SIGIN*DEP/G
            AUX1 = 1.0_JWRU+0.6522_JWRU*SGDLS+0.4622_JWRU*               &
     &             (SGDLS**2)+0.0864_JWRU*(SGDLS**4)+                    &
     &              0.0675_JWRU*(SGDLS**5)
            AUX2 = 1.0_JWRU/(SGDLS+1.0_JWRU/AUX1)
            WVC = SQRT(AUX2*G*REAL(DEP,JWRU))
            WVK = SIGIN/WVC
            WKDEP = WVK*REAL(DEP,JWRU)
            IF (WKDEP > 13.0_JWRU) THEN
               WN = 0.5_JWRU
            ELSE
               WN = 0.5_JWRU*(1.0_JWRU+2.0_JWRU*WKDEP/                    &
     &              SINH(MIN(300.0_JWRU,2.0_JWRU*WKDEP)))
            END IF
            WVCG = WN*WVC
         ELSE
            WVC  = 0.0_JWRU
            WVK  = 10.0_JWRU
            WVCG = 0.0_JWRU
         END IF

!!!!!!debile: make it all deep
!            WVC = G/SIGIN
!            WN = 0.5_JWRU
!            WVK = SIGIN/WVC
!            WVCG = WN*WVC

      END SUBROUTINE WAVEKCG8

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_NODE(IFILE,NKR,NKG)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     XFN header nodes ...

!     Externals. 
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE

      INTEGER(KIND=JWIM), INTENT(IN) :: IFILE
      INTEGER(KIND=JWIM), INTENT(OUT):: NKR, NKG

      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NKR
      READ(IFILE,*)
      READ(IFILE,*) NKG
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)

      END SUBROUTINE RHEADER_NODE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE RHEADER_ELEMENT(IFILE,NELEM)

!     Purpose.: Advects the spectra using non-vectorized RD-schemes
!     --------


!        Explicit arguments :  
!        --------------------   


!        Implicit arguments :     N1.d0
!        --------------------

!     Method.
!     -------
!     XFN header elements

!     Externals.  
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Aron Roland 2011  *ECMWF* based on Aron Roland & Mathieu Dutour 

!     Modifications.
!     --------------
!        Original : Aron Roland, 11,2011
!     --------------------------------------------------------------

      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IFILE
      INTEGER(KIND=JWIM), INTENT(OUT):: NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*) NELEM
      READ(IFILE,*)
      READ(IFILE,*)
      READ(IFILE,*)
      END SUBROUTINE RHEADER_ELEMENT

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE FIND_ELE (Xo,Yo,IE ) !T.s: works good

!     ! Xo, Yo is the point of interest
!     ! IE is the global element id, which contains the point Xo, Yo
!
      use yowpd, only: MNE=>ne_global, INE=>INE_GLOBAL, NODES=>nodes_global
!
      IMPLICIT NONE

! Dummy arguments
!
      REAL(KIND=JWRU), INTENT(IN)       :: Xo , Yo
      INTEGER(KIND=JWIM), INTENT(INOUT)   :: IE 
!
! Local variables
!
      REAL(KIND=JWRU), SAVE   :: xi, xj, xk, yi, yj, yk, dx, dy, f, Xo8, Yo8
      REAL(KIND=JWRU), SAVE   :: xmax , xmin , ymax , ymin
      INTEGER(KIND=JWIM), SAVE  :: if0 , if1 , ijk , k , ki , kj , kk , l
      INTEGER(KIND=JWIM), SAVE  :: i , i0 , idx , ielem
      

      REAL(KIND=JWRU):: THR = TINY(1.0_JWRU)

      DATA idx/0/ , i0/0/ , ielem/1/ , if0/0/ , if1/0/
!
!     Laengste Kannte (DX,DY) bestimmen
!     Dieser Programmabschnitt wird nur beim ersten Aufruf
!     durchlaufen!
!
      Xo8 = REAL(Xo,JWRU)
      Yo8 = REAL(Yo,JWRU)

      IF ( idx/=1 ) THEN
        idx = 1
        DO i = 1 , MNE
          ki = INE(1,i) ! + 1
          kj = INE(2,i) ! + 1
          kk = INE(3,i) ! + 1
          xi = REAL(NODES(ki)%X,JWRU)
          yi = REAL(NODES(ki)%Y,JWRU)
          xj = REAL(NODES(kj)%X,JWRU)
          yj = REAL(NODES(kj)%Y,JWRU)
          xk = REAL(NODES(kk)%X,JWRU)
          yk = REAL(NODES(kk)%Y,JWRU)
          IF ( i==1 ) THEN
            dx = MAX(ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
            dy = MAX(ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
          ELSE
            dx = MAX(dx,ABS(xi-xj),ABS(xi-xk),ABS(xj-xk))
            dy = MAX(dy,ABS(yi-yj),ABS(yi-yk),ABS(yj-yk))
          END IF
        END DO
      END IF
!     ------------------------------------------------------------------
!     TEST, OB DER PUNKT IM ZULETZT ANGESPROCHENEN ELEMENT LIEGT
!     ------------------------------------------------------------------
      IF ( i0==1 .AND. IE/=-1 ) THEN
         IF ( Yo8-ymin > THR ) THEN
            IF ( Yo8-ymax < -THR ) THEN
               IF ( Xo8-xmin > THR ) THEN
                  IF ( Xo8-xmax < -THR ) THEN
                     f = xi*(yj-Yo8) + xj*(Yo8-yi) + Xo8*(yi-yj)
                     IF ( f > THR ) THEN
                        f = xj*(yk-Yo8) + xk*(Yo8-yj) + Xo8*(yj-yk)
                        IF ( f > THR  ) THEN
                           f = xk*(yi-Yo8) + xi*(Yo8-yk) + Xo8*(yk-yi)
                           IF ( f > THR ) THEN
                              IE = ielem ! Element gefunden -->RETURN
                              RETURN
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
      endif
!     ------------------------------------------------------------------
!     Element suchen
!     ------------------------------------------------------------------
      i0 = 1
      i = ielem
      IF ( i<1 ) i = 1
      k = i
      l = i
      ijk = 0

100   DO
         ijk = ijk + 1
!.....   ABFRAGE AUF X-Richtung
         ki = INE(1,i)! + 1
         xi = REAL(NODES(ki)%X,JWRU)
         IF ( DABS(xi-Xo8)<=dx ) THEN
            kj = INE(2,i)! + 1
            kk = INE(3,i)! + 1
            xj = REAL(NODES(kj)%X,JWRU)
            xk = REAL(NODES(kk)%X,JWRU)
!.....    Punkt ausserhalb Element:
            xmin = MIN(xi,xj,xk)
            IF ( Xo8>=xmin ) THEN
               xmax = MAX(xi,xj,xk)
               IF ( Xo8<=xmax ) THEN
!.....        ABFRAGE AUF Y-Richtung
                  yi = REAL(NODES(ki)%Y,JWRU)
                  IF ( DABS(yi-Yo8)<=dy ) THEN
                     yj = REAL(NODES(kj)%Y,JWRU)
                     yk = REAL(NODES(kk)%Y,JWRU)
!.....          Punkt ausserhalb Element:
                     ymin = MIN(yi,yj,yk)
                     IF ( Yo8>=ymin ) THEN
                        ymax = MAX(yi,yj,yk)
                        IF ( Yo8<=ymax ) THEN
!.....              Bis jetzt liegt Punkt innerhalb des das Element
!                   umschlieszenden Rechtecks XMIN/XMAX, YMIN/YMAX
!                   Pruefen, ob Punkt wirklich innerhalb DREIECK-Element
!                   liegt: BERECHNUNG DER TEILFLAECHEN (ohne 0.5)
                           f = xi*(yj-Yo8) + xj*(Yo8-yi) + Xo8*(yi-yj)
                           IF ( f>=0.0_JWRU) THEN
                              f = xj*(yk-Yo8) + xk*(Yo8-yj)+Xo8*(yj-yk)
                              IF ( f>=0.0_JWRU ) THEN
                                 f = xk*(yi-Yo8)+xi*(Yo8-yk)+Xo8*(yk-yi)
                                 IF ( f>=0.0_JWRU ) THEN
                                    IE = i
                                    ielem = IE 
                                    RETURN
                                 endif
                              endif
                           endif
                        endif
                     endif
                  endif
               endif
            endif
         endif
!     SCHLEIFE UEBER ALLE ELEMENTE wird hier folgendermassen hochgezaehlt:
!     beginnend bei IEALT, im Wechsel nach vorn und rueckwaerts suchend
         IF ( k<MNE .AND. if1==0 ) THEN
            if0 = 0
            IF ( l>1 ) if1 = 1
            k = k + 1
            i = k
            IF ( ijk<=MNE ) CYCLE
         endif
         CONTINUE
         EXIT
      ENDDO

      IF ( l>1 .AND. if0==0 ) THEN
         if1 = 0
         IF ( k<MNE ) if0 = 1
         l = l - 1
         i = l
      endif

      IF ( ijk<=MNE ) GOTO 100

      IE = -1

      END SUBROUTINE FIND_ELE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INTELEMENT_IPOL(LLCLST,XYELE,SKALAR,Xo,Yo,IE,PMS8,VAL)

      USE YOWUNPOOL, ONLY : DEGRAD

      IMPLICIT NONE

!     Wrapper for linter interpolation ... 
!     however if LLCLST is true then the value of the closest point is taken.
!     XYELE contain the XY coordinates for the 3 points of the triangle IE
!     Xo, Yo is the point of interest
!     IE is the element number the contains this point
!     SKALAR are the 3 values of interest stored at the nodes of each triangle 
!     VAL is the answer of the interpolation 
!
! Dummy arguments
!
      INTEGER(KIND=JWIM), INTENT(IN)  :: IE

      REAL(KIND=JWRU),  INTENT(IN)  :: XYELE(2,3), SKALAR(3), Xo, Yo, PMS8 
      REAL(KIND=JWRU),  INTENT(OUT) :: VAL

      LOGICAL, INTENT(IN)  :: LLCLST
!
! Local variables

      INTEGER(KIND=JWIM) :: i
      REAL(KIND=JWRU) :: distmin, distmax 
      REAL(KIND=JWRU), DIMENSION(3) :: x, y, dist

      do i=1,3
        x(i) = XYELE(1,i)
        y(i) = XYELE(2,i)
      enddo

      IF (LLCLST) THEN
        do i=1,3
          call SPHERICAL_COORDINATE_DISTANCE(Xo, x(i), Yo, y(i), dist(i))
        enddo
        distmax=maxval(dist)
        do i=1,3
          if(SKALAR(i) == PMS8) dist(i)=distmax
        enddo
        distmin = minval(dist)
        val=SKALAR(1)
        if(dist(2) == distmin) then
          val=SKALAR(2)
        else if(dist(3) == distmin) then
          val=SKALAR(3)
        endif
      ELSE
        call linearInterpolationTriangle(x(1), y(1), SKALAR(1),        &
     &                                   x(2), y(2), SKALAR(2),        &
     &                                   x(3), y(3), SKALAR(3),        &
     &                                   Xo, Yo, val)
      ENDIF

      END SUBROUTINE INTELEMENT_IPOL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine linearInterpolationTriangle(x1, y1, z1, x2, y2,       &
     &                             z2, x3, y3, z3, Xo, Yo, zout)
      !> calc a linear interpolation of a triangle.
      !> describes on this side:
      !> http://www.ems-i.com/smshelp/Data_Module/Interpolation/Linear_Interpolationsms.htm
      !> @param[in] x1, y1, z1 The first Point of the Triangle
      !> @param[in] x2, y2, z2 The second Point of the Triangle
      !> @param[in] x3, y3, z3 The third Point of the Triangle
      !> @param[in] Xo, Yo Point to interpolate
      !> @param[out] zout Interpolated Z-value
        implicit none 
        real(KIND=JWRU), intent(in) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,Xo,Yo
        real(KIND=JWRU), intent(out) :: zout

        real(KIND=JWRU) :: A, B, C, D

        A = y1*(z2 - z3) + y2*(z3 - z1) + y3*(z1 - z2)
        B = z1*(x2 - x3) + z2*(x3 - x1) + z3*(x1 - x2)
        C = x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2)
        D = -A*x1 - B*y1 - C*z1

        zout = - A/C*Xo - B/C*Yo - D/C
      end subroutine linearInterpolationTriangle

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR(FL1,FL3)
      USE MPL_MPIF
      USE YOWFRED  , ONLY : COSTH, SINTH
      USE YOWUNPOOL
      USE yowpd, only: MNE=>ne, MNP=>npa, XP=>x, YP=>y, exchange,Ine, comm
      USE YOWSTAT  , ONLY : IDELPRO
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP
!*--********************************************************************
! calls       AC2      CCON     CG       COSTH    CURTXY   DIFRM
!             EXCHANGE_P2D      EXCHANGE_P3D_WWM  EXCHANGE_P4D_WWM
!             IEN      IE_CELL2 INE      INVSPHTRANS
!             MPI_ALLREDUCE     MY_WTIME POS_CELL2         SI
!             SINTH    SPSIG    WK
! called by   ** NOTHING **
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  CFLXY    COMM     CX       CY       DIFRU    DT4A
!             DT4AI    DTMAX_EXP         DTMAX_GLOBAL_EXP  FL11
!             FL111    FL112    FL12     FL21     FL211    FL212
!             FL22     FL31     FL311    FL312    FL32     FLALL    I
!             I1       I2       I3       ID       IDIFFR   IE
!             IERR     IOBPD    IOBWB    IP       IPOS     IS       IT
!             ITER_MAX IVECTOR  KELEM    KKSUM    KTMP     LAMBDA
!             LCALC    LDIFR    LSECU    LSPHE    LSTCU    MDC
!             MNE      MNP      MPI_MIN  MSC      N        NI
!             NP_RES   ONE      ONEHALF  ONESIXTH REST     RKIND
!             RTYPE    ST       STAT     THR      TIME1    TIME2
!             TMP      TWO      U        U3       USOC     UTILDE
!             UTILDE3  VERYLARGE         VERYSMALL         WILD
!             WILD2D   WVC      ZERO
! uses PARAMs *** NONE ****
!*++********************************************************************
 
      IMPLICIT NONE

      REAL(KIND=JWRB), INTENT(IN) :: FL1(NINF-1:NSUP,NANG,NFRE)
      REAL(KIND=JWRB), INTENT(OUT) :: FL3(NINF-1:NSUP,NANG,NFRE)


! local integer
!
      INTEGER(KIND=JWIM) :: IP, IE, IT, IS, ID
      INTEGER(KIND=JWIM) :: I1, I2, I3, IERR
      INTEGER(KIND=JWIM) :: NI(3), I, IPOS

!
! local double
!

      REAL(KIND=JWRU) :: UTILDE
      REAL(KIND=JWRU) :: DTMAX_GLOBAL_EXP, DTMAX_EXP
      REAL(KIND=JWRU) :: WILD(MNP), WILD2D(MDC,MNP)
      REAL(KIND=JWRU) :: REST, CFLXY
      REAL(KIND=JWRU) :: LAMBDA(2,MSC,MDC), DT4AI
      REAL(KIND=JWRU) :: FL11(MSC,MDC),FL12(MSC,MDC),FL21(MSC,MDC)
      REAL(KIND=JWRU) :: FL22(MSC,MDC),FL31(MSC,MDC),FL32(MSC,MDC)
      REAL(KIND=JWRU) :: KTMP(3,MSC,MDC)
      REAL(KIND=JWRU) :: U3(3)
      REAL(KIND=JWRU) :: KKSUM(MNP,MSC,MDC)
      REAL(KIND=JWRU) :: ST(MNP,MSC,MDC)
      REAL(KIND=JWRU) :: N(MNE,MSC,MDC)
      REAL(KIND=JWRU) :: CX(MSC,MDC,MNP), CY(MSC,MDC,MNP)
      REAL(KIND=JWRU) :: U(MSC,MDC,MNP)
      REAL(KIND=JWRU) :: FLALL(3,MNE,MSC,MDC)
      REAL(KIND=JWRU) :: KELEM(3,MNE,MSC,MDC)
      REAL(KIND=JWRU) :: FL111(MSC,MDC), FL112(MSC,MDC), FL211(MSC,MDC)
      REAL(KIND=JWRU) :: FL212(MSC,MDC), FL311(MSC,MDC), FL312(MSC,MDC)
      REAL(KIND=JWRU) :: UTILDE3(MNE)
      REAL(KIND=JWRU) :: USOC, WVC, DIFRU
      REAL(KIND=JWRU) :: ONEOSIX 
!
! local parameter
!
      REAL(KIND=JWRU) :: TMP(MSC,MDC)
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
         flall = 0.0_JWRU
         kelem = 0.0_JWRU
         kksum = 0.0_JWRU
         st = 0.0_JWRU
         n = 0.0_JWRU
!
! set time step
!
         DT4A = REAL(IDELPRO,JWRU)
!
!
!
         DO IS = 1, MSC
           DO ID = 1, MDC
             DO IP = 1, MNP
               CX(IS,ID,IP) = CG(IP,IS)*SINTH(ID)
               CY(IS,ID,IP) = CG(IP,IS)*COSTH(ID)
               IF (LSPHE) THEN
                 CX(IS,ID,IP) = CX(IS,ID,IP)*                      &
     &                          1.0_JWRU/(DEGRAD*REARTH*COS(YP(IP)*DEGRAD))
                 CY(IS,ID,IP) = CY(IS,ID,IP)*                      &
     &                          1.0_JWRU/(DEGRAD*REARTH)
               END IF
             END DO
           END DO
         END DO
!
!        Calculate K-Values and contour based quantities ...

         ONEOSIX = 1.0_JWRU/6.0_JWRU
!
!!$OMP DO PRIVATE(IE,I1,I2,I3,LAMBDA,KTMP,TMP,FL11,FL12,FL21,FL22,FL31,FL32,FL111,FL112,FL211,FL212,FL311,FL312)
      DO ie = 1 , mne
         i1 = INE(1,ie)
         i2 = INE(2,ie)
         i3 = INE(3,ie)
         lambda(1,:,:) = ONEOSIX*(cx(:,:,i1)+cx(:,:,i2)+cx(:,:,i3))
         lambda(2,:,:) = ONEOSIX*(cy(:,:,i1)+cy(:,:,i2)+cy(:,:,i3))
         kelem(1,ie,:,:) = lambda(1,:,:)*IEN(1,ie) + lambda(2,:,:)    &
     &                     *IEN(2,ie)
         kelem(2,ie,:,:) = lambda(1,:,:)*IEN(3,ie) + lambda(2,:,:)    &
     &                     *IEN(4,ie)
         kelem(3,ie,:,:) = lambda(1,:,:)*IEN(5,ie) + lambda(2,:,:)    &
     &                     *IEN(6,ie)
         ktmp(1,:,:) = kelem(1,ie,:,:)
         ktmp(2,:,:) = kelem(2,ie,:,:)
         ktmp(3,:,:) = kelem(3,ie,:,:)
         tmp(:,:) = SUM(MIN(0.0_JWRU,ktmp(:,:,:)),DIM=1)
         n(ie,:,:) = -1.0_JWRU/MIN(-thr8,tmp(:,:))
         kelem(1,ie,:,:) = MAX(0.0_JWRU,ktmp(1,:,:))
         kelem(2,ie,:,:) = MAX(0.0_JWRU,ktmp(2,:,:))
         kelem(3,ie,:,:) = MAX(0.0_JWRU,ktmp(3,:,:))
!            WRITE(DBG%FHNDL,'(3I10,3F15.4)') IS, ID, IE, KELEM(:,IE)
         fl11 = cx(:,:,i2)*IEN(1,ie) + cy(:,:,i2)*IEN(2,ie)
         fl12 = cx(:,:,i3)*IEN(1,ie) + cy(:,:,i3)*IEN(2,ie)
         fl21 = cx(:,:,i3)*IEN(3,ie) + cy(:,:,i3)*IEN(4,ie)
         fl22 = cx(:,:,i1)*IEN(3,ie) + cy(:,:,i1)*IEN(4,ie)
         fl31 = cx(:,:,i1)*IEN(5,ie) + cy(:,:,i1)*IEN(6,ie)
         fl32 = cx(:,:,i2)*IEN(5,ie) + cy(:,:,i2)*IEN(6,ie)
         fl111 = 2.0_JWRU*fl11 + fl12
         fl112 = 2.0_JWRU*fl12 + fl11
         fl211 = 2.0_JWRU*fl21 + fl22
         fl212 = 2.0_JWRU*fl22 + fl21
         fl311 = 2.0_JWRU*fl31 + fl32
         fl312 = 2.0_JWRU*fl32 + fl31
         flall(1,ie,:,:) = (fl311+fl212)*ONEOSIX + kelem(1,ie,:,:)
         flall(2,ie,:,:) = (fl111+fl312)*ONEOSIX + kelem(2,ie,:,:)
         flall(3,ie,:,:) = (fl211+fl112)*ONEOSIX + kelem(3,ie,:,:)
      ENDDO
 
      IF ( lcalc ) THEN

#ifdef ebug_adv
         write(dbg%fhndl,*) 'starting the ks buisness'
#endif
         kksum = 0.0_JWRU
         DO ie = 1 , mne
            ni = INE(:,ie)
            kksum(ni(1),:,:) = kksum(ni(1),:,:) + kelem(1,ie,:,:)
            kksum(ni(2),:,:) = kksum(ni(2),:,:) + kelem(2,ie,:,:)
            kksum(ni(3),:,:) = kksum(ni(3),:,:) + kelem(3,ie,:,:)
         ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) maxval(kksum), minval(kksum)
         write(dbg%fhndl,*) 'finished the ks buisness'
#endif

         IF ( ivector.EQ.1 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            DO id = 1 , mdc
               DO is = 1 , msc
                  dtmax_global_exp = 10.E10_JWRU 
                  DO ip = 1 , mnp
                     dtmax_exp = SI(ip)/MAX(thr8,kksum(ip,is,id))
                     dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
                  ENDDO
                  dtmax_exp = dtmax_global_exp
                  CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
                  cflxy = dt4a/dtmax_global_exp
                  rest = ABS(MOD(cflxy,1.0_JWRU))
                  IF ( rest.LT.thr8 ) THEN
                     ITER_EXP(is,id) = ABS(NINT(cflxy))
                  ELSEIF ( rest.GT.thr8 .AND. rest.LT.0.5_JWRU ) THEN
                     ITER_EXP(is,id) = ABS(NINT(cflxy)) + 1
                  ELSE
                     ITER_EXP(is,id) = ABS(NINT(cflxy))
                  ENDIF
               ENDDO
            ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', maxval(ITER_EXP)
#endif
         ELSEIF ( ivector.EQ.2 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            dtmax_global_exp = 10.E10
            DO ip = 1 , mnp 
               dtmax_exp = SI(ip)/MAX(thr8,MAXVAL(kksum(ip,:,:)))
               dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
            ENDDO
            dtmax_exp = dtmax_global_exp
            CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
            cflxy = dt4a/dtmax_global_exp
            rest = ABS(MOD(cflxy,1.0_JWRU))
            IF ( rest.LT.thr8 ) THEN
               iter_max = ABS(NINT(cflxy))
            ELSEIF ( rest.GT.thr8 .AND. rest.LT.0.5_JWRU ) THEN
               iter_max = ABS(NINT(cflxy)) + 1
            ELSE
               iter_max = ABS(NINT(cflxy))
            ENDIF
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', iter_max
#endif
         ELSEIF ( ivector.EQ.3 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            DO is = 1 , msc
               dtmax_global_exp = 10.E10_JWRU
               DO ip = 1 , mnp 
                  dtmax_exp = SI(ip)/MAX(thr8,MAXVAL(kksum(ip,is,:)))
                  dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
               ENDDO
               dtmax_exp = dtmax_global_exp
               CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
               cflxy = dt4a/dtmax_global_exp
               rest = ABS(MOD(cflxy,1.0_JWRU))
               IF ( rest.LT.thr8 ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ELSEIF ( rest.GT.thr8 .AND. rest.LT.0.5_JWRU ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy)) + 1
               ELSE
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ENDIF
            ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', maxval(iter_expd)
#endif
         ELSEIF ( ivector.EQ.4 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_max', ivector
#endif
            DO is = 1 , msc
               dtmax_global_exp = 10.E10
               DO ip = 1 , mnp 
                  dtmax_exp = SI(ip)/MAX(thr8,MAXVAL(kksum(ip,is,:)))
                  dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
               ENDDO
               dtmax_exp = dtmax_global_exp
               CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,mpi_real8,mpi_min,comm,ierr)
               cflxy = dt4a/dtmax_global_exp
               rest = ABS(MOD(cflxy,1.0_JWRU))
               IF ( rest.LT.thr8 ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ELSEIF ( rest.GT.thr8 .AND. rest.LT.0.5_JWRU ) THEN
                  ITER_EXPD(is) = ABS(NINT(cflxy)) + 1
               ELSE
                  ITER_EXPD(is) = ABS(NINT(cflxy))
               ENDIF
            ENDDO
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', maxval(iter_expd) 
#endif
         ELSEIF ( ivector.EQ.5 ) THEN
#ifdef ebug_adv
         write(dbg%fhndl,*) 'estimate iter_exp', lvector 
#endif
            dtmax_global_exp = 10.E10_JWRU
            DO ip = 1 , mnp 
               dtmax_exp = SI(ip)/MAX(thr8,MAXVAL(kksum(ip,:,:)))
               dtmax_global_exp = MIN(dtmax_global_exp,dtmax_exp)
            ENDDO
            dtmax_exp = dtmax_global_exp
            CALL MPI_ALLREDUCE(dtmax_exp,dtmax_global_exp,1,MPI_REAL8,mpi_min,comm,ierr)
            cflxy = dt4a/dtmax_global_exp
            rest = ABS(MOD(cflxy,1.0_JWRU))
            IF ( rest.LT.thr8 ) THEN
               iter_max = ABS(NINT(cflxy))
            ELSEIF ( rest.GT.thr8 .AND. rest.LT.0.5_JWRU ) THEN
               iter_max = ABS(NINT(cflxy)) + 1
            ELSE
               iter_max = ABS(NINT(cflxy))
            ENDIF
#ifdef ebug_adv
         write(dbg%fhndl,*) 'iter_max = ', iter_max
#endif
         ENDIF    !IVECTOR
         CALL FLUSH(dbg%fhndl)
      ENDIF     !LCALC


      DO ip = 1 , mnp
         DO is = 1 , msc
            DO id = 1 , mdc
               u(is,id,ip) = FL1(ip,id,is)
            ENDDO
         ENDDO
      ENDDO

#ifdef ebug_adv
      write(dbg%fhndl,*) 'doing time integration'
      call cpu_time(time1)
#endif

      IF ( ivector.EQ.1 ) THEN
         DO id = 1 , mdc
            DO is = 1 , msc
               dt4ai = dt4a/ITER_EXP(is,id)
               DO it = 1 , ITER_EXP(is,id)
                  st(:,is,id) = 0.0_JWRU
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(is,id,ni)
                     utilde = n(ie,is,id)                                  &
     &                        *(flall(1,ie,is,id)*u3(1)+flall(2,ie,is,     &
     &                        id)*u3(2)+flall(3,ie,is,id)*u3(3))
                     st(ni(1),is,id) = st(ni(1),is,id)                     &
     &                                 + kelem(1,ie,is,id)                 &
     &                                 *(u3(1)-utilde)
                     st(ni(2),is,id) = st(ni(2),is,id)                     &
     &                                 + kelem(2,ie,is,id)                 &
     &                                 *(u3(2)-utilde)
                     st(ni(3),is,id) = st(ni(3),is,id)                     &
     &                                 + kelem(3,ie,is,id)                 &
     &                                 *(u3(3)-utilde)
                  ENDDO !IE
                  u(is,id,:) = MAX(0.0_JWRU,u(is,id,:)-dt4ai/SI*st(:,is,id)    &
     &                         *iobwb)*iobpd(id,:)
                  wild = u(is,id,:)
                  CALL exchange(wild)
                  u(is,id,:) = wild
               ENDDO ! IT----> End Iteration
            ENDDO !IS
         ENDDO  !ID
      ELSEIF ( ivector.EQ.2 ) THEN
         dt4ai = dt4a/iter_max
         DO it = 1 , iter_max
            DO id = 1 , mdc
               DO is = 1 , msc
                  st(:,is,id) = 0.0_JWRU
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(is,id,ni)
                     utilde = n(ie,is,id)                                  &
     &                        *(flall(1,ie,is,id)*u3(1)+flall(2,ie,is,     &
     &                        id)*u3(2)+flall(3,ie,is,id)*u3(3))
                     st(ni(1),is,id) = st(ni(1),is,id)                     &
     &                                 + kelem(1,ie,is,id)                 &
     &                                 *(u3(1)-utilde)
                     st(ni(2),is,id) = st(ni(2),is,id)                     &
     &                                 + kelem(2,ie,is,id)                 &
     &                                 *(u3(2)-utilde)
                     st(ni(3),is,id) = st(ni(3),is,id)                     &
     &                                 + kelem(3,ie,is,id)                 &
     &                                 *(u3(3)-utilde)
                  ENDDO!IE
                  u(is,id,:) = MAX(0.0_JWRU,u(is,id,:)-dt4ai/SI*st(:,is,id)    &
     &                         *iobwb)*iobpd(id,:)
               ENDDO !IS
            ENDDO !ID
            CALL exchange(u)
         ENDDO  !IT
      ELSEIF ( ivector.EQ.3 ) THEN
         DO is = 1 , msc
            iter_max = ITER_EXPD(is)
            dt4ai = dt4a/iter_max
            DO it = 1 , iter_max
               DO id = 1 , mdc
                  st(:,is,id) = 0.0_JWRU
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(is,id,ni)
                     utilde = n(ie,is,id)                                  &
     &                        *(flall(1,ie,is,id)*u3(1)+flall(2,ie,is,     &
     &                        id)*u3(2)+flall(3,ie,is,id)*u3(3))
                     st(ni(1),is,id) = st(ni(1),is,id)                     &
     &                                 + kelem(1,ie,is,id)                 &
     &                                 *(u3(1)-utilde)
                     st(ni(2),is,id) = st(ni(2),is,id)                     &
     &                                 + kelem(2,ie,is,id)                 &
     &                                 *(u3(2)-utilde)
                     st(ni(3),is,id) = st(ni(3),is,id)                     &
     &                                 + kelem(3,ie,is,id)                 &
     &                                 *(u3(3)-utilde)
                  ENDDO !IE
                  u(is,id,:) = MAX(0.0_JWRU,u(is,id,:)-dt4ai/SI*st(:,is,id)    &
     &                         *iobwb)*iobpd(id,:)
               ENDDO ! ID
               wild2d = u(is,:,:)
               CALL exchange(wild2d)
               u(is,:,:) = wild2d
            ENDDO !IT
         ENDDO  !IS
      ELSEIF ( ivector.EQ.4 ) THEN
         DO is = 1 , msc
            iter_max = ITER_EXPD(is)
            dt4ai = dt4a/iter_max
            DO it = 1 , iter_max
               DO id = 1 , mdc
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(is,id,ni)
                     utilde3(ie) = n(ie,is,id)                             &
     &                             *(flall(1,ie,is,id)*u3(1)+flall(2,ie,   &
     &                             is,id)*u3(2)+flall(3,ie,is,id)*u3(3))
                  ENDDO
                  st(:,is,id) = 0.0_JWRU
                  DO ip = 1 , mnp
                     DO i = 1 , CCON(ip)
                        ie = IE_CELL2(ip,i)
                        ipos = POS_CELL2(ip,i)
                        st(ip,is,id) = st(ip,is,id)                        &
     &                                 + kelem(ipos,ie,is,id)              &
     &                                 *(u(is,id,ip)-utilde3(ie))
                     ENDDO
                     u(is,id,ip) = MAX(0.0_JWRU,u(is,id,ip)-dt4ai/SI(ip)*st(   &
     &                             ip,is,id)*iobwb(ip))*iobpd(id,ip)
                  ENDDO
               ENDDO
               wild2d = u(is,:,:)
               CALL exchange(wild2d)
               u(is,:,:) = wild2d
            ENDDO !IT
         ENDDO  !IS
      ELSEIF ( ivector.EQ.5 ) THEN
         dt4ai = dt4a/iter_max
         DO it = 1 , iter_max
            DO is = 1 , msc
               DO id = 1 , mdc
                  DO ie = 1 , mne
                     ni = INE(:,ie)
                     u3(:) = u(is,id,ni)
                     utilde3(ie) = n(ie,is,id)                             &
     &                             *(flall(1,ie,is,id)*u3(1)+flall(2,ie,   &
     &                             is,id)*u3(2)+flall(3,ie,is,id)*u3(3))
                  ENDDO
                      !IE
                  st(:,is,id) = 0.0_JWRU
                  DO ip = 1 , mnp
                     DO i = 1 , CCON(ip)
                        ie = IE_CELL2(ip,i)
                        ipos = POS_CELL2(ip,i)
                        st(ip,is,id) = st(ip,is,id)                        &
     &                                 + kelem(ipos,ie,is,id)              &
     &                                 *(u(is,id,ip)-utilde3(ie))
                     ENDDO
                     u(is,id,ip) = MAX(0.0_JWRU,u(is,id,ip)-dt4ai/SI(ip)*st(   &
     &                             ip,is,id)*iobwb(ip))*iobpd(id,ip)
                  ENDDO
                      !IP
               ENDDO
                    !ID
            ENDDO !IS
            CALL exchange(u)
         ENDDO  !IT
      ENDIF
 
      DO ip = 1 , mnp
         DO is = 1 , msc
            DO id = 1 , mdc
               FL3(ip,id,is) = u(is,id,ip)
            ENDDO
         ENDDO
      ENDDO

#ifdef ebug_adv
      call cpu_time(time2)
      write(dbg%fhndl,*) 'end time integration', time2-time1
      call flush(dbg%fhndl)
#endif
 
      END SUBROUTINE EXPLICIT_N_SCHEME_VECTOR

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_ASPAR_BLOCK(ASPAR)
      USE YOWUNPOOL, ONLY : LCUR, LSPHE, CG, CURTXY, DEGRAD, REARTH, &
&                           ONETHIRD, ONESIXTH, IEN, ONE, TWO, THR, TRIA, MDC, &
&                           DT4A, IOBPD, IOBWB, LBCWA, IWBMNP, LINHOM, IWBNDLC, &
&                           SI
      USE yowpd, only: XP=>x, YP=>y
      USE yowpd, only: MNE=>ne, INE, MNP=>npa
      USE YOWFRED  , ONLY : COSTH, SINTH
      USE YOWPARAM , ONLY : NANG, NFRE

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: POS_TRICK(3,2)
      INTEGER(KIND=JWIM) :: I1, I2, I3
      INTEGER(KIND=JWIM) :: IP, ID, IS, IE
      INTEGER(KIND=JWIM) :: I, IPGL1, IPrel

      REAL(KIND=JWRU), INTENT(INOUT) :: ASPAR(NANG, NFRE, NNZ)
      REAL(KIND=JWRU) :: FL11(NANG,NFRE), FL12(NANG,NFRE), FL21(NANG,NFRE), FL22(NANG,NFRE), FL31(NANG,NFRE), FL32(NANG,NFRE)
      REAL(KIND=JWRU) :: CRFS(NANG,NFRE,3), K1(NANG,NFRE), KM(NANG,NFRE,3), K(NANG,NFRE,3), TRIA03
      REAL(KIND=JWRU) :: CXY(2,NANG,NFRE,3)
      REAL(KIND=JWRU) :: DIFRU, USOC, WVC
      REAL(KIND=JWRU) :: DELTAL(NANG,NFRE,3)
      REAL(KIND=JWRU) :: KP(NANG,NFRE,3), NM(NANG,NFRE)
      REAL(KIND=JWRU) :: DTK(NANG,NFRE), TMP3(NANG,NFRE)
      REAL(KIND=JWRU) :: LAMBDA(2,NANG,NFRE)

      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
!
!     Calculate countour integral quantities ...
!
      ASPAR = 0.0_JWRU
      DO IE = 1, MNE
        DO I=1,3
          IP = INE(I,IE)
          DO IS=1,NFRE
            DO ID=1,NANG
              IF (LCUR) THEN
                CXY(1,ID,IS,I) = CG(IP,IS)*COSTH(ID)+CURTXY(IP,1)
                CXY(2,ID,IS,I) = CG(IP,IS)*SINTH(ID)+CURTXY(IP,2)
              ELSE
                CXY(1,ID,IS,I) = CG(IP,IS)*COSTH(ID)
                CXY(2,ID,IS,I) = CG(IP,IS)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,ID,IS,I) = CXY(1,ID,IS,I)/(DEGRAD*REARTH*COS(MIN(ABS(YP(IP)),YPMAX)*DEGRAD))
                CXY(2,ID,IS,I) = CXY(2,ID,IS,I)/(DEGRAD*REARTH)
              END IF
            END DO
          END DO
        END DO

        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(0.0_JWRU,K)
        KP(:,:,:) = MAX(0.0_JWRU,K)
        DELTAL(:,:,:) = CRFS(:,:,:)- KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        DO I=1,3
          IP=INE(I,IE)
          I1=JA_IE(I,1,IE)
          I2=JA_IE(I,2,IE)
          I3=JA_IE(I,3,IE)
          K1(:,:) =  KP(:,:,I)
          DO ID=1,MDC
            DTK(ID,:) =  K1(ID,:) * DT4A * IOBPD(ID,IP) * IOBWB(IP)
          END DO
          TMP3(:,:)  =  DTK(:,:) * NM(:,:)
          ASPAR(:,:,I1) =  TRIA03 + DTK(:,:) - TMP3(:,:) * DELTAL(:,:,I             ) + ASPAR(:,:,I1)
          ASPAR(:,:,I2) =                    - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,1)) + ASPAR(:,:,I2)
          ASPAR(:,:,I3) =                    - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,2)) + ASPAR(:,:,I3)
        END DO
      END DO
      IF (LBCWA) THEN
        DO IP = 1, IWBMNP
          IF (LINHOM) THEN
            IPrel=IP
          ELSE
            IPrel=1
          ENDIF
          IPGL1 = IWBNDLC(IP)
          ASPAR(:,:,I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
        END DO
      END IF
      END SUBROUTINE EIMPS_ASPAR_BLOCK

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR_JAC)
      USE YOWPARAM , ONLY : NANG, NFRE
      IMPLICIT NONE
      REAL(KIND=JWRU), intent(inout) :: ASPAR_JAC(NANG,NFRE,NNZ)
      IF (REFRACTION_IMPL) THEN
        Print *, 'This needs to be written down (Refraction)'
        stop
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        Print *, 'This needs to be written down (Freq Shift)'
        stop
      END IF
      END SUBROUTINE ADD_FREQ_DIR_TO_ASPAR_COMP_CADS

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
      USE YOWUNPOOL
      USE yowpd, only : np_global, NP_RES => np, MNP=>npa
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), INTENT(IN)  :: AC1(NANG,NFRE,MNP)
      REAL(KIND=JWRU), intent(out) :: IMATRA(NANG,NFRE)
      REAL(KIND=JWRU), intent(out) :: IMATDA(NANG,NFRE)
      REAL(KIND=JWRU) eVal
      IF (LNONL) THEN
         !*    Advanced computation, likely call to IMPLSCH
      ELSE
         !*    Put here the fields for IMATRA
      END IF
      eVal=SI(IP)
      IMATRA = IMATRA * eVal
      IMATDA = IMATDA * eVal
      END SUBROUTINE GET_IMATRA_IMATDA

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ACTION_LIMITER_local(IP,eSum,acloc)
      USE YOWPARAM , ONLY : NANG, NFRE
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), intent(INOUT) :: eSum(NANG,NFRE)
      REAL(KIND=JWRU), intent(IN) :: acloc(NANG,NFRE)


      Print *, 'This needs to be written'
      stop
      END SUBROUTINE ACTION_LIMITER_local

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BLOCAL(IP, AC, BLOC)
      USE YOWUNPOOL
      USE yowpd, only : np_global, NP_RES => np, MNP=>npa
      USE YOWPARAM , ONLY : NANG, NFRE
      IMPLICIT NONE
      INTEGER(KIND=JWIM), INTENT(IN) :: IP
      REAL(KIND=JWRU), INTENT(IN)  :: AC(NANG,NFRE,MNP)
      REAL(KIND=JWRU), intent(out) :: BLOC(NANG,NFRE)
      IF (LBCWA) THEN
           !* the boundary if any
      ELSE
        BLOC=AC(:,:,IP)
      END IF
      END SUBROUTINE GET_BLOCAL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IMPLICIT_N_SCHEME_BLOCK(FL1, FL3)
      USE MPL_MPIF
      USE yowpd, only : comm
      USE yowpd, only : np_global, NP_RES => np, MNP=>npa
      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP
      USE yowunpool
      USE yowpd, only: exchange
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: IS, ID, ID1, ID2, IP, J, idx, nbITer, TheVal, is_converged, itmp
      INTEGER(KIND=JWIM) :: I, K, IP_ADJ, IADJ, JDX
      INTEGER(KIND=JWIM) :: ierr

      REAL(KIND=JWRB), INTENT(IN)  :: FL1(NINF-1:NSUP,NANG,NFRE)
      REAL(KIND=JWRB), INTENT(OUT) :: FL3(NINF-1:NSUP,NANG,NFRE)

      REAL(KIND=JWRU) :: MaxNorm, p_is_converged
      REAL(KIND=JWRU) :: eSum(NANG,NFRE)
      REAL(KIND=JWRU) :: IMATRA(NANG,NFRE), IMATDA(NANG,NFRE)
      REAL(KIND=JWRU) :: Norm_L2(NANG,NFRE), Norm_LINF(NANG,NFRE)
      REAL(KIND=JWRU) :: ACLOC(NANG,NFRE)
      REAL(KIND=JWRU) :: CAD(NANG,NFRE), CAS(NANG,NFRE)
      REAL(KIND=JWRU) :: CP_THE(NANG,NFRE), CM_THE(NANG,NFRE)
      REAL(KIND=JWRU) :: CP_SIG(NANG,NFRE), CM_SIG(NANG,NFRE)
      REAL(KIND=JWRU) :: BLOC(NANG,NFRE)
      REAL(KIND=JWRU) :: ASPAR_DIAG(NANG,NFRE)
      REAL(KIND=JWRU) :: A_THE(NANG,NFRE), C_THE(NANG,NFRE)
      REAL(KIND=JWRU) :: A_SIG(NANG,NFRE), C_SIG(NANG,NFRE)
      REAL(KIND=JWRU) :: Norm_L2_gl(NANG,NFRE), Norm_LINF_gl(NANG,NFRE)
      REAL(KIND=JWRU) :: B_SIG(NFRE), eFact, lambda
      REAL(KIND=JWRU) :: Sum_new, Sum_prev, eVal, DiffNew, DiffOld
      REAL(KIND=JWRU) :: AC1(NANG,NFRE,MNP)
      REAL(KIND=JWRU) :: AC2(NANG,NFRE,MNP)
#ifdef DEBUG
      WRITE(740+MyRankGlobal,*) 'JWRU=', JWRU
      WRITE(740+MyRankGlobal,*) 'JWRB=', JWRB
      FLUSH(740+MyRankGlobal)
#endif

      DO IP=1,MNP
        AC2(:,:,IP)=FL1(IP,:,:)
      END DO

      CALL EIMPS_ASPAR_BLOCK(ASPAR_JAC)
      !
      CALL ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR_JAC)
      IF ((.NOT. LNONL) .AND. SOURCE_IMPL) THEN
        DO IP=1,NP_RES
          CALL GET_BLOCAL(IP, AC1, BLOC)
          CALL GET_IMATRA_IMATDA(IP, AC1, IMATRA, IMATDA)
          ASPAR_JAC(:,:,I_DIAG(IP)) = ASPAR_JAC(:,:,I_DIAG(IP)) + IMATDA
          B_JAC(:,:,IP)             = BLOC + IMATRA
        END DO
      END IF
      !
      ! Now the Gauss Seidel iterations
      !
      !SOLVERTHR=10E-8*AVETL!*TLMIN**2
      !
      nbIter=0
      DO
        is_converged = 0
        JDX=0
        DO IP=1,NP_RES
          ACLOC = AC2(:,:,IP)
          Sum_prev = sum(ACLOC)
          ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BLOCAL(IP, AC2, BLOC)
              CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
              ASPAR_DIAG = ASPAR_DIAG + IMATDA
              eSum = BLOC + IMATRA
            ELSE
              eSum = B_JAC(:,:,IP)
            END IF
          ELSE
            CALL GET_BLOCAL(IP, AC2, eSum)
          END IF
          IF (LNANINFCHK) THEN
            IF (SUM(eSum) .ne. SUM(esum)) THEN
              Print *, 'NAN IN SOLVER After assignation of values'
              STOP
            ENDIF
          ENDIF
          DO J=IA(IP),IA(IP+1)-1 
            IF (J .ne. I_DIAG(IP)) eSum = eSum - ASPAR_JAC(:,:,J) * AC2(:,:,JA(J))
          END DO
          IF (LNANINFCHK) THEN
            IF (SUM(eSum) .ne. SUM(esum)) THEN
              Print *, 'NAN IN SOLVER After substraction of entries'
              STOP
            ENDIF
          ENDIF
          IF (REFRACTION_IMPL) THEN
            CAD=CAD_THE(:,:,IP)
            CP_THE = MAX(ZERO,CAD)
            CM_THE = MIN(ZERO,CAD)
            eFact=(DT4D/DDIR)*SI(IP)
            DO ID=1,MDC
              IF (ID .gt. 1) THEN
                ID1=ID-1
              ELSE
                ID1=MDC
              END IF
              IF (ID .lt. MDC) THEN
                ID2=ID + 1
              ELSE
                ID2=1
              END IF
              eSum(ID,:) = eSum(ID,:) + eFact*CP_THE(ID1,:)*ACLOC(ID1,:)
              eSum(ID,:) = eSum(ID,:) - eFact*CM_THE(ID2,:)*ACLOC(ID2,:)
            END DO
          END IF
          IF (FREQ_SHIFT_IMPL) THEN
            CAS=CAS_SIG(:,:,IP)
            CP_SIG = MAX(ZERO,CAS)
            CM_SIG = MIN(ZERO,CAS)
            eFact=DT4F*SI(IP)
            DO ID=1,MDC
              DO IS=2,NFRE
                eSum(ID,IS)=eSum(ID,IS) + eFact*(CP_SIG(IS,IS-1)/DS_INCR(IS-1))*ACLOC(ID,IS-1)
              END DO
              DO IS=1,NFRE-1
                eSum(ID,IS)=eSum(ID,IS) - eFact*(CM_SIG(ID,IS+1)/DS_INCR(IS))*ACLOC(ID,IS+1)
              END DO
            END DO
          END IF
          IF (LNANINFCHK) THEN
            IF (SUM(eSum) .ne. SUM(esum)) THEN
              Print *, 'NAN IN SOLVER before division by ASPAR_DIAG'
              STOP
            ENDIF
          ENDIF
          eSum=eSum/ASPAR_DIAG
          IF (LLIMT) CALL ACTION_LIMITER_LOCAL(IP,eSum,acloc)
          IF (BLOCK_GAUSS_SEIDEL) THEN
            AC2(:,:,IP)=eSum
            IF (LNANINFCHK) THEN
              IF (SUM(eSum) .ne. SUM(esum)) THEN
                Print *, 'NAN IN SOLVER before assignation'
                STOP
              ENDIF
            ENDIF
          ELSE
            U_JACOBI(:,:,IP)=eSum
          END IF
          IF (LCHKCONV) THEN
            Sum_new = sum(eSum)
            IF (Sum_new .gt. thr8) THEN
              DiffNew=sum(abs(ACLOC - eSum))
              p_is_converged = DiffNew/Sum_new
            ELSE
              p_is_converged = zero
            END IF
            IF (p_is_converged .lt. solverthr) is_converged=is_converged+1
          ENDIF
        END DO
        IF (LCHKCONV) THEN
          CALL MPI_ALLREDUCE(is_converged, itmp, 1, MPI_INTEGER, MPI_SUM, COMM, ierr)
          is_converged = itmp
          p_is_converged = (real(np_global) - real(is_converged))/real(np_global) * 100.
        ENDIF
        IF (BLOCK_GAUSS_SEIDEL) THEN
          CALL exchange(AC2)
        ELSE
          CALL exchange(U_JACOBI)
        END IF
        IF (.NOT. BLOCK_GAUSS_SEIDEL) THEN
          AC2 = U_JACOBI
        ENDIF
!
! The termination criterions several can be chosen
!
        WRITE(STAT%FHNDL,'(A10,3I10,E30.20,F10.5)') 'solver', nbiter, is_converged, np_global - is_converged, p_is_converged, pmin
        !
        ! Number of iterations. If too large the exit.
        !
        nbIter=nbIter+1
        IF (nbiter .eq. maxiter) THEN
          EXIT
        ENDIF
        !
        ! Check via number of converged points
        !
        IF (LCHKCONV) THEN
          IF (p_is_converged .le. pmin) EXIT
        ENDIF
        !
        ! Check via the norm
        !
        IF (L_SOLVER_NORM) THEN
          Norm_L2=0
          DO IP=1,NP_RES
            ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
            IF (SOURCE_IMPL) THEN
              IF (LNONL) THEN
                CALL GET_BLOCAL(IP, AC2, BLOC)
                CALL GET_IMATRA_IMATDA(IP, AC2, IMATRA, IMATDA)
                ASPAR_DIAG = ASPAR_DIAG + IMATDA
                eSum = BLOC + IMATRA
              ELSE
                eSum = B_JAC(:,:,IP)
              END IF
            ELSE
              CALL GET_BLOCAL(IP, AC2, eSum)
            END IF
            DO J=IA(IP),IA(IP+1)-1
              idx=JA(J)
              IF (J .eq. I_DIAG(IP)) THEN
                eSum=eSum - ASPAR_DIAG*AC2(:,:,idx)
              ELSE
                eSum=eSum - ASPAR_JAC(:,:,J)*AC2(:,:,idx)
              END IF
            END DO
            IF (REFRACTION_IMPL) THEN
              CAD=CAD_THE(:,:,IP)
              CP_THE = MAX(ZERO,CAD)
              CM_THE = MIN(ZERO,CAD)
              eFact=(DT4D/DDIR)*SI(IP)
              DO ID=1,MDC
                IF (ID .gt. 1) THEN
                  ID1=ID-1
                ELSE
                  ID1=MDC
                END IF
                IF (ID .lt. MDC) THEN
                  ID2=ID + 1
                ELSE
                  ID2=1
                END IF
                eSum(ID,:) = eSum(ID,:) + eFact*CP_THE(ID1,:)*AC2(ID1,:,IP)
                eSum(ID,:) = eSum(ID,:) - eFact*CM_THE(ID2,:)*AC2(ID2,:,IP)
              END DO
            END IF
            IF (FREQ_SHIFT_IMPL) THEN
              CAS=CAS_SIG(:,:,IP)
              CP_SIG = MAX(ZERO,CAS)
              CM_SIG = MIN(ZERO,CAS)
              eFact=DT4F*SI(IP)
              DO ID=1,MDC
                DO IS=2,NFRE
                  eSum(ID,IS)=eSum(ID,IS) + eFact*(CP_SIG(ID,IS-1)/DS_INCR(IS-1))*AC2(ID,IS-1,IP)
                END DO
                DO IS=1,NFRE-1
                  eSum(ID,IS)=eSum(ID,IS) - eFact*(CM_SIG(ID,IS+1)/DS_INCR(IS  ))*AC2(ID,IS+1,IP)
                END DO
              END DO
            END IF
            Norm_L2 = Norm_L2 + (eSum**2)
            Norm_LINF = max(Norm_LINF, abs(eSum))
          END DO
          CALL MPI_ALLREDUCE(Norm_LINF, Norm_LINF_gl, NFRE*MDC,mpi_real8,MPI_MAX,comm,ierr)
          CALL MPI_ALLREDUCE(Norm_L2, Norm_L2_gl, NFRE*MDC, mpi_real8,MPI_SUM,comm,ierr)
          MaxNorm=maxval(Norm_L2_gl)
        END IF
      END DO
      DO IP=1,MNP
        FL3(IP,:,:)=MAX(ZERO,AC2(:,:,IP))
      END DO
      END SUBROUTINE IMPLICIT_N_SCHEME_BLOCK

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXPLICIT_N_SCHEME_VECTOR_HPCF(FL1,FL3)
      USE MPL_MPIF
      USE YOWFRED  , ONLY : COSTH, SINTH
      USE YOWUNPOOL
      USE yowpd, ONLY : XP=>x, YP=>y, comm, ine, mnp => npa, exchange
      USE YOWSTAT  , ONLY : IDELPRO
      USE YOWPARAM , ONLY : NANG, NFRE
      USE yowexchangemodule
      IMPLICIT NONE
!*--********************************************************************
! calls       AC2      CCON     CG       COSTH    CURTXY
!             EXCHANGE_P4D_WWM  IEN      IE_CELL2 INVSPHTRANS
!             IOBPD    IOBWB    MPI_ALLREDUCE     MY_WTIME POS_CELL2
!             SI       SINTH
! called by   ** NOTHING **
! modifies    ** NOTHING **
! uses value  ** NOTHING **
! local vars  CFLXY    COMM     CX       CY       DT4A     DT4AI
!             DTMAX_EXP         DTMAX_GLOBAL_EXP  DTMAX_GLOBAL_EXP_LOC
!             FL11     FL111    FL112    FL12     FL21     FL211
!             FL212    FL22     FL31     FL311    FL312    FL32
!             FLALL    I        I1       I2       I3       ID       IE
!             IERR     INE      IP       IPOS     IS       IT
!             ITER_MAX KELEM    KKSUM    KMAX     KSUM     LAMBDA
!             LCALC    MDC      MNP      MPI_MIN  MSC      N        NI
!             NP_RES   ONE      ONEHALF  ONESIXTH REST     RKIND
!             RTYPE    ST       STAT     THR      TIME1    TIME2
!             TWO      U        U3       UIP      UTILDE3  VERYLARGE
!             ZERO
! uses PARAMs *** NONE ****
!*++********************************************************************
!
! local integer
!
      INTEGER(KIND=JWIM) :: IP, IE, IT, IS, ID
      INTEGER(KIND=JWIM) :: I1, I2, I3
      INTEGER(KIND=JWIM) :: NI(3), I, IPOS, IERR
!
! local double
!
      REAL(KIND=JWRB), INTENT(IN) :: FL1(:,:,:)
      REAL(KIND=JWRB), INTENT(OUT) :: FL3(:,:,:)
      REAL(KIND=JWRU)  :: DTMAX_GLOBAL_EXP, DTMAX_EXP
      REAL(KIND=JWRU)  :: DTMAX_GLOBAL_EXP_LOC
      REAL(KIND=JWRU)  :: REST, CFLXY
      REAL(KIND=JWRU)  :: LAMBDA(2), DT4AI
      REAL(KIND=JWRU)  :: FL11,FL12,FL21,FL22,FL31,FL32
      REAL(KIND=JWRU)  :: U3(3), UIP(MNP)
      REAL(KIND=JWRU)  :: KKSUM, ST, N
      REAL(KIND=JWRU)  :: CX(3), CY(3)
      REAL(KIND=JWRU)  :: U(NFRE,MDC,MNP)
      REAL(KIND=JWRU)  :: FLALL(3)
      REAL(KIND=JWRU)  :: KELEM(3)
      REAL(KIND=JWRU)  :: FL111, FL112, FL211, FL212, FL311, FL312
      REAL(KIND=JWRU)  :: UTILDE3
      REAL(KIND=JWRU)  :: KSUM(MNP), KMAX(MNP)
      REAL(KIND=JWRU)  :: TIME1, TIME2
!
!        Calculate phase speeds for the certain spectral comp1.d0nt ...
!
      flall = 0.0_JWRU
      kelem = 0.0_JWRU
      kksum = 0.0_JWRU
      st = 0.0_JWRU
      n = 0.0_JWRU
 
      IF ( lcalc ) THEN
         DO id = 1 , mdc
            DO is = 1 , msc
               kmax = 0.0_JWRU
               ksum = 0.0_JWRU
               DO ip = 1 , mnp
                  DO i = 1 , CCON(ip)
                     ie = IE_CELL2(ip,i)
                     ipos = POS_CELL2(ip,i)
! get node indices from the element table ...
                     ni = ine(:,ie)
! estimate speed in WAE
                     cx = (CG(ni,is)*SINTH(id)+CURTXY(ni,1))              &
     &                    *1.0_JWRU/(DEGRAD*REARTH*COS(YP(IP)*DEGRAD))
                     cy = (CG(ni,is)*COSTH(id)+CURTXY(ni,2))              &
     &                    *1.0_JWRU/(DEGRAD*REARTH)

! upwind indicators
                     lambda(1) = 1.0_JWRU/6.0_JWRU*SUM(cx)
                     lambda(2) = 1.0_JWRU/6.0_JWRU*SUM(cy)
! flux jacobians
                     kelem(1) = MAX(0.0_JWRU,lambda(1)*IEN(1,ie)+lambda(2)    &
     &                          *IEN(2,ie))                                             ! K
                     kelem(2) = MAX(0.0_JWRU,lambda(1)*IEN(3,ie)+lambda(2)    &
     &                          *IEN(4,ie))
                     kelem(3) = MAX(0.0_JWRU,lambda(1)*IEN(5,ie)+lambda(2)    &
     &                          *IEN(6,ie))
 
                     ksum(ip) = ksum(ip) + kelem(ipos)
!2do check if also stable when abs removed
                     IF ( kelem(ipos).GT.kmax(ip) ) kmax(ip)              &
     &                    = kelem(ipos)
                  ENDDO
               ENDDO
               dtmax_global_exp = 10.E10
               dtmax_global_exp_loc = 10.E10
               DO ip = 1 , mnp 
                  dtmax_exp = SI(ip)/MAX(thr8,ksum(ip))
                  dtmax_global_exp_loc = MIN(dtmax_global_exp_loc,        &
     &               dtmax_exp)
               ENDDO
               CALL MPI_ALLREDUCE(dtmax_global_exp_loc,dtmax_global_exp,  &
     &                            1,MPI_REAL8,mpi_min,comm,ierr)
               cflxy = dt4a/dtmax_global_exp
               rest = ABS(MOD(cflxy,1.0_JWRU))
               IF ( rest.LT.thr8 ) THEN
                  ITER_EXP(is,id) = ABS(NINT(cflxy))
               ELSEIF ( rest.GT.thr8 .AND. rest.LT.0.5_JWRU ) THEN
                  ITER_EXP(is,id) = ABS(NINT(cflxy)) + 1
               ELSE
                  ITER_EXP(is,id) = ABS(NINT(cflxy))
               ENDIF
            ENDDO   ! IS
         ENDDO    ! ID
!         WRITE (dbg%FHNDL,*) 'MAX. ITERATIONS USED IN ADV. SCHEME' ,
!     &                        iter_max , MAXVAL(ITER_EXP)
           FLUSH(STAT%FHNDL)
      ENDIF     ! LCALC
 
      DO ip = 1 , mnp
         DO is = 1 , msc
            DO id = 1 , mdc
               u(is,id,ip) = fl1(ip,id,is)
            ENDDO
         ENDDO
      ENDDO
 
      CALL exchange(u)
 
      iter_max = MAXVAL(ITER_EXP)
      dt4ai = dt4a/iter_max
 
      DO it = 1 , iter_max
         DO id = 1 , mdc
            DO is = 1 , msc
               uip = u(is,id,:)
!!$OMP DO PRIVATE(IP,I,IE,IPOS)
               DO ip = 1 , mnp
                  st = 0.0_JWRU
                  DO i = 1 , CCON(ip)
! get element and the position of IP in the element index
                     ie = IE_CELL2(ip,i)
                     ipos = POS_CELL2(ip,i)
! get node indices from the element table ...
                     ni = ine(:,ie)
                     i1 = ni(1)
                     i2 = ni(2)
                     i3 = ni(3)
! estimate speed in WAE
                     cx = (CG(ni,is)*SINTH(id)+CURTXY(ni,1))             &
     &                    *1.0_JWRU/(DEGRAD*REARTH*COS(YP(IP)*DEGRAD))
                     cy = (CG(ni,is)*COSTH(id)+CURTXY(ni,2))             &
     &                    *1.0_JWRU/(DEGRAD*REARTH)

! flux integration using simpson rule ...
                     fl11 = cx(2)*IEN(1,ie) + cy(2)*IEN(2,ie)
                     fl12 = cx(3)*IEN(1,ie) + cy(3)*IEN(2,ie)
                     fl21 = cx(3)*IEN(3,ie) + cy(3)*IEN(4,ie)
                     fl22 = cx(1)*IEN(3,ie) + cy(1)*IEN(4,ie)
                     fl31 = cx(1)*IEN(5,ie) + cy(1)*IEN(6,ie)
                     fl32 = cx(2)*IEN(5,ie) + cy(2)*IEN(6,ie)
                     fl111 = 2.0_JWRU*fl11 + fl12
                     fl112 = 2.0_JWRU*fl12 + fl11
                     fl211 = 2.0_JWRU*fl21 + fl22
                     fl212 = 2.0_JWRU*fl22 + fl21
                     fl311 = 2.0_JWRU*fl31 + fl32
                     fl312 = 2.0_JWRU*fl32 + fl31
! upwind indicators
                     lambda(1) = 1.0_JWRU/6.0_JWRU*SUM(cx)
                     lambda(2) = 1.0_JWRU/6.0_JWRU*SUM(cy)
 
! flux jacobians
                     kelem(1) = lambda(1)*IEN(1,ie) + lambda(2)       &
     &                          *IEN(2,ie)                                   ! K
                     kelem(2) = lambda(1)*IEN(3,ie) + lambda(2)       &
     &                          *IEN(4,ie)
                     kelem(3) = lambda(1)*IEN(5,ie) + lambda(2)       &
     &                          *IEN(6,ie)
 
! inverse of the positive sum ...
                     n = -1.0_JWRU/MIN(-thr8,SUM(MIN(0.0_JWRU,kelem)))       ! N
 
! positive flux jacobians
                     kelem(1) = MAX(0.0_JWRU,kelem(1))
                     kelem(2) = MAX(0.0_JWRU,kelem(2))
                     kelem(3) = MAX(0.0_JWRU,kelem(3))
 
! simposon integration last step ...
                     flall(1) = (fl311+fl212)*1.0_JWRU/6.0_JWRU + kelem(1)
                     flall(2) = (fl111+fl312)*1.0_JWRU/6.0_JWRU + kelem(2)
                     flall(3) = (fl211+fl112)*1.0_JWRU/6.0_JWRU + kelem(3)
 
! flux conserving upwind contribution
                     u3 = uip(ni)
 
                     utilde3 = n*(flall(1)*u3(1)+flall(2)*u3(2)+flall(3)    &
     &                         *u3(3))
 
! coefficient for the integration in time
                     st = st + kelem(ipos)*(uip(ip)-utilde3)
 
                  ENDDO
! time stepping ...
                  uip(ip) = MAX(0.0_JWRU,uip(ip)-dt4ai/SI(ip)*st*IOBWB(ip))    &
     &                      *IOBPD(id,ip)
               ENDDO  !IP
               u(is,id,:) = uip
            ENDDO   !ID
         ENDDO    !IS
         CALL exchange(u)
      ENDDO     !IT
 
      DO ip = 1 , mnp
         DO is = 1 , msc
            DO id = 1 , mdc
               fl3(ip,id,is) = u(is,id,ip)
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE EXPLICIT_N_SCHEME_VECTOR_HPCF

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COHERENCY_ERROR_KERNEL(ACin, SumError, MaxError,       &
     &                          TotalSum)
      USE MPL_MPIF
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE yowDatapool, only: comm
      USE yowpd, only : np_global, MNP=>npa
      USE yownodepool, ONLY : NP, iplg

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      implicit none

      REAL(KIND=JWRB), intent(in) :: ACin(:)
      REAL(KIND=JWRB), intent(out) :: SumError, MaxError, TotalSum
      integer(KIND=JWIM) :: IP, iProc, IPglob, IS, ID
      integer(KIND=JWIM) :: MNPloc
      integer(KIND=JWIM) :: istat
      integer(KIND=JWIM) :: istatus(MPI_STATUS_SIZE)
      integer(KIND=JWIM) :: eStatus(np_global)
      REAL(KIND=JWRB)    :: ACtotal(np_global)
      REAL(KIND=JWRB)    :: NewVal
      REAL(KIND=JWRB)    :: eReal(3), TheDiff
      integer(KIND=JWIM), allocatable :: eInt(:)
      REAL(KIND=JWRB),    allocatable :: ACloc(:)
      integer(KIND=JWIM) ierr

      IF (IRANK == 1) THEN
        eStatus=0
        DO IP=1,MNP
          IPglob=iplg(IP)
          ACtotal(IPglob)=ACin(IP)
          eStatus(IPglob)=1
        END DO
        SumError=0
        MaxError=0
        DO iProc=2,nproc
          allocate(eInt(1), stat=istat)
          CALL MPI_RECV(eInt,1,MPI_INTEGER,iProc-1, 53, comm, istatus, ierr)
          MNPloc=eInt(1)
          deallocate(eInt)
          !
          allocate(eInt(MNPloc), stat=istat)
          allocate(ACloc(MNPloc), stat=istat)
          CALL MPI_RECV(eInt,MNPloc,MPI_INTEGER,iProc-1, 52, comm, istatus, ierr)
          CALL MPI_RECV(ACloc,MNPloc,MPI_REAL4,iProc-1, 51, comm, istatus, ierr)
          DO IP=1,MNPloc
            IPglob=eInt(IP)
            NewVal=ACloc(IP)
            IF (eStatus(IPglob) .eq. 1) THEN
              TheDiff=abs(NewVal - ACtotal(IPglob))
              SumError=SumError + TheDiff
              MaxError=MAX(MaxError, TheDiff)
            ELSE
              eStatus(IPglob)=1
              ACtotal(IPglob)=NewVal
            END IF
          END DO
          deallocate(eInt, stat=istat)
          deallocate(ACloc, stat=istat)
        END DO
        TotalSum=sum(ACtotal)
        !
        eReal(1)=SumError
        eReal(2)=MaxError
        eReal(3)=TotalSum
        DO iProc=2,nproc
          CALL MPI_SEND(eReal,3,MPI_REAL4, iProc-1, 50, comm, ierr)
        END DO
      ELSE
        allocate(eInt(1))
        eInt(1)=MNP
        CALL MPI_SEND(eInt,1,MPI_INTEGER, 0, 53, comm, ierr)
        deallocate(eInt)
        !
        CALL MPI_SEND(iplg,MNP,MPI_INTEGER, 0, 52, comm, ierr)
        !
        allocate(ACloc(MNP))
        DO IP=1,MNP
          ACloc(IP)=ACin(IP)
        END DO
        CALL MPI_SEND(ACloc,MNP,MPI_REAL4, 0, 51, comm, ierr)
        deallocate(ACloc)
        CALL MPI_RECV(eReal,3,MPI_REAL4,0, 50, comm, istatus, ierr)
        SumError=eReal(1)
        MaxError=eReal(2)
        TotalSum=eReal(3)
      END IF
      END SUBROUTINE COHERENCY_ERROR_KERNEL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COHERENCY_ERROR_KERNEL_THRESHOLD(ACin, Threshold)
      USE MPL_MPIF
      USE YOWMPP   , ONLY : IRANK, NPROC
      USE yowDatapool, only: comm
      USE yowpd, only : np_global, MNP=>npa
      USE yownodepool, ONLY : NP, iplg

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      implicit none

      REAL(KIND=JWRB), intent(in) :: ACin(:)
      REAL(KIND=JWRB), intent(in) :: Threshold
      integer(KIND=JWIM) IP, iProc, IPglob, IS, ID
      integer(KIND=JWIM) :: MNPloc
      integer(KIND=JWIM) :: istat
      integer(KIND=JWIM) :: istatus(MPI_STATUS_SIZE)
      integer(KIND=JWIM) :: eStatus(np_global)
      integer(KIND=JWIM) :: eOrigin(np_global)
      integer(KIND=JWIM) :: eOriginProc(np_global)
      REAL(KIND=JWRB) :: ACtotal(np_global)
      REAL(KIND=JWRB) :: NewVal, SumError
      REAL(KIND=JWRB) :: eReal(2), TheDiff
      integer(KIND=JWIM), allocatable :: eInt(:)
      REAL(KIND=JWRB),    allocatable :: ACloc(:)
      integer(KIND=JWIM) :: ierr

      IF (IRANK == 1) THEN
        eStatus=0
        DO IP=1,MNP
          IPglob=iplg(IP)
          ACtotal(IPglob)=ACin(IP)
          eStatus(IPglob)=1
          eOrigin(IPglob)=IP
          eOriginProc(IPglob)=1
        END DO
        DO iProc=2,nproc
          allocate(eInt(1), stat=istat)
          CALL MPI_RECV(eInt,1,MPI_INTEGER,iProc-1, 53, comm, istatus, ierr)
          MNPloc=eInt(1)
          deallocate(eInt)
          !
          allocate(eInt(MNPloc), stat=istat)
          allocate(ACloc(MNPloc), stat=istat)
          CALL MPI_RECV(eInt,MNPloc,MPI_INTEGER,iProc-1, 52, comm, istatus, ierr)
          CALL MPI_RECV(ACloc,MNPloc,MPI_REAL4,iProc-1, 51, comm, istatus, ierr)
          DO IP=1,MNPloc
            IPglob=eInt(IP)
            NewVal=ACloc(IP)
            IF (eStatus(IPglob) .eq. 1) THEN
              TheDiff=abs(NewVal - ACtotal(IPglob))
              IF (TheDiff .gt. Threshold) THEN
                WRITE(740+MyRankGlobal,*) 'Found divergence at'
                WRITE(740+MyRankGlobal,*) 'TheDiff/IPglob=', TheDiff, IPglob
                WRITE(740+MyRankGlobal,*) 'IP/eOrig=', IP, eOrigin(IPglob)
                WRITE(740+MyRankGlobal,*) 'iProc/iOrig=', iProc, eOriginProc(IPglob)
                FLUSH(740+MyRankGlobal)
              END IF
            ELSE
              ACtotal(IPglob)=NewVal
              eStatus(IPglob)=1
              eOrigin(IPglob)=IP
              eOriginProc(IPglob)=iProc
            END IF
          END DO
          deallocate(eInt, stat=istat)
          deallocate(ACloc, stat=istat)
        END DO
        !
      ELSE
        allocate(eInt(1))
        eInt(1)=MNP
        CALL MPI_SEND(eInt,1,MPI_INTEGER, 0, 53, comm, ierr)
        deallocate(eInt)
        !
        CALL MPI_SEND(iplg,MNP,MPI_INTEGER, 0, 52, comm, ierr)
        !
        allocate(ACloc(MNP))
        DO IP=1,MNP
          ACloc(IP)=ACin(IP)
        END DO
        CALL MPI_SEND(ACloc,MNP,MPI_REAL4, 0, 51, comm, ierr)
        deallocate(ACloc)
      END IF
      END SUBROUTINE COHERENCY_ERROR_KERNEL_THRESHOLD

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COHERENCY_ERROR_1D(ACin, string)

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      IMPLICIT NONE
      REAL(KIND=JWRB), intent(in) :: ACin(:)
      character(*), intent(in) :: string
      REAL(KIND=JWRB) :: SumError, MaxError, TotalSum
      CALL COHERENCY_ERROR_KERNEL(ACin, SumError, MaxError, TotalSum)
      WRITE(740+MyRankGlobal,*) 'COHERENCY ERROR ', TRIM(string)
      WRITE(740+MyRankGlobal,*) 'Error(Sum/Max)=', SumError, MaxError
      FLUSH(740+MyRankGlobal)
      END SUBROUTINE COHERENCY_ERROR_1D

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COHERENCY_ERROR_3D(ACin, string)

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP
      USE yowpd, only: MNP=>npa
      USE YOWUNPOOL, only : MSC, MDC

      IMPLICIT NONE

      INTEGER(KIND=JWIM) :: is, id, IP

      REAL(KIND=JWRB), intent(in) :: ACin(NINF-1:NSUP,NANG,NFRE)
      REAL(KIND=JWRB) :: SumErrorLoc, MaxErrorLoc, TotalSumLoc
      REAL(KIND=JWRB) :: SumError, MaxError, TotalSum
      REAL(KIND=JWRB) :: ACinred(MNP)

      character(*), intent(in) :: string

      SumError=0
      MaxError=0
      TotalSum=0
      DO IS = 1 , MSC
        DO ID = 1 , MDC
          DO IP=1,MNP
            ACinred(IP)=ACin(ip,id,is)
          END DO
          CALL COHERENCY_ERROR_KERNEL(ACinred, SumErrorLoc, MaxErrorLoc, TotalSumLoc)
          SumError = SumError + SumErrorLoc
          TotalSum = TotalSum + TotalSumLoc
          MaxError = MAX(MaxError, MaxErrorLoc)
        ENDDO
      ENDDO
      WRITE(740+MyRankGlobal,*) 'COHERENCY ERROR 3D ', TRIM(string)
      WRITE(740+MyRankGlobal,*) 'Error(Sum/Max),Sum=', SumError, MaxError, TotalSum
      FLUSH(740+MyRankGlobal)
      END SUBROUTINE COHERENCY_ERROR_3D

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECKS_FL1_FL3_SL

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP
      USE yowpd, only: MNP=>npa
      USE YOWSPEC, ONLY   : FL1, FL3, SL
      USE YOWUNPOOL ,ONLY : LLUNSTR
      IMPLICIT NONE
!      character(*), intent(in) :: string
!      WRITE(740+MyRankGlobal,*) 'CHECK ', TRIM(string)
      IF (LLUNSTR) THEN
        IF (ALLOCATED(FL1)) THEN
          CALL COHERENCY_ERROR_3D(FL1, "testing FL1")
        END IF
        IF (ALLOCATED(FL3)) THEN
          CALL COHERENCY_ERROR_3D(FL3, "testing FL3")
        END IF
        IF (ALLOCATED(SL)) THEN
          CALL COHERENCY_ERROR_3D(SL, "testing SL")
        END IF
      END IF
      WRITE(740+MyRankGlobal,*) 'allocated(FL3)=', allocated(FL3)
      WRITE(740+MyRankGlobal,*) 'allocated(SL)=', allocated(SL)
      IF (allocated(FL1)) THEN
        WRITE(740+MyRankGlobal,*) 'sum(FL1)=', sum(FL1)
      END IF
      IF (allocated(FL3)) THEN
        WRITE(740+MyRankGlobal,*) 'sum(FL3)=', sum(FL3)
      END IF
      IF (allocated(SL)) THEN
        WRITE(740+MyRankGlobal,*) 'sum(SL)=', sum(SL)
      END IF
      FLUSH(740+MyRankGlobal)
      END SUBROUTINE CHECKS_FL1_FL3_SL

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECKS_GROWTH

      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal

      USE YOWPARAM , ONLY : NANG, NFRE
      USE YOWMPP   , ONLY : NINF, NSUP
      USE yowpd, only: MNP=>npa
      USE YOWSPEC, ONLY   : FL1, FL3, SL
      USE YOWUNPOOL ,ONLY : LLUNSTR
      IMPLICIT NONE

      IF (allocated(FL1) .and. allocated(FL3)) THEN
        WRITE(740+MyRankGlobal,*) 'sum(FL1/FL3)=', sum(FL1), sum(FL3)
      END IF
      IF (.NOT.(allocated(FL1)) .and. allocated(FL3)) THEN
        WRITE(740+MyRankGlobal,*) 'sum(FL3)=', sum(FL3)
      END IF
      IF (allocated(FL1) .and. .NOT.(allocated(FL3))) THEN
        WRITE(740+MyRankGlobal,*) 'sum(FL1)=', sum(FL1)
      END IF
      FLUSH(740+MyRankGlobal)
      END SUBROUTINE CHECKS_GROWTH
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ALL_CHECKS
      IMPLICIT NONE
      CALL STAT_WHGTTG
      CALL CHECKS_FL1_FL3_SL
      END SUBROUTINE ALL_CHECKS

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PRINT_WHGTTG
      USE YOWINTP  , ONLY : WHGTTG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: I, J
      WRITE(740+MyRankGlobal,*) 'STAT_WHGTTG, alloc=', allocated(WHGTTG)
      IF (ALLOCATED(WHGTTG)) THEN
        DO I=1,NGX
          DO J=1,NGY
            WRITE(740+MyRankGlobal,*) 'I/J/Hs=', I, J, WHGTTG(I,J)
          END DO
        END DO
      END IF
      END SUBROUTINE PRINT_WHGTTG

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STAT_WHGTTG
      USE YOWINTP  , ONLY : WHGTTG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOW_RANK_GLOLOC, ONLY : MyRankGlobal
      IMPLICIT NONE
      INTEGER(KIND=JWIM) :: I, J, nbNAN, nbSea
      REAL(KIND=JWRB) :: sumPt, sumWHGTTG, avgWHGTTG
      REAL(KIND=JWRB) :: ZMISS, eVal
      ZMISS=-999.
      WRITE(740+MyRankGlobal,*) 'STAT_WHGTTG, alloc=', allocated(WHGTTG)
      IF (ALLOCATED(WHGTTG)) THEN
        nbNaN=0
        nbSea=0
        sumPt=0
        sumWHGTTG=0.
        DO I=1,NGX
          DO J=1,NGY
            eVal=WHGTTG(I,J)
            IF (eVal .ne. ZMISS) THEN
              nbSea=nbSea+1
              IF (eVal .ne. eVal) THEN
                nbNaN=nbNan+1
              ELSE
                sumPt=sumPt+1
                sumWHGTTG=sumWHGTTG + eVal
              END IF
            END IF
          END DO
        END DO
        avgWHGTTG=sumWHGTTG/sumPt
        WRITE(740+MyRankGlobal,*) 'nbSea/NaN/Pt/avg=', nbSea, nbNaN, sumPt, avgWHGTTG
      END IF
      END SUBROUTINE STAT_WHGTTG
!**********************************************************************
!*                                                                    *
!**********************************************************************
      END MODULE UNWAM 
