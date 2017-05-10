      SUBROUTINE OUTUNWAMTEST (F, LLAK)

! ----------------------------------------------------------------------

!**** *FEMEAN* - COMPUTATION OF MEAN FREQUENCY AT EACH GRID POINT
!                AND MEAN WAVE NUMBER .
!                THE COMPUTATION OF THE MEAN WAVE ENERGU WAS ALSO
!                ADDED SUCH THAT A CALL TO FEMEAN DOES NOT NEED
!                TO BE PRECEDED BY A CALL TO SEMEAN.

!     S.D. HASSELMANN
!     MODIFIED : P.JANSSEN (INTEGRATION OF F**-4 TAIL)
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER


!*    PURPOSE.
!     --------

!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.

!**   INTERFACE.
!     ----------

!       *CALL* *FEMEAN (F, IJS, IJL, EM, FM, LLAK)*
!              *F*   - SPECTRUM.
!              *IJS* - INDEX OF FIRST GRIDPOINT
!              *IJL* - INDEX OF LAST GRIDPOINT
!              *EM*  - MEAN WAVE ENERGY 
!              *FM*  - MEAN WAVE FREQUENCY 
!              *LLAK*- TRUE IF MEAN WAVE NUMBER IS COMPUTED

!     METHOD.
!     -------

!       NONE.

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------
      USE MPL_MPIF
      USE YOWFRED  , ONLY : FR       ,DFIM     ,DFIMOFR  ,DELTH    ,
     &                WETAIL    ,FRTAIL
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWPCONS , ONLY : G        ,ZPI      ,EPSMIN
      USE YOWSTAT  , ONLY : ISHALLO
      USE YOWSHAL  , ONLY : TFAK     ,INDEP
      USE YOWPD, ONLY     : MNP => npa, iplg, comm, myrank, np, 
     & np_global
      USE YOWUNPOOL, ONLY : FILEDEF, XFN_HS, XFN_TM, DBG
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IJ,M,K
      REAL :: DELT25, DELT2, DEL2
      REAL :: F(MNP,NANG,NFRE)
      REAL,DIMENSION(np) :: TEMP1, TEMP2, EM, FM
      REAL,DIMENSION(np_global) :: EM_GLOBAL, EM_GLO_red
      LOGICAL :: LLAK
      REAL, SAVE :: XFNOUTTIME
 
      DATA XFNOUTTIME /0./

      integer :: i, ierr

! ----------------------------------------------------------------------

!*    1. INITIALISE MEAN FREQUENCY ARRAY AND TAIL FACTOR.
!        ------------------------------------------------

      DO IJ=1,np
! THOMAS init EM with zero. otherwise MPISUM added EPSMIN nProc times
!         EM(IJ) = EPSMIN
        EM(IJ) = 0
        FM(IJ) = EPSMIN
      ENDDO

      DELT25 = WETAIL*FR(NFRE)*DELTH
      DELT2 = FRTAIL*DELTH

!*    2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

!*    2.1 DEEP WATER INTEGRATION.
!         -----------------------

        DO M=1,NFRE
          K=1
          DO IJ=1,np
            TEMP2(IJ) = F(IJ,K,M)
          ENDDO
          DO K=2,NANG
            DO IJ=1,np
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=1,np
            EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
            FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
          ENDDO
        ENDDO

      ELSE

!*    2.2 SHALLOW WATER INTEGRATION.
!         --------------------------

        DO M=1,NFRE
          K=1
          DO IJ=1,np
            TEMP1(IJ) = DFIM(M)/SQRT(TFAK(INDEP(IJ),M))
            TEMP2(IJ) = F(IJ,K,M) 
          ENDDO
          DO K=2,NANG
            DO IJ=1,np
              TEMP2(IJ) = TEMP2(IJ)+F(IJ,K,M)
            ENDDO
          ENDDO
          DO IJ=1,np
            EM(IJ) = EM(IJ)+TEMP2(IJ)*DFIM(M)
            FM(IJ) = FM(IJ)+DFIMOFR(M)*TEMP2(IJ)
          ENDDO
        ENDDO

      ENDIF

!*    3. ADD TAIL CORRECTION TO MEAN FREQUENCY AND
!*       NORMALIZE WITH TOTAL ENERGY.
!        ------------------------------------------

      IF (ISHALLO.EQ.1 .OR. .NOT.LLAK) THEN

        DO IJ=1,np
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
          FM(IJ) = EM(IJ)/FM(IJ)
        ENDDO

      ELSE

        DEL2 = DELT2*SQRT(G)/ZPI
        DO IJ=1,np
          EM(IJ) = EM(IJ)+DELT25*TEMP2(IJ)
          FM(IJ) = FM(IJ)+DELT2*TEMP2(IJ)
          FM(IJ) = EM(IJ)/FM(IJ)
        ENDDO

      ENDIF

      XFNOUTTIME = XFNOUTTIME + 1.

      WRITE(DBG%FHNDL,*) 'WRITING TEST OUTPUT TO BINARY FILE'
      WRITE(DBG%FHNDL,*) 'FIELD RANG'
      WRITE(DBG%FHNDL,*)  MNP
      WRITE(DBG%FHNDL,*) 'TIME STEPPING'
      WRITE(DBG%FHNDL,*) XFNOUTTIME

      WRITE(XFN_TM%FHNDL) XFNOUTTIME
      WRITE(XFN_TM%FHNDL) (0., 0., FM(IJ), IJ =1, np)

      ! parallel output of the Significant wave height

      EM = 4*SQRT(EM)

      EM_GLOBAL = 0
      do IJ=1, np
        EM_GLOBAL(iplg(IJ)) = EM(IJ)
      end do

      EM_GLO_red = 0

#ifdef C2A 
      call mpi_reduce(EM_GLOBAL, EM_GLO_red, np_global, MPI_REAL8,
     &                 MPI_SUM, 0, COMM, IERR)
#else
      call mpi_reduce(EM_GLOBAL, EM_GLO_red, np_global, MPI_REAL4,
     &                 MPI_SUM, 0, COMM, IERR)
#endif

      if(myrank == 0) then
        WRITE(XFN_HS%FHNDL) XFNOUTTIME
        WRITE(XFN_HS%FHNDL) (0., 0., EM_GLO_red(IJ),IJ = 1, np_global)
      endif

      RETURN
      END SUBROUTINE OUTUNWAMTEST 
