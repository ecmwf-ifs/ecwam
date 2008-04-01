      SUBROUTINE CURRENT2WAM (FILNM,IREAD,CDATEIN)

!--------------------------------------------------------------------

!**** *CURRENT2WAM* 

!     J. BIDLOT    ECMWF  NOVEMBER 2002.

!     PURPOSE.
!     --------

!     TO READ THE OCEAN CURRENT IN AND TO TRANSFER THEM ONTO THE
!     WAVE MODEL GRID.

!**   INTERFACE
!     ---------
!     *CALL* *CURRENT2WAM(FILNM,IREAD,CDATEIN)*

!     *FILNM*     DATA INPUT FILENAME.
!     *IREAD*     RANK OF THE PROCESS WHICH INPUTS THE DATA. 
!     *CDATEIN*    DATE OF THE DECODED DATA. 


!     METHOD.
!     -------


!    EXTERNAL.
!    ---------

!---------------------------------------------------------------------

      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP     ,KTAG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IREAD, IUNIT1
      INTEGER :: NBIT, KRET, ISIZE, IVAR, KPLENG, KLEN, KLENG, IDM
      INTEGER :: IFORP, IPARAM, KZLEV, IG, IJ, IC
      INTEGER, ALLOCATABLE :: KGRIB(:), IBUF(:)

      REAL :: ZDUM
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FIELD

      CHARACTER*14 :: CDATEIN, CDATEIN_OLD 
      CHARACTER*24 :: FILNM

      LOGICAL :: FRSTIME

      DATA FRSTIME / .TRUE. /
      SAVE IUNIT1
      SAVE FRSTIME

! --------------------------------------------------------------------  


      IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &    .NOT.LMESSPASS) THEN

        NBIT=175000

        IF (FRSTIME) THEN
          CALL PBOPEN(IUNIT1,FILNM,'r',KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBOPEN      *'
            WRITE (IU06,*) '*  IN CURRENT2WAM                      *'
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILNM
            IF(KRET.EQ.-2) WRITE (IU06,*)'INVALID FILENAME ',FILNM
            IF(KRET.EQ.-3) WRITE (IU06,*)'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
          ENDIF

          FRSTIME = .FALSE.
        ENDIF  

        ISIZE=NBIT
      ENDIF


      ALLOCATE(FIELD(NGX,NGY)) 

      READCURRENT: DO IVAR=1,2

        IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &      .NOT.LMESSPASS) THEN

2001      KPLENG=ISIZE*KIND(ISIZE)
          IF(.NOT.ALLOCATED(KGRIB)) ALLOCATE(KGRIB(ISIZE))

          CALL PBGRIB(IUNIT1,KGRIB,KPLENG,KLEN,KRET)

          IF (KRET.EQ.-1) THEN
            WRITE(IU06,*) ' *****************************************'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *    FATAL ERROR IN SUB. CURRENT2WAM    *'
            WRITE(IU06,*) ' *      ============================     *'
            WRITE(IU06,*) ' * PBGRIB EOF ON FILE ',FILNM 
            WRITE(IU06,*) ' *       NO MORE CURRENTS                *'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *****************************************'
            CALL ABORT1
          ELSE IF (KRET.EQ.-2) THEN
            WRITE(IU06,*) ' *****************************************'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *    FATAL ERROR IN SUB. CURRENT2WAM    *'
            WRITE(IU06,*) ' *      ============================     *'
            WRITE(IU06,*) ' * PBGRIB ERROR ON FILE ',FILNM
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *****************************************'
            CALL ABORT1
          ELSE IF (KRET.EQ.-3) THEN
            DEALLOCATE(KGRIB)
            CALL KGRIBSIZE(IU06, IUNIT1, KLEN, ISIZE, 'CURRENT2WAM')
            GOTO 2001
          ENDIF

          IF (ITEST.GT.1) THEN
            WRITE(IU06,*) ' SUB. CURRENT2WAM - READ FROM ',FILNM
            WRITE(IU06,*) ' IVAR = ',IVAR 
          ENDIF

! ----------------------------------------------------------------------

!*        UNPACK MARS FIELDS.
!         -------------------
          KLENG = KLEN 
          ZDUM=0.
          CALL GRB2WGRD (IU06, ITEST, NKSEK0,NKSEK1,NKSEK2,NKSEK3,
     &                   NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                   KGRIB, KLENG ,
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                   ZDUM, ZDUM,
     &                   CDATEIN, IFORP, IPARAM,KZLEV,IDM,IDM, FIELD)

          IF(IVAR.EQ.2 .AND. CDATEIN.NE.CDATEIN_OLD) THEN
            WRITE(IU06,*) ' *****************************************'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *    FATAL ERROR IN SUB. CURRENT2WAM    *'
            WRITE(IU06,*) ' *      ============================     *'
            WRITE(IU06,*) ' * DATE FOR PARAMETER ',IPARAM 
            WRITE(IU06,*) ' * IS ',CDATEIN
            WRITE(IU06,*) ' * IT SHOULD BE ',CDATEIN_OLD
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
            WRITE(IU06,*) ' *                                       *'
            WRITE(IU06,*) ' *****************************************'
            CALL FLUSH(IU06)
            CALL ABORT1
          ENDIF

          CDATEIN_OLD=CDATEIN

        ENDIF

        IF (LMESSPASS) THEN
!         BROADCAST FIELD TO OTHER PE'S

          CALL MPBCASTGRDFLD(IREAD,KTAG,NGX,NGY,FIELD)
          KTAG=KTAG+1


!         BROADCAST IPARAM AND CDATEIN TO OTHER PE'S

          ALLOCATE(IBUF(15))
          IF(IREAD.EQ.IRANK) THEN
            DO IC=1,14
              READ(CDATEIN(IC:IC),'(I1)') IBUF(IC)
            ENDDO
            IBUF(15)=IPARAM
          ENDIF

          CALL MPBCASTINTFLD(IREAD,KTAG,15,1,IBUF)
          KTAG=KTAG+1

          IF(IREAD.NE.IRANK) THEN
            DO IC=1,14
              WRITE(CDATEIN(IC:IC),'(I1)') IBUF(IC)
            ENDDO
            IPARAM=IBUF(15)
          ENDIF

          DEALLOCATE(IBUF)

        ENDIF

!       it does not work if multi block !!!
        IG=1
        IF(IPARAM.EQ.131) THEN
          CALL MAKEBLO(U(NSTART(IRANK),1),FIELD,
     &                 NSTART(IRANK),NEND(IRANK),IG)
!         some wam model grid points may have a missing data from
!         ocean model. They are set to 0.
          DO IJ = NSTART(IRANK),NEND(IRANK)
            IF(U(IJ,IG).LE.ZMISS) U(IJ,IG)=0.
          ENDDO
!         ECHANGE THE CURRENT FOR HALO POINTS SINCE THEY ARE NEEDED
!         FOR THE CALCULATION OF THE GRADIENTS (SEE GRADI).
          IF (LMESSPASS) THEN
            CALL MPEXCHNG(U,1,1)
          ENDIF
          U(NINF-1,IG)=0.
        ELSEIF(IPARAM.EQ.132) THEN
          CALL MAKEBLO(V(NSTART(IRANK),1),FIELD,
     &                 NSTART(IRANK),NEND(IRANK),IG)
          DO IJ = NSTART(IRANK),NEND(IRANK)
            IF(V(IJ,IG).LE.ZMISS) V(IJ,IG)=0.
          ENDDO
!         ECHANGE THE CURRENT FOR HALO POINTS SINCE THEY ARE NEEDED
!         FOR THE CALCULATION OF THE GRADIENTS (SEE GRADI).
          IF (LMESSPASS) THEN
            CALL MPEXCHNG(V,1,1)
          ENDIF
          V(NINF-1,IG)=0.
        ELSE
          WRITE(IU06,*) ' *****************************************'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. CURRENT2WAM    *'
          WRITE(IU06,*) ' *      ============================     *'
          WRITE(IU06,*) ' * THE INPUT PARAMETER ',IPARAM 
          WRITE(IU06,*) ' * IS NOT RECOGNISED AS BEING A CURRENT !*' 
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *****************************************'
          CALL ABORT1
        ENDIF

!!! need to put message paasing to collect the U and V in the halo.
!!! the index ij is only local for all halo points whereas it is global
!!! for other points


      ENDDO READCURRENT 

      DEALLOCATE(FIELD)

      RETURN
      END SUBROUTINE CURRENT2WAM
