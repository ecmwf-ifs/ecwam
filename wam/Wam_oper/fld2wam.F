      SUBROUTINE FLD2WAM (IPARAM, FILNM, IJS, IJL, IG, VALMISS,
     &                    BLOCK, CDATEIN)

!--------------------------------------------------------------------

!**** *FLD2WAM* 

!     PURPOSE.
!     --------

!     TO READ !ONE! GRIB FIELD IN AND TO TRANSFER IT ONTO THE
!     WAVE MODEL LOCAL SEA POINTS.

!**   INTERFACE
!     ---------
!     *CALL* *FLD2WAM(IPARAM, FILNM, IJS, IJL, IG, VALMISS,
!    &                    BLOCK, CDATEIN)


!     *IPARAM*    GRIB PARAMETER EXPETED.
!     *FILNM*     DATA INPUT FILENAME.
!            !!!  FILNM WILL NOT BE KEPT OPENED !!!!
!     *IJS*       INDEX OF FIRST GRIDPOINT
!     *IJL*       INDEX OF LAST GRIDPOINT
!     *IG*        BLOCK NUMBER
!     *VALMISS*   VALUE GIVEN TO BLOCK IF MISSING DATA ARE DECODED
!     *BLOCK*     LOCAL MODEL SEA POINT VALUES OF INPUT FIELD
!     *CDATEIN*   DATE OF THE DECODED DATA. 


!     METHOD.
!     -------


!    EXTERNAL.
!    ---------

!---------------------------------------------------------------------

      USE YOWCOUT  , ONLY : LGRIBAPI_I
      USE YOWCURR  , ONLY : U        ,V
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3 ,
     &            PPEPS    ,PPREC
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,  
     &            IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECI
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE
      USE GRIB_API_INTERFACE

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NBIT = 550000

      INTEGER :: IPARAM, IREAD, IUNIT1, IJS, IJL, IG
      INTEGER :: KFILE_HANDLE1
      INTEGER :: LFILE, KGRIB_HANDLE
      INTEGER :: KRET, KPLENG, KLEN, KLENG, IDM
      INTEGER :: IRET, KBYTES, ISIZE
      INTEGER :: KK, MM
      INTEGER :: IFORP, IPARAM_IN, KZLEV, IJ, IC, IX, JY
      INTEGER :: IBUF(2)
      INTEGER, ALLOCATABLE :: INGRIB(:)

      REAL :: VALMISS 
      REAL :: ZDUM
      REAL :: ZHOOK_HANDLE
      REAL :: BLOCK(IJS:IJL)
      REAL :: FIELD(NGX,NGY)

      CHARACTER(LEN=14) :: CDATEIN
      CHARACTER(LEN=24) :: FILNM

! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('FLD2WAM',0,ZHOOK_HANDLE)

      IREAD=1

      LFILE = LEN_TRIM(FILNM)

      IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &    .NOT.LMESSPASS) THEN

        IF(LGRIBAPI_I) THEN
          KFILE_HANDLE1=-99
          CALL IGRIB_OPEN_FILE(KFILE_HANDLE1,FILNM(1:LFILE),'r')
        ELSE
          CALL PBOPEN(IUNIT1,FILNM,'r',KRET)
          IF(KRET.LT.0) THEN
            WRITE (IU06,*) '****************************************'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '*  ERROR FOLLOWING CALL TO PBOPEN      *'
            WRITE (IU06,*) '*  IN FLD2WAM                          *'
            IF(KRET.EQ.-1)
     &       WRITE (IU06,*) ' COULD NOT OPEN FILE ',FILNM
            IF(KRET.EQ.-2) WRITE (IU06,*)'INVALID FILENAME ',FILNM
            IF(KRET.EQ.-3) WRITE (IU06,*)'INVALID OPEN MODE SPECIFIED'
            WRITE (IU06,*) '*                                      *'
            WRITE (IU06,*) '****************************************'
            CALL ABORT1
          ENDIF
        ENDIF

        ISIZE=NBIT
      ENDIF


        IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &      .NOT.LMESSPASS) THEN

          IF(LGRIBAPI_I) THEN
!           USING GRIBAPI
1021        ISIZE=NBIT
            KBYTES=ISIZE*NPRECI
            IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
            CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE1,INGRIB,KBYTES,IRET)
            IF(IRET.EQ.JPGRIB_BUFR_TOO_SHORT) THEN
!!!           *IGRIB_READ_FROM_FILE* does not read through the file if
!!!            the size is too small, so figure out the size and read again.
              IUNIT1=-1
              CALL KGRIBSIZE(IU06, IUNIT1, KBYTES, NBIT, 'FLD2WAM')
              DEALLOCATE(INGRIB)
              GOTO 1021
            ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
              WRITE(IU06,*) '**************************************'
              WRITE(IU06,*) '* FLD2WAM: END OF FILE ENCOUNTED *'
              WRITE(IU06,*) '**************************************'
              CALL ABORT1
            ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
               WRITE(IU06,*) '************************************'
               WRITE(IU06,*) '* FLD2WAM: FILE HANDLING ERROR *'
               WRITE(IU06,*) '************************************'
              CALL ABORT1
            ENDIF

          ELSE
!           USING GRIBEX
2001        KPLENG=ISIZE*KIND(ISIZE)
            IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
            CALL GSTATS(1740,0)
            CALL PBGRIB(IUNIT1,INGRIB,KPLENG,KLEN,KRET)
            CALL GSTATS(1740,1)

            IF (KRET.EQ.-1) THEN
              WRITE(IU06,*) ' ****************************************'
              WRITE(IU06,*) ' *                                      *'
              WRITE(IU06,*) ' *    FATAL ERROR IN SUB. FLD2WAM   *'
              WRITE(IU06,*) ' *      ============================    *'
              WRITE(IU06,*) ' * PBGRIB EOF ON FILE ',FILNM 
              WRITE(IU06,*) ' *       NO MORE CURRENTS               *'
              WRITE(IU06,*) ' *                                      *'
              WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS   *'
              WRITE(IU06,*) ' *                                      *'
              WRITE(IU06,*) ' ****************************************'
              CALL ABORT1
            ELSE IF (KRET.EQ.-2) THEN
              WRITE(IU06,*) ' ****************************************'
              WRITE(IU06,*) ' *                                      *'
              WRITE(IU06,*) ' *    FATAL ERROR IN SUB. FLD2WAM   *'
              WRITE(IU06,*) ' *      ============================    *'
              WRITE(IU06,*) ' * PBGRIB ERROR ON FILE ',FILNM
              WRITE(IU06,*) ' *                                      *'
              WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS   *'
              WRITE(IU06,*) ' *                                      *'
              WRITE(IU06,*) ' ****************************************'
              CALL ABORT1
            ELSE IF (KRET.EQ.-3) THEN
              DEALLOCATE(INGRIB)
              CALL KGRIBSIZE(IU06, IUNIT1, KLEN, ISIZE, 'FLD2WAM')
              GOTO 2001
            ENDIF
          ENDIF

          WRITE(IU06,*) ' SUB. FLD2WAM - READ FROM ',FILNM
        ENDIF


!       BROADCAST GRIB DATA TO OTHER PE'S
        CALL GSTATS(622,0)
        IF(LMESSPASS .AND. NPROC.GT.1) THEN

          CALL MPL_BARRIER(CDSTRING='FLD2WAM: INGRIB ')

          IF(IRANK.EQ.IREAD) THEN
            IBUF(1)=ISIZE
            IBUF(2)=KLEN
          ENDIF
          CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=1,
     &                       CDSTRING='FLD2WAM IBUF:')
          IF(IRANK.NE.IREAD) THEN
            ISIZE=IBUF(1)
            KLEN=IBUF(2)
            ALLOCATE(INGRIB(ISIZE))
          ENDIF

          CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=2,
     &                       CDSTRING='FLD2WAM INGRIB:')

        ENDIF
        CALL GSTATS(622,1)
! ----------------------------------------------------------------------

!*        UNPACK GRIB FIELDS.
!         -------------------
          IF(LGRIBAPI_I) THEN
           CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

           CALL GRIB2WGRID (IU06, ITEST, NWAM_BLKS, 
     &                      KGRIB_HANDLE, INGRIB, ISIZE,
     &                      NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                      AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                      ZMISS, PPREC, PPEPS,
     &                      CDATEIN, IFORP, IPARAM_IN,KZLEV,KK,MM,FIELD)

            CALL IGRIB_RELEASE(KGRIB_HANDLE)
          ELSE
            KLENG = KLEN 
            ZDUM=0.
            CALL GRB2WGRD (IU06, ITEST, NKSEK0,NKSEK1,NKSEK2,NKSEK3,
     &                     NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &                     INGRIB, KLENG ,
     &                     NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                     AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                     ZDUM, ZDUM,
     &                     CDATEIN, IFORP,IPARAM_IN,KZLEV,IDM,IDM,FIELD)
          ENDIF


        IF(IPARAM_IN.EQ.IPARAM) THEN
           DO IJ = IJS, IJL 
             IX = IFROMIJ(IJ,IG)
             JY = JFROMIJ(IJ,IG)
             BLOCK(IJ) = FIELD(IX,JY)
!            SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!            OCEAN MODEL. THEY ARE SET TO VALMISS 
             IF(BLOCK(IJ).LE.ZMISS) BLOCK(IJ)=VALMISS
           ENDDO
        ELSE
          WRITE(IU06,*) ' *****************************************'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *    FATAL ERROR IN SUB. FLD2WAM        *'
          WRITE(IU06,*) ' *      ============================     *'
          WRITE(IU06,*) ' * THE INPUT PARAMETER ',IPARAM_IN 
          WRITE(IU06,*) ' * IS NOT THE REQUESTED ONE: ',IPARAM 
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *  PROGRAM ABORTS     PROGRAM ABORTS    *'
          WRITE(IU06,*) ' *                                       *'
          WRITE(IU06,*) ' *****************************************'
          CALL ABORT1
        ENDIF

        DEALLOCATE(INGRIB)


      IF (LHOOK) CALL DR_HOOK('FLD2WAM',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE FLD2WAM
