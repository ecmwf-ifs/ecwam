      SUBROUTINE CURRENT2WAM (FILNM,IREAD,CDATEIN)

!--------------------------------------------------------------------

!**** *CURRENT2WAM* 

!     J. BIDLOT    ECMWF  NOVEMBER 2002.
!     J. BIDLOT    ECMWF  AUGUST 2008: USE *WAMCUR*
!     J. BIDLOT    ECMWF  APRIL 2010: USE GRIBAPI 

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
      USE YOWGRIBHD, ONLY : PPEPS    ,PPREC
      USE YOWGRID  , ONLY : IGL      ,NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,  
     &            IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP     ,
     &            NPRECI
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS    ,EPSMIN
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE
      USE GRIB_API_INTERFACE

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: NBIT = 550000

      INTEGER :: IREAD
      INTEGER :: KFILE_HANDLE1
      INTEGER :: LFILE, KGRIB_HANDLE
      INTEGER :: KRET, IVAR, KPLENG, KLEN, KLENG, IDM, KDM, MDM
      INTEGER :: IRET, ISIZE
      INTEGER :: KK, MM
      INTEGER :: IFORP, IPARAM, KZLEV, IG, IJ, IC, IX, JY
      INTEGER :: IBUF(2)
      INTEGER, ALLOCATABLE :: INGRIB(:)
      INTEGER(KIND=JPKSIZE_T) :: KBYTES

      REAL, PARAMETER :: WLOWEST=0.0001
      REAL :: ZDUM
      REAL :: ZHOOK_HANDLE
      REAL :: FIELD(NGX,NGY)

      CHARACTER(LEN=14) :: CDATEIN, CDATEIN_OLD 
      CHARACTER(LEN=24) :: FILNM

      LOGICAL :: FRSTIME

      DATA FRSTIME / .TRUE. /
      SAVE KFILE_HANDLE1
      SAVE FRSTIME

! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('CURRENT2WAM',0,ZHOOK_HANDLE)


      LFILE = LEN_TRIM(FILNM)

      IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &    .NOT.LMESSPASS) THEN

        IF (FRSTIME) THEN
          KFILE_HANDLE1=-99
          CALL IGRIB_OPEN_FILE(KFILE_HANDLE1,FILNM(1:LFILE),'r')
          FRSTIME = .FALSE.
        ENDIF  

        ISIZE=NBIT
      ENDIF


      READCURRENT: DO IVAR=1,2

        IF ((LMESSPASS .AND. (IRANK.EQ.IREAD)) .OR.
     &      .NOT.LMESSPASS) THEN
1021        ISIZE=NBIT
            KBYTES=ISIZE*NPRECI
            IF(.NOT.ALLOCATED(INGRIB)) ALLOCATE(INGRIB(ISIZE))
            CALL IGRIB_READ_FROM_FILE(KFILE_HANDLE1,INGRIB,KBYTES,IRET)
            IF(IRET.EQ.JPGRIB_BUFFER_TOO_SMALL) THEN
!!!           *IGRIB_READ_FROM_FILE* does not read through the file if
!!!            the size is too small, so figure out the size and read again.
              CALL KGRIBSIZE(IU06, KBYTES, NBIT, 'CURRENT2WAM')
              DEALLOCATE(INGRIB)
              GOTO 1021
            ELSEIF(IRET.EQ.JPGRIB_END_OF_FILE) THEN
              WRITE(IU06,*) '**************************************'
              WRITE(IU06,*) '* CURRENT2WAM: END OF FILE ENCOUNTED *'
              WRITE(IU06,*) '**************************************'
              CALL ABORT1
            ELSEIF(IRET.NE.JPGRIB_SUCCESS) THEN
               WRITE(IU06,*) '************************************'
               WRITE(IU06,*) '* CURRENT2WAM: FILE HANDLING ERROR *'
               WRITE(IU06,*) '************************************'
              CALL ABORT1
            ENDIF

          WRITE(IU06,*) ' SUB. CURRENT2WAM - READ FROM ',FILNM
          WRITE(IU06,*) ' IVAR = ',IVAR 

        ENDIF

!       BROADCAST GRIB DATA TO OTHER PE'S
        CALL GSTATS(622,0)
        IF(LMESSPASS .AND. NPROC.GT.1) THEN

          CALL MPL_BARRIER(CDSTRING='CURRENT2WAM: INGRIB ')

          IF(IRANK.EQ.IREAD) THEN
            IBUF(1)=ISIZE
            IBUF(2)=KLEN
          ENDIF
          CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=IVAR,
     &                       CDSTRING='CURRENT2WAM IBUF:')
          IF(IRANK.NE.IREAD) THEN
            ISIZE=IBUF(1)
            KLEN=IBUF(2)
            ALLOCATE(INGRIB(ISIZE))
          ENDIF

          CALL MPL_BROADCAST(INGRIB(1:ISIZE),KROOT=IREAD,KTAG=IVAR,
     &                       CDSTRING='CURRENT2WAM INGRIB:')

        ENDIF
        CALL GSTATS(622,1)
! ----------------------------------------------------------------------

!*      UNPACK GRIB FIELDS.
!       -------------------
        CALL IGRIB_NEW_FROM_MESSAGE(KGRIB_HANDLE,INGRIB)

        KK=0
        MM=0
        CALL GRIB2WGRID (IU06, ITEST, NWAM_BLKS, 
     &                   KGRIB_HANDLE, INGRIB, ISIZE,
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &                   ZMISS, PPREC, PPEPS,
     &                   CDATEIN, IFORP, IPARAM,KZLEV,KK,MM, FIELD)

        CALL IGRIB_RELEASE(KGRIB_HANDLE)

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


        IF(IPARAM.EQ.131) THEN
           DO IG = 1,IGL
             DO IJ =  NINF, NSUP
!!!            VALUES FOR THE HALO ARE ALSO NEEDED !!!
!!!            HENCE NINF TO NSUP.
               IX = IFROMIJ(IJ,IG)
               JY = JFROMIJ(IJ,IG)
               U(IJ,IG) = FIELD(IX,JY)
!              SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!              OCEAN MODEL. THEY ARE SET TO 0.
!              0. WILL BE USED TO DETECT THE INABILITY TO COMPUTE THE GRADIANT
               IF(ABS(U(IJ,IG)).LE.WLOWEST) U(IJ,IG)=0.0
               IF(U(IJ,IG).LE.ZMISS) U(IJ,IG)=0.0
             ENDDO
             U(NINF-1,IG)=0.0
           ENDDO
        ELSEIF(IPARAM.EQ.132) THEN
           DO IG = 1,IGL
             DO IJ =  NINF, NSUP
!!!            VALUES FOR THE HALO ARE ALSO NEEDED !!!
!!!            HENCE NINF TO NSUP.
               IX = IFROMIJ(IJ,IG)
               JY = JFROMIJ(IJ,IG)
               V(IJ,IG) = FIELD(IX,JY)
!              SOME WAM MODEL GRID POINTS MAY HAVE A MISSING DATA FROM
!              OCEAN MODEL. THEY ARE SET TO 0.
!              0. WILL BE USED TO DETECT THE INABILITY TO COMPUTE THE GRADIANT
               IF(ABS(V(IJ,IG)).LE.WLOWEST) V(IJ,IG)=0.0
               IF(V(IJ,IG).LE.ZMISS) V(IJ,IG)=0.0
             ENDDO
             V(NINF-1,IG)=0.
           ENDDO
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

        DEALLOCATE(INGRIB)

      ENDDO READCURRENT 

      IF (LHOOK) CALL DR_HOOK('CURRENT2WAM',1,ZHOOK_HANDLE)

      RETURN
      END SUBROUTINE CURRENT2WAM
