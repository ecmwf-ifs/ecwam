      SUBROUTINE CURRENT2WAM (FILNM,IREAD,CDATEIN)

!--------------------------------------------------------------------

!**** *CURRENT2WAM* 

!     J. BIDLOT    ECMWF  NOVEMBER 2002.
!     J. BIDLOT    ECMWF  AUGUST 2008: USE *WAMCUR*

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
      USE YOWGRID  , ONLY : IGL      ,NLONRGG
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,  
     &            IFROMIJ  ,JFROMIJ
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : ZMISS
      USE YOWSPEC  , ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST

      USE MPL_MODULE

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: IREAD, IUNIT1
      INTEGER :: NBIT, KRET, ISIZE, IVAR, KPLENG, KLEN, KLENG, IDM
      INTEGER :: IFORP, IPARAM, KZLEV, IG, IJ, IC, IX, JY
      INTEGER :: IBUF(2)
      INTEGER, ALLOCATABLE :: KGRIB(:)

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

        NBIT=350000

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

        ENDIF


!       BROADCAST GRIB DATA TO OTHER PE'S
        IF(LMESSPASS .AND. NPROC.GT.1) THEN

          CALL MPL_BARRIER(CDSTRING='CURRENT2WAM: KGRIB ')

          CALL GSTATS(622,0)
          IF(IRANK.EQ.IREAD) THEN
            IBUF(1)=ISIZE
            IBUF(2)=KLEN
          ENDIF
          CALL MPL_BROADCAST(IBUF(1:2),KROOT=IREAD,KTAG=IVAR,
     &                       CDSTRING='CURRENT2WAM IBUF:')
          IF(IRANK.NE.IREAD) THEN
            ISIZE=IBUF(1)
            KLEN=IBUF(2)
            ALLOCATE(KGRIB(ISIZE))
          ENDIF

          CALL MPL_BROADCAST(KGRIB(1:ISIZE),KROOT=IREAD,KTAG=IVAR,
     &                       CDSTRING='CURRENT2WAM KGRIB:')

          CALL GSTATS(622,1)
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
               IF(U(IJ,IG).LE.ZMISS) U(IJ,IG)=0.
             ENDDO
             U(NINF-1,IG)=0.
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
               IF(V(IJ,IG).LE.ZMISS) V(IJ,IG)=0.
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

        DEALLOCATE(KGRIB)

      ENDDO READCURRENT 

      DEALLOCATE(FIELD)

      RETURN
      END SUBROUTINE CURRENT2WAM
