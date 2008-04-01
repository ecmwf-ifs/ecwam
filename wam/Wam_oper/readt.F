!------------------------------------------------------------------

      SUBROUTINE READT(IU06,FILNM,ITEST,FIELD,IPARAM)

!----------------------------------------------------------------------

!**** *READT*  READS TEMPERATURE OR SEA ICE FACTION FIELD.

!     R. PORTZ      MPI HAMBURG OCTOBER 1992
!     J. BIDLOT     ECMWF THE FIELD IS CONVERTED TO
!                         THE WAVE MODEL GRID. 

!*    PURPOSE.
!     --------

!       INPUT OF GLOBAL TEMPARATURE OR SEA ICE FACTION.

!**   INTERFACE.
!     ----------

!       *CALL* *READT* *(IU06,FILNM,ITEST,FIELD,IPARAM)*

!*     VARIABLE.   TYPE.     PURPOSE.
!      ---------   -------   --------
!      *IU06*      INTEGER   OUTPUT UNIT FOR STANDARD OUTPUT.
!      *FILNM*     CHARACTER*24 DATA INPUT FILENAME.
!      *ITEST*     TEST MESSAGE LEVEL.
!      *FIELD*     REAL      GLOBAL FIELD USED TO DETERMINE THE ICE MASK
!      *IPARAM*    INTEGER   GRIB PARAMETER OF FIELD

!     METHOD.
!     -------

!       FORMATTED READ FROM UNIT 1.

!     EXTERNALS.
!     ----------

!       INTMARS 
!       OUTPP

!     REFERENCE.
!     ----------

!       NONE.

!-------------------------------------------------------------------

!*    *PARAMETER*  FOR ARRAY DIMENSIONS.

      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : NLONRGG
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO

!-----------------------------------------------------------------------

      REAL,DIMENSION(NGX,NGY) :: FIELD
      REAL,ALLOCATABLE :: WORK(:,:)

      CHARACTER*14 CTDATE
      CHARACTER*24 FILNM
      CHARACTER TITL*100

!-----------------------------------------------------------------------

!*    1. INPUT OF TEMPERATURE DATA.
!     -----------------------------

      CALL GSTATS(1770,0)
      CALL INTMARS ( IU06, ITEST, FILNM, NKSEK0, NKSEK1, NKSEK2,
     &               NKSEK3, NKSEK4, NPSEK2, NPSEK3, KSEK3, PSEK3,
     &               NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &               AMOWEP, AMOSOP, AMOEAP, AMONOP,
     &               CTDATE, IFORP, IPARAM, KZLEV, FIELD )



      WRITE(IU06,*)' SUB. READT: ICE BOUNDARY DATE=',CTDATE
      CALL FLUSH(IU06)

      IF (ITEST.GT.4) THEN
        ALLOCATE(WORK(NGX,NGY))
        IF(IPARAM.EQ.139) THEN
          WRITE(IU06,*) 'TEMPERATURE GRID:'
          WORK=FIELD-273.15
          TITL= '   TEMPERATURES IN CELCIUS         '
        ELSEIF(IPARAM.EQ.31) THEN
          WRITE(IU06,*) 'SEA ICE FRACTION GRID:'
          TITL=  'SEA ICE FRACTION '
          WORK=FIELD*100
        ELSE
          WRITE(IU06,*) 'PARAMETER ',IPARAM
          TITL=  'UNKNOW FIELD '
        ENDIF
        WRITE(IU06,*) '-----------------'
        WRITE(IU06,*)
        WRITE(IU06,*) '  RESOLUTION LAT-LON  ',XDELLA, XDELLO
        WRITE(IU06,*) '   SOUTHERN LAT: ',AMOSOP,
     &                '   NORTHERN LAT: ',AMONOP
        WRITE(IU06,*) '   WESTERN LONG: ',AMOWEP,
     &                '   EASTERN LONG: ',AMOEAP
        WRITE(IU06,*)
        CALL OUTPP (CTDATE, IU06, NGX, NGX, NGX, NGY,
     &              TITL, 1.,
     &              WORK, AMOWEP, AMOSOP, AMOEAP, AMONOP)
        DEALLOCATE(WORK)
        CALL FLUSH(IU06)
      ENDIF
      CALL GSTATS(1770,1)

      RETURN
      END SUBROUTINE READT
