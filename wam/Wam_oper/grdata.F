      SUBROUTINE GRDATA (WHME, USME, NSTART, NEND)
C
C--------------------------------------------------------------------
C
C**** *GRDATA* - TRANSFERS MEASUREMENTS TO MODEL GLOBAL GRID,
C****            DISCARDING UNRELIABLE DATA.
C
C     P.LIONELLO     ECMWF       APRIL 1990
C
C     PURPOSE.
C     --------
C
C       PREPARE THE DATA FOR OPTIMAL INTERPOLATION.
C
C**   INTERFACE.
C     ----------
C
C       *CALL* *GRDATA (WHME, USME)*
C
C        *WHME*   REAL    OUTPUT FIELD OF MEASURED WAVE HEIGHTS.
C        *USME*   REAL    OUTPUT FIELD OF MEASURED WIND SPEEDS.
C        *NSTART* INTEGER INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C        *NEND*   INTEGER INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C
C
C     METHOD.
C     -------
C
C        A BOX OF SIZE DELLA*DELLO IS CENTERED IN EACH GRID POINT.
C        THE AVERAGE VALUE OF THE SATELLITE MEASUREMENTS INSIDE THE
C        BOX IS TAKEN TO UPDATE THE WAM MODEL. THE VALUE IS DISCARDED
C        IF THE ROOT MEAN SQUARE ERROR IS TOO BIG.
C        MOREOVER THE VALUE IS DISCARDED IF THERE ARE LESS THEN FOUR
C        MEASUREMENTS INSIDE THE BOX.
C
C      EXTERNALS.
C      ----------
C
C         READSAT  - READ A MEASUREMENT.
C
C      MODIFICATONS.
C      -------------
C         FEB 92  Bjoern Hansen :
C                 Interface to readsat extended in order to transfer
C                 the flags added to the satellite information by the
C                 quality control preura.
C                 The as unreliable flagged satellite measurements
C                 are rejected.
C                 The input unit of the measurements is automatically
C                 looked up each time grdata is called. It is assumed
C                 that the needed input file is present in the working
C                 directory, contains data of the same period of time
C                 as the output time step and is named RFLyyyymmddhhmm,
C                 where yyyymmddhhmm is a DATE/TIME group.
C
C          Aug 97 Jean Bidlot : in coupled mode IDELWO is not the
C                 appropriate time span for the selection of the data
C                 to be used but rather IDELALT, A NEW VARIABLE WHICH
C                 IS SUPPLIED IN WAMINPUT
C
C          Oct 97 Jean Bidlot : use of the new ice mask instead of 
C                 SST field to determine whether an observation which falls
C                 on a model sea point.
C
C-----------------------------------------------------------------------
C
#include "param.h"
C
#include "comice.h"
C
#include "comintp.h"
C
#include "commap.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comgrid.h"
C
#include "comstat.h"
C
#include "comtest.h"
C
C ----------------------------------------------------------------------
C
      LOGICAL EOFD
      CHARACTER*12 CBEGINDT, CENDDT
      CHARACTER*14 CDATE,CDATEO
      DIMENSION WHME(NGX,NGY), USME(NGX,NGY),
     1          WHSE(NGX,NGY), USSE(NGX,NGY),
     2          NUMBWH(NGX,NGY), NUMBUS(NGX,NGY),NR(NGX,NGY)
C
C*     VARIABLE     TYPE     PURPOSE.
C      --------     ----     --------
C
C     *CBEGINDT*    CHAR*12  BEGINNING DATE OF THE TIME WINDOW
C                            DURING WHICH DATA ARE USED FOR ANALYSIS.
C     *CENDDT*      CHAR*12  END DATE OF THE TIME WINDOW
C                            DURING WHICH DATA ARE SAVE FOR ANALYSIS.
C     *WHME*        REAL     GRID WAVE HEIGHT.
C     *USME*        REAL     GRID WIND SPEED.
C     *WHSE*        REAL     WAVE HEIGHT ROOT MEAN SQUARE ERROR.
C     *USSE*        REAL     WIND SPEED ROOT MEAN SQUARE ERROR.
C     *NUMBWH*      INTEGER  NUMBER OF HS MEASUREMENTS INSIDE A BOX.
C     *NUMBUS*      INTEGER  NUMBER OF US MEASUREMENTS INSIDE A BOX.
C
C ----------------------------------------------------------------------
C
      INTEGER,DIMENSION(NPROC) :: NSTART, NEND
C ----------------------------------------------------------------------
      INTEGER,ALLOCATABLE :: ICEMASKG(:,:)

      PARAMETER (jpice=1000)

      REAL,ALLOCATABLE :: RLAICE(:),RLOICE(:),SWHICE(:),WSICE(:)
      CHARACTER*12,ALLOCATABLE :: CDAICE(:)
      LOGICAL flagged, opened
C
C ----------------------------------------------------------------------
C
C
C*    1. FIX TIME WINDOW EQUAL TO ONE WIND TIME STEP.
C        --------------------------------------------
C
      CBEGINDT = CDTPRO
      CALL INCDATE (CBEGINDT,-IDELALT/2)
      CENDDT = CBEGINDT
      CALL INCDATE (CENDDT,IDELALT)

      WRITE(IU06,*) ' GRDATA : SELECTING ALTIMETER MEASUREMENTS '
      WRITE(IU06,*) '          PERIOD FROM ', CBEGINDT,' TO ', CENDDT
C
C     FIX THE MINIMUM NUMBER OF OBSERVATIONS PER GRID BOX REQUIRED
C     IN ORDER TO ACCEPT THE DATA AT THE GRIB POINTS
C
      IF(XDELLA.GE.0.5) THEN
        NNMIN=4
      ELSE
        NNMIN=2
      ENDIF
      WRITE(IU06,*) '   GRIDDED DATA WILL BE KEPT WHEN THE NUMBER OF'
      WRITE(IU06,*) '   OBSERVATIONS PER GRID BOX IS LARGER THAN ',
     &  NNMIN-1
C
C ----------------------------------------------------------------------
C
C*    2. INITIALIZE ARRAYS, counters and units.
C        --------------------------------------
C
      ALLOCATE(RLAICE(JPICE),RLOICE(JPICE),SWHICE(JPICE),WSICE(JPICE))
      ALLOCATE(CDAICE(JPICE))
      ICOUNT=0
      ICOUNTT=0
      ICOUNTD=0
      icountf=0
      icountb=0
      icountice=0
      icountqcf=0
      DO 2001 J=1,NGY
         DO 2002 I=1,NGX
            WHME(I,J) = 0.
            USME(I,J) = 0.
            WHSE(I,J) = 0.
            USSE(I,J) = 0.
            NUMBWH(I,J) = 0
            NUMBUS(I,J) = 0
            NR(I,J)=0
 2002    CONTINUE
 2001 CONTINUE
C
C     INITIALISE THE ICE MASK ARRAY
C
      ALLOCATE(ICEMASKG(NGX,NGY))

      ICEMASKG=1

      IF(LMESSPASS) THEN
         DO IJ = NSTART(1), NEND(NPROC)
            IF (ICEMASK(IJ,1).EQ.0) THEN
               IX = IXLG(IJ,1)
               IY = NGY- KXLT(IJ,1) +1
               ICEMASKG(IX,IY)=0
            ENDIF
         ENDDO
      ELSE
         DO IG = 1, IGL 
            DO IJ = IJS(IG),IJL(IG) 
               IF (ICEMASK(IJ,IG).EQ.0) THEN
                  IX = IXLG(IJ,IG)
                  IY = NGY- KXLT(IJ,IG) +1
                  ICEMASKG(IX,IY)=0
               ENDIF
            ENDDO
         ENDDO
      ENDIF
C
      IOLD = 0
      JOLD = 0
      CDATEO = CBEGINDT(1:12)//'00'

      DO 2003 iume=33,99,1
        INQUIRE ( UNIT=iume, OPENED=opened)
        IF ( .NOT. opened) THEN
           WRITE(iu06,*)
     .     ' GRDATA : TAKES UNIT ',iume,' AS ALTIMETER DATA INPUT UNIT'
           GOTO 3000
        ENDIF
 2003 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    3. READ THE MEASUREMENT TAPE AND CUMULATE THE DATA ON GRID POINTS.
C        ---------------------------------------------------------------
C
 3000 CONTINUE
      flagged = .TRUE.
      CALL READSAT (IU06, iume, CDTPRO,
     .              CDATE, RLAT, RLON, SWH, WS, flagged, EOFD)
      IF (EOFD) GOTO 4000
      ICOUNTD = ICOUNTD+1

C     CHECK DATE.

      IF (CDATE(1:12).LT.CBEGINDT)  then
         icountb = icountb + 1
         GOTO 3000
      endif

      IF (CDATE(1:12).LE.CENDDT) THEN

C     CHECK QUALITY.

         IF (flagged) THEN
           icountf=icountf+1
           GOTO 3000
         ENDIF

C     DATE IS INSIDE THE TIME WINDOW and record flagged as reliable.
C     THE INDECES ARE COMPUTED

         ICOUNTT = ICOUNTT+1
         J       = NGY-NINT((RLAT-AMOSOP)/XDELLA)
         XI      = MOD(RLON-AMOWEP+720.,360.)

         JSN     = NGY-J+1
         IF((J.LE.0).OR.(J.GT.NGY)) JSN=1
         NRGG    = NLONRGG(JSN)
         I       = NINT(XI/ZDELLO(JSN)) + 1
         IF (IPER.EQ.1.AND.I.EQ.NRGG+1) I = 1
C
C     IF THE INDICES ARE INSIDE THE GRID THE MEASUREMENT IS CUMULATED
C
         IF((I.LE.NRGG).AND.(I.GE.1) .AND. (J.LE.NGY).AND.(J.GE.1)) THEN
C
C     IF THE GRID POINT IS SET TO BE ICE TAKE NEXT MEASUREMENT.
C
            IF (ICEMASKG(I,J).EQ.0) THEN
              icountice = icountice + 1
              if (icountice .le. jpice) then
                cdaice(icountice)=cdate(1:12)
                rlaice(icountice)=rlat
                rloice(icountice)=rlon
                swhice(icountice)=swh
                wsice(icountice)=ws
              endif
              GOTO 3000
            ENDIF
C
C     SOME MEASURENTS ARE VERY CLOSE TO LAND, PASSED THE QC AS
C     BEING ON SEA BUT REGARDED AS BEING ON LAND BY GRDATA.
C     THESE ARE ELIMINATED QUIETLY HERE.
C     NOTE: CDG CAN BE = 0.0 FOR SEA POINTS WHEN GRDATA IS USED BY *RFL4WAM*
C
            IF (CDG(I,J).LT.0.000) THEN
              icountQCF = icountQCF + 1
              GOTO 3000
            ENDIF


            ICOUNT=ICOUNT+1
            IF (I.NE.IOLD .OR. J.NE.JOLD) THEN
C
C               NEW BOX
C
                CALL DIFDATE(CDATEO,CDATE,ISHIFT)
                IF (ISHIFT.LT.3.AND.ICOUNT.GT.1) THEN
                   IF (SWH-SWHO.gt.2.0) THEN
                      NR(I,J)=NR(I,J)+1
                      GO TO 3000
                   ENDIF
                END IF
                IOLD = I
                JOLD = J
                SWHO = SWH
                WSO  = WS
            ELSE
C
C     SUDDEN POSITVE STEPS ARE REFUSED.
C
               IF (SWH-SWHO.GT.1.0) THEN
                  NR(I,J) = NR(I,J)+1
                  GO TO 3000
               ENDIF
            ENDIF
C
            IF (SWH.GE.0.) THEN
               WHME(I,J) = WHME(I,J) + SWH
               WHSE (I,J) = WHSE(I,J) + SWH**2
               NUMBWH(I,J) = NUMBWH(I,J) + 1
               IF (WS.GE.0.) THEN
                   USME(I,J) = USME(I,J) + WS
                   USSE (I,J) = USSE(I,J) + WS**2
                   NUMBUS(I,J) = NUMBUS(I,J) + 1
               ENDIF
            ENDIF
C
C     A PREVIOUS SUDDEN NEGATIVE STEP IS COMPENSATED.
C
            IF (SWHO-SWH.GT.2.0) then
               NR(I,J) = NR(I,J)+1
               WHME(I,J) = WHME(I,J) - SWHO
               WHSE (I,J) = WHSE(I,J) - SWHO**2
               NUMBWH(I,J) = NUMBWH(I,J) - 1
               IF (WSO.GE.0.) THEN
                  USME(I,J) = USME(I,J) - WSO
                  USSE (I,J) = USSE(I,J) - WSO**2
                  NUMBUS(I,J) = NUMBUS(I,J) - 1
               ENDIF
            ENDIF
C
            SWHO = SWH
            WSO = WS
            CDATEO = CDATE
         ENDIF
C*

         GOTO 3000
      ENDIF

C
C ----------------------------------------------------------------------
C
C*    4. DETERMINE AVERAGE VALUE AND ERROR OF THE MEASUREMENTS IN THE
C*       BOX AROUND GRIDPOINT I,J.
C        ------------------------------------------------------------
C
 4000 CONTINUE
      nuwh=0
      nuus=0
      DO j2=1,ngy
        j2sn = ngy-j2+1     
        DO j1=1,nlonrgg(j2sn)
            nuwh=nuwh+numbwh(j1,j2)
            nuus=nuus+numbus(j1,j2)
        ENDDO
      ENDDO

      WRITE (IU06,'(" GRDATA: INPUT STATISTICS:")')
      WRITE (IU06,'("         NUMBER OF RECORDS READ FROM DATA FILE IS",
     1              "......:", I8)') ICOUNTD
      WRITE (IU06,'("         THEREOF NUMBER OF RECORDS WITH FLAGS IS ",
     1              "......:", I8)') ICOUNTF
      WRITE (IU06,'("         NUMBER OF DATA IN SELECTED PERIOD IS ...",
     1              "......:", I8)') ICOUNTT
      WRITE (IU06,'("         NUMBER OF DATA BEFORE SELECTED PERIOD IS",
     1              "......:", I8)') ICOUNTB
      WRITE (IU06,'("         NUMBER OF DATA IN BOXES AROUND SEA POINT",
     1              "S IS..:", I8)') ICOUNT
      WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE LAND IS.......",
     1              "......:", I8)') ICOUNTQCF
      WRITE (IU06,'("         NUMBER OF GOOD DATA ABOVE ICE IS........",
     1              "......:", I8)') ICOUNTICE
      WRITE (iu06,'("         NUMBER OF DATA ACCEPTED BY ''GRDATA'' IS",
     1              " FOR wh :", I8, /, 
     2"                                            AND FOR us :", I8)')
     3    nuwh, nuus

      DO 4150 iubl=30,99,1
        INQUIRE ( UNIT=iubl, OPENED=opened)
        IF ( .NOT. opened.and.ITEST.GT.0.and.ITEST.GT.0) THEN
           WRITE(iu06,*)
     .     ' GRDATA: TAKES UNIT ',iubl,' AS BLACKLIST OUTPUT UNIT'
           GOTO 4160
        ENDIF
 4150 CONTINUE
 4160 CONTINUE
      OPEN (IUBL,FILE='BLL'//CDTPRO)
      WRITE(IUBL,'(A12)') CDTPRO
      DO 4200 j=1,MIN(icountice,jpice)
        WRITE(IUBL,'(a12,4f8.2)') cdaice(j),rlaice(j),rloice(j),
     .                            swhice(j),wsice(j)
 4200 continue
C
      DO 4400 J=1,NGY
         JSN = NGY-J+1
         DO 4300 I=1,NLONRGG(JSN)
            NN = NUMBWH(I,J)
            IF (NN.GE.NNMIN) THEN
C
C           MEAN VALUES
C
               WHME(I,J) = WHME(I,J)/FLOAT(NN)
C
C           COMPUTE THE STANDARD DEVIATION.
C
               WHSE2 =  (WHSE(I,J)-WHME(I,J)**2*FLOAT(NN))
     1                              / (FLOAT(NN)-1.)
               WHSE(I,J) = SQRT( MAX(WHSE2,0.))
            ENDIF
 4300    CONTINUE
 4400 CONTINUE
C
      DO 4600 J=1,NGY
         JSN = NGY-J+1
         DO 4500 I=1,NLONRGG(JSN)
            NN = NUMBUS(I,J)
            IF (NN.GE.NNMIN) THEN
C
C           MEAN VALUES
C
               USME(I,J) = USME(I,J)/FLOAT(NN)
C
C           COMPUTE THE STANDARD DEVIATION.
C
               USSE2 =  (USSE(I,J)-USME(I,J)**2*FLOAT(NN))
     1                              / (FLOAT(NN)-1.)
               USSE(I,J) = SQRT( MAX(USSE2,0.))
            ENDIF
 4500    CONTINUE
 4600 CONTINUE
C
C ----------------------------------------------------------------------
C
C*    5. FLAG GRID POINTS, IF THERE ARE FEW OR STRONGLY SCATTERED DATA.
C        --------------------------------------------------------------
C
 5000 CONTINUE
      NWH = 0
      DO 5001 J=1,NGY
         JSN = NGY-J+1
         DO 5002 I=1,NLONRGG(JSN)
            NN = NUMBWH(I,J)
            IF (NN.EQ.0) THEN

C           IF THERE ARE NO MEASUREMENTS

               WHME(I,J)=-3.
            ELSE IF (NN.GT.0. AND. NN.LT.NNMIN) THEN

C           IF THERE ARE FEW MEASUREMENTS.

               WHME(I,J)=-2.
            ELSE IF (NN.GE.NNMIN) THEN

C           IF THE VARIANCE IS TOO LARGE.
C           OR SWH IS TOO SMALL  OR THERE ARE TOO MANY SPIKES

               XX = FLOAT(NR(I,J))/FLOAT(NN)
               RFWM = MAX(0.5,0.25*WHME(I,J))
               IF (WHSE(I,J).GT.RFWM .OR. WHME(I,J).LT.0.50 .OR.
     1             XX.GT.0.1) THEN
                  WHME(I,J) = -WHME(i,j)
               ELSE
                  NWH = NWH+1
               ENDIF
            ENDIF
 5002    CONTINUE
 5001 CONTINUE

      NUS = 0
      DO 5003 J=1,NGY
         JSN = NGY-J+1
         DO 5004 I=1,NLONRGG(JSN)
            NN = NUMBUS(I,J)
            IF (NN.EQ.0) THEN

C           IF THERE ARE NO MEASUREMENTS

               USME(I,J)=-3.
            ELSE IF (NN.GT.0. AND. NN.LT.NNMIN) THEN

C           IF THERE ARE FEW MEASUREMENTS.

               USME(I,J)=-2.
            ELSE IF (NN.GE.NNMIN)THEN

C           IF THE VARIANCE IS TOO LARGE.
C           OR THERE IS NO CORRESPONDING SWH

               RFU = 0.0
               IF (USME(I,J).GT.0.0001) RFU = USSE(I,J)/USME(I,J)
               IF (RFU.GE.0.5.OR.WHME(I,J).LE.0.) THEN
                  USME(I,J) = -USME(I,J)
               ELSE
                  NUS = NUS + 1
               ENDIF
            ENDIF
 5004    CONTINUE
 5003 CONTINUE
C
C*
C
      WRITE (IU06,*) ' GRDATA: GRIDDED MEASUREMENTS ARE SELECTED:'
      WRITE (IU06,*) '         WAVE HEIGHT AT NWH = ', NWH,' GRIDPOINTS'
      WRITE (IU06,*) '         WIND SPEED  AT NUS = ', NUS,' GRIDPOINTS'

      CLOSE (iume,IOSTAT=ioszx)
      IF(ITEST.GT.0)
     & WRITE(iu06,*) ' UNIT ',iume,' CLOSED WITH IOSTAT= ',ioszx
      CLOSE (iubl,IOSTAT=ioszx)
      IF(ITEST.GT.0)
     & WRITE(iu06,*) ' UNIT ',iubl,' CLOSED WITH IOSTAT= ',ioszx

      IF(ALLOCATED(ICEMASKG)) DEALLOCATE(ICEMASKG)
      DEALLOCATE(RLAICE,RLOICE,SWHICE,WSICE,CDAICE)

      RETURN
      END
