      SUBROUTINE OUTINT (IU06, IU20, IU21, IU30, ITEST, Kfdb)

C ----------------------------------------------------------------------
C
C**** *OUTINT* - OUTPUT OF INTEGRATED FIELDS IS WRITTEN.
C
C     H. GUNTHER       ECMWF  JUNE     1990
C     A. GUILLAUME     ECMWF  JUNE     1992 HANDLING OF MISSING VALUES,
C                                           AND FDB OUTPUT.
C     A. GUILLAUME     ECMWF  JULY     1992 DIRECTIONS.
C
C     P.A.E.M. JANSSEN ECMWF  AUGUST   1995
C
C     J. BIDLOT        ECMWF  APRIL    1996 MESSAGE PASSING.
C     J. BIDLOT        ECMWF  OCTOBER  1997 ADD OUTPUT OF U10G.
C     B. HANSEN        ECMWF  FEBRUARY 1998 INTERFACE TO FDB.
C
C*    PURPOSE.
C     --------
C
C       OUTPUT OF WAVE AND WIND FIELDS IN GRIDDED FORMAT
C       TO DISK ,PRINTER AND FDB.
C
C
C**   INTERFACE.
C     ----------
C       *CALL* *OUTINT (IU06, IU20, IU21, IU30, ITEST, Kfdb)
C          *IU06*   - PRINTER OUTPUT UNIT.
C          *IU20*   - OUTPUT UNIT OF PARAMETERS OF TOTAL SPECTRUM.
C          *IU21    - OUTPUT UNIT OF SEA AND SWELL PARAMETER.
C          *IU30    - OUTPUT UNIT OF GRIBBED PARAMETERS OF TOTAL
C                     SPECTRUM.
C          *ITEST*  - TESTOUTPUT GIVEN IF ITEST LT 2.
C          *Kfdb*   - REFERENCE TO THE FDB DATABASE.
C
C     EXTERNALS.
C     ----------
C
C       *OUTPP*     - PRINTER OUTPUT OF AN ARRAY.
C       *GRIBPAC*   -
C       *MPGATHERGRDFLD* - GATHER GRID FIELD ONTO ONE PROCESS
C
C     METHOD.
C     -------
C
C       NONE.
C
C     REFERENCE.
C     ----------
C
C       NONE.
C
C ----------------------------------------------------------------------
C
#include "param.h"
C
#include "parcons.h"
C
#include "comcoup.h"
C
#include "txtcoup.h"
C
#include "comcout.h"
C
#include "comgrid.h"
C
#include "comintp.h"
C
#include "comints.h"
C
#include "commap.h"
C
#include "commesspass.h"
C
#include "commpp.h"
C
#include "comstat.h"
C
C ----------------------------------------------------------------------
C
      CHARACTER TITL*100, CLDOMAIN*1
      LOGICAL LGRIBPAC
C
C ----------------------------------------------------------------------
C
      LGRIBPAC = .FALSE.

C     0. COLLECT INTEGRATED PARAMETERS FOR OUTPUT ON A SPECIFIC PROCESS
C        --------------------------------------------------------------
C
      IF(LMESSPASS.AND.(FFLAG20.OR.PFLAG20.OR.GFLAG20)) THEN
C     0.1 COLLECT TOTAL SEA PARAMETER ARRAYS TO PROCESS NPROC 
C         ---------------------------------------------------
        IRECV=NPROC

        IF (FFLAG(1)) CALL MPGATHERGRDFLD(IPFGTBL(1),IRECV,1,WHGTTG)
        IF (FFLAG(2)) CALL MPGATHERGRDFLD(IPFGTBL(2),IRECV,2,WDIRTG)
        IF (FFLAG(3)) CALL MPGATHERGRDFLD(IPFGTBL(3),IRECV,3,WMNFTG)
        IF (FFLAG(4)) CALL MPGATHERGRDFLD(IPFGTBL(4),IRECV,4,USTARG)
        IF (FFLAG(5)) CALL MPGATHERGRDFLD(IPFGTBL(5),IRECV,5,UDIRG)
        IF (FFLAG(6)) CALL MPGATHERGRDFLD(IPFGTBL(6),IRECV,6,WPKFTG)
        IF (FFLAG(7)) CALL MPGATHERGRDFLD(IPFGTBL(7),IRECV,7,CDG)
        IF (FFLAG(8)) CALL MPGATHERGRDFLD(IPFGTBL(8),IRECV,8,TAUWG)
        IF (FFLAG(9)) CALL MPGATHERGRDFLD(IPFGTBL(9),IRECV,9,SMEANG)
        IF (FFLAG(10))CALL MPGATHERGRDFLD(IPFGTBL(10),IRECV,10,U10G)

      ENDIF

      IF(LMESSPASS) THEN
C     0.2 COLLECT INTEGRATED PARAMETER ARRAYS TO PROCESS 1
C         ONLY FOR PRINT OUTPUT OR GRIB IF OUTPUT TO FILE
C         ------------------------------------------------
        IRECV=1

        IF (PFLAG(1).OR.(.NOT.LFDB.AND.GFLAG(1)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(1),IRECV,101,WHGTTG)
        IF (PFLAG(2).OR.(.NOT.LFDB.AND.GFLAG(2)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(2),IRECV,102,WDIRTG)
        IF (PFLAG(3).OR.(.NOT.LFDB.AND.GFLAG(3)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(3),IRECV,103,WMNFTG)
        IF (PFLAG(4).OR.(.NOT.LFDB.AND.GFLAG(4)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(4),IRECV,104,USTARG)
        IF (PFLAG(5).OR.(.NOT.LFDB.AND.GFLAG(5)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(5),IRECV,105,UDIRG)
        IF (PFLAG(6).OR.(.NOT.LFDB.AND.GFLAG(6)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(6),IRECV,106,WPKFTG)
        IF (PFLAG(7).OR.(.NOT.LFDB.AND.GFLAG(7)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(7),IRECV,107,CDG)
        IF (PFLAG(8).OR.(.NOT.LFDB.AND.GFLAG(8)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(8),IRECV,108,TAUWG)
        IF (PFLAG(9).OR.(.NOT.LFDB.AND.GFLAG(9)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(9),IRECV,109,SMEANG)
        IF (PFLAG(10).OR.(.NOT.LFDB.AND.GFLAG(10)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(10),IRECV,110,U10G)
        IF (PFLAG(11).OR.(.NOT.LFDB.AND.GFLAG(11)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(11),IRECV,111,WHGTWG)
        IF (PFLAG(12).OR.(.NOT.LFDB.AND.GFLAG(12)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(12),IRECV,112,WHGTSG)
        IF (PFLAG(13).OR.(.NOT.LFDB.AND.GFLAG(13)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(13),IRECV,113,WDIRWG)
        IF (PFLAG(14).OR.(.NOT.LFDB.AND.GFLAG(14)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(14),IRECV,114,WDIRSG)
        IF (PFLAG(15).OR.(.NOT.LFDB.AND.GFLAG(15)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(15),IRECV,115,WMNFWG)
        IF (PFLAG(16).OR.(.NOT.LFDB.AND.GFLAG(16)))
     1      CALL MPGATHERGRDFLD(IPFGTBL(16),IRECV,116,WMNFSG)

      ENDIF
C
      IF(LMESSPASS.AND.(FFLAG21.OR.PFLAG21.OR.GFLAG21)) THEN
C     0.3 COLLECT SWELL AND WIND SEA PARAMETER ARRAYS TO PROCESS 2 or 1
C         -------------------------------------------------------------
        IRECV=MIN(2,NPROC)

        IF (FFLAG(11)) CALL MPGATHERGRDFLD(IPFGTBL(11),IRECV,11,WHGTWG)
        IF (FFLAG(12)) CALL MPGATHERGRDFLD(IPFGTBL(12),IRECV,12,WHGTSG)
        IF (FFLAG(13)) CALL MPGATHERGRDFLD(IPFGTBL(13),IRECV,13,WDIRWG)
        IF (FFLAG(14)) CALL MPGATHERGRDFLD(IPFGTBL(14),IRECV,14,WDIRSG)
        IF (FFLAG(15)) CALL MPGATHERGRDFLD(IPFGTBL(15),IRECV,15,WMNFWG)
        IF (FFLAG(16)) CALL MPGATHERGRDFLD(IPFGTBL(16),IRECV,16,WMNFSG)

      ENDIF


C*    1. OUTPUT OF INTEGRATED PARAMETERS OF TOTAL SEA.
C        ---------------------------------------------
C
      DNGX = FLOAT(NGX)
      DNGY = FLOAT(NGY)
      IF (CDTINTT.EQ.CDTPRO) THEN
C
C*    1.1 CONVERT TO WAVE HEIGHT AND DEGREES.
C         -----------------------------------
C
C
           DO J=1,NGY
             DO I=1,NGX
               IF(WHGTTG(I,J).NE.ZMISS) THEN
                 WHGTTG(I,J) = 4.0 * SQRT(WHGTTG(I,J))
               ENDIF
               IF(WDIRTG(I,J).NE.ZMISS) THEN
                 WDIRTG(I,J) = MOD(DEG * WDIRTG(I,J) +180.,360.)
               ENDIF
               IF(UDIRG(I,J).NE.ZMISS) THEN
                 UDIRG(I,J) = MOD(DEG * UDIRG(I,J) +180.,360.)
               ENDIF
             ENDDO
           ENDDO
C
C*    1.2 WRITE PARAMETERS OF TOTAL SEA AND WIND TO FILE.
C         -----------------------------------------------
C
         
         IF(FFLAG20) THEN
         IF((LMESSPASS.AND.(IRANK.EQ.NPROC)).OR.(.NOT.LMESSPASS)) THEN
            WRITE (IU20) CDTPRO, DNGX, DNGY, AMOWEP,AMOSOP,AMOEAP,
     1                   AMONOP, ZDELLO
            WRITE (IU20) (FFLAG(J),J=1,9)
            IF (FFLAG(1)) WRITE (IU20) WHGTTG
            IF (FFLAG(2)) WRITE (IU20) WDIRTG
            IF (FFLAG(3)) WRITE (IU20) WMNFTG
            IF (FFLAG(4)) WRITE (IU20) USTARG
            IF (FFLAG(5)) WRITE (IU20) UDIRG
            IF (FFLAG(6)) WRITE (IU20) WPKFTG
            IF (FFLAG(7)) WRITE (IU20) CDG
            IF (FFLAG(8)) WRITE (IU20) TAUWG
            IF (FFLAG(9)) WRITE (IU20) SMEANG
            IF (FFLAG(10))WRITE (IU20) U10G
         ENDIF
         ENDIF
C
C*    1.3 PRINT PARAMETERS OF TOTAL SEA AND WIND.
C         ---------------------------------------
C
         IF (PFLAG(1)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' SIGNIFICANT WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     1                 WHGTTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(2)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WAVE DIRECTION ( DEGREE FROM NORTH TO )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     1                 WDIRTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(3)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WAVE FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 WMNFTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(4)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' FRICTION VELOCITY ( 0.01 METRES/SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 USTARG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(5)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' FRICTION DIRECTION ( DEGREE FROM NORTH TO )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     1                 UDIRG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(6)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' PEAK WAVE FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 WPKFTG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(7)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' DRAG COEFFICIENT (IN  0.1 *PROMILLE )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10000.,
     1                 CDG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(8)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' NORMALISED WAVE STRESS (IN % )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 TAUWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(9)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SQUARE SLOPE (IN % )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 SMEANG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(10)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' WIND SPEED ( METRES/SECOND )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     1                 U10G, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    2. OUTPUT OF INTEGRATED PARAMETERS OF SEA AND SWELL.
C        -------------------------------------------------
C
      IF (CDTINTS.EQ.CDTPRO) THEN
C
C*    2.1 CONVERT TO WAVE HEIGHT AND DEGREES.
C         -----------------------------------
C
C
           DO J=1,NGY
             DO I=1,NGX
               IF(WHGTSG(I,J).NE.ZMISS) THEN
                 WHGTSG(I,J) = 4.0 * SQRT(WHGTSG(I,J))
               ENDIF
               IF(WDIRSG(I,J).NE.ZMISS) THEN
                 WDIRSG(I,J) = MOD(DEG * WDIRSG(I,J) + 180.,360.)
               ENDIF
               IF(WHGTWG(I,J).NE.ZMISS) THEN
                 WHGTWG(I,J) = 4.0 * SQRT(WHGTWG(I,J))
               ENDIF
               IF(WDIRWG(I,J).NE.ZMISS) THEN
                 WDIRWG(I,J) = MOD(DEG * WDIRWG(I,J) + 180.,360.)
               ENDIF
             ENDDO
           ENDDO
C
C
C*    2.2 WRITE PARAMETERS OF SWELL AND WIND DIRECTION TO FILE.
C         -----------------------------------------------------
C
         IF(FFLAG21) THEN
         IF((LMESSPASS.AND.(IRANK.EQ.MIN(2,NPROC)))
     $      .OR.(.NOT.LMESSPASS)) THEN
            WRITE (IU21) CDTPRO, DNGX, DNGY, AMOWEP,AMOSOP,AMOEAP,
     1                   AMONOP, ZDELLO
            WRITE (IU21) (FFLAG(J),J=11,16)
            IF (FFLAG(11)) WRITE (IU21) WHGTWG
            IF (FFLAG(13)) WRITE (IU21) WDIRWG
            IF (FFLAG(15)) WRITE (IU21) WMNFWG
            IF (FFLAG(12)) WRITE (IU21) WHGTSG
            IF (FFLAG(14)) WRITE (IU21) WDIRSG
            IF (FFLAG(16)) WRITE (IU21) WMNFSG
         ENDIF
         ENDIF
C
C
C*    2.3 PRINT PARAMETERS OF SWELL AND WIND WAVE DIRECTION.
C         --------------------------------------------------
         IF (PFLAG(11)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' WIND SEA WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     1                 WHGTWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(13)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SEA DIRECTION ( DEGREE FROM NORTH FROM )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     1                 WDIRWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(15)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN WIND SEA FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 WMNFWG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(12)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' SWELL WAVE HEIGHT ( 0.1 METRES )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 10.,
     1                 WHGTSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(14)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL DIRECTION ( DEGREE FROM NORTH FROM )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 1.,
     1                 WDIRSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
         IF (PFLAG(16)) THEN
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            TITL = ' MEAN SWELL FREQUENCY ( 0.01 HERTZ )'
            CALL OUTPP (CDTPRO, IU06, NGX, NGY, NGX, NGY, TITL, 100.,
     1                 WMNFSG, AMOWEP, AMOSOP, AMOEAP, AMONOP)
          ENDIF
         ENDIF
      ENDIF
C
C ----------------------------------------------------------------------
C
C*    3.  PACK INTEGRATED PARAMETERS INTO GRIB FM94.
C         ------------------------------------------
C
      IF(IRGG.EQ.1) THEN
        NTGRIB=0
        DO J=1,NGY
          NTGRIB=NTGRIB+NLONRGG(J)
        ENDDO
      ELSE
#if region=='g'
        NTGRIB=(NINT( 180./((AMONOP-AMOSOP)/REAL(NGY-1)) )+1)*NGX
#elif region=='m'
        NTGRIB=NGY*NGX
#endif
      ENDIF

      IF (GFLAG20) THEN
        IF(CDTPRO.LE.CDATEF) THEN
C
C*    3.1.  THIS IS AN ANALYSIS DATE
C           ------------------------
C
           READ(CDTPRO,'(I4,4I2)') IY, IM, ID, IH, IMN
           IFCST = 0
        ELSE
C
C*    3.2.  THIS IS A  FORECAST DATE
C           ------------------------
C
          READ(CDATEF,'(I4,4I2)')IY,IM,ID,IH,IMN
          CALL DIFDATE(CDATEF,CDTPRO,IFCST)
          IFCST=IFCST/3600
        ENDIF
C
C*       WAVE HEIGHTS.
C        -------------
C
        IF (GFLAG(1)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(1)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=229
          CALL GRIBPAC(IU06,IU30,NTGRIB,WHGTTG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*          WAVE DIRECTION.
C           ---------------
C
        IF (GFLAG(2)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(2)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=230
          CALL GRIBPAC(IU06,IU30,NTGRIB,WDIRTG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*          WAVE MEAN PERIOD.
C           -----------------
C
        IF (GFLAG(3)) THEN
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(3)))
     1      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2      .OR.(.NOT.LMESSPASS)) THEN
            IPARAM=232
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WMNFTG(JX,JY).NE.0.)THEN
                  IF(WMNFTG(JX,JY).NE.ZMISS)THEN
                    WMNFTG(JX,JY)=1./WMNFTG(JX,JY)
                  ELSE
                    WMNFTG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WMNFTG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTGRIB,WMNFTG,NGX,NGY,IPARAM,IY,IM,
     1                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     2                   YEXPVER,CFDBSF,NENSFNB,NTOTENS,Kfdb)
            LGRIBPAC = .TRUE.
          ENDIF
        ENDIF
C
C*          WAVE PEAK PERIOD.
C           -----------------
C
        IF (GFLAG(6)) THEN
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(6)))
     1      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2      .OR.(.NOT.LMESSPASS)) THEN
            IPARAM=231
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WPKFTG(JX,JY).NE.0.)THEN
                  IF(WPKFTG(JX,JY).NE.ZMISS)THEN
                    WPKFTG(JX,JY)=1./WPKFTG(JX,JY)
                  ELSE
                    WPKFTG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WPKFTG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTGRIB,WPKFTG,NGX,NGY,IPARAM,IY,IM,
     1                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     2                   YEXPVER,CFDBSF,NENSFNB,NTOTENS,Kfdb)
            LGRIBPAC = .TRUE.
          ENDIF
        ENDIF
C
C*          DRAG COEFFICIENT
C           ----------------
C
        IF (GFLAG(7)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(7)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=233
          IZLEV=NINT(XNLEV(1))
          CALL GRIBPAC(IU06,IU30,NTGRIB,CDG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,IZLEV,LFDB,ITEST,MARSTYPE,YCLASS,
     2                 YEXPVER,CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*          MEAN SQUARE SLOPE.
C           ------------------
C
        IF (GFLAG(9)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(9)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=244
          CALL GRIBPAC(IU06,IU30,NTGRIB,SMEANG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*          WIND SPEED .
C           -----------
C
        IF (GFLAG(10)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(10)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=245
          IZLEV=NINT(XNLEV(1))
          CALL GRIBPAC(IU06,IU30,NTGRIB,U10G,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,IZLEV,LFDB,ITEST,MARSTYPE,YCLASS,
     2                 YEXPVER,CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF

      ENDIF
C
C ----------------------------------------------------------------------
C
C*    4.  PACK WIND SEA AND SWELL PARAMETERS INTO GRIB FM94.
C         -------------------------------------------------
C
      IF (GFLAG21) THEN
         IF(CDTPRO.LE.CDATEF) THEN
C
C*    4.1.  THIS IS AN ANALYSIS DATE
C           ------------------------
C
           READ(CDTPRO,'(I4,4I2)') IY, IM, ID, IH, IMN
           IFCST = 0

         ELSE
C
C*    4.2.  THIS IS A  FORECAST DATE
C           ------------------------
C
          READ(CDATEF,'(I4,4I2)')IY,IM,ID,IH,IMN
          CALL DIFDATE(CDATEF,CDTPRO,IFCST)
          IFCST=IFCST/3600.
        ENDIF
C
C*       WIND SEA WAVE HEIGHTS.
C        ---------------------
C
        IF (GFLAG(11)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(11)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=234
          CALL GRIBPAC(IU06,IU30,NTGRIB,WHGTWG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*       WIND SEA WAVE DIRECTION.
C           ---------------
C
C
        IF (GFLAG(13)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(13)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=235
          CALL GRIBPAC(IU06,IU30,NTGRIB,WDIRWG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*       WIND SEA MEAN PERIOD.
C        --------------------
C
        IF (GFLAG(15)) THEN
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(15)))
     1      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2      .OR.(.NOT.LMESSPASS)) THEN
            IPARAM=236
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WMNFWG(JX,JY).NE.0.)THEN
                  IF(WMNFWG(JX,JY).NE.ZMISS)THEN
                    WMNFWG(JX,JY)=1./WMNFWG(JX,JY)
                  ELSE
                    WMNFWG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WMNFWG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTGRIB,WMNFWG,NGX,NGY,IPARAM,IY,IM,
     1                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     2                   YEXPVER,CFDBSF,NENSFNB,NTOTENS,Kfdb)
            LGRIBPAC = .TRUE.
          ENDIF
        ENDIF
C
C*       SWELL WAVE HEIGHTS.
C        -------------------
C
        IF (GFLAG(12)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(12)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=237
          CALL GRIBPAC(IU06,IU30,NTGRIB,WHGTSG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*       SWELL WAVE DIRECTION.
C        --------------------
C
C
        IF (GFLAG(14)) THEN
         IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(14)))
     1   .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2   .OR.(.NOT.LMESSPASS)) THEN
          IPARAM=238
          CALL GRIBPAC(IU06,IU30,NTGRIB,WDIRSG,NGX,NGY,IPARAM,IY,IM,
     1                 ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,YEXPVER,
     2                 CFDBSF,NENSFNB,NTOTENS,Kfdb)
          LGRIBPAC = .TRUE.
         ENDIF
        ENDIF
C
C*       SWELL SEA MEAN PERIOD.
C        ---------------------
C
        IF (GFLAG(16)) THEN
          IF ((LMESSPASS.AND.LFDB.AND.(IRANK.EQ.IPFGTBL(16)))
     1      .OR.(LMESSPASS.AND.(.NOT.LFDB).AND.(IRANK.EQ.1))
     2      .OR.(.NOT.LMESSPASS)) THEN
            IPARAM=239
            DO JY=1,NGY
              DO JX=1,NGX
                IF(WMNFSG(JX,JY).NE.0.)THEN
                  IF(WMNFSG(JX,JY).NE.ZMISS)THEN
                    WMNFSG(JX,JY)=1./WMNFSG(JX,JY)
                  ELSE
                    WMNFSG(JX,JY)=ZMISS
                  ENDIF
                ELSE
                  WMNFSG(JX,JY)=0.
                ENDIF
              ENDDO
            ENDDO
            CALL GRIBPAC(IU06,IU30,NTGRIB,WMNFSG,NGX,NGY,IPARAM,IY,IM,
     1                   ID,IH,IFCST,0,LFDB,ITEST,MARSTYPE,YCLASS,
     2                   YEXPVER,CFDBSF,NENSFNB,NTOTENS,Kfdb)
            LGRIBPAC = .TRUE.
          ENDIF
        ENDIF
      ENDIF  
C
      RETURN
      END
