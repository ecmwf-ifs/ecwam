      SUBROUTINE MPDECOMP(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)

!****  *MPDECOMP* - EVEN DECOMPOSITION OF THE GRID DOMAIN
!****               AMONG PROCESSES AFTER INPUT OF PREPROC GRID
!                   AND UBUF FILES.

!     J. BIDLOT    ECMWF   MARCH 1996  MESSAGE PASSING
!     J. BIDLOT    ECMWF   JANUARY 1998 introduce NPR 
!     J. BIDLOT    ECMWF   OCTOBER 1998 COMPLETE READING OF
!                                       IU07 AND IU08 

!     PURPOSE.
!     --------
!     IT WILL READ IU07 AND IU08 AND IF MESSAGE PASSASING IT WILL
!     DETERMINE AN EVEN DECOMPOSITION OF THE GRID ARRAYS FOR USE ON
!     A DISTRIBUTED MEMORY COMPUTER USING A MESSAGE PASSING PROTOCOL
!     FOR THE EXCHANGE OF INFORMATION ACROSS THE DIFFERENT PROCESSES
!     THE LENGTH OF THE MESSAGES BETWEEN DIFFERENT PROCESS IS ALSO
!     COMPUTED.
!     IT IS ASSUMED THAT EACH GRID SUB DOMAIN WILL ONLY SEND OR
!     RECEIVE MESSAGES FROM THE PREVIOUS AND NEXT PROCESS 

!*    INTERFACE.
!     ----------
!     CALL *MPDECOMP*(NPR,NSTART,NEND,KLENTOP,KLENBOT,MAXLEN)

!      *NPR*       NUMBER OF SUBDOMAINS (USUALLY THE NUMBER OF PE'S )
!      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!      *KLENTOP*   LENGTH OF THE MESSAGE SENT BY THE NEXT
!                  SUB GRID DOMAIN (IN GRID POINT UNIT)
!      *KLENBOT*   LENGTH OF THE MESSAGE SENT BY THE PREVIOUS
!                  SUB GRID DOMAIN (IN GRID POINT UNIT)
!      *MAXLEN*    MAXIMUM NUMBER OF POINTS IN ANY SUB DOMAIN

!     METHOD.
!     -------
!     SINCE THE GLOBAL GRID DOMAIN IS MAPPED ONTO A 1 DIMENSIONAL ARRAY
!     FOLLOWING INCREASING LATITUDE LINES, SUB GRID DOMAINS ARE SEGMENTS
!     OF THAT LONG ARRAY. THE LENGTH OF THOSE SEGMENTS IS DETERMINED BY
!     DISTRIBUTING AS EVENLY AS POSSIBLE THE TOTAL NUMBER OF POINTS.

!     THE ADVECTION SCHEME IS SUCH THAT IT ONLY USES
!     NEIGHBORING GRID POINTS IN THE 2D GRID. THEREFORE,(see figure) 
!     FOR A GIVEN SEGMENT (PROCESS) INFORMATION FROM
!     PREVIOUS GRID POINTS ARE ONLY NEEDED FROM FEW GRID POINTS
!     STARTING FROM THE GRID POINT LOCATED, IN THE 2-D GRID,
!     JUST BELOW (ON THE PREVIOUS LATITUDE LINE) THE STARTING POINT
!     OF THAT PARTICULAR SEGMENT PROVIDED THIS PREVIOUS GRID POINT
!     IS A VALID SEA POINT. IF NOT A VALID GRID POINT, FIND THE NEXT ONE
!     WHICH IS A VALID SEA POINT BUT DOES NOT BELONG TO THE SEGEMENT
!     IN QUESTION.
!     NATURALLY, THE FIRST SEGMENT DOES NOT REQUIRE ANY INFORMATION
!     FROM PREVIOUS SEGMENT. 



!     FIGURE: 2-D GRID
!     -------

!       +    +    +    +    +    +    +    +
!       +    +    +    +    +    +    +    +
!       +    +    +    I    I    I    I    I
!       I    I    I    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    *    *
!       *    *    *    *    *    *    I    I
!       I    I    I    I    I    I    +    +
!       +    +    +    +    +    +    +    +
!       +    +    +    +    +    +    +    +

!     FIGURE LEGEND :
!     --------------
!     *  GRID POINTS OF THE SELECTED SEGMENT (PROCESS)
!     I  GRID POINTS OF THE PREVIOUS/FOLLOWING SEGEMENT FOR WHICH
!        INFORMATION IS NECESSARY IF THE ADVECTION SCHEME IS ONLY 
!        COMPUTED AT THE * POINTS

!     +  OTHER GRID POINTS THAT DO NOT AFFECT THE * POINTS 


!     SIMILARLY, INFORMATION FROM FOLLOWING GRID POINTS ARE ONLY NEEDED
!     UP TO THE GRID POINT JUST ABOVE (IN THE 2D GRID) THE LAST POINT
!     OF THE SELECTED SEGMENT, PROVIDED THIS FOLLOWING GRID POINT IS A
!     VALID SEA POINT, OTHERWISE, FIND THE FIRST PREVIOUS ONE WHICH IS
!     A VALID ONE BUT DOES NOT BELONG TO THE SEGMENT IN QUESTION.

!     EXTERNALS.
!     ----------
!          NONE

!     REFERENCES.
!     -----------
!          NONE

! -------------------------------------------------------------------

      USE YOWCOUP  , ONLY : LWCOU
      USE YOWGRID  , ONLY : IGL      ,IJS      ,IJL      ,IJLT
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP     ,KTAG 
      USE YOWPARAM , ONLY : NBLO     ,NIBLO
      USE YOWSTAT  , ONLY : IREFRA 
      USE YOWTEST  , ONLY : IU06
      USE YOWUBUF  , ONLY : KLAT     ,KLON
      USE YOWUNIT  , ONLY : IU07     ,IU08

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPR) :: NSTART,NEND,KLENTOP,KLENBOT
      CHARACTER*80 LOGFILENAME

! ----------------------------------------------------------------------

      INTEGER,ALLOCATABLE :: ICOMBUF(:), KDUM(:,:) 

!     STANDARD OUTPUT UNIT
!     --------------------
      IF(LMESSPASS.AND.LWCOU) THEN
        IU06=20
      ELSEIF (LMESSPASS.AND..NOT.LWCOU) THEN
        IU06=66
        LOGFILENAME='logfile.%p'
        CALL EXPAND_STRING(IRANK,NPR,0,0,LOGFILENAME,1)
        OPEN(IU06,FILE=LOGFILENAME,STATUS='UNKNOWN')
        WRITE(*,*) ' STDOUT FILE ', LOGFILENAME, ' OPENED ON PE ', IRANK
        WRITE(IU06,*) ' STDOUT OF PE ', IRANK
        WRITE(IU06,*) ' ============ '
        CALL FLUSH (IU06)
      ELSE
        IU06=6
      ENDIF

!     GRID AND UBUF FILES UNITS
!     -------------------------

      IF(LWCOU) THEN
        IU07=29
        IU08=31
      ELSE
        IU07=7
        IU08=8
      ENDIF

!     0. READ GRID INPUT FROM PREPROC 
!        ----------------------------

      IREAD=1
      KTAG=1

      IF((LMESSPASS.AND.IRANK.EQ.IREAD).OR.(.NOT.LMESSPASS)) THEN
        CALL READPRE (IU07, IREFRA)
        CLOSE (UNIT=IU07, STATUS='KEEP')
      ENDIF
!     SEND INFORMATION FROM READPRE TO ALL PE's
      IF(LMESSPASS) CALL MPBCASTGRID(IU06,IREAD,KTAG,IREFRA)

      WRITE(IU06,*) '  PREPROC GRID INFORMATION READ IN '


!*    1. INPUT NEIGHBOURING GRID POINT INDICES (UBUF) 
!        --------------------------------------------

      ALLOCATE(KLAT(NIBLO,2))  ! THE SIZE OF KLAT IS READJUSTED SEE BELOW
      ALLOCATE(KLON(NIBLO,2))  ! THE SIZE OF KLON IS READJUSTED SEE BELOW

      IF(LMESSPASS) THEN

!       READ IOU8 ON PE IREAD
        IF(IRANK.EQ.IREAD) THEN
          READ (IU08) KLAT, KLON
        ENDIF

!       SEND KLAT AND KLON TO OTHER PE'S
        IF(NPR.GT.1) THEN
          ITAG=KTAG
          MPLENGTH=4*NIBLO
          ALLOCATE(ICOMBUF(MPLENGTH))

          IF(IRANK.EQ.IREAD) THEN
            KCOUNT=0
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                ICOMBUF(KCOUNT)=KLAT(IJ,IC)
              ENDDO
            ENDDO
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                ICOMBUF(KCOUNT)=KLON(IJ,IC)
              ENDDO
            ENDDO

            CALL MPE_BROADCAST(ICOMBUF,MPLENGTH,1,IRANK,
     &       ITAG,0,0,0,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     &       ('MPE_BROADCAST ERROR in MPDECOMP' )
          ELSE

            CALL MPE_RECV(ICOMBUF,MPLENGTH,1,IREAD,ITAG,0,0,0,
     &       KRCOUNT,KRFROM,KRTAG,IERR)
            IF(IERR.LT.0) CALL MPEI_ABORT
     &       ('MPE_RECV ERROR in MPDECOMP' )
            IF(KRCOUNT.NE.MPLENGTH) CALL MPEI_ABORT
     &       ('MPE_RECV ERROR in MPDECOMP:MISMATCHED MSG LENGTH')
            IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     &       ('MPE_RECV ERROR in MPRECOMP MISMATCHED TAGS' )


            KCOUNT=0
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                KLAT(IJ,IC)=ICOMBUF(KCOUNT)
              ENDDO
            ENDDO
            DO IC=1,2
              DO IJ=1,NIBLO
                KCOUNT=KCOUNT+1
                KLON(IJ,IC)=ICOMBUF(KCOUNT)
              ENDDO
            ENDDO

          ENDIF

          DEALLOCATE(ICOMBUF)

        ENDIF

      ELSEIF(.NOT.LMESSPASS) THEN
        READ (IU08) KLAT, KLON
      ENDIF


!*    2. FIND THE NUMBER OF POINTS PER PROCESS, THE START AND END INDEX
!        --------------------------------------------------------------
      NMEAN=IJL(NBLO)/NPR
      NREST=IJL(NBLO)-NMEAN*NPR
      
      NSTART(1)=1
      IF(NREST.GT.0) THEN
        NPTS=NMEAN+1
        NREST=NREST-1 
      ELSE
        NPTS=NMEAN
      ENDIF
      NEND(1)=NSTART(1)+NPTS-1

      DO IP=2,NPR
        NSTART(IP)=NSTART(IP-1)+NPTS
        IF(NREST.GT.0) THEN
          NPTS=NMEAN+1
          NREST=NREST-1 
        ELSE
          NPTS=NMEAN
        ENDIF
        NEND(IP)=NSTART(IP)+NPTS-1
      ENDDO


!     3. DETERMINE THE LENGTH OF THE MESSAGE THAT WILL BE EXCHANGED 
!        BETWEEN NEIGHBORING SUB GRID DOMAIN
!        -----------------------------------------------------------

      KLENBOT(1)=0
      DO IP=2,NPR
        I=-1
20      I=I+1
        ISB=KLAT(NSTART(IP)+I,1)
        IF(ISB.LE.0) GOTO 20

        KLENBOT(IP)=MAX(0,NSTART(IP)-ISB)
      ENDDO

      KLENTOP(NPR)=0
      DO IP=1,NPR-1
        I=-1
40      I=I+1
        IET=KLAT(NEND(IP)-I,2)
        IF(IET.LE.0) GOTO 40

        KLENTOP(IP)=MAX(0,IET-NEND(IP))
      ENDDO

      MAXLEN=0
      DO IP=1,NPR
        MAXLEN=MAX(MAXLEN,NEND(IP)-NSTART(IP)+1)
      ENDDO


!     4. KEEP THE PART OF KLAT AND KLON WHICH IS NECESSARY
!        -------------------------------------------------

      IF(LMESSPASS.AND.NPR.GT.1) THEN

        ALLOCATE(KDUM(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM(IJ,IC)=KLAT(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLAT)
        ALLOCATE(KLAT(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLAT(IJ,IC)=KDUM(IJ,IC)
          ENDDO
        ENDDO

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KDUM(IJ,IC)=KLON(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KLON)
        ALLOCATE(KLON(NSTART(IRANK):NEND(IRANK),2))

        DO IC=1,2
          DO IJ=NSTART(IRANK),NEND(IRANK)
            KLON(IJ,IC)=KDUM(IJ,IC)
          ENDDO
        ENDDO

        DEALLOCATE(KDUM)

      ENDIF

!     5. MODIFY KLAT AND KLON SUCH THAT POINT INDICES FOR LAND IS
!        NINF-1.
!        ---------------------------------------------------------

      NINF=NSTART(IRANK)-KLENBOT(IRANK)
      NSUP=NEND(IRANK)+KLENTOP(IRANK)
      IF(.NOT.LMESSPASS) NSUP=NIBLO

      IF(LMESSPASS) THEN
        DO IC=1,2
          DO IJ = NSTART(IRANK),NEND(IRANK)
            IF(KLAT(IJ,IC).EQ.0) KLAT(IJ,IC) = NINF-1
            IF(KLON(IJ,IC).EQ.0) KLON(IJ,IC) = NINF-1
          ENDDO
        ENDDO
      ENDIF

!     6. TEST AND RESET
!        -------------- 
      IF(LMESSPASS) THEN
        IF(IGL.GT.1) THEN
          WRITE (IU06,*) '******************************************'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '*    FATAL ERROR                         *'
          WRITE (IU06,*) '*    ===========                         *'
          WRITE (IU06,*) '* A MULTI BLOCK VERSION RUN WAS REQUESTED*'
          WRITE (IU06,*) '* THIS OPTION IS NOT IMPLEMENTED FOR     *'
          WRITE (IU06,*) '* MESSAGE PASSING PROTOCOL               *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '* PROGRAM ABORTS.   PROGRAM ABORTS.      *'
          WRITE (IU06,*) '*                                        *'
          WRITE (IU06,*) '******************************************'
          CALL ABORT1
        ELSE
! SPECIFY THE LIMITS OF THE SEGMENT OVER WHICH THE PE HAS DIRECT ACCESS
          IJS(1)=NSTART(IRANK)
          IJL(1)=NEND(IRANK)
          IJLT(1)=NSUP
        ENDIF
      ENDIF

      RETURN
      END SUBROUTINE MPDECOMP
