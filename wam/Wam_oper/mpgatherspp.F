      SUBROUTINE MPGATHERSPP(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD,
     &                      FL3,FL1,FLPTS,
     &                      EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,THSWELL,
     &                      THWISEA,U10,THW,US)

!****  *MPGATHERSSPP* - GATHERS        SPECTRUM AND SCALAR FIELDS ONTO
!****                   A SINGLE PROCESS FOR OUTPUT IN OUTSPP AT
!****                   SELECTED GRID POINTS

!     J. BIDLOT    ECMWF   APRIL 1996  MESSAGE PASSING

!     PURPOSE.
!     --------

!     GATHER NSPFLD SPECTRUM AND NSCFLD SCALAR FIELDS DISTRIBUTED 
!     ACROSS THE DIFFERENT PROCESSES ONTO THE SINGLE PROCESS IRECV
!     FOR SELECTED GRID POINTS.

!*    INTERFACE.
!     ----------

!     CALL *MPGATHERSPP*(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD, 
!    1                   FL1,FL3,FLPTS,
!    2                   EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,THSWELL, 
!    3                   THWISEA,U10,THW,US)


!     *IRECV*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED 
!     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS
!                 THE USE OF A SEPARATE TAG FOR EACH SPECIFIC CALL
!                 ENSURES THAT NO CONFLICT ARISES BETWEEN DIFFERENT
!                 MESSAGES.
!     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!     *NSPFLD*    NUMBER OF SPECTRUM FIELDS
!     *NSCFLD*    NUMBER OF SCALAR FILEDS
!     *FL3*       BLOCK OF SPECTRA.
!     *FL1*       SWELL SPECTRUM.
!     *FLPTS*     SPECTRA AT THE SELECTED GRID POINTS (OUTPUT) 
!     *EMEAN*     TOTAL ENERGY.
!     *FMEAN*     MEAN FREQUENCY (1./MEAN PERIOD).
!     *THQ*       MEAN WAVE DIRECTION IN RADIANS.
!     *SMEAN*     MEAN SQUARE SLOPE.
!     *ESWELL*    TOTAL SWELL ENERGY
!     *FSWELL*    MEAN SWELL FREQUENCY
!     *THSWELL*   MEAN SWELL DIRECTION
!     *THWISEA*   MEAN WIND INTERACTING SEA DIRECTION
!     *U10*    NEW WIND SPEED IN M/S.
!     *THW*    WIND DIRECTION IN RADIANS
!     *US*     NEW FRICTION VELOCITY IN M/S.

!     METHOD.
!     -------
!     MPL SEND OF A MESSAGE BUFFER CONTAINING ALL THE FIELDS AT ALL 
!     SELECTED GRID POINTS HELD BY THAT PROCESS TO THE PROCESS
!     CORRESPONDING TO IRECV
!     FOR ALL PROCESS EXCEPT FOR THE PROCESS CORRESPONDING TO IRECV
!     WHERE THE MESSAGE BUFFERS ARE RECEIVED AND REORDERED.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_SEND
!         MPL_RECV
!         MPL_ABORT 

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWCOUT  , ONLY : NGOUT    ,IJAR
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE MPL_MODULE
      USE YOMGSTATS, ONLY : LSYNCSTATS

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL1, FL3
      REAL,DIMENSION(NSPFLD,NGOUT,NANG,NFRE) :: FLPTS
      REAL,DIMENSION(NIBLO) :: EMEAN,FMEAN,THQ,SMEAN,ESWELL,FSWELL,
     &                         THSWELL,THWISEA
      REAL,DIMENSION(NIBLO) :: U10,THW,US  
!----------------------------------------------------------------------

      REAL,ALLOCATABLE :: ZCOMBUF(:)

      MAXLENGTH=(NSCFLD+NSPFLD*NANG*NFRE)*NGOUT
      ALLOCATE(ZCOMBUF(MAXLENGTH))

!     1.0 DEFAULT ACTION IF NO FIELD GATHERING 
!         ------------------------------------
      IF((IRECV.EQ.0).OR.(NPROC.EQ.1)) THEN
        DO NGOU=1,NGOUT
          IJ=IJAR(NGOU)
          DO  M = 1, NFRE
            DO K = 1, NANG
              FLPTS(1,NGOU,K,M)=FL3(IJ,K,M)
            ENDDO
          ENDDO
          DO  M = 1, NFRE
            DO K = 1, NANG
              FLPTS(2,NGOU,K,M)=FL1(IJ,K,M)
            ENDDO
          ENDDO
        ENDDO

      ELSEIF(IRANK.NE.IRECV) THEN
!     1.1 SEND TO THE PROCESS THAT GATHERS THE WHOLE FIELD
!         ------------------------------------------------

!     LOAD COMMUNICATION BUFFER
        KCOUNT=0
        DO NGOU=1,NGOUT

          IJ=IJAR(NGOU)
          IF((IJ.GE.NSTART(IRANK)).AND.(IJ.LE.NEND(IRANK))) THEN
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=EMEAN(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FMEAN(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=SMEAN(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=THQ(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=ESWELL(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FSWELL(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=THSWELL(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=THWISEA(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=U10(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=THW(IJ)
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=US(IJ)

            DO  M = 1, NFRE
              DO K = 1, NANG
                KCOUNT=KCOUNT+1
                ZCOMBUF(KCOUNT)=FL3(IJ,K,M)
              ENDDO
            ENDDO

            DO  M = 1, NFRE
              DO K = 1, NANG
                KCOUNT=KCOUNT+1
                ZCOMBUF(KCOUNT)=FL1(IJ,K,M)
              ENDDO
            ENDDO
          ENDIF

        ENDDO

!     SEND BUFFER
        IF (.NOT.LSYNCSTATS) CALL GSTATS(673,0)
        CALL MPL_SEND(ZCOMBUF(1:KCOUNT),KDEST=IRECV,KTAG=ITAG,
     &   CDSTRING='MPGATHERSPP:')
        IF (.NOT.LSYNCSTATS) CALL GSTATS(673,1)


      ELSE
!     1.2.1  RECEIVE CONTRIBUTION TO THE FIELD FROM THE OTHER PROCESSES
!            ---------------------------------------------------------- 
        DO IP=1,NPROC-1
          IF (.NOT.LSYNCSTATS) CALL GSTATS(673,0)
          CALL MPL_RECV(ZCOMBUF(1:MAXLENGTH),KTAG=ITAG,KFROM=KRFROM,
     &     KOUNT=KRCOUNT,KRECVTAG=KRTAG,CDSTRING='MPGATHERSPP:')
          IF(KRTAG.NE.ITAG) CALL MPL_ABORT
     &     ('MPL_RECV ERROR in MPGATHERSPP:  MISMATCHED TAGS' )
          IF (.NOT.LSYNCSTATS) CALL GSTATS(673,1)

!  DECODE BUFFER
          KCOUNT=0

          DO NGOU=1,NGOUT
            IJ=IJAR(NGOU)
            IF((IJ.GE.NSTART(KRFROM)).AND.(IJ.LE.NEND(KRFROM))) THEN
              KCOUNT=KCOUNT+1
              EMEAN(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              FMEAN(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              SMEAN(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              THQ(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              ESWELL(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              FSWELL(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              THSWELL(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              THWISEA(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              U10(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              THW(IJ)=ZCOMBUF(KCOUNT)
              KCOUNT=KCOUNT+1
              US(IJ)=ZCOMBUF(KCOUNT)

              DO  M = 1, NFRE
                DO K = 1, NANG
                  KCOUNT=KCOUNT+1
                  FLPTS(1,NGOU,K,M)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO

              DO  M = 1, NFRE
                DO K = 1, NANG
                  KCOUNT=KCOUNT+1
                  FLPTS(2,NGOU,K,M)=ZCOMBUF(KCOUNT)
                ENDDO
              ENDDO

            ENDIF
          ENDDO
          IF(KRCOUNT.NE.KCOUNT) CALL MPL_ABORT
     &     ('MPL_RECV ERROR in MPGATHERSPP: WRONG LENGTH OF MESSAGE')

        ENDDO
!     1.2.2  CONTRIBUTION FROM RECEIVING PROCESS 
!            ----------------------------------- 

        DO NGOU=1,NGOUT
          IJ=IJAR(NGOU)
          IF((IJ.GE.NSTART(IRANK)).AND.(IJ.LE.NEND(IRANK)))THEN
            DO  M = 1, NFRE
              DO K = 1, NANG
                FLPTS(1,NGOU,K,M)=FL3(IJ,K,M)
              ENDDO
            ENDDO
            DO  M = 1, NFRE
              DO K = 1, NANG
                FLPTS(2,NGOU,K,M)=FL1(IJ,K,M)
              ENDDO
            ENDDO
          ENDIF
        ENDDO

      ENDIF

      DEALLOCATE(ZCOMBUF)

      RETURN
      END SUBROUTINE MPGATHERSPP
