      SUBROUTINE MPBCASTSARIN(ISEND, ITAG, LNODATA) 

!****  *MPBCASTSARIN* - BROADCASTS INPUT SAR DATA 

!     J. BIDLOT    ECMWF   NOVEMBER 1999

!     PURPOSE.
!     --------

!     BROADCASTS SAR INPUT FROM PE ISEND TO ALL OTHER PE'S 

!*    INTERFACE.
!     ----------

!     CALL *MPBCASTSARIN*(ISEND, ITAG)

!     *ISEND*     PE ON WHICH THE SAR INPUT WAS READ IN.
!     *ITAG*      TAG ASSOCIATED WITH A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS
!                 THE USE OF A SEPARATE TAG FOR EACH SPECIFIC CALL
!                 SHOULD ENSURE THAT NO CONFLICT ARISES BETWEEN 
!                 DIFFERENT MESSAGES

!     METHOD.
!     -------

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE YOWICE   , ONLY : FLMINFR
      USE YOWMPP   , ONLY : IRANK    ,NPROC
      USE YOWPARAM , ONLY : NANG     ,NFRE
      USE YOWSARAS , ONLY : NSPEC    ,SPEC     ,LONG     ,LAT     ,
     &              XKLOW   ,XKUP    ,XKCUTOFF
      USE MPL_MODULE

!----------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER :: ISEND, ITAG, MZC, KCOUNT
      INTEGER :: ISPEC, IFRE, IANG
!      INTEGER :: BUFFER_SIZE_BCSARIN
      INTEGER :: IBUF(2)
      REAL,ALLOCATABLE :: ZCOMBUF(:)
      LOGICAL :: LNODATA


      IF(NPROC.EQ.1) THEN
        RETURN
      ELSE
        IF(IRANK.EQ.ISEND) THEN
          IF(LNODATA) THEN
            IBUF(1) = 1
          ELSE
            IBUF(1) = 0
          ENDIF
          IBUF(2) = NSPEC
        ENDIF

        CALL MPL_BROADCAST(IBUF,KROOT=ISEND,KTAG=ITAG,
     &    CDSTRING='MPBCASTSARIN 0:')
        ITAG=ITAG+1

        IF(IRANK.NE.ISEND) THEN
          IF(IBUF(1).EQ.1) THEN
            LNODATA = .TRUE.
          ELSE
            LNODATA = .FALSE.
          ENDIF
          NSPEC=IBUF(2)
        ENDIF

        IF(LNODATA) GOTO 5000

        MZC=BUFFER_SIZE_BCSARIN(NSPEC,NANG,NFRE)
        ALLOCATE(ZCOMBUF(MZC))

        IF(IRANK.EQ.ISEND) THEN
          KCOUNT=0
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=LONG(ISPEC,2) 
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=LAT(ISPEC,2) 
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=XKLOW(ISPEC) 
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=XKUP(ISPEC) 
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=XKCUTOFF(ISPEC) 
          ENDDO 
!!!  always check that the size of ZCOMBUF is in agreement with the
!!!  number of entries (see BUFFER_SIZE_BCSARIN) !!!
          DO IFRE = 1,NFRE
            DO IANG = 1,NANG
              DO ISPEC = 0,NSPEC
                KCOUNT=KCOUNT+1
                ZCOMBUF(KCOUNT)=SPEC(ISPEC,IANG,IFRE,2)
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF

        CALL MPL_BROADCAST(ZCOMBUF(1:MZC),KROOT=ISEND,KTAG=ITAG,
     &    CDSTRING='MPBCASTSARIN 1:')
        ITAG=ITAG+1

        IF(IRANK.NE.ISEND) THEN
          ALLOCATE(SPEC(0:NSPEC,NANG,NFRE,2))
          DO IFRE = 1,NFRE
            DO IANG = 1,NANG
              DO ISPEC = 0,NSPEC
                KCOUNT=KCOUNT+1
                SPEC(ISPEC,IANG,IFRE,1)=0.25*FLMINFR(1,IFRE)
              ENDDO 
            ENDDO 
          ENDDO 
          ALLOCATE(LONG(NSPEC,2))
          ALLOCATE(LAT(NSPEC,2))
          ALLOCATE(XKLOW(NSPEC))
          ALLOCATE(XKUP(NSPEC))
          ALLOCATE(XKCUTOFF(NSPEC))


          KCOUNT=0
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            LONG(ISPEC,2)=ZCOMBUF(KCOUNT)
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            LAT(ISPEC,2)=ZCOMBUF(KCOUNT)
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            XKLOW(ISPEC)=ZCOMBUF(KCOUNT) 
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            XKUP(ISPEC)=ZCOMBUF(KCOUNT) 
          ENDDO 
          DO ISPEC = 1,NSPEC
            KCOUNT=KCOUNT+1
            XKCUTOFF(ISPEC)=ZCOMBUF(KCOUNT) 
          ENDDO 
          DO IFRE = 1,NFRE
            DO IANG = 1,NANG
              DO ISPEC = 0,NSPEC
                KCOUNT=KCOUNT+1
                SPEC(ISPEC,IANG,IFRE,2)=ZCOMBUF(KCOUNT)
              ENDDO 
            ENDDO 
          ENDDO 
        ENDIF

        DEALLOCATE(ZCOMBUF)
      ENDIF


5000  CONTINUE
      RETURN

      CONTAINS

      INTEGER FUNCTION BUFFER_SIZE_BCSARIN (NSPEC,NANG,NFRE)
        IMPLICIT NONE

        INTEGER :: NSPEC,NANG,NFRE
        BUFFER_SIZE_BCSARIN = 1+5*NSPEC+(NSPEC+1)*NANG*NFRE
      END FUNCTION BUFFER_SIZE_BCSARIN
      END SUBROUTINE MPBCASTSARIN
