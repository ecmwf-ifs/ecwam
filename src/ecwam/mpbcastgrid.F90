! (C) Copyright 1989- ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!

SUBROUTINE MPBCASTGRID(IU06, ISEND, ITAG)

! ----------------------------------------------------------------------
!**** *MPBCASTGRID* - BROADCAST THE CONTENT OF THE GRID FILE READ BY
!                     READPRE TO THE OTHER PE'S. IT WILL ALSO ALLOCATE
!                     THE NECESSARY ARRAYS ON THE RECEIVING PE'S AS THEY
!                     WERE NOT ON THOSE PE'S SINCE READPRE WAS NOT
!                     CALLED FOR THEM.

!     J. BIDLOT    ECMWF   OCTOBER 1997

!     PURPOSE.
!     --------
!     BROADCAST THE CONTENT OF THE GRID FILE READ BY READPRE
!     TO THE OTHER PE'S
!*    INTERFACE.
!     ----------

!     CALL *MPBCASTGRID*(IU06,ISEND,ITAG)

!     *IU06*      UNIT FOR PRINTER MESSAGES.
!     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED
!     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
!                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS

!     METHOD.
!     -------
!     MPL_BROADCAST FROM PROCESSOR CORRESPONDING TO ISEND TO
!     ALL OTHER PROCESSORS.

!     EXTERNALS.
!     ----------
!     MPL PACKAGE :
!         MPL_BROADCAST

!     REFERENCES.
!     -----------
!         NONE
! -------------------------------------------------------------------

      USE PARKIND_WAVE, ONLY : JWIM, JWRB, JWRU

      USE YOWFRED  , ONLY : FR       ,DFIM     ,GOM      ,C        ,    &
     &            DELTH    ,DELTR    ,TH       ,COSTH    ,SINTH
      USE YOWGRID  , ONLY : DELPHI   ,DELLAM   ,SINPH    ,COSPH    ,    &
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : NX       ,NY       ,              &
     &            IPER     ,IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,    &
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,NLONRGG  ,    &
     &            IQGAUSS
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NPRECR   ,NPRECI
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NFRE_RED ,              &
     &            NGX      ,NGY      ,                                  &
     &            NOVER    ,NIBL1    ,CLDOMAIN
      USE YOWSHAL  , ONLY : NDEPTH   ,BATHY    ,DEPTHA   ,DEPTHD   ,    &
     &            TCGOND   ,TFAK     ,TSIHKD   ,TFAC_ST

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
      USE MPL_MODULE, ONLY : MPL_BROADCAST

!----------------------------------------------------------------------

      IMPLICIT NONE
#include "abort1.intfb.h"

      INTEGER(KIND=JWIM), INTENT(IN) :: IU06, ISEND
      INTEGER(KIND=JWIM), INTENT(INOUT) :: ITAG
      INTEGER(KIND=JWIM), PARAMETER :: MFIRST=11
      INTEGER(KIND=JWIM) :: I, J, IJ, K, K1, K2, M, M1, M2, IC, L, KDEPTH, NGOU 
      INTEGER(KIND=JWIM) :: IKCOUNT, KCOUNT
      INTEGER(KIND=JWIM) :: MIC, MZC 
      INTEGER(KIND=JWIM),ALLOCATABLE :: ICOMBUF(:)

      REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
      REAL(KIND=JWRB),ALLOCATABLE :: ZCOMBUF(:)

!----------------------------------------------------------------------

      IF (LHOOK) CALL DR_HOOK('MPBCASTGRID',0,ZHOOK_HANDLE)

      IF (ISEND == 0 .OR. NPROC == 1) THEN
         WRITE (IU06,*) ''
!     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
!         ------------------------------------------------
      ELSE

!       BUFFER SIZE MESSAGE AND THE FEW DIMENSIONS NEEDED TO 
!       ALLOCATE ALL ARRAYS
        ALLOCATE(ICOMBUF(MFIRST))

        IF (IRANK == ISEND) THEN
          IKCOUNT=0
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NANG
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NFRE
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NFRE_RED
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGX
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NGY
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NOVER
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NIBL1
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=ICHAR(CLDOMAIN)
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NDEPTH
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NANGH
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NMAX
          IF (IKCOUNT /= MFIRST) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* IKCOUNT .NE. MFIRST !!!*' 
            WRITE (IU06,*) '* ON IRANK = ',IRANK
            WRITE (IU06,*) '* IKCOUNT = ',IKCOUNT
            WRITE (IU06,*) '* MFIRST  = ',MFIRST
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF,KROOT=ISEND,KTAG=ITAG,CDSTRING='MPBCASTGRID:')
        ITAG=ITAG+1

        IF (IRANK /= ISEND) THEN
          IKCOUNT=0
          IKCOUNT=IKCOUNT+1
          NANG=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NFRE=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NFRE_RED=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NGX=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NGY=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NOVER=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NIBL1=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          CLDOMAIN=CHAR(ICOMBUF(IKCOUNT))
          IKCOUNT=IKCOUNT+1
          NDEPTH=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NANGH=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NMAX=ICOMBUF(IKCOUNT)
          IF (IKCOUNT /= MFIRST) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* IKCOUNT .NE. MFIRST !!!*' 
            WRITE (IU06,*) '* ON IRANK = ',IRANK
            WRITE (IU06,*) '* IKCOUNT = ',IKCOUNT
            WRITE (IU06,*) '* MFIRST  - ',MFIRST
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF
        ENDIF
        DEALLOCATE(ICOMBUF)

        MIC=7+NGY
        MZC=11+(4+4*NDEPTH)*NFRE+3*NANG+4*NGY+NGX*NGY

!       ENCODE MAIN MESSAGE BUFFERS (ON PE=ISEND) AND
!       ALLOCATE ALL ARRAYS NEEDED TO KEEP THE BUFFERS ON THE OTHER PE'S

        ALLOCATE(ICOMBUF(MIC))
        ALLOCATE(ZCOMBUF(MZC))

        IF (IRANK /= ISEND) THEN

          IF (.NOT.ALLOCATED(FR)) ALLOCATE(FR(NFRE))
          IF (.NOT.ALLOCATED(DFIM)) ALLOCATE(DFIM(NFRE))
          IF (.NOT.ALLOCATED(GOM)) ALLOCATE(GOM(NFRE))
          IF (.NOT.ALLOCATED(C)) ALLOCATE(C(NFRE))
          IF (.NOT.ALLOCATED(TH)) ALLOCATE(TH(NANG))
          IF (.NOT.ALLOCATED(COSTH)) ALLOCATE(COSTH(NANG))
          IF (.NOT.ALLOCATED(SINTH)) ALLOCATE(SINTH(NANG))
          IF (.NOT.ALLOCATED(DELLAM)) ALLOCATE(DELLAM(NGY))
          IF (.NOT.ALLOCATED(NLONRGG)) ALLOCATE(NLONRGG(NGY))
          IF (.NOT.ALLOCATED(SINPH)) ALLOCATE(SINPH(NGY))
          IF (.NOT.ALLOCATED(COSPH)) ALLOCATE(COSPH(NGY))
          IF (.NOT.ALLOCATED(ZDELLO)) ALLOCATE(ZDELLO(NGY))

          IF (ALLOCATED(BATHY)) DEALLOCATE(BATHY)
          ALLOCATE(BATHY(NGX,NGY))

          IF (.NOT.ALLOCATED(TCGOND)) ALLOCATE(TCGOND(NDEPTH,NFRE))
          IF (.NOT.ALLOCATED(TFAK)) ALLOCATE(TFAK(NDEPTH,NFRE))
          IF (.NOT.ALLOCATED(TSIHKD)) ALLOCATE(TSIHKD(NDEPTH,NFRE))
          IF (.NOT.ALLOCATED(TFAC_ST)) ALLOCATE(TFAC_ST(NDEPTH,NFRE))

        ELSE 
          KCOUNT=0
          IKCOUNT=0
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FR(M)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DFIM(M)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=GOM(M)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=C(M)
          ENDDO
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DELTH
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DELTR
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=TH(K)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=COSTH(K)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=SINTH(K)
          ENDDO

          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DELPHI
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=DELLAM(J)
          ENDDO
          DO J=1,NGY
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=NLONRGG(J)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=SINPH(J)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=COSPH(J)
          ENDDO
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IJS
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IJL

          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NX
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=NY
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IPER
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOWEP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOSOP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMOEAP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=AMONOP
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XDELLA
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=XDELLO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=ZDELLO(J)
          ENDDO
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IRGG
          IKCOUNT=IKCOUNT+1
          ICOMBUF(IKCOUNT)=IQGAUSS

          DO K=1,NGY
            DO I=1,NGX
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=BATHY(I,K)
            ENDDO
          ENDDO

          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DEPTHA
          KCOUNT=KCOUNT+1
          ZCOMBUF(KCOUNT)=DEPTHD
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TCGOND(KDEPTH,M)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TFAK(KDEPTH,M)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TSIHKD(KDEPTH,M)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              ZCOMBUF(KCOUNT)=TFAC_ST(KDEPTH,M)
            ENDDO
          ENDDO

          IF (IKCOUNT /= MIC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT NE MIC PRIOR   *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
            WRITE (IU06,*) '* MIC =',MIC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
          IF (KCOUNT /= MZC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* KCOUNT NE MZC PRIOR    *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* KCOUNT =',KCOUNT
            WRITE (IU06,*) '* MZC =',MZC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
        ENDIF

        CALL MPL_BROADCAST(ICOMBUF,KROOT=ISEND,KTAG=ITAG,               &
     &                     CDSTRING='MPBCASTGRID 1:')
        ITAG=ITAG+1

        CALL MPL_BROADCAST(ZCOMBUF,KROOT=ISEND,KTAG=ITAG,               &
     &                     CDSTRING='MPBCASTGRID 2:')
        ITAG=ITAG+1

        IF (IRANK /= ISEND) THEN
          KCOUNT=0
          IKCOUNT=0
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            FR(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            DFIM(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            GOM(M)=ZCOMBUF(KCOUNT)
          ENDDO
          DO M=1,NFRE
            KCOUNT=KCOUNT+1
            C(M)=ZCOMBUF(KCOUNT)
          ENDDO
          KCOUNT=KCOUNT+1
          DELTH=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          DELTR=ZCOMBUF(KCOUNT)
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            TH(K)=ZCOMBUF(KCOUNT)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            COSTH(K)=ZCOMBUF(KCOUNT)
          ENDDO
          DO K=1,NANG
            KCOUNT=KCOUNT+1
            SINTH(K)=ZCOMBUF(KCOUNT)
          ENDDO

          KCOUNT=KCOUNT+1
          DELPHI=ZCOMBUF(KCOUNT)
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            DELLAM(J)=ZCOMBUF(KCOUNT)
          ENDDO
          DO J=1,NGY
            IKCOUNT=IKCOUNT+1
            NLONRGG(J)=ICOMBUF(IKCOUNT)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            SINPH(J)=ZCOMBUF(KCOUNT)
          ENDDO
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            COSPH(J)=ZCOMBUF(KCOUNT)
          ENDDO
          IKCOUNT=IKCOUNT+1
          IJS=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IJL=ICOMBUF(IKCOUNT)

          IKCOUNT=IKCOUNT+1
          NX=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          NY=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IPER=ICOMBUF(IKCOUNT)
          KCOUNT=KCOUNT+1
          AMOWEP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMOSOP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMOEAP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          AMONOP=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XDELLA=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          XDELLO=ZCOMBUF(KCOUNT)
          DO J=1,NGY
            KCOUNT=KCOUNT+1
            ZDELLO(J)=ZCOMBUF(KCOUNT)
          ENDDO
          IKCOUNT=IKCOUNT+1
          IRGG=ICOMBUF(IKCOUNT)
          IKCOUNT=IKCOUNT+1
          IQGAUSS=ICOMBUF(IKCOUNT)

          DO K=1,NGY
            DO I=1,NGX
              KCOUNT=KCOUNT+1
              BATHY(I,K)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO

          KCOUNT=KCOUNT+1
          DEPTHA=ZCOMBUF(KCOUNT)
          KCOUNT=KCOUNT+1
          DEPTHD=ZCOMBUF(KCOUNT)
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TCGOND(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TFAK(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TSIHKD(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO
          DO M=1,NFRE
            DO KDEPTH=1,NDEPTH
              KCOUNT=KCOUNT+1
              TFAC_ST(KDEPTH,M)=ZCOMBUF(KCOUNT)
            ENDDO
          ENDDO

          IF (IKCOUNT /= MIC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* IKCOUNT NE MIC AFTER   *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
            WRITE (IU06,*) '* MIC =',MIC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
          IF (KCOUNT /= MZC) THEN
            WRITE (IU06,*) '**************************'
            WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
            WRITE (IU06,*) '* KCOUNT NE MZC AFTER    *'
            WRITE (IU06,*) '* CALL TO MPL_BROADCAST  *'
            WRITE (IU06,*) '* KCOUNT =',KCOUNT
            WRITE (IU06,*) '* MZC =',MZC
            WRITE (IU06,*) '**************************'
            CALL ABORT1
          ENDIF 
        ENDIF

        DEALLOCATE(ICOMBUF,ZCOMBUF)

      ENDIF

      IF (LHOOK) CALL DR_HOOK('MPBCASTGRID',1,ZHOOK_HANDLE)

END SUBROUTINE MPBCASTGRID
