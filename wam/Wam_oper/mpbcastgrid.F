      SUBROUTINE MPBCASTGRID(IU06,ISEND,ITAG,NSTART,NEND,IREFRA)
C
C ----------------------------------------------------------------------
C
C**** *MPBCASTGRID* - BROADCAST THE CONTENT OF THE GRID FILE READ BY READPRE
C                     TO THE OTHER PE'S
C
C     J. BIDLOT    ECMWF   OCTOBER 1997
C
C     PURPOSE.
C     --------
C     BROADCAST THE CONTENT OF THE GRID FILE READ BY READPRE
C     TO THE OTHER PE'S
C*    INTERFACE.
C     ----------
C
C     CALL *MPBCASTGRID*(IU06,ISEND,ITAG,NSTART,NEND,IREFRA)
C
C     *IU06*      UNIT FOR PRINTER MESSAGES.
C     *ISEND*     RANK OF THE PROCESS ONTO WHICH FIELD IS COLLECTED
C     *ITAG*      TAG ASSOCIATED WITH AS A PARTICULAR CALL TO SUBROUTINE
C                 THIS IS NECESSARY TO DIFFERENTIATE THE DIFFERENT CALLS
C     *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C     *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C     *IREFRA*    REFRACTION OPTION.
C
C     METHOD.
C     -------
C     MPE_BROADCAST FROM PROCESSOR CORRESPONDING TO ISEND TO
C     ALL OTHER PROCESSORS.
C
C     EXTERNALS.
C     ----------
C     MPE PACKAGE :
C         MPE_BROADCAST
C         MPE_RECV
C         MPEI_ABORT
C
C     REFERENCES.
C     -----------
C         NONE
C -------------------------------------------------------------------
C
#include "param.h"
C
#include "comcoup.h"
C
#include "comcurr.h"
C
#include "comfred.h"
C
#include "comgrid.h"
C
#include "comindn.h"
C
#include "commap.h"
C
#include "commpp.h"
C
#include "comcout.h"
C
#include "comshal.h"
C
#include "comtabl.h"
C
C----------------------------------------------------------------------
C     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS
      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
C----------------------------------------------------------------------

      INTEGER,ALLOCATABLE :: ICOMBUF(:)
      REAL,ALLOCATABLE :: ZCOMBUF(:)

C
      MIC=6+5*NBLO+NGY+(2*NBLO+1)*NIBLO+4*(NFRE+4)+8*NANG+2*MOUTP

      MZC=25+(4+3*NDEPTH)*NFRE+5*(NFRE+4)+3*NANG+4*NGY+KFRH+JPLEVC+
     $    JPLEVT*(1+JUMAX)*(1+ITAUMAX)+(1+IUSTAR)*(1+IALPHA)+NBLO*NIBLO
      IF(IREFRA.EQ.2) MZC=MZC+2*NBLC*(1+NIBLC)

      IF((ISEND.EQ.0).OR.(NPROC.EQ.1)) THEN
         RETURN
C     1.1 SEND TO ALL PROCESSORS OTHER THAN ISEND
C         ------------------------------------------------
      ELSE IF(IRANK.EQ.ISEND) THEN
         ALLOCATE(ICOMBUF(MIC))
         ALLOCATE(ZCOMBUF(MZC))

C  ENCODE MESSAGE BUFFERS
C
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
         ICOMBUF(IKCOUNT)=IGL
         DO IG=1,NBLO
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJS(IG)
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJL2(IG)
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJLS(IG)
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJL(IG)
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJLT(IG)
         ENDDO

         DO IG=1,NBLO
            DO IJ=1,NIBLO
               IKCOUNT=IKCOUNT+1
               ICOMBUF(IKCOUNT)=IXLG(IJ,IG)
            ENDDO
         ENDDO
         DO IG=1,NBLO
            DO IJ=1,NIBLO
               IKCOUNT=IKCOUNT+1
               ICOMBUF(IKCOUNT)=KXLT(IJ,IG)
            ENDDO
         ENDDO
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

         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKP(M)
         ENDDO
         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKP1(M)
         ENDDO
         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKM(M)
         ENDDO
         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IKM1(M)
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K1W(K,IC)
           ENDDO
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K2W(K,IC)
           ENDDO
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K11W(K,IC)
           ENDDO
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              ICOMBUF(IKCOUNT)=K21W(K,IC)
           ENDDO
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=AF11(M)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAP(M)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAP1(M)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAM(M)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FKLAM1(M)
         ENDDO
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=ACL1
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=ACL2
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=CL11
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=CL21
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=DAL1
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=DAL2
         DO IC=1,KFRH
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=FRH(IC)
         ENDDO


         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=BETAMAX
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=ZALP
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=ALPHA
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=XKAPPA
         DO IC=1,JPLEVC
            KCOUNT=KCOUNT+1
            ZCOMBUF(KCOUNT)=XNLEV(IC)
         ENDDO


         DO ILEV=1,JPLEVT
            DO J=0,JUMAX
               DO I=0,ITAUMAX
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=TAUT(I,J,ILEV)
               ENDDO
            ENDDO
         ENDDO
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=DELTAUW
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=DELU
         DO J=0,IUSTAR
            DO I=0,IALPHA
               KCOUNT=KCOUNT+1
               ZCOMBUF(KCOUNT)=TAUHFT(I,J)
            ENDDO
         ENDDO
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=DELUST
         KCOUNT=KCOUNT+1
         ZCOMBUF(KCOUNT)=DELALP

         IKCOUNT=IKCOUNT+1
         ICOMBUF(IKCOUNT)=NGOUT
         DO NGOU=1,MOUTP
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IGAR(NGOU)
         ENDDO
         DO NGOU=1,MOUTP
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=IJAR(NGOU)
         ENDDO

         DO IG=1,NBLO
            DO IJ=1,NIBLO
               KCOUNT=KCOUNT+1
               ZCOMBUF(KCOUNT)=DEPTH(IJ,IG)
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
         DO IJ=1,NIBLO
            IKCOUNT=IKCOUNT+1
            ICOMBUF(IKCOUNT)=INDEP(IJ)
         ENDDO

         IF (IREFRA.EQ.2) THEN
            DO IG=1,NBLC
               DO IJ=0,NIBLC
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=U(IJ,IG)
               ENDDO
            ENDDO
            DO IG=1,NBLC
               DO IJ=0,NIBLC
                  KCOUNT=KCOUNT+1
                  ZCOMBUF(KCOUNT)=V(IJ,IG)
               ENDDO
            ENDDO
         ENDIF

         IF(IKCOUNT.NE.MIC) THEN
             WRITE (IU06,*) '**************************'
             WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
             WRITE (IU06,*) '* IKCOUNT NE MIC PRIOR   *'
             WRITE (IU06,*) '* CALL TO MPE_BROADCAST  *'
             WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
             WRITE (IU06,*) '* MIC =',MIC
             WRITE (IU06,*) '**************************'
             CALL ABORT1
         ENDIF 
         IF(KCOUNT.NE.MZC) THEN
             WRITE (IU06,*) '**************************'
             WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
             WRITE (IU06,*) '* KCOUNT NE MZC PRIOR    *'
             WRITE (IU06,*) '* CALL TO MPE_BROADCAST  *'
             WRITE (IU06,*) '* KCOUNT =',KCOUNT
             WRITE (IU06,*) '* MZC =',MZC
             WRITE (IU06,*) '**************************'
             CALL ABORT1
         ENDIF 
C
         CALL MPE_BROADCAST(ICOMBUF,MIC,1,IRANK,ITAG,0,0,0,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                 ('MPE_BROADCAST ERROR 1 in MPBCASTGRID' )
         ITAG=ITAG+1
         CALL MPE_BROADCAST(ZCOMBUF,MZC,2,IRANK,ITAG,0,0,0,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                 ('MPE_BROADCAST ERROR 2 in MPBCASTGRID' )
         ITAG=ITAG+1
C
C
C     1.2 RECEIVE ON ALL OTHER PE THAN ISEND
C         ----------------------------------
      ELSE
         ALLOCATE(ICOMBUF(MIC))
         ALLOCATE(ZCOMBUF(MZC))

         CALL MPE_RECV(ICOMBUF,MIC,1,ISEND,ITAG,0,0,0,
     $                 KRCOUNT,KRFROM,KRTAG,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_RECV ERROR 1 in MPBCASTGRID' )
         IF(KRCOUNT.NE.MIC) CALL MPEI_ABORT
     1   ('MPE_RECV ERROR 1 in MPBCASTGRID:MISMATCHED MSG LENGTH')
         IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1   ('MPE_RECV ERROR 1 in MPBCASTGRID MISMATCHED TAGS' )

         ITAG=ITAG+1

         CALL MPE_RECV(ZCOMBUF,MZC,2,ISEND,ITAG,0,0,0,
     $                 KRCOUNT,KRFROM,KRTAG,IERR)
         IF(IERR.LT.0) CALL MPEI_ABORT
     1                    ('MPE_RECV ERROR 2 in MPBCASTGRID' )
         IF(KRCOUNT.NE.MZC) CALL MPEI_ABORT
     1   ('MPE_RECV ERROR 2 in MPBCASTGRID:MISMATCHED MSG LENGTH')
         IF(KRTAG.NE.ITAG) CALL MPEI_ABORT
     1   ('MPE_RECV ERROR 2 in MPBCASTGRID MISMATCHED TAGS' )

         ITAG=ITAG+1

C DECODE MESSAGE BUFFERS
C
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
         IGL=ICOMBUF(IKCOUNT)
         DO IG=1,NBLO
            IKCOUNT=IKCOUNT+1
            IJS(IG)=ICOMBUF(IKCOUNT)
            IKCOUNT=IKCOUNT+1
            IJL2(IG)=ICOMBUF(IKCOUNT)
            IKCOUNT=IKCOUNT+1
            IJLS(IG)=ICOMBUF(IKCOUNT)
            IKCOUNT=IKCOUNT+1
            IJL(IG)=ICOMBUF(IKCOUNT)
            IKCOUNT=IKCOUNT+1
            IJLT(IG)=ICOMBUF(IKCOUNT)
         ENDDO

         DO IG=1,NBLO
            DO IJ=1,NIBLO
               IKCOUNT=IKCOUNT+1
               IXLG(IJ,IG)=ICOMBUF(IKCOUNT)
            ENDDO
         ENDDO
         DO IG=1,NBLO
            DO IJ=1,NIBLO
               IKCOUNT=IKCOUNT+1
               KXLT(IJ,IG)=ICOMBUF(IKCOUNT)
            ENDDO
         ENDDO
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

         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            IKP(M)=ICOMBUF(IKCOUNT)
         ENDDO
         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            IKP1(M)=ICOMBUF(IKCOUNT)
         ENDDO
         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            IKM(M)=ICOMBUF(IKCOUNT)
         ENDDO
         DO M=1,NFRE+4
            IKCOUNT=IKCOUNT+1
            IKM1(M)=ICOMBUF(IKCOUNT)
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K1W(K,IC)=ICOMBUF(IKCOUNT)
           ENDDO
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K2W(K,IC)=ICOMBUF(IKCOUNT)
           ENDDO
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K11W(K,IC)=ICOMBUF(IKCOUNT)
           ENDDO
         ENDDO
         DO IC=1,2
           DO K=1,NANG
              IKCOUNT=IKCOUNT+1
              K21W(K,IC)=ICOMBUF(IKCOUNT)
           ENDDO
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            AF11(M)=ZCOMBUF(KCOUNT)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            FKLAP(M)=ZCOMBUF(KCOUNT)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            FKLAP1(M)=ZCOMBUF(KCOUNT)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            FKLAM(M)=ZCOMBUF(KCOUNT)
         ENDDO
         DO M=1,NFRE+4
            KCOUNT=KCOUNT+1
            FKLAM1(M)=ZCOMBUF(KCOUNT)
         ENDDO
         KCOUNT=KCOUNT+1
         ACL1=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         ACL2=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         CL11=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         CL21=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         DAL1=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         DAL2=ZCOMBUF(KCOUNT)
         DO IC=1,KFRH
            KCOUNT=KCOUNT+1
            FRH(IC)=ZCOMBUF(KCOUNT)
         ENDDO


         KCOUNT=KCOUNT+1
         BETAMAX=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         ZALP=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         ALPHA=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         XKAPPA=ZCOMBUF(KCOUNT)
         DO IC=1,JPLEVC
            KCOUNT=KCOUNT+1
            XNLEV(IC)=ZCOMBUF(KCOUNT)
         ENDDO


         DO ILEV=1,JPLEVT
            DO J=0,JUMAX
               DO I=0,ITAUMAX
                  KCOUNT=KCOUNT+1
                  TAUT(I,J,ILEV)=ZCOMBUF(KCOUNT)
               ENDDO
            ENDDO
         ENDDO
         KCOUNT=KCOUNT+1
         DELTAUW=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         DELU=ZCOMBUF(KCOUNT)
         DO J=0,IUSTAR
            DO I=0,IALPHA
               KCOUNT=KCOUNT+1
               TAUHFT(I,J)=ZCOMBUF(KCOUNT)
            ENDDO
         ENDDO
         KCOUNT=KCOUNT+1
         DELUST=ZCOMBUF(KCOUNT)
         KCOUNT=KCOUNT+1
         DELALP=ZCOMBUF(KCOUNT)

         IKCOUNT=IKCOUNT+1
         NGOUT=ICOMBUF(IKCOUNT)
         DO NGOU=1,MOUTP
            IKCOUNT=IKCOUNT+1
            IGAR(NGOU)=ICOMBUF(IKCOUNT)
         ENDDO
         DO NGOU=1,MOUTP
            IKCOUNT=IKCOUNT+1
            IJAR(NGOU)=ICOMBUF(IKCOUNT)
         ENDDO

         DO IG=1,NBLO
            DO IJ=1,NIBLO
               KCOUNT=KCOUNT+1
               DEPTH(IJ,IG)=ZCOMBUF(KCOUNT)
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
         DO IJ=1,NIBLO
            IKCOUNT=IKCOUNT+1
            INDEP(IJ)=ICOMBUF(IKCOUNT)
         ENDDO

         IF (IREFRA.EQ.2) THEN
            DO IG=1,NBLC
               DO IJ=0,NIBLC
                  KCOUNT=KCOUNT+1
                  U(IJ,IG)=ZCOMBUF(KCOUNT)
               ENDDO
            ENDDO
            DO IG=1,NBLC
               DO IJ=0,NIBLC
                  KCOUNT=KCOUNT+1
                  V(IJ,IG)=ZCOMBUF(KCOUNT)
               ENDDO
            ENDDO
         ENDIF
         IF(IKCOUNT.NE.MIC) THEN
             WRITE (IU06,*) '**************************'
             WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
             WRITE (IU06,*) '* IKCOUNT NE MIC AFTER   *'
             WRITE (IU06,*) '* CALL TO MPE_RECV       *'
             WRITE (IU06,*) '* IKCOUNT =',IKCOUNT
             WRITE (IU06,*) '* MIC =',MIC
             WRITE (IU06,*) '**************************'
             CALL ABORT1
         ENDIF 
         IF(KCOUNT.NE.MZC) THEN
             WRITE (IU06,*) '**************************'
             WRITE (IU06,*) '* ERROR IN MPBCASTGRID   *'
             WRITE (IU06,*) '* KCOUNT NE MZC AFTER    *'
             WRITE (IU06,*) '* CALL TO MPE_RECV       *'
             WRITE (IU06,*) '* KCOUNT =',KCOUNT
             WRITE (IU06,*) '* MZC =',MZC
             WRITE (IU06,*) '**************************'
             CALL ABORT1
         ENDIF 
      ENDIF

      DEALLOCATE(ICOMBUF,ZCOMBUF)
      RETURN
      END
