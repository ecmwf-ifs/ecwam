      SUBROUTINE MPEXCHNG(FL,NSTART,NEND,KLENBOT,KLENTOP)
C
C****  *MPEXCHNG* - EXHANGES MESSAGE BETWEEN THE PROCESS IRANK 
C****               AND THE PREVIOUS AND NEXT ONE. 
C
!     J. BIDLOT    ECMWF   MARCH 1996  MESSAGE PASSING
C
C     PURPOSE.
C     --------
C     EXHANGE MESSAGE BETWEEN ONE PROCESS AND THE PREVIOUS AND NEXT PROCESS
C     USING MESSAGE PASSING PROTOCOL. FOR EACH PROCESS THERE ARE 2 SEND
C     AND 2 RECEIVE TO BE PERFORMED, EXCEPT FOR THE FIRST AND LAST PROCESS
C     WHICH ONLY NEED 1 OF EACH.
C
C*    INTERFACE.
C     ----------
C     CALL *MPEXCHNG*(FL,NSTART,NEND,KLENTOP,KLENBOT)
C
C      *FL*        FREQUENCY SPECTRUM
C      *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
C      *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
C      *KLENTOP*   LENGTH OF THE MESSAGE SENT BY THE NEXT
C                  SUB GRID DOMAIN (IN GRID POINT UNIT)
C      *KLENBOT*   LENGTH OF THE MESSAGE SENT BY THE PREVIOUS
C                  SUB GRID DOMAIN (IN GRID POINT UNIT)

C     METHOD.
C     -------
C
C     USING NON-BLOCKING SEND, YOU CAN FIRST SEND ALL MESSAGES THEN
C     CLLECT THEM ACCORDINGLY WITH A BLOCKING RECEIVE.
C
C     NB BY DEFAULT THE MPE_SEND IS NON BLOCKING BUT FOR CONSITANCY WE
C     WILL SET THE APPROPRIATE PARAMETER FOR NON BLOCKING SEND
C
C     EXTERNALS.
C     ----------
C     MPE PACKAGE :
C         MPE_SEND
C         MPE_RECV
C         MPEI_ABORT 
C
C     REFERENCES.
C     -----------
C     CHAPTER 4 OF
C     USING MPI, PORTABLE PARALLEL PROGRAMMING WITH THE MESSAGE PASSING 
C     INTERFACE. W.CROPP, E LUSK, A SKJELLUM. MIT PRESS 1995 
C
C -------------------------------------------------------------------
C
#include "param.h"
C
#include "comtest.h"
#include "comubuf.h"
#include "txtubuf.h"
C
#include "commpp.h"
C
#include "txtmpp.h"

!----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND,KLENTOP,KLENBOT
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL

!----------------------------------------------------------------------
      logical ldwait,ldflag
 
      REAL,ALLOCATABLE :: zcombuf(:)

c     0.1 find the maximum length of the message buffer
c         ---------------------------------------------
      if(nproc.eq.1) then
        nbufmax=0
      else if (irank.eq.1) then
        nbufmax=max(klenbot(irank+1),klenbot(irank),klentop(irank))
      else if (irank.eq.nproc) then
        nbufmax=max(klentop(irank-1),klenbot(irank),klentop(irank))
      else
        nbufmax=max(klenbot(irank+1),klenbot(irank),klentop(irank),
     &              klentop(irank-1))
      endif
c     0.2 allocate the message buffer space
c         ---------------------------------
      nbufmax=nbufmax*nang*nfre
      allocate(zcombuf(nbufmax))

c     0.3 find index of the the start or end point of the different messages
c         ------------------------------------------------------------------
      ist=nend(irank)-klenbot(irank+1)+1
      isb=nstart(irank)-klenbot(irank)
      ieb=nstart(irank)+klentop(irank-1)-1
      iet=nend(irank)+klentop(irank)
 


c     1.1.1 send to next process
c           --------------------
      if(nnext.ne.0) then

c       fill the send buffer
        kcount=0
        do M=1,NFRE
           do K=1,NANG
              do IJ=ist,ist+klenbot(nnext)-1
                 kcount=kcount+1
                 zcombuf(kcount)=FL(IJ,K,M)
              end do
           end do
        end do

c       send the buffer to next process
        call MPE_SEND(zcombuf,kcount,2,nnext,ktag,0,0,1,ierr)
        if(ierr.lt.0) call MPEI_ABORT('ERROR in MPE_SEND ')
      endif

c     1.2.1 send to previous process
c           -------------------------
      if(nprevious.ne.0) then

c       fill the send buffer
        kcount=0
        do M=1,NFRE
           do K=1,NANG
              do IJ=nstart(irank),nstart(irank)+klentop(nprevious)-1
                 kcount=kcount+1
                 zcombuf(kcount)=FL(IJ,K,M)
              end do
           end do
        end do

c       send the buffer to previous process
        call MPE_SEND(zcombuf,kcount,2,nprevious,ktag+1,0,0,1,ierr)
        if(ierr.lt.0) call MPEI_ABORT('ERROR in MPE_SEND ')

      endif
c
c     1.1.2 receive from previous process
c           -----------------------------
      if(nprevious.ne.0) then

c       receive the buffer from previous process
        kcount=klenbot(irank)*NANG*NFRE
        call MPE_RECV(zcombuf,kcount,2,nprevious,ktag,0,0,0,
     $                krcount,krfrom,krtag,ierr)
        if(ierr.lt.0) call MPEI_ABORT('ERROR in MPE_RECV ')

c       decode the received buffer
        kcount=0
        do M=1,NFRE
           do K=1,NANG
              do IJ=isb,isb+klenbot(irank)-1
                 kcount=kcount+1
                 FL(IJ,K,M)=zcombuf(kcount)
              end do
           end do
        end do

      endif

c     1.1.2 receive from next process
c           -------------------------
      if(nnext.ne.0) then

c       receive the buffer from next process
        kcount=klentop(irank)*NANG*NFRE
        call MPE_RECV(zcombuf,kcount,2,nnext,ktag+1,0,0,0,
     $                krcount,krfrom,krtag,ierr)
        if(ierr.lt.0) call MPEI_ABORT('ERROR in MPE_RECV ')

c       decode the received buffer
        kcount=0
        do M=1,NFRE
           do K=1,NANG
              do IJ=nend(irank)+1,nend(irank)+klentop(irank)
                 kcount=kcount+1
                 FL(IJ,K,M)=zcombuf(kcount)
              end do
           end do
        end do

      endif
c
c
      ktag=ktag+2

      deallocate(zcombuf)
c
      return
      end

