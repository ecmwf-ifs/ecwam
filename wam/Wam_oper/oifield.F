      SUBROUTINE OIFIELD (XMO, XOI, SIGMOD, DIST, DISTMAX) 

!**** *OIFIELD* - OPTIMUM INTERPOLATION.                                

!     P.LIONELLO     ECMWF       APRIL 1990                             
!     J. BIDLOT      ECMWF       SEPTEMBER 96  MESSAGE PASSING

!     PURPOSE.                                                          
!     --------                                                          

!       TO PRODUCE A MAP OF THE FIELD X, MERGING MEASUREMENT            
!       AND MODEL , BY OPTIMUM INTERPOLATION. THE ARRAY XOI
!       AT THE END OF THE SUBROUTINE CONTAINS THE VALUES                
!       TO BE USED TO ANALYSE THE SPECTRA, HAVING NEGATIVE RETURN       
!       CODES WHERE O.I. PRODUCED NO RESULTS.                           
!       THE FIRST GUESS FIELD XMO IS NOT MODIFIED                       
!       IN THIS SUBROUTINE.                                             

!**   INTERFACE.                                                        
!     ----------                                                        

!       *CALL* *OIFIELD (XMO, SIGMOD, DIST, DISTMAX)*      

!         *XMO*     REAL     FIELD FROM MODEL (FIRST GUESS)
!         *XOI*     REAL     FIELD FROM O.I. ON OUTPUT.
!         *SIGMOD*  REAL     MODEL ERROR ESTIMATE.                             
!         *DIST*    REAL     CORRELATION DISTANCE
!         *DISTMAX* REAL     MAXIMIN SPREADING DISTANCE.

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : NIJALT   ,NALTDT   ,IJALT    ,INTLMAX  ,
     &            NOBSPE   ,ALTDATA
      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,NLONRGG  ,IGL      ,
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,XDELLA   ,
     &            ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      INTEGER :: NOBS(NSTART(IRANK):NEND(IRANK))
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IOBS4IJ

      REAL :: XMO(NIBLO)
      REAL :: XOI(NINF:NSUP)
      REAL :: XNEW(NSTART(IRANK):NEND(IRANK))
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KALTMIN, KALTMAX
      REAL, ALLOCATABLE, DIMENSION(:) :: XLONOBS, SIGRATIO2, DIFFALTFG
      REAL, ALLOCATABLE, DIMENSION(:) :: V
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XM, D, P 
      REAL, ALLOCATABLE, DIMENSION(:,:) :: W

! ----------------------------------------------------------------------


!*    1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.
!     -----------------------------------------------------

      LMAX = NINT(DEG*DISTMAX/XDELLA)

! ??? should be adapted if more than one satellite !!!
      NDIM2 = 5*(2*LMAX+1)

!     THE MEASUREMENTS ARE STORED IN THE O.I. FIELD.                  
!     NEGATIVE VALUES WILL REMAIN WHERE O.I. WILL PRODUCE NO RESULTS  
      XOI = -1.

!       LOOP OVER GRID POINTS.                                          
!       ----------------------                                          

!     FIND THE INDEX OF ALL OBSERVATIONS INFLUENCING A GIVEN GRID POINT

      IF(NOBSPE(IRANK).GT.0) THEN
        ALLOCATE(KALTMIN(NOBSPE(IRANK)))
        ALLOCATE(KALTMAX(NOBSPE(IRANK)))
        ALLOCATE(XLONOBS(NOBSPE(IRANK)))
        ALLOCATE(SIGRATIO2(NOBSPE(IRANK)))
        ALLOCATE(DIFFALTFG(NOBSPE(IRANK)))
      ENDIF

1000  CONTINUE
      ALLOCATE(IOBS4IJ(NSTART(IRANK):NEND(IRANK),NDIM2))
      ALLOCATE(W(NSTART(IRANK):NEND(IRANK),NDIM2))

      IOBS4IJ=0
      NOBS=0
      W=1.
      NOBSMAX=0
      IG = 1
      DO IOBS=1,NOBSPE(IRANK)
        IALT = IXLG(IJALT(IOBS,1),IG)
        KALT = KXLT(IJALT(IOBS,1),IG)
        KALTMIN(IOBS) = KALT-LMAX
        KALTMAX(IOBS) = KALT+LMAX
        XLONOBS(IOBS) = REAL(IALT-1)*ZDELLO(KALT)
        SIGRATIO2(IOBS) = (ALTDATA(IOBS,2)/SIGMOD)**2
      ENDDO

      DO IOBS=1,NOBSPE(IRANK)
        IF(ALTDATA(IOBS,1).GT.0.) THEN
          DIFFALTFG(IOBS) =  ALTDATA(IOBS,1)-XMO(IJALT(IOBS,1))
        ELSE
          DIFFALTFG(IOBS) = 0.
        ENDIF
      ENDDO

      DO IOBS=1,NOBSPE(IRANK)
        IF(ALTDATA(IOBS,1).GT.0.) THEN
          KALT = KXLT(IJALT(IOBS,1),IG)
          DO IJ = IJS(IG),IJL(IG)
            IF(XMO(IJ).GE.0.01) THEN
              I = IXLG(IJ,IG)
              K = KXLT(IJ,IG)
              IF(K.GE.KALTMIN(IOBS) .AND. K.LE.KALTMAX(IOBS)) THEN

                IF(IJ.EQ.IJALT(IOBS,1)) THEN
                  DOBS = 0.
                ELSE
                  DELLON = XLONOBS(IOBS) - REAL(I-1)*ZDELLO(K)
                  COSLON = COS(DELLON*RAD)
                  DOBS2  = COSLON*COSPH(KALT)*COSPH(K) +
     &             SINPH(KALT)*SINPH(K)
                  DOBS = ACOS(MAX(MIN(DOBS2,1.),-1.))
                ENDIF

                IF(DOBS.LE.DISTMAX) THEN
                  NOBS(IJ)=NOBS(IJ)+1
                  IF(NOBS(IJ).LE.NDIM2) THEN
                    IOBS4IJ(IJ,NOBS(IJ))=IOBS
                    W(IJ,NOBS(IJ)) = EXP(-DOBS/DIST)
                  ELSE
                    NOBSMAX=MAX(NOBSMAX,NOBS(IJ))
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO

!     IN CASE THE MAXIMUM NUMBER OF OBSERVATIONS INFLUENCING A POINT IS
!     TOO SMALL, RESET IT AND TRY AGAIN.
      IF(NOBSMAX.GT.0) THEN
        WRITE(*,*) ''
        WRITE(*,*) 'WARNING IN WAVE MODEL :'
        WRITE(*,*) '**********************'
        WRITE(*,*) 'IN OIFIELD : THE MAXIMUM NUMBER OF OBSERVATIONS'
        WRITE(*,*) '             INFLUENCING A POINT WAS RESET !!!'
        WRITE(*,*) 'IF THAT HAPPEN OFTEN, YOU SHOULD CHANGE NDIM2 !!'
        WRITE(*,*) 'IRANK = ',IRANK
        WRITE(*,*) 'NDIM2 = ',NDIM2
        WRITE(*,*) 'LMAX = ',LMAX
        WRITE(*,*) 'NOBSMAX = ',NOBSMAX
        WRITE(*,*) ''
        WRITE(IU06,*) ''
        WRITE(IU06,*) 'IN OIFIELD : THE MAXIMUM NUMBER OF OBSERVATIONS'
        WRITE(IU06,*) '             INFLUENCING A POINT WAS RESET !!!'
        WRITE(IU06,*) 'IF THAT HAPPEN OFTEN, YOU SHOULD CHANGE NDIM2 !!'
        WRITE(IU06,*) 'NDIM2 = ',NDIM2
        WRITE(IU06,*) 'LMAX = ',LMAX
        WRITE(IU06,*) 'NOBSMAX = ',NOBSMAX
        WRITE(IU06,*) ''
        CALL FLUSH(IU06)
        NDIM2=MIN(2*NDIM2,NOBSPE(IRANK))
        IF(ALLOCATED(IOBS4IJ)) DEALLOCATE(IOBS4IJ)
        IF(ALLOCATED(W)) DEALLOCATE(W)
        GOTO 1000
      ENDIF

      IF(ALLOCATED(KALTMIN)) DEALLOCATE(KALTMIN)
      IF(ALLOCATED(KALTMAX)) DEALLOCATE(KALTMAX)

!     LOOP ON EACH GRID POINT

      XNEW = 0.
      IG = 1
      DO IJ = IJS(IG),IJL(IG)

        IF (NOBS(IJ).GT.0) THEN

          ALLOCATE(D(NOBS(IJ),NOBS(IJ)))
          ALLOCATE(P(NOBS(IJ),NOBS(IJ)))
          ALLOCATE(XM(NOBS(IJ),NOBS(IJ)))

!         DETERMINE OBSERVATION CORRELATION MATRIX  D

          D = 0.

          DO IOBS=1,NOBS(IJ)
            D(IOBS,IOBS) = SIGRATIO2(IOBS4IJ(IJ,IOBS))
          ENDDO

!         DETERMINE MODEL CORRELATION MATRIX AT OBSERVATION
!         POINTS.
                                                                        
          DO IOBS=1,NOBS(IJ)
            I1 = IXLG(IJALT(IOBS4IJ(IJ,IOBS),1),IG) 
            K1 = KXLT(IJALT(IOBS4IJ(IJ,IOBS),1),IG) 

            DO JOBS=IOBS+1,NOBS(IJ)

              I2 = IXLG(IJALT(IOBS4IJ(IJ,JOBS),1),IG)
              K2 = KXLT(IJALT(IOBS4IJ(IJ,JOBS),1),IG) 

              IF(IJALT(IOBS4IJ(IJ,IOBS),1).EQ.IJALT(IOBS4IJ(IJ,JOBS),1))
     &        THEN
                DOBS = 0.
              ELSE
                DELLON = XLONOBS(IOBS4IJ(IJ,IOBS))
     &           - XLONOBS(IOBS4IJ(IJ,JOBS))
                COSLON = COS(DELLON*RAD)
                DOBS2  = COSLON*COSPH(K1)*COSPH(K2) +
     &           SINPH(K1)*SINPH(K2)
                DOBS = ACOS(MAX(MIN(DOBS2,1.),-1.))
              ENDIF

              P(IOBS,JOBS) = EXP(-DOBS/DIST)
              P(JOBS,IOBS) = P(IOBS,JOBS)
            ENDDO
            P(IOBS,IOBS) = 1.
          ENDDO

          XM=P+D

          DEALLOCATE(D)
          DEALLOCATE(P)

!         INVERSE MATRIX XM

          COND=0.
          EPS=0. 
          ALLOCATE(V(NOBS(IJ)))            
          V = 0.

          CALL SYMINV(XM,NOBS(IJ),NOBS(IJ),COND,V,EPS)
          DO IOBS=1,NOBS(IJ)-1
            DO JOBS=IOBS+1,NOBS(IJ)
             XM(IOBS,JOBS) = XM(JOBS,IOBS)
            ENDDO
          ENDDO

          DEALLOCATE(V)            

!         PRODUCE FIELD AT POINT I,J ACCORDING TO OPTIMUM INTERPOLATION

          DO IOBS=1,NOBS(IJ)
            DIFF = 0.
            DO JOBS=1,NOBS(IJ)
              KOBS=IOBS4IJ(IJ,JOBS)
              IJOBS=IJALT(KOBS,1)
              DIFF=DIFF+XM(IOBS,JOBS)*DIFFALTFG(KOBS)
            ENDDO
            XNEW(IJ)=XNEW(IJ)+W(IJ,IOBS)*DIFF
          ENDDO

          DEALLOCATE(XM)            

        ENDIF

      ENDDO


      IG = 1
      DO IJ = IJS(IG),IJL(IG)
        IF (NOBS(IJ).GT.0) THEN
          XOI(IJ) = XMO(IJ) + XNEW(IJ)
        ENDIF
      ENDDO
                                                                        
      IF(ALLOCATED(XLONOBS)) DEALLOCATE(XLONOBS)
      IF(ALLOCATED(SIGRATIO2)) DEALLOCATE(SIGRATIO2)
      IF(ALLOCATED(DIFFALTFG)) DEALLOCATE(DIFFALTFG)
      IF(ALLOCATED(IJALT)) DEALLOCATE(IJALT)
      IF(ALLOCATED(ALTDATA)) DEALLOCATE(ALTDATA)
      IF(ALLOCATED(IOBS4IJ)) DEALLOCATE(IOBS4IJ)
      IF(ALLOCATED(W)) DEALLOCATE(W)
                                                                  
      RETURN                                                            
      END SUBROUTINE OIFIELD
