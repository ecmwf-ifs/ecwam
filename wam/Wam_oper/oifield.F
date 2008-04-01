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

      USE YOWALTAS , ONLY : NALTAVLB ,IJALT    ,INTLMAX  ,
     &            NOBSPE   ,ALTDATA  ,LALTPASSIV,
     &            XLONOBS  , SIGRATIO2, DIFFALTFG
      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,NLONRGG  ,IGL      ,
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,XDELLA   ,
     &            ZDELLO
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      INTEGER :: JC
      INTEGER :: NOBS(NSTART(IRANK):NEND(IRANK))
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IOBS4IJ
      INTEGER JKGLO,KIJS,KIJL,NPROMA


      REAL :: XMO(NIBLO)
      REAL :: XOI(NINF:NSUP)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KALTMIN, KALTMAX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: W

! ----------------------------------------------------------------------


!*    1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.
!     -----------------------------------------------------

      LMAX = NINT(DEG*DISTMAX/XDELLA)

      NDIM2 = (5+2*(NALTAVLB-1))*(2*LMAX+1)

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
!           TAKE ONLY POSITIVE ACTIVE HS
            IF(XMO(IJ).GE.0.01 .AND. .NOT.LALTPASSIV(IJALT(IOBS,2)))THEN
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
        WRITE(*,*) 'IF THAT HAPPENS OFTEN, YOU SHOULD CHANGE NDIM2 !!'
        WRITE(*,*) 'IRANK = ',IRANK
        WRITE(*,*) 'NDIM2 = ',NDIM2
        WRITE(*,*) 'LMAX = ',LMAX
        WRITE(*,*) 'NOBSMAX = ',NOBSMAX
        WRITE(*,*) ''
        WRITE(IU06,*) ''
        WRITE(IU06,*) 'IN OIFIELD : THE MAXIMUM NUMBER OF OBSERVATIONS'
        WRITE(IU06,*) '             INFLUENCING A POINT WAS RESET !!!'
        WRITE(IU06,*) 'IF THAT HAPPENS OFTEN, YOU SHOULD CHANGE NDIM2 !'
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

      IG = 1
      NPROMA=(IJL(IG)-IJS(IG)+1)/NWAM_BLKS+1
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
      DO JKGLO=IJS(IG),IJL(IG),NPROMA
        KIJS=JKGLO
        KIJL=MIN(KIJS+NPROMA-1,IJL(IG))
        CALL WAMOI(NOBS(KIJS:KIJL),IOBS4IJ(KIJS:KIJL,:),W(KIJS:KIJL,:),
     &             KIJS,KIJL,NDIM2,
     &             DIST, XMO(KIJS:KIJL), XOI(KIJS:KIJL), IG)
      ENDDO
!$OMP END PARALLEL DO



      IF(ALLOCATED(XLONOBS)) DEALLOCATE(XLONOBS)
      IF(ALLOCATED(SIGRATIO2)) DEALLOCATE(SIGRATIO2)
      IF(ALLOCATED(DIFFALTFG)) DEALLOCATE(DIFFALTFG)
      IF(ALLOCATED(IJALT)) DEALLOCATE(IJALT)
      IF(ALLOCATED(ALTDATA)) DEALLOCATE(ALTDATA)
      IF(ALLOCATED(IOBS4IJ)) DEALLOCATE(IOBS4IJ)
      IF(ALLOCATED(W)) DEALLOCATE(W)
                                                                  
      RETURN                                                            
      END SUBROUTINE OIFIELD
