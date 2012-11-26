      SUBROUTINE OIFIELD (MINIJS,MAXIJL,XMO,CICVR,
     &                    XOI,SIGMOD,DIST,DISTMAX,EXTENDMAX)

!**** *OIFIELD* - OPTIMUM INTERPOLATION.                                

!     P.LIONELLO     ECMWF       APRIL 1990                             
!     J. BIDLOT      ECMWF       SEPTEMBER 96  MESSAGE PASSING
!     J. BIDLOT      ECMWF       OCTOBER 2009 : INTRODUCE VARIABLE DIST

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

!       *CALL* *OIFIELD (MINIJS,MAXIJL,XMO,CICVR,
!    &                    XOI,SIGMOD,DIST,DISTMAX,EXTENDMAX)

!         *MINIJS*  INTEGER  MINUMUM INDEX OF XMO AND CICVR
!         *MAXIJL*  INTEGER  MAXIMUM INDEX OF XMO AND CICVR 
!         *XMO*     REAL     FIELD FROM MODEL (FIRST GUESS)
!         *XOI*     REAL     FIELD FROM O.I. ON OUTPUT.
!         *CICVR*   REAL     MODEL SEA ICE
!         *SIGMOD*  REAL     MODEL ERROR ESTIMATE.                             
!         *DIST*    REAL     CORRELATION DISTANCE (in radian)
!         *DISTMAX* REAL     MAXIMUM OF DIST
!         *EXTENDMAX* REAL   EXTENDMAX*DIST = MAXIMUM SPREADING DISTANCE

! ----------------------------------------------------------------------

      USE YOWALTAS , ONLY : NALTAVLB ,IJALT    ,INTLMAX  ,
     &            NOBSPE   ,ALTDATA  ,LALTPASSIV,
     &            XLONOBS  ,SIGRATIO2,DIFFALTFG
      USE YOWICE   , ONLY : CITHRSH_SAT 
      USE YOWGRID  , ONLY : SINPH    ,COSPH    ,NLONRGG  ,IGL      ,
     &            IJS      ,IJL
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,IPER     ,XDELLA   ,
     &            ZDELLO   ,AMOWEP
      USE YOWMPP   , ONLY : IRANK    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY
      USE YOWPCONS , ONLY : RAD      ,DEG      ,R
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWSTAT  , ONLY : NWAM_BLKS
      USE YOWTEST  , ONLY : IU06

! ----------------------------------------------------------------------

      INTEGER :: MINIJS, MAXIJL 
      INTEGER :: JC
      INTEGER :: NOBS(NSTART(IRANK):NEND(IRANK))
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IOBS4IJ
      INTEGER :: JKGLO,KIJS,KIJL,NPROMA

      REAL, DIMENSION(MINIJS:MAXIJL) :: XMO, CICVR
      REAL :: DIST(NSTART(IRANK):NEND(IRANK))
      REAL :: XOI(NSTART(IRANK):NEND(IRANK))
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KALTMIN, KALTMAX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: W

      LOGICAL :: LLINVIEW
      LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: LLSIMASK

! ----------------------------------------------------------------------


!*    1. INITIALISE SEARCH DISTANCE AND CORRELATION LENGTH.
!     -----------------------------------------------------

      IG=1

      DIMAX = EXTENDMAX*DISTMAX

      LMAX = INT(DEG*DIMAX/XDELLA)+1

      NDIM2 = (5+2*(NALTAVLB-1))*(2*LMAX+1)

!     THE MEASUREMENTS ARE STORED IN THE O.I. FIELD.                  
!     NEGATIVE VALUES WILL REMAIN WHERE O.I. WILL PRODUCE NO RESULTS  
      XOI = -1.

!      LOOP OVER GRID POINTS.                                          
!      ----------------------                                          

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
      KMIN=NGY
      KMAX=1
      DO IOBS=1,NOBSPE(IRANK)
        IALT = IXLG(IJALT(IOBS,1),IG)
        KALT = KXLT(IJALT(IOBS,1),IG)
        KALTMIN(IOBS) = KALT-LMAX
        KMIN=MIN(KMIN,KALTMIN(IOBS))
        KALTMAX(IOBS) = KALT+LMAX
        KMAX=MAX(KMAX,KALTMAX(IOBS))
        XLONOBS(IOBS) = REAL(IALT-1)*ZDELLO(KALT)
        SIGRATIO2(IOBS) = (ALTDATA(IOBS,2)/SIGMOD)**2
      ENDDO
      KMIN=MAX(1,KMIN)
      KMAX=MIN(NGY,KMAX)

      DO IOBS=1,NOBSPE(IRANK)
        IF(ALTDATA(IOBS,1).GT.0.) THEN
          DIFFALTFG(IOBS) =  ALTDATA(IOBS,1)-XMO(IJALT(IOBS,1))
        ELSE
          DIFFALTFG(IOBS) = 0.
        ENDIF
      ENDDO

!     SET LAND SEA ICE MASK (LLSIMASK=.true. if over land or sea ice)
!     BASED ON FIRST GUESS HS AND SEA ICE COVER WITH THRESHOLD CITHRSH_SAT
      IF(NOBSPE(IRANK).GT.0) THEN

        IF(.NOT.ALLOCATED(LLSIMASK))THEN
           ALLOCATE(LLSIMASK(-LMAX:NGX+LMAX+1,KMIN:KMAX))
        ENDIF
        DO K = KMIN,KMAX
          DO I = 1,NLONRGG(K)
            LLSIMASK(I,K) = .TRUE. 
          ENDDO
        ENDDO
        DO IJ = MINIJS, MAXIJL
          K = KXLT(IJ,IG)
          IF( K.GE.KMIN .AND. K.LE.KMAX ) THEN
            I = IXLG(IJ,IG)
            IF(XMO(IJ).GT. 0.0 .AND. CICVR(IJ).LE.CITHRSH_SAT) THEN
              LLSIMASK(I,K) = .FALSE. 
            ENDIF
          ENDIF
        ENDDO
!       periodicity
        IF(IPER.EQ.1) THEN
          DO K = KMIN,KMAX
            DO I = -LMAX,0 
              LLSIMASK(I,K) = LLSIMASK(NLONRGG(K)+I-1,K)
            ENDDO
            DO I = NLONRGG(K)+1,NLONRGG(K)+LMAX+1
              LLSIMASK(I,K) = LLSIMASK(I-NLONRGG(K),K)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

!     FIND HOW MANY OBSERVATIONS CAN INFLUENCE EACH GRID POINT
      DO IOBS=1,NOBSPE(IRANK)
        IF(ALTDATA(IOBS,1).GT.0.) THEN
          IALT = IXLG(IJALT(IOBS,1),IG)
          KALT = KXLT(IJALT(IOBS,1),IG)
!         ALTIMETER LONGITUDE WITH RESPECT TO AMOWEP
          ALTLON=(IALT-1)*ZDELLO(KALT)
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

!               FIND WHETHER THE OBSERVATION IS WITHIN THE MAXIMUM
!               SPREADING DISTANCE FROM THE GRID POINT.

                IF(DOBS.LE.EXTENDMAX*DIST(IJ)) THEN

!                 FIND IF LAND OR ICE IS NOT BLOCKING THE
!                 LINE OF VIEW BETWEEN OBSERVATION AND GRID POINT.

!                 MODEL LONGITUDE WITH RESPECT TO AMOWEP
                  XMODLON=(I-1)*ZDELLO(K)

                  LLINVIEW=.TRUE.

                  IF(K.EQ.KALT) THEN

                    ! periodicity ?
                    IF(IPER.EQ.1 .AND. ABS(IALT-I).GT.NLONRGG(K)/2) THEN
                      IF(IALT.GT.I) THEN
                        IIALT=IALT-NLONRGG(K)
                      ELSE
                        IIALT=IALT+NLONRGG(K)
                      ENDIF
                    ELSE
                      IIALT=IALT
                    ENDIF
                    IIMIN=MAX(MIN(I,IIALT),-LMAX)
                    IIMAX=MIN(MAX(I,IIALT),NGX+LMAX+1)
                    IC=IIMIN
                    DO WHILE (IC.LE.IIMAX .AND. LLINVIEW)
                      IF(LLSIMASK(IC,K)) THEN
                         LLINVIEW=.FALSE.
                         EXIT
                      ENDIF
                      IC=IC+1
                    ENDDO

                  ELSE
                    KKMIN=MIN(K+1,KALT)
                    KKMAX=MAX(K-1,KALT)
                    NBLK=0
                    DO KK=KKMIN,KKMAX
                      XII = XMODLON/ZDELLO(KK)+1.
                      II = MAX(1,INT(XII))
                      XIIALT = ALTLON/ZDELLO(KK)+1.
                      IIALT = MAX(1,INT(XIIALT))
                      ! periodicity ?
                      IF(IPER.EQ.1 .AND.
     &                   ABS(IIALT-II).GT.NLONRGG(KK)/2) THEN
                        IF(IIALT.GT.II) THEN
                          IIALT=IIALT-NLONRGG(KK)
                        ELSE
                          IIALT=IIALT+NLONRGG(KK)
                        ENDIF
                      ENDIF
                      IIMIN=MAX(MIN(II,IIALT),-LMAX)
                      IIMAX=MIN(MAX(II,IIALT),NGX+LMAX+1)
                      IC=IIMIN
                      DO WHILE (IC.LE.IIMAX)
                        IF(LLSIMASK(IC,KK)) NBLK=NBLK+1
                        IC=IC+1
                      ENDDO
                    ENDDO
                    IF(NBLK.GE.(KKMAX-KKMIN+1)) LLINVIEW=.FALSE.

                  ENDIF

                  IF(LLINVIEW) THEN
                    NOBS(IJ)=NOBS(IJ)+1
                    IF(NOBS(IJ).LE.NDIM2) THEN
                      IOBS4IJ(IJ,NOBS(IJ))=IOBS
                      W(IJ,NOBS(IJ)) = EXP(-DOBS/DIST(IJ))
                    ELSE
                      NOBSMAX=MAX(NOBSMAX,NOBS(IJ))
                    ENDIF
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
     &             DIST(KIJS:KIJL), XMO(KIJS:KIJL), XOI(KIJS:KIJL), IG)
      ENDDO
!$OMP END PARALLEL DO



      IF(ALLOCATED(XLONOBS)) DEALLOCATE(XLONOBS)
      IF(ALLOCATED(SIGRATIO2)) DEALLOCATE(SIGRATIO2)
      IF(ALLOCATED(IOBS4IJ)) DEALLOCATE(IOBS4IJ)
      IF(ALLOCATED(W)) DEALLOCATE(W)
      IF(ALLOCATED(LLSIMASK)) DEALLOCATE(LLSIMASK)
                                                                  
      RETURN                                                            
      END SUBROUTINE OIFIELD
