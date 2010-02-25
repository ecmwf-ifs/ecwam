      SUBROUTINE IFSTOWAM (CDTWIS, NFIELDS, NGPTOTG, NCA, NRA, FCRANGE,
     &                     ISEC1, ISEC2, FIELDS, LWCUR, MASK_IN)

!***  *IFSTOWAM* - INTERPOLATES FORCING FIELDS FROM IFS TO WAM GRID. 


!     PURPOSE                                                       
!     -------                                                      

!     *IFSTOWAM* - INTERPOLATES FORCING FIELDS FROM IFS TO WAM GRID. 
!                  ONLY FOR COUPLED RUNS.

!     INTERFACE                                                         
!     ---------                                                         

!     *CALL* *IFSTOWAM (CDTWIS, NFIELDS, NGPTOTG, NCA, NRA, FCRANGE,
!    &                  ISEC1, ISEC2, FIELDS, LWCUR, MASK_IN)
!
!        *CDTWIS* - EXPECTED DATE/TIME OF THE DATA OBTAINED FROM IFS.
!        *NFIELDS*- NUMBER OF FIELDS HOLDING ATMOSPHERIC DATA
!        *NGPTOTG*- NUMBER OF ATMOSPHERIC GRID POINTS
!        *NCA*    - NUMBER OF ATM. COLUMNS OF LONGITUDE NEAR EQUATOR
!        *NRA*    - NUMBER OF ATM. ROWS OF LATITUDES
!        *FCRANGE*- FORECAST RANGE IN SECONDS
!        *ISEC1*  - GRIB ISEC1 (TO EXTRACT DATE/TIME OF ATM. FIELDS)
!        *ISEC2*  - GRIB ISEC2 (TO EXTRACT ATM. GRID CHARACTERISTICS)
!        *FIELDS* - ATMOSPHERIC FIELDS AS FOLLOWS:
!                   FIELDS(:,1) = U COMPONENT OF WIND SPEED
!                   FIELDS(:,2) = V COMPONENT OF WIND SPEED
!                   FIELDS(:,3) = AIR DENSITY
!                   FIELDS(:,4) = w* USED FOR GUSTINESS
!                   FIELDS(:,5) = SEA ICE FRACTION 
!                   FIELDS(:,6) = U COMPONENT OF SURFACE CURRENT 
!                   FIELDS(:,7) = V COMPONENT OF SURFACE CURRENT 
!        *LWCUR* - LOGICAL CONTROLLING THE PRESENCE OF MEANINGFUL
!                   SURFACE CURRENTS (I.E. NOT ALL ZEROS).
!        *MASK_IN*- MASK FOR FIELDS TO ONLY POINTS TO THE PART
!                   OF FIELDS THAT ARE NEEDED ON EACH PE.

!     EXTERNALS                                                         
!     ---------                                                         

!     *ABORT1*         TERMINATES PROCESSING
!     *INCDATE*        INCREMENTS DATE 
!     *INITIALINT*     COMPUTES INTERPOLATION COEFFICIENTS (FOR FLDINTER)
!     *FLDINTER*       INTERPOLATES ATMOSPHERIC FIELDS TO WAM GRID

! --------------------------------------------------------------------- 

      USE YOWCOUP  , ONLY : LWCOU    ,LWCOUNORMS,
     &            I_MASK_IN,N_MASK_IN
      USE YOWGRIBHD, ONLY : NKSEK0   ,NKSEK1   ,NKSEK2   ,NKSEK3,
     &            NKSEK4   ,NPSEK2   ,NPSEK3   ,KSEK3    ,PSEK3
      USE YOWGRID  , ONLY : IGL      ,NLONRGG
      USE YOWICE   , ONLY : IPARAMCI
      USE YOWMAP   , ONLY : IRGG     ,AMOWEP   ,AMOSOP   ,AMOEAP   ,
     &            AMONOP   ,XDELLA   ,XDELLO   ,ZDELLO   ,IFROMIJ  , 
     &            JFROMIJ 
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWPARAM , ONLY : NGX      ,NGY      ,NIBLO    ,LL1D
      USE YOWPCONS , ONLY : ZMISS    ,ROAIR    ,WSTAR0
      USE YOWSTAT  , ONLY : LADEN    ,LGUST    ,NWAM_BLKS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE YOWWNDG  , ONLY : DLAM     ,DPHI     ,RLATS    ,RLATN    ,
     &            RLONL    ,RLONR    ,KCOL     ,KROW     ,ICODE    ,
     &            IWPER    ,ICOORD
      USE YOWWIND  , ONLY : NC       ,NR       ,FIELDG   ,IUNITW   ,
     &            CWDFILE  ,LLNEWCURR,NORMS_IN

      USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
      USE MPL_MODULE


! --------------------------------------------------------------------- 

      INTEGER, ALLOCATABLE :: KGRIB(:)

      INTEGER NFIELDS, NGPTOTG, NCA, NRA
      INTEGER ISEC1(NKSEK1), ISEC2(22+NRA)
      INTEGER, ALLOCATABLE, SAVE :: RLONRGG(:), IJBLOCK(:,:)
! Mod for OPENMP
      INTEGER JKGLO,KIJS,KIJL,NPROMA,I,JF
! End Mod for OPENMP

      REAL FIELDS(NGPTOTG,NFIELDS)
      INTEGER :: IC
      INTEGER MASK_IN(NGPTOTG)
      LOGICAL, SAVE :: LLINTERPOL

      INTEGER :: ISENDREQ, ITAG
      REAL, ALLOCATABLE :: ZCOMBUF(:)

      INTEGER, ALLOCATABLE, SAVE :: JJ(:), II(:,:), IIP(:,:)
      INTEGER, ALLOCATABLE, SAVE :: ISTART(:)
      REAL, ALLOCATABLE, SAVE :: DJ1(:), DII1(:,:), DIIP1(:,:)
      REAL, DIMENSION(NFIELDS,4) :: NORMS ! 1:average, 2:min, 3:max, 4:variance 

      CHARACTER*14 CDTWIR, CDTWIS
      CHARACTER*16 :: AFLABEL(NFIELDS)

      LOGICAL :: LLDUM
      LOGICAL FRSTMATM
      LOGICAL :: LWCUR
      LOGICAL :: LFRSTIME
      LOGICAL :: LFRST

      DATA FRSTMATM / .TRUE. /
      DATA LFRSTIME / .TRUE. /
      DATA LFRST /.TRUE./

      SAVE IPERIODIC

      REAL :: ZHOOK_HANDLE

! --------------------------------------------------------------------  

      IF (LHOOK) CALL DR_HOOK('IFSTOWAM',0,ZHOOK_HANDLE)

!     ONLY ACTIVE IF COUPLED TO IFS.
!     -----------------------------
      IF(LWCOU) THEN
        IF (FRSTMATM) THEN
!        
!         1.1 FIRST TIME ONLY, COMPUTE INTERPOLATION COEFFICIENTS:
!             ----------------------------------------------------
!
          IF (.NOT.ALLOCATED(RLONRGG)) ALLOCATE (RLONRGG(NRA))
          IF (.NOT.ALLOCATED(JJ)) ALLOCATE (JJ(NGY))
          IF (.NOT.ALLOCATED(II)) ALLOCATE (II(NGX,NGY))
          IF (.NOT.ALLOCATED(IIP)) ALLOCATE (IIP(NGX,NGY))
          IF (.NOT.ALLOCATED(DJ1)) ALLOCATE (DJ1(NGY))
          IF (.NOT.ALLOCATED(DII1)) ALLOCATE (DII1(NGX,NGY))
          IF (.NOT.ALLOCATED(DIIP1)) ALLOCATE (DIIP1(NGX,NGY))

          
          CALL INITIALINT(IU06, ITEST, NCA, NRA, ISEC2,
     &                    NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                    AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                    RLONRGG, 
     &                    LLINTERPOL,DJ1, DII1, DIIP1, JJ, II, IIP)

          IF (.NOT.ALLOCATED(IJBLOCK)) ALLOCATE(IJBLOCK(0:NCA+1,NRA))

          L = 0
          DO J=1,NRA
            JSN=NRA-J+1
            DO I=1,RLONRGG(JSN)
              L = L+1
              IJBLOCK(I,J) = L
            ENDDO
          ENDDO
          IF(NGPTOTG.NE.L)THEN
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ####  WARNING IN  IFSTOWAM WARNING  ###'
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ##  NGPTOTG  IS  NOT  EQUAL TO  L  ###'
            WRITE(IU06,*)' ##  NGPTOTG = ', NGPTOTG
            WRITE(IU06,*)' ##        L = ', L
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' ##  NCA = ', NCA, '   NRA = ', NRA
            WRITE(IU06,*)' ##  SUM RLONRGG = ',SUM(RLONRGG(1:NRA))
            WRITE(IU06,*)' ##  RLONRGG = ', RLONRGG(1:NRA)
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' #######################################'
            WRITE(IU06,*)' #######################################'
            CALL FLUSH(IU06)
          ENDIF

          IF(IPERIODIC.EQ.1) THEN
            DO J=1,NRA
            JSN=NRA-J+1
            IJBLOCK(0,J)= IJBLOCK(RLONRGG(JSN),J)
            IJBLOCK(RLONRGG(JSN)+1,J)= IJBLOCK(1,J)
            ENDDO
          ENDIF
     
          FRSTMATM = .FALSE.
          
          IF (ITEST.GT.1)
     &     WRITE(IU06,*) ' WAM_IFSTOWAM - INTERPOLATION COEFFICIENTS'
     &                     //' COMPUTED'
     
          IF (LADEN.AND.NFIELDS.LT.3) THEN 
            LADEN = .FALSE.
            WRITE(IU06,*) ' WAM_IFSTOWAM - WARNING: '
            WRITE(IU06,*) '      AIR DENSITY RUN WAS REQUESTED BUT '
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 3'
            WRITE(IU06,*) '      RESET  LADEN  TO  .FALSE.'
          ENDIF
          
          IF (LADEN) THEN 
            WRITE(IU06,*) ' WAM_IFSTOWAM - VARIABLE AIR DENSITY RUN'
          ELSE
            WRITE(IU06,*) ' WAM_IFSTOWAM - AIR DENSITY IS '
     &                     //' CONSTANT = ', ROAIR
          ENDIF
          
          IF (LGUST.AND.NFIELDS.LT.4) THEN 
            LGUST = .FALSE.
            WRITE(IU06,*) ' SUB. IFSTOWAM - WARNING: '
            WRITE(IU06,*) '      GUSTINESS RUN WAS REQUESTED BUT '
     &                 //'NOT ENOUGH FIELDS PASSED.'
            WRITE(IU06,*) '      NFIELDS = ', NFIELDS
            WRITE(IU06,*) '      NFIELDS MUST BE AT LEAST 4'
            WRITE(IU06,*) '      RESET  LGUST  TO  .FALSE.'
          ENDIF
          
          IF (LGUST) THEN 
            WRITE(IU06,*) ' WAM_IFSTOWAM - GUSTINESS RUN'
          ELSE
            WRITE(IU06,*) ' WAM_IFSTOWAM - NO GUSTINESS EFFECT.'
          ENDIF
        ENDIF
!        
!       COMPUTATION OF THE NORMS OF INPUT FIELDS
!
        AFLABEL(1)='Neutral 10m U   '
        AFLABEL(2)='Neutral 10m V   '
        IF (NFIELDS.GE.3) AFLABEL(3)='Air density     '
        IF (NFIELDS.GE.4) AFLABEL(4)='w*              '
        IF (NFIELDS.GE.5) AFLABEL(5)='Sea ice fraction'
        IF (NFIELDS.GE.6) AFLABEL(6)='U-ocean-current '
        IF (NFIELDS.GE.7) AFLABEL(7)='V-ocean-current '

!       COMPUTATION OF THE NORMS OF INPUT FIELDS
!       ----------------------------------------

!       NORMAS ARE ALWAYS NEEDED AS THEY ARE ALSO USED TO CHECK
!       IF CURRENTS AND SEA ICE FRACTION FIELDS WERE UPDATED !!!

!       *GLOBAL* ONLY AVAILABLE THE FIRST TIME OR
!       IF REQUESTED BY SETTING LWCOUNORMS=TRUE IN INPUT NAMELIST.
        IF( LFRSTIME .OR. LWCOUNORMS )THEN
          WRITE(IU06,*) ' '
          WRITE(IU06,*) ' *GLOBAL* NORMS OF INPUT ATMOSPHERIC FIELDS :'
          CALL GSTATS(1436,0)
!$OMP     PARALLEL DO SCHEDULE(STATIC) PRIVATE(JF,I,IC)
          DO IFLD=1,NFIELDS
            JF=IFLD
            DO IC=1,3
              NORMS(JF,IC)=FIELDS(1,JF)
            ENDDO
            NORMS(JF,4)=FIELDS(1,JF)**2
            DO I=2,NGPTOTG 
              NORMS(JF,1)=NORMS(JF,1)+FIELDS(I,JF)
              NORMS(JF,2)=MIN(NORMS(JF,2),FIELDS(I,JF))
              NORMS(JF,3)=MAX(NORMS(JF,3),FIELDS(I,JF))
              NORMS(JF,4)=NORMS(JF,4)+FIELDS(I,JF)**2
            ENDDO
            NORMS(JF,1)=NORMS(JF,1)/NGPTOTG
            NORMS(JF,4)=NORMS(JF,4)/NGPTOTG
            NORMS(JF,4)=NORMS(JF,4)-NORMS(JF,1)**2
          ENDDO
!$OMP     END PARALLEL DO
          CALL GSTATS(1436,1)
          LFRSTIME=.FALSE.
          IF (.NOT.ALLOCATED(ISTART)) ALLOCATE (ISTART(NFIELDS))
          ISTART=0
          N_MASK_IN=NGPTOTG

          DO IFLD=1,NFIELDS
            WRITE(IU06,*) ' ',AFLABEL(IFLD),
     &       NORMS(IFLD,1),NORMS(IFLD,2),NORMS(IFLD,3),N_MASK_IN
          ENDDO

        ELSE
!       *LOCAL* NORM

          IRECV=1

          IF(IRANK.EQ.IRECV) THEN
          WRITE(IU06,*) ' NORMS OF GLOBAL INPUT ATMOSPHERIC FIELDS :'
          ELSE
          WRITE(IU06,*) ' *Local* NORMS OF INPUT ATMOSPHERIC FIELDS : '
          ENDIF
          WRITE(IU06,*) ' *Can only be compared if the log files'
          WRITE(IU06,*) '  are for processor: ',IRANK 
          WRITE(IU06,*) '  have the same number of processors: ',NPROC 
          IF(LL1D) THEN
          WRITE(IU06,*) '  and have the same model decomposition: ONE D'
          ELSE
          WRITE(IU06,*) '  and have the same model decomposition: TWO D'
          ENDIF

!         INITIALISE I POINTER TO MASK_IN
          IF(LFRST) THEN
            LFRST=.FALSE.
            JF=1
            N_MASK_IN=0
            DO I=1,NGPTOTG 
              IF(MASK_IN(I).EQ.1) THEN
                N_MASK_IN=N_MASK_IN+1
              ENDIF
            ENDDO
            ALLOCATE(I_MASK_IN(MAX(1,N_MASK_IN)))
            I_MASK_IN(1)=1

            IC=0
            DO I=1,NGPTOTG 
              IF(MASK_IN(I).EQ.1) THEN
                IC=IC+1
                I_MASK_IN(IC)=I
              ENDIF
            ENDDO
          ENDIF

!         COMPUTE NORMS
          DO IFLD=1,NFIELDS
            JF=IFLD
            VAL=FIELDS(I_MASK_IN(1),JF)
            NORMS(JF,1)=VAL
            NORMS(JF,2)=VAL
            NORMS(JF,3)=VAL
            NORMS(JF,4)=VAL**2
            DO IC=2,N_MASK_IN
                VAL=FIELDS(I_MASK_IN(IC),JF)
                NORMS(JF,1)=NORMS(JF,1)+VAL
                NORMS(JF,2)=MIN(NORMS(JF,2),VAL)
                NORMS(JF,3)=MAX(NORMS(JF,3),VAL)
                NORMS(JF,4)=NORMS(JF,4)+VAL**2
            ENDDO
            NORMS(JF,1)=NORMS(JF,1)/MAX(1,N_MASK_IN)
            NORMS(JF,4)=NORMS(JF,4)/MAX(1,N_MASK_IN)
            NORMS(JF,4)=NORMS(JF,4)-NORMS(JF,1)**2
          ENDDO

!         FOR PRIMARY PE (IRECV), COLLECT ALL THE NORMS TO PRODUCE A
!         PSEUDO GLOBAL NORM.
 
          NCOMLOC=1+4*NFIELDS
          NCOMBUF=NCOMLOC*NPROC
          ALLOCATE(ZCOMBUF(NCOMBUF))

          IST=1+(IRANK-1)*NCOMLOC
          IED=IST+NCOMLOC-1
          ICOUNT=IST
          ZCOMBUF(ICOUNT)=N_MASK_IN
          DO IFLD=1,NFIELDS
            DO IC=1,4
              ICOUNT=ICOUNT+1
              ZCOMBUF(ICOUNT)=NORMS(IFLD,IC)
            ENDDO
          ENDDO

          ITAG=1
          CALL MPL_SEND(ZCOMBUF(IST:IED),KDEST=IRECV,KTAG=ITAG,
     &                  KMP_TYPE=JP_NON_BLOCKING_STANDARD,
     &                  KREQUEST=ISENDREQ,
     &                  CDSTRING='IFSTOWAM: MPL_SEND' )
          IF(IRANK.EQ.IRECV) THEN
            DO IP=1,NPROC
              IST=1+(IP-1)*NCOMLOC
              IED=IST+NCOMLOC-1
              CALL MPL_RECV(ZCOMBUF(IST:IED),
     &                      KSOURCE=IP,KTAG=ITAG,
     &                      KMP_TYPE=JP_BLOCKING_STANDARD,
     &                      CDSTRING='IFSTOWAM: MPL_RECV' )

            ENDDO
            CALL MPL_WAIT(ZCOMBUF(:),KREQUEST=ISENDREQ,
     &                    CDSTRING='IFSTOWAM: WAIT MPL_SEND')
          ENDIF


!         COMPUTE PSEUDO GLOBAL NORM
          IF(IRANK.EQ.IRECV) THEN
            IST=1+(IRANK-1)*NCOMLOC
            ICOUNT=IST
            NMASK=ZCOMBUF(ICOUNT)
            NTOT=NMASK
            DO JF=1,NFIELDS
              ICOUNT=ICOUNT+1
              NORMS(JF,1)=NMASK*ZCOMBUF(ICOUNT)
              ICOUNT=ICOUNT+3
              NORMS(JF,4)=NMASK*ZCOMBUF(ICOUNT)
            ENDDO
            DO IP=1,NPROC
              IF(IP.NE.IRECV) THEN
                IST=1+(IP-1)*NCOMLOC
                ICOUNT=IST
                NMASK=ZCOMBUF(ICOUNT)
                NTOT=NTOT+NMASK
                DO JF=1,NFIELDS
                  ICOUNT=ICOUNT+1
                  NORMS(JF,1)=NORMS(JF,1)+NMASK*ZCOMBUF(ICOUNT)
                  ICOUNT=ICOUNT+1
                  NORMS(JF,2)=MIN(NORMS(JF,2),ZCOMBUF(ICOUNT))
                  ICOUNT=ICOUNT+1
                  NORMS(JF,3)=MAX(NORMS(JF,3),ZCOMBUF(ICOUNT))
                  ICOUNT=ICOUNT+1
                  NORMS(JF,4)=NORMS(JF,4)+NMASK*ZCOMBUF(ICOUNT)
                ENDDO
              ENDIF
            ENDDO
            DO JF=1,NFIELDS
              NORMS(JF,1)=NORMS(JF,1)/NTOT
              NORMS(JF,4)=NORMS(JF,4)/NTOT
            ENDDO
          ELSE
            NTOT=N_MASK_IN
          ENDIF

          DO IFLD=1,NFIELDS
            WRITE(IU06,*) ' ',AFLABEL(IFLD),
     &         NORMS(IFLD,1),NORMS(IFLD,2),NORMS(IFLD,3),NTOT,
     &                    IRANK, NPROC, LL1D
          ENDDO

          DEALLOCATE(ZCOMBUF)
        ENDIF
        WRITE(IU06,*) ' '
        IF(ITEST.GT.0) CALL FLUSH(IU06)


!       CHECK WHETHER INPUT SURFACE CURRENTS HAVE CHANGED
        IF(.NOT.LWCUR) THEN
          LLNEWCURR=.FALSE.
        ELSE IF(NFIELDS.GE.7) THEN
          LLDUM=.TRUE.
          DO IC=1,4
            DO IFLD=6,7
              LLDUM=(LLDUM .AND. (NORMS(IFLD,IC).EQ.NORMS_IN(IFLD,IC)))
            ENDDO
          ENDDO
          LLNEWCURR=.NOT.LLDUM
        ENDIF

!       SAVE NORMS
        DO IC=1,4
          DO IFLD=1,NFIELDS
            NORMS_IN(IFLD,IC)=NORMS(IFLD,IC)
          ENDDO
        ENDDO

        IWPER = 1
        ICOORD = 1

! Mod for OPENMP
        NPROMA=(NSUP-NINF+1)/NWAM_BLKS+1
        CALL GSTATS(1437,0)
!$OMP   PARALLEL DO SCHEDULE(DYNAMIC,1) PRIVATE(JKGLO,KIJS,KIJL)
        DO JKGLO=NINF, NSUP,NPROMA
          KIJS=JKGLO
          KIJL=MIN(KIJS+NPROMA-1,NSUP)

          CALL FLDINTER (IU06, ITEST, NGPTOTG, NCA, NRA, NFIELDS,FIELDS,
     &                   NGX, NGY, IRGG, NLONRGG, XDELLA, ZDELLO,
     &                   IFROMIJ, JFROMIJ, NINF, NSUP, IGL, KIJS, KIJL,
     &                   AMOWEP, AMOSOP, AMOEAP, AMONOP, IPERIODIC,
     &                   RLONRGG, IJBLOCK, ZMISS, 
     &                   LADEN, ROAIR, LGUST, WSTAR0,LWCUR,LLINTERPOL,
     &                   DJ1, DII1, DIIP1, JJ, II, IIP, MASK_IN)
        ENDDO
!$OMP   END PARALLEL DO
        CALL GSTATS(1437,1)

        IF (LADEN) THEN 
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. IFSTOWAM - AIR DENSITY COMPUTED'
     &                     //' AND INTERPOLATED'
        ELSE
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. IFSTOWAM - AIR DENSITY USED AS '
     &                    //' CONSTANT = ', ROAIR
        ENDIF

        IF (LGUST) THEN 
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. IFSTOWAM - w* INTERPOLATED'
        ENDIF

        IF (LWCUR) THEN 
          IF (ITEST.GT.1)
     &      WRITE(IU06,*) ' SUB. IFSTOWAM - CURRENTS INTERPOLATED'
        ENDIF

! Decode the date of the fields:
!       START DATE.

        IYYYY=(ISEC1(21)-1)*100+ISEC1(10)
        WRITE(CDTWIR(1:12),'(I4.4,4I2.2)') IYYYY,ISEC1(11),ISEC1(12),
     &                             ISEC1(13),ISEC1(14)
        CDTWIR(13:14)='00'


!       DETERMINE FORECAST PERIOD

        IFORP = NINT(FCRANGE)

!       GET DATE OF THE FORECAST INSTEAD OF STARTING DATE

        CALL INCDATE (CDTWIR,IFORP)

      
        IF (CDTWIR.LT.CDTWIS) THEN
!         DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
!         ------------------------------------------------
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. IFSTOWAM         *'
          WRITE (IU06,*) ' *      ==========================        *'
          WRITE (IU06,*) ' * WIND DATE IS EARLIER THAN EXPECTED     *'
          WRITE (IU06,*) ' * DATE FROM IFS IS CDTWIR = ', CDTWIR
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ELSEIF (CDTWIR.GT.CDTWIS) THEN
!         DATE OF INPUT FIELD IS LATER THAN REQUESTED DATE
!         ------------------------------------------------
          WRITE (IU06,*) ' ******************************************'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *      FATAL ERROR SUB. IFSTOWAM*'
          WRITE (IU06,*) ' *      ==========================        *'
          WRITE (IU06,*) ' * WIND DATE IS LATER THAN EXPECTED       *'
          WRITE (IU06,*) ' * DATE FROM IFS IS CDTWIR = ', CDTWIR
          WRITE (IU06,*) ' * DATE EXPECTED IS CDTWIS = ', CDTWIS
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' *   PROGRAM ABORTS  PROGRAM ABORTS       *'
          WRITE (IU06,*) ' *                                        *'
          WRITE (IU06,*) ' ******************************************'
          CALL ABORT1
        ENDIF

        IPARAMU=165
        IPARAMV=166
        IPARAMCI=31

!*        1.1.2 DETERMINE CODE FOR WIND FIELD TYPE:
!               165 = U10 FIELD , 180 = STRESS FIELD.
!               ------------------------------------

        IF (IPARAMU.EQ.165 .AND. IPARAMV.EQ.166) THEN
          ICODE = 3
        ELSEIF (IPARAMU.EQ.33 .AND. IPARAMV.EQ.34) THEN
          ICODE = 3
        ELSEIF (IPARAMU.EQ.131 .AND. IPARAMV.EQ.132) THEN
          ICODE = 3
        ELSEIF (IPARAMU.EQ.180 .AND. IPARAMV.EQ.181) THEN
          ICODE = 2
        ELSE
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +      WARNING ERROR IN SUB. IFSTOWAM     +'
          WRITE(IU06,*) ' +      ==============================     +'
          WRITE(IU06,*) ' + SUSPICIOUS WIND FIELD TYPE IN MARS BLOCK+'
          WRITE(IU06,*) ' + U TYPE IS = ', IPARAMU
          WRITE(IU06,*) ' + V TYPE IS = ', IPARAMV
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +    PROGRAM CONTINUES                    +'
          WRITE(IU06,*) ' +                                         +'
          WRITE(IU06,*) ' +++++++++++++++++++++++++++++++++++++++++++'
          CALL ABORT1
        ENDIF

!*      1.2 GRID DEFINITION BLOCK.
!           ---------------------

!*        1.2.1 DETERMINE GRID PARAMETERS.
!               --------------------------

        KCOL = NGX
        KROW = NGY
        RLATN = AMONOP
        RLATS = AMOSOP
        RLONL = AMOWEP
        RLONR = AMOEAP
        IF(KCOL.NE.1) THEN
          DLAM  = (RLONR-RLONL)/REAL(KCOL-1)
        ELSE
          DLAM  = 0.
        ENDIF
        IF(KROW.NE.1) THEN
          DPHI  = (RLATN-RLATS)/REAL(KROW-1)
        ELSE
          DPHI  = 0.
        ENDIF

      ENDIF

      IF (LHOOK) CALL DR_HOOK('IFSTOWAM',1,ZHOOK_HANDLE)

      RETURN

      END SUBROUTINE IFSTOWAM
