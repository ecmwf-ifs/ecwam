      SUBROUTINE OUTERS (FL3, IG, IGL, CDTPRO, NSTART, NEND,
     &                   U10NEW,THWNEW,USNEW)
! -----------------------------------------------------------------     

!**** *OUTERS* -  OUTPUT OF SATELLITE COLOCATION SPECTRA.

!     H. GUNTHER     GKSS/ECMWF            JULY 1991                    
!     J. BIDLOT      ECMWF           FEBRARY 1996  MESSAGE PASSING
!     J. BIDLOT      ECMWF           JANUARY 1997  MODIFY DEFINITION
!                                                  OF TH0 

!*** INTERFACE.                                                         
!    ----------                                                         

!       *CALL  OUTERS (FL3, IG, IGL, CDTPRO, NSTART, NEND,
!                      U10NEW,THWNEW,USNEW)* 
!          *FL3*    -   SPECTRUM.                                       
!          *IG*     -   BLOCK NUMBER.                                   
!          *IGL*    -   NUMBER OF BLOCKS.                               
!          *CDTPRO* -   MODEL PROPAGATION TIME.                         
!          *NSTART*    INDEX OF THE FIRST POINT OF THE SUB GRID DOMAIN
!          *NEND*      INDEX OF THE LAST POINT OF THE SUB GRID DOMAIN
!          *U10NEW*    NEW WIND SPEED IN M/S.
!          *THWNEW*    WIND DIRECTION IN RADIANS IN OCEANOGRAPHIC
!                      NOTATION (POINTING ANGLE OF WIND VECTOR,
!                      CLOCKWISE FROM NORTH).
!          *USNEW*     NEW FRICTION VELOCITY IN M/S.

!    EXTERNALS.                                                         
!    ----------                                                         

!      ERSFILE   - FETCHES COLOCATION INFORMATION DONE BY PREUWA.       

!    METHOD.                                                            
!    -------                                                            

!      NONE.                                                            

!    REFERENCES.                                                        
!    -----------                                                        

!      NONE.                                                            

!---------------------------------------------------------------------- 

      USE YOWCOER  , ONLY : NERS     ,CDTERS   ,IERS     ,IJERS    ,
     &            IGERS
      USE YOWFRED  , ONLY : FR       ,TH       ,FRATIO
      USE YOWMAP   , ONLY : IXLG     ,KXLT     ,AMOWEP   ,AMOSOP   ,
     &            XDELLA   ,ZDELLO
      USE YOWMEAN  , ONLY : EMEAN    ,FMEAN    ,THQ
      USE YOWMPP   , ONLY : IRANK    ,NPROC    ,NINF     ,NSUP
      USE YOWMESPAS, ONLY : LMESSPASS
      USE YOWPARAM , ONLY : NANG     ,NFRE     ,NIBLO
      USE YOWPCONS , ONLY : DEG      ,ZMISS
      USE YOWTEST  , ONLY : IU06     ,ITEST
      USE MPL_MODULE

! ----------------------------------------------------------------------
!     ALLOCATABLE ARRAYS THAT ARE PASSED AS SUBROUTINE ARGUMENTS 

      INTEGER,DIMENSION(NPROC) :: NSTART,NEND
      REAL,DIMENSION(NINF-1:NSUP,NANG,NFRE) :: FL3
      REAL,DIMENSION(NINF:NSUP) :: U10NEW,THWNEW,USNEW
      CHARACTER CDTPRO*14

! ----------------------------------------------------------------------

!     WORK ARRAYS
      INTEGER,ALLOCATABLE :: IBUFF(:)
      REAL,ALLOCATABLE :: FLPTS(:,:,:,:)
      REAL,ALLOCATABLE,DIMENSION(:) :: EMPTS,FMPTS,THQPTS,U10PTS,THWPTS,
     &                                 USPTS

      IU91 = 91                                                         
      IU92 = 92                                                         

      IF (IG.EQ.1) THEN
        CALL ERSFILE (IU91, IU92, IG, IGL)

        IF(LMESSPASS.AND.NPROC.GT.1.AND.CDTPRO.EQ.CDTERS) THEN

!*  EXCHANGE THE VALUE OF NERS, IERS,IGERS, AND IJERS WITH THE OTHER PE

          ITAG=92
          ALLOCATE(IBUFF(2))

          IF(IRANK.EQ.1) THEN
            IBUFF(1) = NERS
            IBUFF(2) = IERS
          ENDIF
          CALL MPL_BROADCAST(IBUFF,KROOT=1,KTAG=ITAG,CDSTRING='OUTERS:')
          ITAG=ITAG+1
          IF(IRANK.NE.1) THEN
            NERS = IBUFF(1)
            IERS = IBUFF(2)
          ENDIF
          DEALLOCATE(IBUFF)

          IF(IERS.GT.0.) THEN
            KBUFRLENGTH=2*NERS
            ALLOCATE(IBUFF(KBUFRLENGTH))

            IF(IRANK.EQ.1) THEN
              KCOUNT=0
              DO NGOU=1,NERS
                KCOUNT=KCOUNT+1
                IBUFF(KCOUNT)=IGERS(NGOU)
              ENDDO
              DO NGOU=1,NERS
                KCOUNT=KCOUNT+1
                IBUFF(KCOUNT)=IJERS(NGOU)
              ENDDO
            ENDIF

            CALL MPL_BROADCAST(IBUFF(1:KBUFRLENGTH),KROOT=1,KTAG=ITAG,
     &       CDSTRING='OUTERS 1:')

            IF(IRANK.NE.1) THEN
              ALLOCATE(IGERS(NERS))
              ALLOCATE(IJERS(NERS))

              KCOUNT=0
              DO NGOU=1,NERS
                KCOUNT=KCOUNT+1
                IGERS(NGOU)=IBUFF(KCOUNT)
              ENDDO
              DO NGOU=1,NERS
                KCOUNT=KCOUNT+1
                IJERS(NGOU)=IBUFF(KCOUNT)
              ENDDO
            ENDIF

            DEALLOCATE(IBUFF)
          ENDIF

        ENDIF
      ENDIF

!*   1. LOOP OVER OUTPUT POINTS.                                        
!       ------------------------                                        

      IF (CDTPRO.EQ.CDTERS.AND.IERS.GT.0) THEN                          

        IF(LMESSPASS) THEN
!       COLLECT NECESSARY FIELDS TO PROCESS 1

          IRECV=1
          ITAG=ITAG+1
          NSPFLD=1
          NSCFLD=6

          ALLOCATE(EMPTS(IERS))
          ALLOCATE(FMPTS(IERS))
          ALLOCATE(THQPTS(IERS))
          ALLOCATE(U10PTS(IERS))
          ALLOCATE(THWPTS(IERS))
          ALLOCATE(USPTS(IERS))

          ALLOCATE(FLPTS(NSPFLD,IERS,NANG,NFRE))

          CALL MPGATHERERSFILE(IRECV,ITAG,NSTART,NEND,NSPFLD,NSCFLD,
     &                         FL3,FLPTS,
     &                         U10NEW,THWNEW,USNEW,
     &                         EMPTS,FMPTS,THQPTS,U10PTS,THWPTS,USPTS)
        ENDIF

        IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
          XANG = REAL(NANG)
          XFRE = REAL(NFRE)
          DO NGOU=1,IERS
            IF (IG.EQ.IGERS(NGOU)) THEN
              IJ = IJERS(NGOU)                                         
              XLON = AMOWEP + REAL(IXLG(IJ,IG)-1)*ZDELLO(KXLT(IJ,IG))
              XLAT = AMOSOP + REAL(KXLT(IJ,IG)-1)*XDELLA               

!*    1.1 WRITE INFORMATION TO FILE IU92.                               
!         -------------------------------                               

              IF (ITEST.GE.3) THEN
                WRITE(IU06,*)' SUB : outers  LOOP 1002 ngou=', ngou
                CALL FLUSH (IU06)
              ENDIF
              WRITE(IU92) XLON, XLAT, CDTPRO, XANG, XFRE,               
     &                    TH(1), FR(1), FRATIO
              IF(LMESSPASS) THEN
                IF (EMPTS(NGOU) .GT. 0.0 ) THEN
                  WRITE(IU92) 4.*SQRT(EMPTS(NGOU)), DEG*THQPTS(NGOU),
     &                      FMPTS(NGOU), USPTS(NGOU), DEG*THWPTS(NGOU),
     &                      U10PTS(NGOU)
                ELSE
                  WRITE(IU92) ZMISS, ZMISS, ZMISS, ZMISS, ZMISS, ZMISS
                ENDIF
                WRITE(IU92) ((FLPTS(1,NGOU,K,M),K=1,NANG),M=1,NFRE)
              ELSE
                IF (EMEAN(IJ) .GT. 0.0 ) THEN
                  WRITE(IU92) 4.*SQRT(EMEAN(IJ)), DEG*THQ(IJ),
     &                      FMEAN(IJ), USNEW(IJ), DEG*THWNEW(IJ),
     &                      U10NEW(IJ)
                ELSE
                  WRITE(IU92) ZMISS, ZMISS, ZMISS, ZMISS, ZMISS, ZMISS
                ENDIF
                WRITE(IU92) ((FL3(IJ,K,M),K=1,NANG),M=1,NFRE)
              ENDIF
            ENDIF                                                       
          ENDDO
        ENDIF
      ENDIF                                                             

      IF (CDTPRO.EQ.CDTERS) THEN                          
        IF (IG.EQ.IGL) THEN                                               
          IF((LMESSPASS.AND.(IRANK.EQ.1)).OR.(.NOT.LMESSPASS)) THEN
            CALL ERSFILE (IU91, IU92, -9, -9)
          ENDIF
        ENDIF                                                             
      ENDIF                                                             

      IF(ALLOCATED(IJERS))  DEALLOCATE (IJERS)
      IF(ALLOCATED(IGERS))  DEALLOCATE (IGERS)
      IF(ALLOCATED(EMPTS)) DEALLOCATE (EMPTS)
      IF(ALLOCATED(FMPTS)) DEALLOCATE (FMPTS)
      IF(ALLOCATED(THQPTS)) DEALLOCATE (THQPTS)
      IF(ALLOCATED(U10PTS)) DEALLOCATE (U10PTS)
      IF(ALLOCATED(THWPTS)) DEALLOCATE (THWPTS)
      IF(ALLOCATED(USPTS)) DEALLOCATE (USPTS)
      IF(ALLOCATED (FLPTS)) DEALLOCATE (FLPTS)

      RETURN                                                            
      END SUBROUTINE OUTERS
