      SUBROUTINE MPPEWITHINDIST(DISTMAX,INTLMAX,KMINLMAX,KMAXLMAX)

! ----------------------------------------------------------------------

!****  *MPPEWITHINDIST* -  DETERMINES A TABLE WHICH TELLS WHETHER
!                          A PE IS WITHIN A DISTANCE DISTMAX FROM 
!                          THE OTHERS.

!      J. BIDLOT       ECMWF     MAY 2002 

!     PURPOSE.                                                          
!     --------                                                          

!     DETERMINES A TABLE WHICH TELLS WHETHER A PE IS WITHIN A DISTANCE
!     DISTMAX FROM ANOTHER PE.

!*    INTERFACE.                                                        
!     ----------                                                        

!     *CALL* *MPPEWITHINDIST(DISTMAX,INTLMAX,KMINLMAX,KMAXLMAX)

!     *DISTMAX*  REAL     DISTANCE IN RADIAN !!!!
!     *INTLMAX*  INTEGER  TABLE INDICATING WHETHER A PE IS WITHIN 
!                         DISTMAX OF ANOTHER PE.
!     *KMINLMAX* INTEGER  SOUTHERN LATITUDE INDEX OF THE MOST SOUTHERN
!                         LATITUDE DISTMAX AWAY FROM EACH PE
!     *KMAXLMAX* INTEGER  NORTHERN LATITUDE INDEX OF THE MOST NORTHERN
!                         LATITUDE DISTMAX AWAY FROM EACH PE

!     METHOD.                                                           
!     -------                                                           

!     SEARCH THE DISTANCE BETWEEN GRID POINTS FROM ONE PE AND THE OTHER
!     ONES UNTIL THE DISTANCE IS LESS THE DISTMAX OR THERE IS NO POINT
!     LEFT.

!     EXTERNALS.                                                        
!     ----------                                                        

!     REFERENCES.                                                       
!     -----------                                                       

! ----------------------------------------------------------------------

      USE YOWGRID  , ONLY : SINPH    ,COSPH
      USE YOWMAP   , ONLY : XDELLA   ,IXLG    ,KXLT    ,ZDELLO   ,
     &              LLBOUND
      USE YOWMPP   , ONLY : NPROC
      USE YOWPARAM , ONLY : NGY      ,NIBLO   ,LL1D
      USE YOWPCONS , ONLY : DEG      ,RAD
      USE YOWSPEC,   ONLY : NSTART   ,NEND
      USE YOWTEST  , ONLY : IU06     ,ITEST

!-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: INTLMAX(NPROC,NPROC)
      INTEGER :: KMINLMAX(NPROC),KMAXLMAX(NPROC)
      INTEGER :: IJ, IJJ, IR, JR, KI, KJ, IC
      INTEGER :: LMAX

      REAL :: DISTMAX, COSDISTMAX
      REAL :: COSLON, DIS2 
      REAL, ALLOCATABLE :: XLONG(:)
      REAL, ALLOCATABLE,DIMENSION(:,:) :: CPH2, SPH2

      LMAX = NINT(DEG*DISTMAX/XDELLA)
!     add a small number to DISTMAX as a safety precaution
!     i.e each PE will get slightly more data that they should
!     strictly have.
      COSDISTMAX=COS(DISTMAX+0.1)

      IF (ITEST.GE.2) THEN
        WRITE(IU06,*) '   SUB. MPPEWITHINDIST: COMPUTING INTLMAX'
      ENDIF


      IF(LL1D) THEN
!       FOR THE SPECIAL CASE OF 1-D DECOMPOSITION
        DO IR = 1, NPROC
          IJ=NSTART(IR)
          KMINLMAX(IR) = MAX(1,KXLT(IJ,1)-LMAX)

          IJ=NEND(IR)
          KMAXLMAX(IR) = MIN(NGY,KXLT(IJ,1)+LMAX)

          DO JR = 1, IR-1
            INTLMAX(JR,IR)=0

            IJ=NSTART(JR)
            IF(KXLT(IJ,1) .GE. KMINLMAX(IR) .AND.
     &         KXLT(IJ,1) .LE. KMAXLMAX(IR) ) THEN
              INTLMAX(JR,IR)=1
            ENDIF

            IJ=NEND(JR)
            IF(KXLT(IJ,1) .GE. KMINLMAX(IR) .AND.
     &         KXLT(IJ,1) .LE. KMAXLMAX(IR) ) THEN
              INTLMAX(JR,IR)=1
            ENDIF

            INTLMAX(IR,JR) = INTLMAX(JR,IR)
          ENDDO
          INTLMAX(IR,IR)=1
        ENDDO

      ELSE
!       FOR THE MORE GENERAL CASE OF 2-D DECOMPOSITION
        DO IR = 1, NPROC
          KMINLMAX(IR) = NGY
          KMAXLMAX(IR) = 1
          DO IJ=NSTART(IR),NEND(IR)
            KMINLMAX(IR) = MIN(KMINLMAX(IR),KXLT(IJ,1)-LMAX)
            KMAXLMAX(IR) = MAX(KMAXLMAX(IR),KXLT(IJ,1)+LMAX)
          ENDDO
          KMINLMAX(IR) = MAX(1,KMINLMAX(IR))
          KMAXLMAX(IR) = MIN(NGY,KMAXLMAX(IR))
        ENDDO

        ALLOCATE(XLONG(NSTART(1):NEND(NPROC)))
        DO IR = 1, NPROC
          DO IJ=NSTART(IR),NEND(IR)
             XLONG(IJ) = REAL(IXLG(IJ,1)-1)*ZDELLO(KXLT(IJ,1))*RAD
          ENDDO
        ENDDO

        ALLOCATE(CPH2(NGY,NGY))
        ALLOCATE(SPH2(NGY,NGY))
        DO KJ=1,NGY
          DO KI=1,NGY
            CPH2(KI,KJ) = COSPH(KJ)*COSPH(KI)
            SPH2(KI,KJ) = SINPH(KJ)*SINPH(KI)
          ENDDO
        ENDDO


        DO IR = 1, NPROC
          DO JR = 1, NPROC
            INTLMAX(JR,IR)=0
          ENDDO
        ENDDO

        DO IR = 1, NPROC
          DO JR = 1, IR-1

            OUTER : DO IJJ=NSTART(JR),NEND(JR)
!             limit search to boundary points
              IF(LLBOUND(IJJ)) THEN
                KJ = KXLT(IJJ,1)
                IF(KJ.GE.KMINLMAX(IR) .AND. KJ.LE.KMAXLMAX(IR) ) THEN
                  DO IJ=NEND(IR), NSTART(IR), -1
!                   limit search to boundary points
                    IF(LLBOUND(IJ)) THEN
                      KI = KXLT(IJ,1)
                      COSLON = COS(XLONG(IJJ)-XLONG(IJ))
                      DIS2 = COSLON*CPH2(KJ,KI)+SPH2(KJ,KI)

                      IF(DIS2.GE.COSDISTMAX) THEN
                        INTLMAX(JR,IR)=1
                        EXIT OUTER 
                      ENDIF
                    ENDIF
                  ENDDO
                ENDIF 
              ENDIF 
            ENDDO OUTER 

            INTLMAX(IR,JR) = INTLMAX(JR,IR)

          ENDDO

          INTLMAX(IR,IR)=1

        ENDDO

        DEALLOCATE(XLONG)
        DEALLOCATE(CPH2)
        DEALLOCATE(SPH2)

        IF(ITEST.GT.0) THEN
          WRITE(IU06,*) '  IN MPPEWITHINDIST: INTLMAX '
          DO IR = 1, NPROC
            DO JR = 1, NPROC
              WRITE(IU06,*) IR,JR,INTLMAX(JR,IR)
            ENDDO
          ENDDO
          CALL FLUSH(IU06)
        ENDIF

      ENDIF
          
      RETURN 
      END SUBROUTINE MPPEWITHINDIST
