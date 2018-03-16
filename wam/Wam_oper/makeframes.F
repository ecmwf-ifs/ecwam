      SUBROUTINE MAKEFRAMES( 
!                  DIMENSIONS
     &               NSPEC ,ISPECS ,ISPECE,  NGAP ,  NANG , NFRE,
!                  INPUT
     &               SPEC )

!-------------------------------------------------------------------

!     PURPOSE:
!     --------

!     BUILDS FRAMES AROUND GAPS IN THE TRANSFORMED SPECTRA.

!     AUTHOR.
!     -------
!     S.HASSELMANN, MPI HAMBURG, 1993.
!     J.WASZKEWITZ, MPI HAMBURG, 1993.

!     EXTERNALS:
!     ----------
!     *RESIZE_GAP_ARRAY  - RESIZE ARRAYS IF MGAP BECOMES TOO SMALL.

!-------------------------------------------------------------------

!     MODULE :
!     --------

      USE YOWGAP   , ONLY : GAPVAL   ,MGAP     ,LANG     ,HANG    ,
     &            LFRE     ,HFRE     ,EQUIGAP

!-------------------------------------------------------------------

      IMPLICIT NONE
#include "resize_gap_array.intfb.h"

!     INTERFACE:
!     ----------

      INTEGER :: NGAP(ISPECS:ISPECE)
      INTEGER NSPEC, ISPECS, ISPECE, NANG, NFRE
!                  ARRAY DIMENSIONS.
      REAL SPEC(0:NSPEC,NANG,NFRE)
!                   2D SPECTRUM.

!     LOCAL VARIABLES:
!     ----------------

      INTEGER, ALLOCATABLE ::  GAP(:,:,:)
!                   ASSIGNMENT OF A GAP NUMBER TO A 
!                   FREQUENCY-DIRECTIONAL BIN.
      INTEGER ISPEC, IGAP, JGAP, KGAP, IANG, IFRE
!                    LOOP INDICES.
      INTEGER IANGP1, IANGM1, IFREM1
!                    "IANG PLUS 1" , "IANG MINUS 1" , "IFRE PLUS 1".
      INTEGER LANG0, HANG0, LFRE0, HFRE0
!                    FOR TEMPORARILY STORAGE.
      INTEGER MAXNGAP
!                    THE MAXIMAL NUMBER OF GAPS OCCURING IN ONE 
!                    SPECTRUM.
      INTEGER GAP0
!                    TEMPORARILY USE.
      LOGICAL :: NEWGAP(ISPECS:ISPECE)
      LOGICAL WRAPI, WRAPJ
!                    FLAG INDICATING IF GAPI RESP. GAPJ WRAPP 
!                    AROUND 0 DEGREE ?

!-------------------------------------------------------------------

!     1. INITIALIZATION OF ARRAYS:
!     ----------------------------
      ALLOCATE(GAP(ISPECS:ISPECE,NANG,NFRE))

      DO IANG = 1,NANG
        DO IFRE = 1,NFRE
          DO ISPEC = ISPECS,ISPECE
            GAP(ISPEC,IANG,IFRE) = 0
          END DO
        END DO
      END DO
      DO ISPEC = ISPECS,ISPECE
        NGAP(ISPEC) = 0
      END DO
      DO IGAP = 1,MGAP
        DO JGAP = 1,MGAP
          DO ISPEC = ISPECS,ISPECE 
            EQUIGAP(ISPEC,IGAP,JGAP) = .FALSE.
          END DO
        END DO
      END DO

!------------------------------------------------------------------

!     2. FIND CORNERS OF FRAMES:
!     ---------------------------

      DO IANG = 1,NANG
        IF( IANG.LT.NANG )THEN
          IANGP1 = IANG + 1
        ELSE
          IANGP1 = 1
        END IF
        IF( IANG.GT.1 )THEN
          IANGM1 = IANG - 1
        ELSE
          IANGM1 = NANG
        END IF
        DO IFRE = 1,NFRE
          NEWGAP = .FALSE.
          DO ISPEC = ISPECS,ISPECE 
            IF( SPEC(ISPEC,IANG,IFRE).EQ.GAPVAL )THEN
              NEWGAP(ISPEC) = .TRUE.
              IF( IFRE.GT.1 )THEN
                IFREM1 = IFRE - 1
                IF( GAP(ISPEC,IANG,IFREM1).NE.0 )THEN
                  NEWGAP(ISPEC) = .FALSE.
                  GAP0 = GAP(ISPEC,IANG,IFREM1)
                  GAP(ISPEC,IANG,IFRE) = GAP0
                  IF( HFRE(ISPEC,GAP0).EQ.IFREM1 )
     &                HFRE(ISPEC,GAP0) = IFRE
                END IF
              END IF
              IF( GAP(ISPEC,IANGM1,IFRE).NE.0 )THEN
                NEWGAP(ISPEC) = .FALSE.
                IF( GAP(ISPEC,IANG,IFRE).EQ.0 )THEN
                  GAP0 = GAP(ISPEC,IANGM1,IFRE)
                  GAP(ISPEC,IANG,IFRE) = GAP0
                  IF( HANG(ISPEC,GAP0).EQ.IANGM1
     &                .AND.LANG(ISPEC,GAP0).NE.IANG )
     &                HANG(ISPEC,GAP0) = IANG
                ELSE
                  EQUIGAP(ISPEC,
     &                GAP(ISPEC,IANGM1,IFRE),GAP(ISPEC,IANG,IFRE))
     &                = .TRUE.
                END IF
              END IF
              IF( GAP(ISPEC,IANGP1,IFRE).NE.0 )THEN
                NEWGAP(ISPEC) = .FALSE.
                IF( GAP(ISPEC,IANG,IFRE).EQ.0 )THEN
                  GAP0 = GAP(ISPEC,IANGP1,IFRE)
                  GAP(ISPEC,IANG,IFRE) = GAP0
                  IF( LANG(ISPEC,GAP0).EQ.IANGP1
     &                .AND.HANG(ISPEC,GAP0).NE.IANG )
     &                LANG(ISPEC,GAP0) = IANG
                ELSE
                  EQUIGAP(ISPEC,
     &                GAP(ISPEC,IANGP1,IFRE),GAP(ISPEC,IANG,IFRE))
     &                = .TRUE.
                END IF
              END IF
            END IF
          END DO

          DO ISPEC = ISPECS,ISPECE 
            IF( NEWGAP(ISPEC) )THEN
              NGAP(ISPEC) = NGAP(ISPEC) + 1
              IF( NGAP(ISPEC).GT.MGAP ) THEN
                CALL RESIZE_GAP_ARRAY(ISPECS,ISPECE)
              ENDIF
            END IF
          END DO
          DO ISPEC = ISPECS,ISPECE 
            IF( NEWGAP(ISPEC) )THEN
              GAP(ISPEC,IANG,IFRE) = NGAP(ISPEC)
            END IF
          END DO
          DO ISPEC = ISPECS,ISPECE 
            IF( NEWGAP(ISPEC) )THEN
              GAP0 = NGAP(ISPEC)
              LANG(ISPEC,GAP0) = IANG
              HANG(ISPEC,GAP0) = IANG
              LFRE(ISPEC,GAP0) = IFRE
              HFRE(ISPEC,GAP0) = IFRE
            END IF
          END DO

        END DO
      END DO

      DEALLOCATE(GAP)

      MAXNGAP = 0
      DO ISPEC = ISPECS,ISPECE 
        MAXNGAP = MAX(MAXNGAP,NGAP(ISPEC))
      END DO

!------------------------------------------------------------------

!     2.1 COMBINE FRAMES THAT OVERLAP OR ARE CONTAINED IN ONE 
!         ANOTHER.
!     --------------------------------------------------------------

      DO IGAP = 1,MAXNGAP
        DO JGAP = 1,MAXNGAP

          DO ISPEC = ISPECS,ISPECE 
            IF( EQUIGAP(ISPEC,IGAP,JGAP).AND.IGAP.NE.JGAP )THEN
              IF( HANG(ISPEC,IGAP)+1.EQ.LANG(ISPEC,IGAP) )THEN
                LANG(ISPEC,IGAP) = 1
                HANG(ISPEC,IGAP) = NANG
              END IF
              IF( HANG(ISPEC,JGAP)+1.EQ.LANG(ISPEC,JGAP) )THEN
                LANG(ISPEC,JGAP) = 1
                HANG(ISPEC,JGAP) = NANG
              END IF
              WRAPI = LANG(ISPEC,IGAP).GT.HANG(ISPEC,IGAP)
              WRAPJ = LANG(ISPEC,JGAP).GT.HANG(ISPEC,JGAP)
              IF( .NOT.WRAPI.AND..NOT.WRAPJ )THEN
!               CASES:
!               |III     |  |IIIII   |  | IIII   |  |  IIII  |
!               |     JJJ|  |   JJJJJ|  |   JJJJ |  | JJJJJJ |

!               | IIIIII |  |   IIII |  |   IIIII|  |     III|
!               |  JJJJ  |  | JJJJ   |  |JJJJJ   |  |JJJ     |
                IF( MIN(HANG(ISPEC,IGAP),HANG(ISPEC,JGAP))+1 .LT.
     &              MAX(LANG(ISPEC,IGAP),LANG(ISPEC,JGAP)) )THEN
                  LANG0 = MAX( LANG(ISPEC,IGAP),LANG(ISPEC,JGAP) )
                  HANG0 = MIN( HANG(ISPEC,IGAP),HANG(ISPEC,JGAP) )
                ELSE
                  LANG0 = MIN(LANG(ISPEC,IGAP),LANG(ISPEC,JGAP))
                  HANG0 = MAX(HANG(ISPEC,IGAP),HANG(ISPEC,JGAP))
                END IF
              ELSE IF( WRAPI.AND..NOT.WRAPJ )THEN
!               CASES:
!               |IIIII  I|  |III    I|  |III  III|  |I    III|  |I   IIIII|
!               | JJJ    |  | JJJJ   |  | JJJJJJ |  |   JJJJ |  |     JJJ |
                IF( HANG(ISPEC,JGAP)+1.LT.LANG(ISPEC,IGAP) )THEN
                  LANG0 = LANG(ISPEC,IGAP)
                ELSE
                  LANG0 = MIN(LANG(ISPEC,IGAP),LANG(ISPEC,JGAP))
                END IF
                IF( LANG(ISPEC,JGAP)-1.GT.HANG(ISPEC,IGAP) )THEN
                  HANG0 = HANG(ISPEC,IGAP)
                ELSE
                  HANG0 = MAX(HANG(ISPEC,IGAP),HANG(ISPEC,JGAP))
                END IF
                IF( LANG0.LE.HANG0 )THEN
                  LANG0 = 1
                  HANG0 = NANG
                END IF
              ELSE IF( WRAPJ.AND..NOT.WRAPI )THEN
!               CASES:
!               |JJJJJ  J|  |JJJ    J|  |JJJ  JJJ|  |J    JJJ|  |J   JJJJJ|
!               | III    |  | IIII   |  | IIIIII |  |   IIII |  |     III |
                IF( HANG(ISPEC,IGAP)+1.LT.LANG(ISPEC,JGAP) )THEN
                  LANG0 = LANG(ISPEC,JGAP)
                ELSE
                  LANG0 = MIN(LANG(ISPEC,IGAP),LANG(ISPEC,JGAP))
                END IF
                IF( LANG(ISPEC,IGAP)-1.GT.HANG(ISPEC,JGAP) )THEN
                  HANG0 = HANG(ISPEC,JGAP)
                ELSE
                  HANG0 = MAX(HANG(ISPEC,IGAP),HANG(ISPEC,JGAP))
                END IF
                IF( LANG0.LE.HANG0 )THEN
                  LANG0 = 1
                  HANG0 = NANG
                END IF
              ELSE
!               CASES:
!               |IIIII  I|  |III    I|  |I      I|
!               |J  JJJJJ|  |J    JJJ|  |JJJ  JJJ|

!               |III  III|  |I    III|  |I  IIIII|
!               |J      J|  |JJJ    J|  |JJJJJ  J|
                LANG0 = MIN(LANG(ISPEC,IGAP),LANG(ISPEC,JGAP))
                HANG0 = MAX(HANG(ISPEC,IGAP),HANG(ISPEC,JGAP))
                IF( LANG0.LE.HANG0 )THEN
                  LANG0 = 1
                  HANG0 = NANG
                END IF
              END IF
              LFRE0 = MIN(LFRE(ISPEC,IGAP),LFRE(ISPEC,JGAP))
              HFRE0 = MAX(HFRE(ISPEC,IGAP),HFRE(ISPEC,JGAP))
              LANG(ISPEC,IGAP) = LANG0
              LANG(ISPEC,JGAP) = LANG0
              HANG(ISPEC,IGAP) = HANG0
              HANG(ISPEC,JGAP) = HANG0
              LFRE(ISPEC,IGAP) = LFRE0
              LFRE(ISPEC,JGAP) = LFRE0
              HFRE(ISPEC,IGAP) = HFRE0
              HFRE(ISPEC,JGAP) = HFRE0
            END IF
          END DO
        END DO
      END DO

!     2.2 REDUCE INDICES OF FRAMES :
!     ------------------------------

      DO ISPEC = ISPECS,ISPECE 
        DO IGAP = 1,NGAP(ISPEC)-1
          JGAP = IGAP+1
          DO WHILE( JGAP.LE.NGAP(ISPEC) )
            DO WHILE( (EQUIGAP(ISPEC,IGAP,JGAP)
     &          .OR.EQUIGAP(ISPEC,JGAP,IGAP))
     &          .AND. JGAP.LE.NGAP(ISPEC) )
!             COMBINE IGAP AND JGAP TO IGAB AND TRANSFER NGAP(ISPEC) TO JGAP
              LANG(ISPEC,JGAP) = LANG(ISPEC,NGAP(ISPEC))
              HANG(ISPEC,JGAP) = HANG(ISPEC,NGAP(ISPEC))
              LFRE(ISPEC,JGAP) = LFRE(ISPEC,NGAP(ISPEC))
              HFRE(ISPEC,JGAP) = HFRE(ISPEC,NGAP(ISPEC))
              DO KGAP = 1,NGAP(ISPEC)
                EQUIGAP(ISPEC,IGAP,KGAP) = EQUIGAP(ISPEC,IGAP,KGAP)
     &              .OR.EQUIGAP(ISPEC,JGAP,KGAP)
                EQUIGAP(ISPEC,JGAP,KGAP)
     &              = EQUIGAP(ISPEC,NGAP(ISPEC),KGAP)
              END DO
              DO KGAP = 1,NGAP(ISPEC)
                EQUIGAP(ISPEC,KGAP,IGAP) = EQUIGAP(ISPEC,KGAP,IGAP)
     &              .OR.EQUIGAP(ISPEC,KGAP,JGAP)
                EQUIGAP(ISPEC,KGAP,JGAP)
     &              = EQUIGAP(ISPEC,KGAP,NGAP(ISPEC))
              END DO
              NGAP(ISPEC) = NGAP(ISPEC) - 1
            END DO
            JGAP = JGAP + 1
          END DO
        END DO

!       2.3 COMBINE GAPS IF CORNERS COINCIDE.
!       -------------------------------------

        IGAP=1
        DO WHILE(IGAP.LE.NGAP(ISPEC))
         IF(LANG(ISPEC,IGAP).EQ.HANG(ISPEC,IGAP).OR.
     1      LFRE(ISPEC,IGAP).EQ.HFRE(ISPEC,IGAP)) THEN
          DO JGAP=igap,NGAP(ISPEC)-1
           LANG(ISPEC,JGAP) = LANG(ISPEC,JGAP+1)
           HANG(ISPEC,JGAP) = HANG(ISPEC,JGAP+1)
           LFRE(ISPEC,JGAP) = LFRE(ISPEC,JGAP+1)
           HFRE(ISPEC,JGAP) = HFRE(ISPEC,JGAP+1)
          END DO
           NGAP(ISPEC)=NGAP(ISPEC)-1
         ELSE
          IGAP=IGAP+1
         END IF
        END DO
      END DO

      RETURN

      END SUBROUTINE MAKEFRAMES
