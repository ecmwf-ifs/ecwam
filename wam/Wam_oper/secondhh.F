      SUBROUTINE SECONDHH

!----------------------------------------------------------------

!**** *SECONDHH* - COMPUTATION OF SECOND ORDER HARMONICS AND
!                  RELEVANT TABLES FOR THE ALTIMETER CORRECTIONS.

!     P.A.E.M. JANSSEN

!     PURPOSE.
!     ---------

!          COMPUTE SECOND HARMONICS

!**   INTERFACE.
!     ----------

!          *CALL* *SECONDHH*

!     METHOD.
!     -------

!          SEE REFERENCE.

!     EXTERNALS.
!     ----------

!         VMIN_D
!         VPLUS_D          

!     REFERENCES.
!     -----------

!          V E ZAKHAROV(1967)

!-------------------------------------------------------------------

      USE YOWFRED  , ONLY : FR       ,DELTH    ,TH       ,FRATIO   ,
     &            COSTH    ,SINTH 
      USE YOWTABL  , ONLY : NFREHF   ,FAC0     ,FAC1     ,FAC2     ,
     &            FAC3     ,FAK      ,FRHF     ,DFIMHF
      USE YOWPARAM , ONLY : NANG
      USE YOWPCONS , ONLY : G        ,ZPI

!-------------------------------------------------------------------
      REAL, PARAMETER :: DEL1=1.0E-8

      REAL, ALLOCATABLE :: B(:,:,:,:)

!-----------------------------------------------------------------------

!*    1. INITIALISE RELEVANT QUANTITIES.

      ALLOCATE(B(NANG,NANG,NFREHF,NFREHF))

      ALLOCATE(FAC0(NANG,NANG,NFREHF,NFREHF))
      ALLOCATE(FAC1(NANG,NANG,NFREHF,NFREHF))
      ALLOCATE(FAC2(NANG,NANG,NFREHF,NFREHF))
      ALLOCATE(FAC3(NANG,NANG,NFREHF,NFREHF))

      ALLOCATE(FAK(NFREHF))
      ALLOCATE(FRHF(NFREHF))
      ALLOCATE(DFIMHF(NFREHF))

      FRHF(1)  = FR(1)
      DO M=2,NFREHF
        FRHF(M) = FRATIO*FRHF(M-1)
      ENDDO

      DO M=1,NFREHF
         FAK(M) = (ZPI*FRHF(M))**2/G
      ENDDO

      CO1 = 0.5*(FRATIO-1.)*DELTH
      DFIMHF(1) = CO1*FRHF(1)
      DO M=2,NFREHF-1
         DFIMHF(M)=CO1*(FRHF(M)+FRHF(M-1))
      ENDDO
      DFIMHF(NFREHF)=CO1*FRHF(NFREHF-1)


      DO M2=1,NFREHF
        XK2 = FAK(M2)
        XK2SQ = FAK(M2)**2
        DO  M1=1,NFREHF
          XK1 = FAK(M1)
          XK1SQ = FAK(M1)**2
          DO K1=1,NANG
            DO K2=1,NANG
              COSDIFF = COS(TH(K1)-TH(K2))
              X12 = XK1*XK2*COSDIFF
              XK3 = XK1SQ + XK2SQ +2.*X12 +DEL1
              XK3 = SQRT(XK3)
              X13 = XK1SQ+X12
              X32 = X12+XK2SQ
              OM1 = SQRT(G*XK1)
              OM2 = SQRT(G*XK2)
              OM3 = SQRT(G*XK3)
              F1 = SQRT(XK1/(2.*OM1))
              F2 = SQRT(XK2/(2.*OM2))
              F3 = SQRT(XK3/(2.*OM3))
              VM = ZPI*VMIN_D(XK3,XK1,XK2,X13,X32,X12,OM3,OM1,OM2)
              VP = ZPI*VPLUS_D(-XK3,XK1,XK2,-X13,-X32,X12,OM3,OM1,OM2)
              DELOM1 = OM3-OM1-OM2+DEL1
              DELOM2 = OM3+OM1+OM2+DEL1
              FAC0(K1,K2,M1,M2) = -F3/(F1*F2)*(VM/(DELOM1)+
     &                            VP/(DELOM2))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO M2=1,NFREHF
        XK2 = FAK(M2)
        XK2SQ = FAK(M2)**2
        DO  M1=1,NFREHF
          XK1 = FAK(M1)
          XK1SQ = FAK(M1)**2
          DO K1=1,NANG
            DO K2=1,NANG
              COSDIFF = COS(TH(K1)-TH(K2))
              X12 = XK1*XK2*COSDIFF
              XK3 = XK1SQ + XK2SQ - 2.*X12 + DEL1
              XK3 = SQRT(XK3)
              X13 = XK1SQ-X12
              X32 = X12-XK2SQ
              OM1 = SQRT(G*XK1)
              OM2 = SQRT(G*XK2)
              OM3 = SQRT(G*XK3)+DEL1
              F1 = SQRT(XK1/(2.*OM1))
              F2 = SQRT(XK2/(2.*OM2))
              F3 = SQRT(ABS(XK3)/(2.*OM3))
              VM = ZPI*VMIN_D(XK1,XK3,XK2,X13,X12,X32,OM1,OM3,OM2)
              VP = ZPI*VMIN_D(XK2,-XK3,XK1,-X32,X12,-X13,OM2,OM3,OM1)
              DELOM321 = OM3+OM2-OM1+DEL1
              DELOM312 = OM3+OM1-OM2+DEL1
              B(K1,K2,M1,M2) = -F3/(F1*F2)*(VM/(DELOM321)+
     &                         VP/(DELOM312))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO M2=1,NFREHF
        XK2SQ = FAK(M2)**2
        DO M1=1,NFREHF
          XK1SQ = FAK(M1)**2
          DO K2=1,NANG
            DO K1=1,NANG
              C22 = FAC0(K1,K2,M1,M2)+B(K1,K2,M1,M2)
              S22 = B(K1,K2,M1,M2)-FAC0(K1,K2,M1,M2)
              FAC1(K1,K2,M1,M2) =
     &             (XK1SQ*COSTH(K1)**2 + XK2SQ*COSTH(K2)**2)*C22
     &             -FAK(M1)*FAK(M2)*COSTH(K1)*COSTH(K2)*S22
              FAC2(K1,K2,M1,M2) =
     &             (XK1SQ*SINTH(K1)**2 + XK2SQ*SINTH(K2)**2)*C22
     &             -FAK(M1)*FAK(M2)*SINTH(K1)*SINTH(K2)*S22
              FAC3(K1,K2,M1,M2) =
     &             (XK1SQ*SINTH(K1)*COSTH(K1) +
     &              XK2SQ*SINTH(K2)*COSTH(K2))*C22
     &             -FAK(M1)*FAK(M2)*COSTH(K1)*SINTH(K2)*S22
              FAC0(K1,K2,M1,M2) = C22
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(B)

!     -----------------------------------------------------------------

      RETURN
      END SUBROUTINE SECONDHH 
