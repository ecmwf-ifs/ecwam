      SUBROUTINE SECONDHH(NFRE,NANG,TH,COSTH,SINTH,FAK,A,B,FAC1,
     1                    FAC2,FAC3 )

!----------------------------------------------------------------

!**** *SECONDHH* - COMPUTATION OF SECOND ORDER HARMONICS.

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

!          NONE

!     REFERENCES.
!     -----------

!          V E ZAKHAROV(1967)

!-------------------------------------------------------------------

      DIMENSION TH(NANG),COSTH(NANG),SINTH(NANG),FAK(NFRE)
      DIMENSION A(NANG,NANG,NFRE,NFRE),B(NANG,NANG,NFRE,NFRE),
     &          FAC1(NANG,NANG,NFRE,NFRE),
     &          FAC2(NANG,NANG,NFRE,NFRE),
     &          FAC3(NANG,NANG,NFRE,NFRE)

!-----------------------------------------------------------------------

!*    1. INITIALISE RELEVANT QUANTITIES.

      DEL1 = 10.**(-8)
      PI   = 4.*ATAN(1.)
      ZPI  = 2.*PI
      G    = 9.806
      DELTH = ZPI/FLOAT(NANG)

      DO K1=1,NANG
        DO K2=1,NANG
          COSDIFF = COS(TH(K1)-TH(K2))
          DO  M1=1,NFRE
            XK1 = FAK(M1)
            DO M2=1,NFRE
              XK2 = FAK(M2)
              X12 = XK1*XK2*COSDIFF
              XK3 = XK1**2 + XK2**2 +2.*X12 +DEL1
              XK3 = SQRT(XK3)
              X13 = XK1**2+X12
              X32 = X12+XK2**2
              OM1 = SQRT(G*XK1)
              OM2 = SQRT(G*XK2)
              OM3 = SQRT(G*XK3)
              F1 = SQRT(XK1/(2.*OM1))
              F2 = SQRT(XK2/(2.*OM2))
              F3 = SQRT(XK3/(2.*OM3))
              VM = ZPI*VMIN(XK3,XK1,XK2,X13,X32,X12,OM3,OM1,OM2)
              VP = ZPI*VPLUS(-XK3,XK1,XK2,-X13,-X32,X12,OM3,OM1,OM2)
              DELOM1 = OM3-OM1-OM2+DEL1
              DELOM2 = OM3+OM1+OM2+DEL1
              A(K1,K2,M1,M2) = -F3/(F1*F2)*(VM/(DELOM1)+
     &                          VP/(DELOM2))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO K1=1,NANG
        DO K2=1,NANG
          COSDIFF = COS(TH(K1)-TH(K2))
          DO M1=1,NFRE
            XK1 = FAK(M1)
            DO M2=1,NFRE
              XK2 = FAK(M2)
              X12 = XK1*XK2*COSDIFF
              XK3 = XK1**2 + XK2**2 - 2.*X12 + DEL1
              XK3 = SQRT(XK3)
              X13 = XK1**2-X12
              X32 = X12-XK2**2
              OM1 = SQRT(G*XK1)
              OM2 = SQRT(G*XK2)
              OM3 = SQRT(G*XK3)+DEL1
              F1 = SQRT(XK1/(2.*OM1))
              F2 = SQRT(XK2/(2.*OM2))
              F3 = SQRT(ABS(XK3)/(2.*OM3))
              VM = ZPI*VMIN(XK1,XK3,XK2,X13,X12,X32,OM1,OM3,OM2)
              VP = ZPI*VMIN(XK2,-XK3,XK1,-X32,X12,-X13,OM2,OM3,OM1)
              DELOM321 = OM3+OM2-OM1+DEL1
              DELOM312 = OM3+OM1-OM2+DEL1
              B(K1,K2,M1,M2) = -F3/(F1*F2)*(VM/(DELOM321)+
     &                         VP/(DELOM312))
            ENDDO
          ENDDO
        ENDDO
      ENDDO

      DO K1=1,NANG
        DO K2=1,NANG
          DO M1=1,NFRE
            XK1 = FAK(M1)
            DO M2=1,NFRE
              XK2 = FAK(M2)
              C22 = A(K1,K2,M1,M2)+B(K1,K2,M1,M2)
              S22 = B(K1,K2,M1,M2)-A(K1,K2,M1,M2)
              FAC1(K1,K2,M1,M2) =
     &             (XK1**2*COSTH(K1)**2 + XK2**2*COSTH(K2)**2)*C22
     &             -XK1*XK2*COSTH(K1)*COSTH(K2)*S22
              FAC2(K1,K2,M1,M2) =
     &             (XK1**2*SINTH(K1)**2 + XK2**2*SINTH(K2)**2)*C22
     &             -XK1*XK2*SINTH(K1)*SINTH(K2)*S22
              FAC3(K1,K2,M1,M2) =
     &             (XK1**2*SINTH(K1)*COSTH(K1) +
     &              XK2**2*SINTH(K2)*COSTH(K2))*C22
     &             -XK1*XK2*COSTH(K1)*SINTH(K2)*S22
              A(K1,K2,M1,M2) = C22
              B(K1,K2,M1,M2) = S22
            ENDDO
          ENDDO
        ENDDO
      ENDDO


!     -----------------------------------------------------------------

      RETURN
      END SUBROUTINE SECONDHH 
