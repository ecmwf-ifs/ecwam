      PROGRAM cmpbls 

! ======================================================================

!**** cmp_bls   compares two blsp files.

!        B. HANSEN   *ECMWF*   17/05/93

!     Interface.
!     ----------
!         INPUT.
!         ------
!           *jpuin*     standard user input to read name of blsp 
!                       and grs file.
!           *jpubls1*    input unit opened for 1st blsp file.
!           *jpubls2*    input unit opened for 2nd blsp file.
!           *jpupar*    input unit for WAModel parameter file.

!         OUTPUT.
!         -------
!           *jpuso*     output unit opened for standard output. 

!     Method.
!     -------
!       NONE.

!     Modification.
!     -------------
!       NONE.

! ----------------------------------------------------------------------

      CHARACTER* 4 cpugrd
      CHARACTER* 6 cpupar
      PARAMETER (
     & jpuin  =         1, ! input unit for user information.
     & jpuso  =         6, ! output unit for printer output.
     & cpugrd =    'GRID', ! input unit for file GRID.
     & cpupar =  'PARWAM'  ! input unit for WAModel parameter file.
     &          )

      CHARACTER*256 clfile(2) ! For filenames.
      CHARACTER*  3 yostat    ! For the status of a file to be opened.

      INTEGER iubls(2)

      LOGICAL llexist

! ----------------------------------------------------------------------

!*       PARAMETER NAMELIST AS IN COMMON PARAM OF THE WAM-MODEL.
!        -------------------------------------------------------

      LOGICAL llper          !  the land sea mask is periodic (T)
      NAMELIST /PARWAM/
     &  iang, ifre, igx, igy, iblo, iiblo, iover, ioutp, ioutt,
     &  imaxc, imaxf, ibinp, iibl1, iibld, ibld, iiblc, iblc ,
     &  itaumax, iumax, iustar, ialpha, idepth, irefra, llper


! ----------------------------------------------------------------------

#include <outctrl.h>
      DATA csubna /"cmpbls "/

! ----------------------------------------------------------------------

!*    1. INITIALIZATIONS.
!        ----------------

      msub = 1
      mtest = 0
      msuvi = 0
      csuvi = ' '
      ctest = ' '

      INQUIRE (FILE='ts', EXIST=llexist)
      IF (llexist) THEN
        iuts   = i_get_unit (-1, 'ts',       's', 'f', 0)
        READ (iuts, ts)
      ENDIF

      CALL uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(jpuso,*)' ',csubna,': '

      iu_grd = i_get_unit (jpuso, cpugrd, 's', 'u', 0)
      iu_par = i_get_unit (jpuso, cpupar, 's', 'f', 0)
      READ(iu_par,parwam)
      IF ( lltest ) THEN
        WRITE(jpuso,'("READ PARWAM:")')
        WRITE(jpuso, PARWAM)
      ENDIF

! ----------------------------------------------------------------------

!* 2.0   READ FILENMAES AND OPEN.
!        ------------------------

      WRITE(jpuso,*)'  CMP_BLS:  '
      WRITE(jpuso,*)'  READS WITH TWO READ STATEMENTS FILENAMES'
      WRITE(jpuso,*)'  1. NAME OF THE FILE FOR BLOCKED SPECTRA '
      WRITE(jpuso,*)'  2. NAME OF THE FILE FOR BLOCKED SPECTRA '
      WRITE(jpuso,*)'  BOTH FILES MUST EXIST '
      WRITE(jpuso,*)'  CMP_BLS WILL COMPARE BOTH DATA FILES '

      DO ju = 1, 2
        READ (jpuin,'(a30)') clfile(ju)
        il1 = LEN_TRIM(clfile(ju))
        INQUIRE(FILE=clfile(ju)(:il1),EXIST=llexist)
        IF (llexist) THEN
          iubls(ju) = i_get_unit (jpuso, clfile(ju)(:il1), 's', 'u', 0)
        ELSE
          WRITE(jpuso,'(A50)')clfile(ju)(:il1)
          STOP ' nofile'
        ENDIF
      ENDDO

! ----------------------------------------------------------------------

!*    3.0 READ GRID INFORMATION FILE AND COMPARE BLS FILES.
!         -------------------------------------------------

      CALL BLSCOMP ( jpuso, iu_grd, iubls, clfile,
     &  iang, ifre, igx, igy, iblo, iiblo, iover, ioutp, ioutt,
     &  imaxc, imaxf, ibinp, iibl1, iibld, ibld, iiblc, iblc,
     &  itaumax, iumax, iustar, ialpha, idepth, irefra)
      END
      SUBROUTINE blscomp( kuso, kugrd, kubls, cdfile,
     &  NANG, NFRE, NGX, NGY, NBLO, NIBLO, NOVER, NGOUT, NOUTT,
     &  NMAXC, NMAXF, NBINP, NIBL1, NIBLD, NBLD, NIBLC, NBLC,
     &  ITAUMAX, JUMAX, IUSTAR, IALPHA, NDEPTH, IREFRA)

! ----------------------------------------------------------------------

!**** *BLSCOMP* CALLS MODIFIED READPRE IN ORDER TO INPUT THE GRID
!               INFORMATION FILE OF THE WAM MODEL AND INPUTS THEN
!               BLOCKED WAM MODEL SPECTRA FROM TWO FILES AND COMPARES 
!               THEM.

!      B. HANSEN         ECMWF          MAY      1993

!     PURPOSE.
!     --------

!          S.O.

!**   INTERFACE.
!     ----------
!          VARIABLES.   PURPOSE.
!          ----------   --------
!         I/  * kuso *  unit number for printer output.
!         I/  * kugrd*  input unit number for grid information file.
!         I/  * kubls*  input unit numbers for restart files blsp.
!         I/  * cdfile* For filenames.
!         the rest is input to this subroutine descritpion in the
!         wam model.

!     LIBRARIES.
!     ----------
!          ALTLIB.

!     EXTERNALS.
!     ----------
!          READPRE.

!     METHOD.
!     -------
!          NONE

!     REFERENCE.
!     ----------
!          NONE

! ----------------------------------------------------------------------

!* 0.    DEFINITIONS.
!  ------------------

      CHARACTER*  4 nn(3)
      CHARACTER* 20 yoname
      CHARACTER* 40 MSG
      CHARACTER*256 cdfile(2) ! For filenames.
      CHARACTER TITL*100

      INTEGER kubls(2)

      REAL DIFF(0:NIBLO), GDIFF(ngx,ngy)

! ----------------------------------------------------------------------

      PARAMETER (jplev=1)
      DIMENSION
     &     FR(NFRE), DFIM(NFRE), GOM(NFRE), C(NFRE),
     &     TH(NANG), COSTH(NANG), SINTH(NANG),
     &     SINPH(NGY), COSPH(NGY),
     &     IJS(NBLO), IJL2(NBLO), IJLS(NBLO),
     &     IJL(NBLO), IJLT(NBLO),
     &     IXLG(NIBLO,NBLO), KXLT(NIBLO,NBLO),
     &     IKP(NFRE+4), IKP1(NFRE+4),
     &     IKM(NFRE+4), IKM1(NFRE+4),
     &     K1W(NANG,2), K2W(NANG,2), K11W(NANG,2),
     &     K21W(NANG,2), AF11(NFRE+4), FKLAP(NFRE+4),
     &     FKLAP1(NFRE+4), FKLAM(NFRE+4), FKLAM1(NFRE+4)
      REAL
     &     FRH(30), IGAR(NGOUT), IJAR(NGOUT),
     &     DEPTH(NIBLO, NBLO),
     &     TCGOND(NDEPTH,NFRE), TFAK(NDEPTH,NFRE),
     &     TSIHKD(NDEPTH,NFRE), TAUT(0:ITAUMAX,0:JUMAX,jplev),
     &     TAUHFT(0:IUSTAR,0:IALPHA),
     &     U(0:NIBLC,NBLC), V(0:NIBLC,NBLC),
     &     FL3(0:NIBLO,NANG,NFRE), FL2(0:NIBLO,NANG,NFRE),
     &     FL(NANG,NFRE), NLONRGG(NGY), zdello(NGY), dellam(NGY),
     &     EMEAN3(NIBLO),EMEAN2(NIBLO),TEMP3(NIBLO),TEMP2(NIBLO)

      REAL X4(2)
      DATA X4 / 1.0, 1.0 /

! ----------------------------------------------------------------------

#include <outctrl.h>
      DATA csubna /"blscomp"/

! ----------------------------------------------------------------------

!*    1. INITIALIZATIONS.
!        ----------------

      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      IF (llsuvi .OR. lltest) WRITE(kuso,*)' ',csubna,': '

      CTEST(1)='       '
      ZG=-9999.
      NPRECR = KIND(X4)

! ----------------------------------------------------------------------

!* 2.    READ GRID INFORMATION FILE.
!  ---------------------------------


      CALL READPREB (kuso, kugrd, IREFRA,
     &     NANG, NFRE, NGY, NBLO, NIBLO, NGOUT, NIBLC, NBLC,
     &     NDEPTH, ITAUMAX, JUMAX, IUSTAR, IALPHA, jplev,
     &     FR, DFIM, GOM, C, DELTH, DELTR, TH, COSTH, SINTH,
     &     DELPHI, DELLAM, NLONRGG, SINPH, COSPH, IGL, IJS, IJL2, IJLS,
     &     IJL, IJLT, IXLG, KXLT, NX, NY, IPER,
     &     AMOWEP, AMOSOP, AMOEAP, AMONOP, XDELLA, XDELLO, ZDELLO,
     &     IRGG, IKP, IKP1, IKM, IKM1, K1W, K2W, K11W, K21W, AF11,
     &     FKLAP, FKLAP1, FKLAM, FKLAM1, ACL1, ACL2,  CL11,
     &     CL21, DAL1, DAL2, FRH,
     &     BETAMAX, ZALP, ALPHA, XKAPPA, XNLEV, TAUT, DELTAUW,
     &     DELU, TAUHFT, DELUST, DELALP, NGOUT, IGAR, IJAR,
     &     DEPTH, DEPTHA, DEPTHD, TCGOND, TFAK, TSIHKD, U, V)

! ----------------------------------------------------------------------

!* 3.    LOOP OVER BLOCKS.
!  -----------------------


      IOS = 0
      DIFF=0.0
      WRITE (kuso,*)
      nn(1)(1:2)='d1'  
      nn(2)(1:2)='d2'  
      nn(3)(1:2)='dd'  
      DO IG=1,NBLO
        maxblock=-1
        maxang  =-1
        maxfreq =-1
        WRITE (kuso,*) ' IJS / IJL ', IJS(IG), IJL(IG)
        WRITE (kuso,*) ' IJLS / IJL2 ', IJLS(IG), IJL2(IG)
        WRITE (kuso,*) ' IJLT        ', IJLT(IG)
        IF (IOS .EQ. 0) THEN
          READ(kubls(1),ERR=9000, IOSTAT=IOS)  FL3
          READ(kubls(2),ERR=9000, IOSTAT=IOS)  FL2
        ENDIF
 9000   IF (ios .NE. 0 ) THEN
          IF ( IG.EQ.1) THEN
            DO ju=1,2
              CLOSE (kubls(ju))
              ifl = LEN_TRIM(cdfile(ju))
              WRITE(kuso,'("blscomp: open ", a50)') cdfile(ju)(:ifl)
              CALL PBOPEN(kubls(ju),cdfile(ju)(:ifl),'r',IRET)
              IF (IRET.LT.0) THEN
                WRITE(kuso,'("blscomp: pbopen failed iret=", i3)') iret
                CALL abort
              ENDIF
            ENDDO
            IOUNT = (NIBLO+1)*NANG*NFRE*NPRECR
            WRITE(kuso,'("blscomp: IOUNT =", i8, " NIBLO=", i6, 
     &       " NANG=", i3, "NFRE=", i3, "NPRECR", i3)')
     &       IOUNT, NIBLO, NANG, NFRE, NPRECR
          ENDIF
          CALL PBREAD(kubls(1),FL3,IOUNT,IRET3)
          CALL PBREAD(kubls(2),FL2,IOUNT,IRET2)
          IF(IRET2.LT.0.OR.IRET3.LT.0) THEN
            WRITE (kuso,*) '*****************************************'
            WRITE (kuso,*) '*                                       *'
            WRITE (kuso,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
            WRITE (kuso,*) '*    IN READFL                          *'
            WRITE (kuso,*) '*    FILE 1 ',cdfile(1)
            WRITE (kuso,*) '*    FILE 2 ',cdfile(2)
            WRITE (kuso,*) '*                                       *'
            WRITE (kuso,*) '*    RET2 = ', IRET2
            WRITE (kuso,*) '*    RET3 = ', IRET3
            WRITE (kuso,*) '*                                       *'
            WRITE (kuso,*) '*****************************************'
            CALL ABORT
          ELSEIF(IRET2.NE.IOUNT .OR. IRET3.NE.IOUNT) THEN
            WRITE (kuso,*) '*****************************************'
            WRITE (kuso,*) '*                                       *'
            WRITE (kuso,*) '*    ERROR FOLLOWING CALL TO PBREAD     *'
            WRITE (kuso,*) '*    IN READFL                          *'
            WRITE (kuso,*) '*    FILE 1 ',cdfile(1)
            WRITE (kuso,*) '*    FILE 2 ',cdfile(2)
            WRITE (kuso,*) '*    NUMBER OF BYTES READ : ',IRET2, IRET3
            WRITE (kuso,*) '*    IS NOT EQUAL TO THE NUMBER OF BYTES*'
            WRITE (kuso,*) '*    REQUESTED : ',IOUNT
            WRITE (kuso,*) '*                                       *'
            WRITE (kuso,*) '*****************************************'
            CALL ABORT
          ELSE
            WRITE (kuso,*) '*READ FILE 1 ',cdfile(1)(:40), IRET2
            WRITE (kuso,*) '*READ FILE 2 ',cdfile(2)(:40), IRET3
          ENDIF
        ENDIF
        WRITE(nn(1)(3:4),'(i2.2)')IG
        WRITE(nn(2)(3:4),'(i2.2)')IG
        WRITE(nn(3)(3:4),'(i2.2)')IG


!          compute wave height

        WRITE (kuso,*) ' compute wave height'
        IF (llsuvi.OR.lltest) CALL flush (kuso)
        EMEAN3=0.
        EMEAN2=0.
        DO J3=1,NFRE
          TEMP3=0.
          TEMP2=0.
          DO J2=1,NANG
            DO J1=IJS(IG),IJL(IG)
              TEMP3(J1) = TEMP3(j1)+FL3(j1,j2,j3)
              TEMP2(J1) = TEMP2(j1)+FL2(j1,j2,j3)
            ENDDO
          ENDDO
          DO J1=IJS(IG),IJL(IG)
            EMEAN3(j1) = EMEAN3(j1)+DFIM(j3)*TEMP3(J1)
            EMEAN2(j1) = EMEAN2(j1)+DFIM(j3)*TEMP2(j1)
          ENDDO
        ENDDO
!        add tail
        DELT25 = FR(NFRE)/4.*DELTH
        DO J1=IJS(IG),IJL(IG)
          EMEAN3(j1) = EMEAN3(j1)+DELT25*TEMP3(j1)
          EMEAN3(j1) = 4.*SQRT(EMEAN3(j1))
          EMEAN2(j1) = EMEAN2(j1)+DELT25*TEMP2(j1)
          EMEAN2(j1) = 4.*SQRT(EMEAN2(j1))
        ENDDO

! find the largest difference 

        WRITE (kuso,*) ' find the largest difference' 
        delmax=-999999.
        delmin= 999999.
        DO J1=IJS(IG),IJL(IG)
          DIFF(J1) = EMEAN3(j1)-EMEAN2(j1)
          del=abs(DIFF(J1))
          IF ( delmax .LT. del ) maxblock=J1 
          IF ( delmin .GT. del ) minblock=J1 
          delmax=MAX(delmax,del)
          delmin=MIN(delmin,del)
        ENDDO

        j1=maxblock

        WRITE (kuso,*) ' at maxblock = ',maxblock
        WRITE (kuso,*) ' largest wave height difference' 
        WRITE (kuso,*) ' diff = ', EMEAN3(j1)-EMEAN2(j1) 
        WRITE (kuso,*) ' Hs1 = ', EMEAN3(j1)
        WRITE (kuso,*) ' Hs2 = ', EMEAN2(j1)
 
        IF (maxblock .EQ. 1 .AND. EMEAN3(j1)-EMEAN2(j1) .EQ. 0.0 ) THEN
          WRITE (kuso, *) 'FILES ARE IDENTICAL '
          RETURN
        ENDIF

        inc = MIN ( 6, maxblock)
        Do j1= maxblock-inc, maxblock+inc, 2

          WRITE (kuso,*) ' THIS BLOCK  = ', j1
          WRITE (kuso,*) ' wave height difference' 
          WRITE (kuso,*) ' diff = ', EMEAN3(j1)-EMEAN2(j1) 
          WRITE (kuso,*) ' Hs1 = ', EMEAN3(j1)
          WRITE (kuso,*) ' Hs2 = ', EMEAN2(j1)

          IOLG=IXLG(j1,ig)
          IOLT=KXLT(j1,ig)
          LAT= REAL(IOLT-1)*XDELLA+AMOSOP
          LON= MOD( REAL(IOLG-1)*ZDELLO(IOLT)+AMOWEP+720.,360.)
          WRITE (kuso,*) ' THIS IS AT NGY ', KXLT(j1,ig), 'LAT = ',LAT
          WRITE (kuso,*) ' THIS IS AT NGX ', IXLG(j1,ig), 'LON = ',LON
          WRITE (kuso,*) ' '
          WRITE (kuso,*) ' with spectral components : '
          WRITE (kuso,*) ' '
          WRITE (kuso,*) ' IG J1  J2  J3  FL3   FL2'
          DO J3=1,NFRE
            DO J2=1,NANG
              WRITE(kuso,'(i3, 1x, i6, i3,1x, i5, 2f17.10)')
     &         IG, J1, J2, J3, FL3(j1,j2,j3),  FL2(j1,j2,j3)
            ENDDO
          ENDDO

          DO J2=1,NANG
            DO J3=1,NFRE
              FL(j2,j3)=FL3(j1,j2,j3)-FL2(j1,j2,j3)
            ENDDO
          ENDDO
        ENDDO

        OPEN(50,FILE=nn(1))
        OPEN(51,FILE=nn(2))
        OPEN(52,FILE=nn(3))
        WRITE(50,21) ((FL2(maxblock,j2,j3),j2=1,NANG),j3=1,NFRE)
        WRITE(51,21) ((FL3(maxblock,j2,j3),j2=1,NANG),j3=1,NFRE)
        WRITE(52,21) ((FL(j2,j3),j2=1,NANG),j3=1,NFRE)
        CLOSE(50)
        CLOSE(51)
        CLOSE(52)

        WRITE (kuso,*)
        
        CALL makegridb (DIFF, GDIFF, IJS, IJL, IG,
     &   NGX, NGY, NIBLO, NBLO,
     &   IXLG, KXLT)

      ENDDO

      TITL=" !!!!! wave height difference in mm !!!!!"
      CALL OUTPP ('xxxxxxxxxx', kuso, ngx, ngy, ngx, ngy, TITL,1000.0,
     &             GDIFF, AMOWEP, AMOSOP, AMOEAP, AMONOP)

!  12 FORMAT("1", i5,2i4,2e20.10)
   21 FORMAT("1", e20.15)
      END SUBROUTINE blscomp


