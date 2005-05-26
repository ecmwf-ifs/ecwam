! ======================================================================

      SUBROUTINE blspcon( kuso, kugrd, kubls, kugrs, cdtask,
     &  kang, kfre, kgx, kgy, kblo, kiblo, kover, koutp, koutt,
     &  kmaxc, kmaxf, kbinp, kibl1, kibld, kbld, kiblc, kblc,
     &  ktaumax, kumax, kustar, kalpha, kdepth, krefra,
     &  kyear, kmonth, kday, khour, ldoutspp, kspecs, ldlog, ldfdb,
     &  kfcst, cdclass, cdexpver, cddomain,  ksect1_3 )

! ======================================================================

!**** *BLSPCON* CALLS MODIFIED READPRE IN ORDER TO INPUT THE GRID
!               INFORMATION FILE OF THE WAModel AND INPUTS THEN
!               BLOCKED WAModel SPECTRA DISTRIBUTE THEM ON A
!               GRID AND WRITES THIS GRID TO THE OUTPUT FILE
!               ON UNIT KUGRS..
!CCC
!               it's coded the old way, ie for parameter 250
!CCCC

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
!         I/  * kubls*  input unit number for restart file blsp.
!         I/  * kugrs*  output unit number for the gridded spctra.
!         the rest is input to this subroutine descritpion in the
!         WAModel.

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

      USE YOWPCONS , ONLY : ZMISS


!* 0.    DEFINITIONS.
!  ------------------

      REAL, PARAMETER :: PPMISS = -3.0
!     REAL, PARAMETER :: PPMISS = -10.0

!     PPMISS: ALL SPECTRAL VALUES LESS OR EQUAL PPMISS ARE REPLACED
!             BY THE MISSING DATA INDICATOR FOR THE SECONDARY BIT 
!             MASK AND HAVE TO BE SET TO ZERE AFTER RETORING FROM
!             GRIB.

      INTEGER ix, iy
      REAL    zx, zy
      REAL zlatitude, zlongitude
      NAMELIST /location/  zlatitude, zlongitude

      CHARACTER*  1 clfdbsf, cldomain
      CHARACTER*  2 cdclass, cltype
      CHARACTER*  4 cdexpver
      CHARACTER*  5 cdtask
      CHARACTER*  6 cddomain
      CHARACTER* 12 cldate
      CHARACTER* 20 clname
      CHARACTER* 40 clmsg
      CHARACTER*100 cltitl
      LOGICAL ldlog, ldoutspp, ldfdb

      CHARACTER* 28 clct
      DATA clct / '__COMPILED_ON__' /

! ----------------------------------------------------------------------

      INTEGER, PARAMETER :: jplev = 2
      DIMENSION
     &     fr(kfre), dfim(kfre), gom(kfre), c(kfre),
     &     th(kang), costh(kang), sinth(kang),
     &     sinph(kgy), cosph(kgy),
     &     ijs(kblo), ijl2(kblo), ijls(kblo),
     &     ijl(kblo), ijlt(kblo),
     &     ixlg(kiblo,kblo), kxlt(kiblo,kblo),
     &     ikp(kfre+4), ikp1(kfre+4),
     &     ikm(kfre+4), ikm1(kfre+4),
     &     k1w(kang,2), k2w(kang,2), k11w(kang,2),
     &     k21w(kang,2), af11(kfre+4), fklap(kfre+4),
     &     fklap1(kfre+4), fklam(kfre+4), fklam1(kfre+4),
     &     frh(30), igar(koutp), ijar(koutp),
     &     depth(kiblo, kblo),
     &     tcgond(kdepth,kfre), tfak(kdepth,kfre),
     &     tsihkd(kdepth,kfre),
     &     taut(0:ktaumax,0:kumax,jplev),
     &     tauhft(0:kustar,0:kalpha),
     &     u(0:kiblc,kblc), v(0:kiblc,kblc),
     &     pdello(kgy), dellam(kgy), ilonrgg(kgy), kdomrgg(kgy),
     &     fl3(0:kiblo,kang,kfre)

      INTEGER ihlp(kgy)

      INTEGER ilon(kgy)

      REAL zspectrm(kang,kfre)
!     REAL zsec4(kgx*kgy*kang*kfre)
      REAL, ALLOCATABLE :: zsec4(:)
      REAL X4(2)

      CHARACTER*1 cl_bitmap(kgx,kgy)

!     FOR TIMING:
      REAL wam_user_clock,
     &     time0,           ! start time (reusable).
     &     time_m           ! time used to put fl3 into psec4.
     &     time_r           ! time used for reading.

! ----------------------------------------------------------------------

!* 0.    COMMON OUTPUT CONTROL.
!  ----------------------------


#include <outctrl.h>
      DATA csubna /"BLSPCON"/

! ----------------------------------------------------------------------

!*    1. INITIALIZATIONS.
!        ----------------
      X4=1.
      jpidl = KIND(X4)
      msub = msub + 1
      CALL  uraopc (csubna, msub, llsuvi, lltest)
      WRITE(kuso,*)' ',csubna,'        ', clct

      IF (lltest)  THEN
        WRITE(kuso,*)'  llsuvi: ',llsuvi,' lltest: ',lltest,
     &               ' ldoutspp: ',ldoutspp,' ldlog: ',ldlog,
     &               ' ldfdb: ',ldfdb
        WRITE(kuso,'(10a7)') ctest
      ENDIF

      cl_bitmap = 'l'
      clfdbsf = ' '
      itest = 0

! ----------------------------------------------------------------------

!* 2.    READ GRID INFORMATION FILE.
!  ---------------------------------

      CALL READPREB (kuso, kugrd, krefra,
     &     kang, kfre, kgy, kblo, kiblo, koutp, kiblc, kblc,
     &     kdepth, ktaumax, kumax, kustar, kalpha, jplev,
     &     fr, dfim, gom, c, delth, deltr, th, costh, sinth,
     &     delphi, dellam, ilonrgg, sinph, cosph, igl, ijs, ijl2, ijls,
     &     ijl, ijlt, ixlg, kxlt, nx, ny, iper,
     &     zmowep, zmosop, zmoeap, zmonop, zdella, zdello, pdello,
     &     irgg, ikp, ikp1, ikm, ikm1, k1w, k2w, k11w, k21w, af11,
     &     fklap, fklap1, fklam, fklam1, acl1, acl2,  cl11,
     &     cl21, dal1, dal2, frh,
     &     betamax, zalp, alpha, xkappa, xnlev, taut, deltauw,
     &     delu, tauhft, delust, delalp, ngout, igar, ijar,
     &     depth, deptha, depthd, tcgond, tfak, tsihkd, u, v)

      IF (irgg.EQ.1) THEN
        idl=0
        DO jgy = 1, kgy
          idl = ilonrgg(jgy) + idl
        ENDDO
        idl = idl * kang * kfre + 1
      ELSE
        idl = kgy * kgx * kfre * kang
      ENDIF

      ALLOCATE ( zsec4(idl) )

! ----------------------------------------------------------------------

!* 3.    LOOP OVER BLOCKS.
!  -----------------------

      pprec=0.0  ! -3.
      ppeps=1.0e-10
      waveht=   2.
      wavedir=  3.
      wavemfr=  4.
      windspd=  5.
      winddir=  6.
      wind10=   7.
      clname= ' LOC '

      zsec4=zmiss

      IF (llsuvi) WRITE(kuso,*)' WHICH TASK (task):',cdtask
      IF (ldoutspp) iu_locs =  i_get_unit (kuso, 'locs', 's', 'f', 0)

      IF ( cdtask .EQ. 'GRID ') THEN
        WRITE (cldate,'(I4.4, 4I2.2)') kyear, kmonth, kday, khour, 0
        DO jblock=1,kblo
          IF (llsuvi) time0=-wam_user_clock()
          READ (kubls) (((fl3(j1,j2,j3),
     &                                  j1=0,kiblo),
     &                                  j2=1,kang),
     &                                  j3=1,kfre)


          IF (llsuvi) THEN
            time_r=time0+wam_user_clock()
            time_r=time_r*1E-06
            WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
            WRITE (kuso,*) ' + TIME USED FOR READING FL3  +'
            WRITE (kuso,*) ' +                            +'
            WRITE (kuso,*) ' + ', time_r
            WRITE (kuso,*) ' +                            +'
            WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
          ENDIF

          IF (lltest) THEN
            DO jbl=100, 110
              WRITE(kuso,*)' SPECTRUM ',jbl
              WRITE(kuso,'(12(x,f7.3))')
     &        ((FL3(jbl,js1,js2),js1=1,kang),js2=1,kfre)
              WRITE(kuso,*)' '
            ENDDO
          ENDIF
          IF (ldlog) THEN
            IF (llsuvi) WRITE(kuso,*)' SUB: ',csubna,' DO LOG10'
            DO j3=1,kfre
              DO j2=1,kang
                DO j1=1,kiblo
                  fl3(j1,j2,j3)=log10(fl3(j1,j2,j3)+ppeps)+abs(pprec)
                  IF ( fl3(j1,j2,j3) .LE. ppmiss )
     &              fl3(j1,j2,j3)=zmiss
!---------->        fl3(j1,j2,j3)=MAX(fl3(j1,j2,j3),0.)
                ENDDO
              ENDDO
            ENDDO
          ENDIF
          IF (llsuvi) time0=-wam_user_clock()
          CALL MAKEGRIDS (kuso, fl3, zsec4, idl,
     &                    ijs(jblock), ijl(jblock), jblock,
     &                    ilonrgg, kgx, kgy, kang, kfre, ixlg, kxlt,
     &                    kblo, kiblo )
          IF (llsuvi) THEN
            time_m=time0+wam_user_clock()
            time_m=time_m*1E-06
            WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
            WRITE (kuso,*) ' + TIME USED FOR MAKEGRIDS    +'
            WRITE (kuso,*) ' +                            +'
            WRITE (kuso,*) ' + ', time_m
            WRITE (kuso,*) ' +                            +'
            WRITE (kuso,*) ' ++++++++++++++++++++++++++++++'
          ENDIF
        ENDDO

!***********************************************************************

        IF (lltest) THEN
          nx=1
          ny=1
          zem=ZMISS*kfre*kang
          DO J1=1, kgx*kgy   ! *kfre*kang
            sum=0.0
            js = ((j1-1)*kfre*kang)+1
            je =  (j1-1)*kfre*kang +kfre*kang
            DO j2=js,je
              sum = sum + zsec4(j2)
            ENDDO
            IF ( sum.GT.zem ) cl_bitmap(nx,ny)='.'
            nx=nx+1
            IF (nx.GT.ilonrgg(ny)) THEN
              ny=ny+1
              IF (ny .GT. kgy) EXIT
              nx=1
            ENDIF
          ENDDO
   
          DO j1=1,kgy
            WRITE(kuso,'(132a1)') 
     &       (cl_bitmap(j2,j1),j2=1,MIN(132,ilonrgg(j1)))
          ENDDO
          IF (kgx.gt.132) THEN
            WRITE(kuso,'(/132a1)') '-'
            DO j1=1,kgy
              WRITE(kuso,'(132a1)') (cl_bitmap(j2,j1),j2=133,kgx)
            ENDDO
          ENDIF
        ENDIF

!***********************************************************************

        IF (llsuvi ) THEN
          WRITE(kuso,*)' ALL BLOCKS DONE, jbl=',jblock
          WRITE(kuso,*)' ANY SPECTRA TO PRINT, ldoutspp=',ldoutspp      
        ENDIF
!       IF (lltest) WRITE(99) zsec4

        iparam = 250
 
        IF (llsuvi .OR. ldoutspp) THEN
          WRITE(kuso,*)' NUMBER OF SPECTRA TO PRINT, kspecs=',kspecs

          DO j1=1,kspecs
            READ ( iu_locs, location,END=3800)
            zy=zlatitude
            zx=zlongitude
            WRITE(kuso,location)

            iy = NINT((zlatitude-zmonop)/(-zdella)) + 1

            IF ( iy.GE.1 .AND. iy.LE.kgy ) THEN
              zx = MOD(zlongitude-zmowep+720.,360.)
              zx = zx / pdello(iy)
              ix = NINT(zx) +1
            ENDIF

            IF ( ix.LT.1 .OR. ix.GT.ilonrgg(iy) .OR.
     &           iy.LT.1 .OR. iy.GT.kgy ) THEN
              WRITE(kuso,*)zlatitude,'X',zlongitude,' IS OUTSIDE GRID'
              WRITE(kuso,*) ' IX=',ix,' KGX=',ilonrgg(iy)
              WRITE(kuso,*) ' IY=',iy,' KGY=',kgy
              WRITE(kuso,*) ' ZMOWEP=',zmowep
              WRITE(kuso,*) ' ZDELLO=',zdello
              WRITE(kuso,*) ' ZMONOP=',zmonop
              WRITE(kuso,*) ' ZDELLA=',zdella
            ELSE
              WRITE( kuso,*)' PRINT SPECTRUM AT IX=',ix,
     &          ' AND AT IY=',iy
              i_off = ispoff(ilonrgg, kgy, kang, kfre, iy-1, ix)
              IF (llsuvi) WRITE(kuso,*) ' OFFSET= ',i_off
              DO j3=1,kfre
                DO j2=1,kang
                  jx=(j3-1)*kang+1+(j2-1)+i_off

                  IF (ZSEC4(jx) .NE. ZMISS) THEN
                    zspectrm(j2,j3) = 10.**(zsec4(jx)-ABS(pprec))- ppeps
                  ELSE
                    zspectrm(j2,j3) = 0.0
                  ENDIF
                ENDDO
              ENDDO
              CALL PRSPP (kuso, cldate, zlongitude, zlatitude, clname,
     &                    waveht, wavedir, wavemfr, windspd, winddir,
     &                    wind10,
     &                    fr, th, zspectrm, kang, kfre, kang, kfre)
            ENDIF
          ENDDO
        ENDIF
 3800   CONTINUE
        IF ( ldoutspp) THEN
          WRITE(kuso,*) 'END OF FILE "LOCS" REACHED'
          WRITE(kuso,*) 'LAST RECORD READ:'
          WRITE(kuso,location)
        ENDIF

!* PACK FULL FIELD.

        istart=1
        iend  =kgy

        IF (cddomain .EQ. 'medite' .OR. cddomain .EQ. 'global' ) THEN

          itgrib=0
          DO ji = 1, kgy
            itgrib = itgrib+ilonrgg(ji)
          ENDDO
          cldomain = cddomain(1:1)
          CALL gribpacs
     &       ( kuso, kugrs, zsec4, itgrib, kgx, kgy, kang, kfre,
     &         irgg, th, fr, iparam, kyear, kmonth, kday, khour, kfcst,
     &         zmowep, zmosop, zmoeap, zmonop, ilonrgg, istart, iend,
     &         zdella, zdello, zmiss, cdclass, cdexpver, cldomain,ldfdb,
     &         ksect1_3, ifdb, jpidl)

        ELSE
!         IF ( ksect1_3 .GT. 300 .AND. ksect1_3 .LT. 400 .AND.
!    .         zmosop .LT. 0.0 .AND. zmonop .GE. 0.0 ) THEN
          IF ((cddomain(1:5) .EQ. 'north' .OR.
     &         cddomain(1:5) .EQ. 'south' .OR.
     &         cddomain(1:4) .EQ. 'both' ) .AND.
     &         zmosop .LT. 0.0 .AND. zmonop .GE. 0.0 ) THEN
 
            zlatitude=0.0
            zlongitude=zmoeap
            iy = NINT((zlatitude-zmonop)/(-zdella)) + 1
               
            zx = MOD(zlongitude-zmowep+720.,360.)
            zx = zx / pdello(iy)
            ix = NINT(zx) +1

            itgrib=0
            DO ji = 1, iy
              itgrib = itgrib+ilonrgg(ji)
            ENDDO
            itsp = ispoff(ilonrgg, kgy, kang, kfre, iy-1, ix) +
     &             (kang*kfre)

            IF (llsuvi) WRITE(*,*) "  ITSP  ", itsp

!* PACK NORTHERN HEMISPHERE.

            IF (cddomain(1:5) .EQ. 'north' .OR.
     &          cddomain(1:4) .EQ. 'both' ) THEN
!             isect1_3 = 106 ! 304
              zsouth = 0.0
              iend = iy
              cldomain = 'n'
              CALL gribpacs
     &       ( kuso, kugrs, zsec4, itgrib, kgx, kgy, kang, kfre,
     &         irgg, th, fr, iparam, kyear, kmonth, kday, khour, kfcst,
     &         zmowep, zsouth, zmoeap, zmonop, ilonrgg, istart, iend,
     &         zdella, zdello, zmiss, cdclass, cdexpver,cldomain, ldfdb,
     &         ksect1_3, ifdb, jpidl)
            ENDIF
 
!* PACK SOUTHERN HEMISPHERE.

            IF (cddomain(1:5) .EQ. 'south' .OR.
     &          cddomain(1:4) .EQ. 'both' ) THEN
!             isect1_3 =  106 ! 305
              znorth = 0.0 - zdella
              istrt = 1+itsp
              itsp  = itgrib-itsp
              istart = iy + 1
              iend = kgy

              itgrib=0
              DO ji = iy+1, kgy
                itgrib = itgrib+ilonrgg(ji)
              ENDDO
              IF (llsuvi) WRITE(kuso,*) " SOUTH: istrt=", istrt
              cldomain = 's'
              CALL gribpacs
     &       ( kuso, kugrs, zsec4(istrt), itgrib, kgx, kgy, kang, kfre,
     &         irgg, th, fr, iparam, kyear, kmonth, kday, khour, kfcst,
     &         zmowep, zmosop, zmoeap, znorth, ilonrgg, istart, iend,
     &         zdella, zdello, zmiss, cdclass, cdexpver,cldomain, ldfdb,
     &         ksect1_3, ifdb, jpidl)
            ENDIF
          ELSE
            WRITE(kuso,*) ' TRY TO SPLIT A NONGLOBAL FIELD '
            CALL abort
          ENDIF
        ENDIF
        IF (ldfdb) THEN
          IF ( kfcst .EQ. 0 ) THEN
            cltype='an'
          ELSE
            cltype='fc'
          ENDIF
          idir = 0
          ifre = 0
          IF (llsuvi.OR.lltest) itest=1
          CALL FLD2FDB (kuso, itest, idummy, idummy, idummy, 'close',
     &                  ierr,
     &                 cdclass, cdexpver, clfdbsf, cldomain,
     &                 idir, ifre,  ifdb)
          IF (llsuvi) WRITE(kuso,*) ' FDB ', ifdb, " CLOSED:", ierr
          IF (ierr .NE. 0) THEN
            WRITE(kuso,*) ' ------------------------'
            WRITE(kuso,*) ' ERROR CLOSING FDB '
            WRITE(kuso,*) ' FDB ERROR CODE IS ',ierr
            WRITE(kuso,*) ' ------------------------'
            CALL abort
          ENDIF
        ENDIF


      ELSEIF (cdtask .EQ. 'BLOCK') THEN

! ----------------------------------------------------------------------

!* 4.    DECODING AND BLOCKING.
!  ----------------------------

        IF (llsuvi) WRITE(kuso,*)' SUB: ',csubna,': CALL INMARSS',
     &    ' - DOMAIN = ', cddomain
        TH=-99.
        FR=-99.

        CALL inmarss (kuso, kugrs, kgy, idl, kang, kfre, zmiss,
     &                cldate, iforp, iparam, krgg, kdomrgg, ingx, ingy,
     &                zamowep, zamosop, zamoeap, zamonop,
     &                inang, infre, zsec4, th, fr, kret, clmsg)
         
        IF (llsuvi.OR.lltest) THEN
          WRITE(kuso,*) ' SUB: ',csubna,': INMARSS DONE '
          WRITE(kuso,*) csubna,': GRIBEX RETURN CODE = ',kret, clmsg
        ENDIF  
        IF (lltest) THEN
          WRITE(kuso,*)' SUB: ',csubna,': INMARSS DONE'
          WRITE(kuso,*)' CLMSG : ', clmsg
          WRITE(kuso,*)' KRET  : ', kret
          WRITE(kuso,*)' TH    : '
          WRITE(kuso,'(10f8.2)') (th(j1),j1=1,inang)
          WRITE(kuso,*)' FR    : '
          WRITE(kuso,'(10f8.2)') (fr(j1),j1=1,infre)
          WRITE(kuso,*)' INANG : ', inang          
          WRITE(kuso,*)' INFRE : ', infre              
          WRITE(kuso,*)' NORTH : ', zamonop
          WRITE(kuso,*)' WEST  : ', zamowep
          WRITE(kuso,*)' EAST  : ', zamoeap
          WRITE(kuso,*)' SOUTH : ', zamosop
          WRITE(kuso,*)' INGX  : ', ingx
          WRITE(kuso,*)' INGY  : ', ingy
          WRITE(kuso,*)' IPARAM: ',iparam
        ENDIF
        IF (cddomain(1:5) .EQ. 'both' ) THEN
          zlatitude=0.0
          zlongitude=zmoeap
          iy = NINT((zlatitude-zmonop)/(-zdella)) + 1
          zx = MOD(zlongitude-zmowep+720.,360.)
          zx = zx / pdello(kgy-iy+1)
          ix = NINT(zx) +1
          itsp = ispoff(ilonrgg, kgy, kang, kfre, iy-1, ix) +
     &           (kang*kfre)
          istrt = 1+itsp

          IF (llsuvi) WRITE(kuso,*) " SOUTH: istrt=", istrt
          CALL inmarss (kuso, kugrs, kgy, idl-istrt+1, kang, kfre,zmiss,
     &                  cldate, iforp, iparam, krgg,kdomrgg, ingx, ingy,
     &                  zamowep, zamosop, zamoeap, zamonop,
     &                  inang, infre, zsec4(istrt), th, fr, kret, clmsg)
          IF (lltest) THEN
            WRITE(kuso,*)' SUB: ',csubna,': INMARSS DONE'
            WRITE(kuso,*)' CLMSG : ', clmsg
            WRITE(kuso,*)' KRET  : ', kret
            WRITE(kuso,*)' TH    : '
            WRITE(kuso,'(10f8.2)') (th(j1),j1=1,inang)
            WRITE(kuso,*)' FR    : '
            WRITE(kuso,'(10f8.2)') (fr(j1),j1=1,infre)
            WRITE(kuso,*)' INANG : ', inang          
            WRITE(kuso,*)' INFRE : ', infre              
            WRITE(kuso,*)' NORTH : ', zamonop
            WRITE(kuso,*)' WEST  : ', zamowep
            WRITE(kuso,*)' EAST  : ', zamoeap
            WRITE(kuso,*)' SOUTH : ', zamosop
            WRITE(kuso,*)' INGX  : ', ingx
            WRITE(kuso,*)' INGY  : ', ingy
            WRITE(kuso,*)' IPARAM: ',iparam
          ENDIF
        ENDIF
        IF (ldlog) THEN
          IF (llsuvi) WRITE(kuso,*)' SUB: ',csubna,' : DO LOG10'
          DO jx=1,idl
            IF (zsec4(jx) .NE. zmiss) THEN
              zsec4(jx) = 10.**(zsec4(jx)-abs(pprec))- ppeps
            ELSE
              zsec4(jx) = 0.0
            ENDIF
          ENDDO
        ENDIF
        IF (ldoutspp) THEN
          DO j1=1,kspecs
            READ ( iu_locs, location,END=4500)
            zy=zlatitude
            zx=zlongitude
            IF (llsuvi) WRITE(kuso,location)

            iy = NINT((zlatitude-zmonop)/(-zdella)) + 1

            IF ( iy.GE.1 .AND. iy.LE.kgy ) THEN
              zx = MOD(zlongitude-zmowep+720.,360.)
              zx = zx / pdello(iy)
              ix = NINT(zx) +1
            ENDIF

            IF ( ix.LT.1 .OR. ix.GT.ilonrgg(iy) .OR.
     &           iy.LT.1 .OR. iy.GT.kgy) THEN
              WRITE(kuso,*)zlatitude,'X',zlongitude,' IS OUTSIDE GRID'
              WRITE(kuso,*) ' IX=',ix,' KGX=',ilonrgg(iy)
              WRITE(kuso,*) ' IY=',iy,' KGY=',kgy
              WRITE(kuso,*) ' ZMOWEP=',zmowep
              WRITE(kuso,*) ' ZDELLO=',zdello
              WRITE(kuso,*) ' ZMONOP=',zmonop
              WRITE(kuso,*) ' ZDELLA=',zdella
            ELSE
              IF (llsuvi) 
     &        WRITE( kuso,*)' PRINT SPECTRUM AT IX=',ix,
     &        ' AND AT IY=',iy
              i_off = ispoff(ilonrgg, kgy, kang, kfre, iy-1, ix)
              DO j3=1,kfre
                DO j2=1,kang
                  jx=(j3-1)*kang+1+(j2-1)+i_off

                  zspectrm(j2,j3) = zsec4(jx)
                ENDDO
              ENDDO
              CALL PRSPP (kuso, cldate, zlongitude,  zlatitude, clname,
     &                    waveht, wavedir, wavemfr, windspd, winddir,
     &                    wind10,
     &                    fr, th, zspectrm, kang, kfre, kang, kfre)
            ENDIF
          ENDDO
        ENDIF
 4500   CONTINUE
        fl3=0.
        DO jblock=1,kblo
          CALL makeblos (fl3, zsec4, idl, 
     &                   ijs(jblock), ijl(jblock), jblock,
     &                   ilonrgg, kgx, kgy, kang, kfre, ixlg, kxlt,
     &                   kblo, kiblo)
          WRITE(kubls) FL3
        ENDDO
      ENDIF
      IF (ALLOCATED(ZSEC4)) DEALLOCATE(ZSEC4)

      IF (ldoutspp) CLOSE (iu_locs)

! ----------------------------------------------------------------------

!* 5.    TERMINATION.
!  ------------------
      msub = msub - 1
      RETURN
      END SUBROUTINE blspcon
