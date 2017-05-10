
#.=====================================================================
#.  Define user input (needed for preset and wamodel 
#.=====================================================================
#.
# This example runs 6 hrs with analysed winds (but no data assimilation)
typeset -Z10 start_date=2003072912     # starting date (YYYYMMDDHH)
typeset -i anlength=6                  # Length of analysis in hours
typeset -i fclength=0                  # Length of forcast in hours
#.
#.
typeset -i anlen=anlength*3600
typeset -Z7 antime=$anlen
typeset -i fclen=fclength*3600
typeset -Z7 fctime=$fclen
beginfcdt=$(newdate $start_date $anlength)
enddt=$(newdate $beginfcdt $fclength)
typeset -Z12 begofrn=${start_date}00   # BEGin date OF RUn
typeset -Z12 endofrn=${enddt}00        # END   date OF RUn
typeset -Z12 begoffo=${beginfcdt}00    # BEGin date OF FOrcast.
#                                        This date must equal to endofrn
#                                        when analysis is only required
typeset -Z12 outofrf=${endofrn}        # Date to output restart binary files
                                       # or the gribbed spectra if grib output
                                       # was requested. This date is on top of
                                       # all other dates that might be set by
                                       # cycling from the initial date with
                                       # step given by IDELRES.
                                       # Set it to 0000000000 if determined
                                       # by namelist NAOS (which is not used in 
                                       # this example)
#.
CLASS=RD                               # user class (used when coding the
#.                                                   output data in grib)
#.! expver is the experiment id. It is used when gribbing the output data. 
typeset -l expver=wave
#.
#if region=='s'
iassi=0
laltas=F
lsarinv=F
lsaras=F
#else
iassi=0                                # iassi=0 no data assimilation
#. Specify what type of data will be used for the assimilation
laltas=F                               # altimeter
lsarinv=F                              # SAR inversion
lsaras=F                               # SAR
#endif
#.
#.define file storage directory, as needed for output destination
##################################################################
#.
export STORAGE_PATH=$WDIR
#.
#.=====================================================================
#.  end user input (needed for preset and wamodel 
#.=====================================================================
