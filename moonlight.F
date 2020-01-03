      program moonlight
      implicit real*8(a-h,o-z)
c23456789012345678901234567890123456789012345678901234567890123456789012
c ***
c *** Program to calculate the relative amount of illumination at each
c *** hour on specified nights at a given location.
c ***
c ***  Version 2.2   (January 2020)
c ***  For unix systems running Fortran 90.
c ***
c ***  Reference:
c ***   Austin, R.H., Phillips, B.F. and Webb, D.J. (1976)  A method
c ***       for calculation moonlight illuminance at the Earth's
c ***       surface.  J. Appl. Ecol., 13(3), 741-748.
c ***
c ***  Format for input file
c ***
c ***   Unless stated fields should be right justified and without
c ***   a decimal point.  Dates and times should be zero filled
c ***   (i.e. use zeros instead of blanks).
c ***
c ***   Column   Contents
c ***   1-2      Latitude - degrees.
c ***   4-5      Latitude - minutes.
c ***   7        N or S for North or South latitude.
c ***   11-13    Longitude - degrees.
c ***   15-16    Longitude - minutes.
c ***   18       E or W for East or West longitude
c ***   21-29    Time correction.  The number of hours that the
c ***            time is ahead of or behind greenwich mean time.
c ***            Positive if ahead, negative if behind.  If no
c ***            decimal point, decimal point assumed to lie to
c ***            the right of column 29.  Otherwise punch left
c ***            justified with decimal point.
c ***   31-35    Number of days for which the user requires the
c ***            calculation to be performed.
c ***   37-44    Date (DDMMYYYY) of first day of period for which
c ***            calculation is to be performed
c ***   46-49    Time (HHMM) of day of the first instance for
c ***            which the illumination is to be calculated.
c ***   51-54    Time (HHMM) of day of the last instance for
c ***            which the illumination is to be calculated.
c ***            If the period between the first and second
c ***            times includes midnight, add 2400 to the time
c ***            (i.e. use 3000 for 6 am the following day).
c ***   56-59    Time period (HHMM) between calculations.
c ***            (Thus 0130 calculates every 1 hour and 30
c ***            minutes).
c ***
c ***   An end-of-file or blank line ends the input data.
c ***
c ***  Example of the use of this program:
c ***
c ***    To obtain the relative amount of illuminance for 10 nights
c ***    at each hour from 1800 hours to 0600 hours the following day
c ***    from 26/09/77 to 05/10/77 at Perth, Western Australia
c ***    (Latitude 31 degrees, 58 minutes south, longitude 115 degrees
c ***    51 minutes east).  Local time is 8 hours ahead of Greenwich
c ***    mean time.  The control cards are therefore:
c *** Column 12345678901234567890123456789012345678901234567890123456789
c ***        31 58 S   115 51 E  +8.0         10 26091977 1800 3000 0100
c ***        Blank line
c ***
c *** The authors woul;d appreciate being informed of any errors
c *** discovered in this routine, or any problems in its use.  Please
c *** contact: Dr B.F. Phillips
c ***          CSIRO, P.O. Box 20, North Beach, Western Australia 6020
c ***          Australia.
c ***
c ***          Dr D.J. Webb
c ***          National Oceanography Centre, Southampton SO14 3ZH,
c ***          United Kingdom.
c ***
c *** Notes on the conversion of this program to other computer systems
c ***
c *** 1.  The program was originally written for a CDC Cyber 76 computer
c ***     in Fortran IV.
c *** 2.  The version of Fortran used is designed to comply with American
c ***     National Standards institute (ANSI) Fortran language as
c ***     described in X3.9-1966, but has been extended by CDC.
c *** 3.  The programme listing contains details of non-ansi coding.
c *** 4.  The word size of the Cyber is 60 bits with-
c ***      (i)   Integers ranging from -2**59 - 1 to 2**59 - 1
c ***      (ii)  Real numbers ranging from 10**(-293) to 10**(322) with
c ***            precision of approximately 14 decimal digits.  For
c ***            double precision the precision is approximately 29
c ***            decimal digits.
c ***      (iii) Each word is capable of holding 10 alphanumeric
c ***            characters
c *** 5.  Some subroutines are non-standard CDC library routines
c ***     e.g. 'data' called by newpage returns an alphanumeric date
c ***     in the format 10HBYY/MM/DDB, where B=Blank, YY=Year, MM=Month,
c ***     DD=day.
c ***     Also acos called by routines esunmoon and psunmo.
c ***
c ***
c ***  Changes and corrections
c ***
c ***  A number of minor alterations have been made to this programme
c ***  since its inception, in order to tidy up the code and to extend
c ***  its versitility.  The method presented by Austin et al. (1976)
c ***  has not been modified, however the original programme contained
c ***  errors.  The following changes have been made to correct these
c ***  and users should ensure that earlier versions of the programme
c ***  are corrected in the same way:
c ***
c ***  1.  The follwing statements originally found in subroutine
c ***      bright
c ***        VE=EV2(IZ)+(Z-(IZ-1)*2.)*(EV2(IZ1)-EV2(IZ))/2.
c ***        VE=EV3(IZ)+(Z-(IZ-1))*(EV3(IZ1)-EV3(IZ))
c ***     have been replaced by
c ***        VE=EV2(IZ)+(Z-50.-(IZ-1)*2.)*(EV2(IZ1)-EV2(IZ))/2.
c ***        VE=EV3(IZ)+(Z-70.-(IZ-1))*(EV3(IZ1)-EV3(IZ))
c ***     The error was found by Dr John T. Collins, Fisheries Research
c ***     Station, Tehkummah, Ontario, Canada.  It exaggerated the effect
c ***     of the angular distance of the moon above the horizon.
c ***                                                         (17/8/78)
c ***
c ***  2.  The code has been modified to work with Fortran 90 running
c ***      on a Linux workstation.
c ***
c ***  Version 2.1
c ***
c ***  3.  Dates outside the range 1900-2000 can be used.
c ***
c ***  4.  The name of the data file can be provided on the command line.
c ***      For example:
c ***      ./moonlight PhillipIsland.dat
c ***
c ***  Version 2.2
c ***
c ***  5.  File README gives further information on changes of accuracy
c ***      with time.
c ***
c

      character*1 latdrn,lngdrn,isnmn*2
      common latdeg,latmin,lngdeg,lngmin,timcor,lincnt,latdrn,lngdrn
      real*8 long
      character*4 cj
      integer ::  narg
      character(len=256) :: input_string, input_file
c
      pi = 4d0*atan(1d0)
!
!  Look for input file on command line
!
       narg = command_argument_count()
       if(narg.eq.0)then
         input_file = "input.dat"
       else
         call get_command_argument(1,input_string)
         input_file = adjustl(input_string)
         print *," Using input file = ",trim(input_file)
       endif
c
c  choose a new date
c
      open(5,file=trim(input_file))
   5  read(5,60,end=300)latdeg,latmin,latdrn,lngdeg,lngmin,lngdrn,
     &              timcor,ndays,iday,imth,iyr,ihr,imn,jhr,jmn,khr,kmn
      if(ndays.eq.0)goto 300
  60  format(i2,1x,i2,1x,a1,3x,i3,1x,i2,1x,a1,2x,f9.0,1x,i5,1x,2i2,i4,
     &       1x,2i2,1x,2i2,1x,2i2)
c
      lincnt = 100
!     iyr = 1900+iyr
c
c carry out checks
c
      if(latdrn.ne.'N'.and.latdrn.ne.'n'.and.latdrn.ne.'S'
     &                .and.latdrn.ne.'s')goto 5
      if(lngdrn.ne.'E'.and.lngdrn.ne.'e'.and.lngdrn.ne.'W'
     &                .and.lngdrn.ne.'w')goto 5
      if(latdrn.eq.'N'.or.latdrn.eq.'n')then
        plat = 90d0-(latdeg+latmin/60d0)
      else
        plat = 90d0+(latdeg+latmin/60d0)
      endif
      if(lngdrn.eq.'E'.or.lngdrn.eq.'e')then
        long = lngdeg+lngmin/60d0
      else
        long = 360d0-(lngdeg+lngmin/60d0)
      endif
      gmt = timcor
      if(khr.lt.0.or.kmn.lt.0)then
        print *,khr,kmn
        print *,'negative time increment - check'
        goto 5
      endif
      if(ihr.lt.0.or.imn.lt.0.or.jhr.lt.0.or.jmn.lt.0)then
        print *,'negative start or end time - check'
        goto 5
      endif
      if(khr.eq.0.and.kmn.eq.0)khr=1
      if(jhr.lt.ihr.or.(jhr.eq.ihr.and.jmn.lt.imn))then
        jhr = ihr
        jmn = imn
      endif
c
c check start time and time difference
c
      if(ihr.gt.23)then
        print *,' starting hour larger than 23 - check'
        goto 5
      endif
      idifmn = (jhr-ihr)*60d0+jmn-imn
      if(idifmn.gt.1440)then
        print *,' Difference of start and end times is larger ',
     &          'than 24 hours - check'
        goto 5
      endif
c
c  main loop over days
c
      do 200 kd = 1,ndays
c  obtain date of next day
      call addone(iday,imth,iyr,iiday,iimth,iiyr)
c  set code to detect first printout for the period following midnight
      lprcod = 1
c
c   starting from the first time instant, go through the day
c   incrementing the time until the last time instant is exceeded
c
      lhr = ihr
      lmn = imn
      goto 201
 202  lhr = lhr+khr
      lmn=lmn+kmn
 204  if(lmn.gt.60)then
        lmn = lmn-60
        lhr = lhr+1
        goto 204
      endif
 203  if(lhr.lt.jhr)goto 201
      if(lhr.gt.jhr)goto 100
      if(lmn.gt.jmn)goto 100
 201  hih = lhr+lmn/60d0
      if(hih.gt.24.)goto 206
      call esunmoon(as,am,el,rm,hih,iday,imth,iyr,gmt,plat,long)
      goto 207
 206  hih=hih-24d0
      call esunmoon(as,am,el,rm,hih,iiday,iimth,iiyr,gmt,plat,long)
c
 207  llhr = lhr*100d0+lmn
      if(lhr.gt.24)llhr = (lhr-24)*100+lmn
      write(cj,"(i4)")llhr
c
c  test if the sun is above the horizon or less than 12 degrees
c  below the horizon
c
      if(as.le.-12)goto 213
      if(am.lt.0.)goto 214
      isnmn = '*'
      goto 101
 214  isnmn = '*$'
      goto 101
 213  if(am.ge.0.)goto 215
      isnmn = '$'
      goto 101
 101  if(lhr.eq.ihr.and.lmn.eq.imn) got o 102
      if(lhr.ge.24.and.lprcod.eq.1) goto 103
      call newpag
      print 104,cj,as,am,el,rm,isnmn
 104  format(15x,a4,1x,3(5x,g12.5),3x,g12.5,6x,a2)
      goto 110
 102  call newpag
      print 105, iday,imth,iyr,cj,as,am,el,rm,isnmn
 105  format(2x,2(i2,2h- ),i4,1x,a4,1x,3(5x,g12.5),3x,g12.5,6x,a2)
      goto 110
 103  call newpag
      lprcod = 2
      print 105,iiday,iimth,iiyr,cj,as,am,el,rm,isnmn
      goto 110


c
 215  if(lhr.eq.ihr.and.lmn.eq.imn)ell=el
      call bright(rm,el,ell,am,brt)
      brt = brt*100d0
      siglev = 0.001
      if(brt.lt.siglev)then
 106    if(lhr.eq.ihr.and.lmn.eq.imn) goto 107
        if(lhr.ge.24.and.lprcod.eq.1)goto 108
        call newpag
        print 109, cj,as,am,el,rm,siglev
 109  format(15x,a4,1x,3(5x,g12.5),3x,g12.5,6x,'less than ',f5.3)
        goto 110
 107    call newpag
        print 111, iday,imth,iyr,cj,as,am,el,rm,siglev
 111  format(2x,2(i2,2h- ),i4,1x,a4,1x,3(5x,g12.5),3x,g12.5,6x,
     & 'less than ',f5.3)
        goto 110
 108    call newpag
        lprcod = 2
        print 111, iiday,iimth,iiyr,cj,as,am,el,rm,siglev
        goto 110
      endif

      if(lhr.eq.ihr.and.lmn.eq.imn)goto 112
      if(lhr.gt.24.and.lprcod.eq.1) goto 113
      call newpag
      print 114, cj,as,am,el,rm,brt
 114  format(15x,a4,1x,3(5x,g12.5),3x,g12.5,6x,f7.3)
      goto 110
 112  call newpag
      print 115,iday,imth,iyr,cj,as,am,el,rm,brt
 115  format(2x,2(i2,2h- ),i4,1x,a4,1x,3(5x,g12.5),3x,g12.5,6x,f7.3)
      goto 110
 113  call newpag
      lprcod = 2
      print 115,iiday,iimth,iiyr,cj,as,am,el,rm,brt
 110  continue
c
c store the elongation for use at the next iteration
c
      ell = el
c
c go back for next iteration
c
      goto 202
c
c  no more iterations this day.  Print footnotes and go onto next day
c
 100  continue
      call newpag
      print 120
 120  format('     $ Moon below horizon    * Sun less than 12 ',
     &           'degrees below horizon')
      iday = iiday
      imth = iimth
      iyr = iiyr
 200  continue
      goto 5
 300  continue
      stop
      end

      subroutine bright(rm,el,ell,am,brt)
      implicit real*8(a-h,o-z)
c
c  routine to find brightness of moon
c
c  Input:
c
c  am - real*8 - angle in degrees of moon above horizon
c  rm - real*8 - distance to moon in terms of 1/(mean equatorial
c                parallax)
c  el - real*8 - elongation in degrees
c ell - real*8 - elongation at alst instant (used to determine if
c                the moon is waxing or waining.
c
      dimension en(11),ex(11),ev1(11),ev2(11),ev3(18)
      data pi/3.1415926532/
      data en/0d0,0.1d0,1d0,3.2d0,7.5d0,15d0,26d0,41d0,65d0,100d0,100d0/
      data ex/0d0,0.1d0,2d0,4.8d0,9.5d0,19d0,30d0,46d0,72d0,100d0,100d0/
      data ev1/0d0,0d0,0d0,0d0,0.01d0,0.02d0,0.03d0,0.04d0,0.06d0,
     &         0.09d0,0.12d0/
      data ev2/0.12d0,0.14d0,0.16d0,0.18d0,0.20d0,0.23d0,0.26d0,
     &         0.30d0,0.34d0,0.39d0,0.45d0/
      data ev3/0.45d0,0.48d0,0.52d0,0.56d0,0.60d0,0.65d0,0.71d0,
     &         0.77d0,0.83d0,0.91d0,0.99d0,1.08d0,1.19d0,1.33d0,
     &         1.52d0,1.77d0,2.12d0,2.61d0/
c
c      print *,' bright entered with - rm, el, ell, am, brt:'
c      print 11, rm,el,ell,am
c  11  format(4f12.6)
c
      radian = pi/180d0
      d=1/rm**2
      le = (el/20)+1.
      le1 = le+1
      if(el.lt.126..and.el.le.ell)goto 30
c  waxing moon
      ep = ex(le)+(el-(le-1)*20d0)*(ex(le1)-ex(le))/20d0
      goto 35
  30  continue
c  waning moon
      ep = en(le)+(el-(le-1)*20d0)*(en(le1)-en(le))/20d0
  35  continue
      ep=ep/100d0
      z=90d0-am
      if(z.ge.50d0)goto 40
      iz=(z/5d0)+1d0
      iz1=iz+1
      ve = ev1(iz)+(z-(iz-1)*5d0)*(ev1(iz1)-ev1(iz))/5d0
      goto 100
  40  if(z.ge.70d0)goto 60
      iz = (z/2d0)-24d0
      iz1=iz+1
      ve=ev2(iz)+(z-50d0-(iz-1)*2d0)*(ev2(iz1)-ev2(iz))/2d0
      goto 100
  60  if(z.ge.87)goto 80
      iz = z-69.
      iz1=iz+1
      ve=ev3(iz)+(z-70d0-(iz-1))*(ev3(iz1)-ev3(iz))
      goto 100
  80  ve = 2.61
 100  ve = -0.4d0*ve
      at = 10d0**ve
      brt=d*ep*at
      amr=am*radian
      ami = sin(amr)
      brt=brt*ami
      return
      end
#if 0
      subroutine moonp(as,am,el,rm,hour,iday,month,iyear,plsgmt,
     &                 pole, rlon)
      implicit real*8(a-h,o-z)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c  subroutine to calculate angles of the moon and sun above the horizon
c
c  Input:
c
c    hour   -  real*8  - hours since midnight
c    iday   -  integer - day of month (1=1st)
c    month  -  integer - number of month (1=January)
c    iyear  -  integer - year in range 1900 to 2000
c    plsgmt -  real*8  - time zone relative to Greenwich
c    pole   -  real*8  - polar angle of the place of interest in degrees
c                        (equals 90 - latitude north of equator)
c    rlon   -  real*8  - longitude of place of interest in degrees
c
c  Output:
c
c    as     -  real*8  - angle in degrees of sun above the horizon
c    am     -  real*8  - angle in degrees of moon above the horizon
c    el     -  real*8  - elongation in degrees
c    rm     -  real*8  - distance of the moon in units of 1/(mean
c                                   equatorial parallax)
c
c***********************************************************************
c
      data pi/3.1415926532d0/
c
      radian = pi/180d0
      call skytime(hour,iday,month,iyear,plsgmt,hrgmt,centj)
      call psunmo(hrgmt,gra,ps,rlonsu,rs,pm,rlonmo,rm,t15)
      c1 = cos(pole*radian)
      c2 = cos(rlon*radian)
      s1 = sin(pole*radian)
      s2 = sin(rlon*radian)
      zs = acos(sin(ps)*s1*(cos(rlonsu)*c2+sin(rlonsu)*s2)+cos(ps)*c1)
      zm = acos(sin(pm)*s1*(cos(rlonmo)*c2+sin(rlonmo)*s2)+cos(pm)*c1)
      as = 90d0-zs/radian
      am = 90d0-zm/radian
      el = acos(sin(ps)*sin(pm)*(cos(rlonsu)*cos(rlonmo)
     &             +sin(rlonsu)*sin(rlonmo))+cos(ps)*cos(pm))/radian
      return
      end

      subroutine psunmo(hrgmt,gra,ps,rlonsu,rs,pm,rlonmo,rm,t15)
      implicit real*16(a-h,o-z)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c***********************************************************************
c
c  subroutine to calculate position of the moon and sun using
c  astronomical constants for epoch 1900
c
c  Input:
c
c    hrgmt  -  real*8  - hours since 0000 1st January 1900
c
c  Output:
c
c    gra    -  real*8  - right ascension of greenwich
c    ps     -  real*8  - polar angle of the sun
c    rlonsu -  real*8  - longitude of sun east of Greenwich
c    rs     -  real*8  - relative distance of the sun
c    pm     -  real*8  - polar angle of the moon
c    rlonmo -  real*8  - longitude of moon east of Greenwich
c    rm     -  real*8  - distance of the moon in units of 1/(mean
c                                   equatorial parallax)
c
c  Note:  1.  all angles in radians
c         2.  All internal calculations use real*16 for accuracy
c
c***********************************************************************
c
      real*8  hrgmt,gra,ps,rlonsu,rs,pm,rlonmo,rm,t15
      real*16 c(19),ls,lm,nm,t
c     data in/0/
      save in
      in = 0
c
      if (in.eq.0)goto 200
c
 100  t = (hrgmt+12d0)/(36525d0*24d0)
c
c  calculate sun's position
c
      hs = mod((52d-7*t+628.3319509d0)*t+4.8816280d0,pi2)
      ps = mod((79d-7*t+0.0300053d0  )*t+4.9082295d0,pi2)
      es = (-126e-9*t-418e-7)*t+0.01675104
      cos1 = cos(hs-ps)
      sin1 = sin(hs-ps)
      ls = (2.5d0*es*cos1+2.)*es*sin1+hs
      rs = ((2d0*cos1**2-1d0)*es+cos1)*es+1d0
      ps = acos(sin(ls)*c(18))
c
      tval = (-0.5d0*pi+ls+ps)*0.5d0
      fis = 2d0*atan(c(1)*sin(tval)/cos(tval))
      t15 = mod(15d0*pi*hrgmt/180d0,pi2)
      gra = mod(t15-12d0*15d0*radian+hs,pi2)
      rlonsu = mod(fis-gra,pi2)
      if(t15.lt.0d0)t15 = t15+pi2
      if(gra.lt.0d0)gra = gra+pi2
      if(rlonsu.lt.0d0)rlonsu=rlonsu+pi2
c
c      print 11,hrgmt,t
c      print 12,hs*degree,ps*degree,es*degree
c      print 12,cos1,sin1,c(18)
c      print 12,es*degree,cos1,sin1,hs*degree
c      print 12,ls*degree,ps*degree
c      print 12,tval,fis,t15*degree
c      print 12,gra*degree,rlonsu*degree
c      print 12,90d0-ps*degree,rlonsu*degree
c  11  format(3f12.4)
c  12  format(4f12.6)
      in=in+1
      if(in.ge.3)stop
c
c calculate moon's position
c
      hm = mod((346d-7*t+8399.7092745d0)*t+4.7200089d0,pi2)
      pm = mod((-1801d-7*t+71.0180412d0)*t+5.8351526d0,pi2)
      nm = mod((363d-7*t-33.7571463d0)*t+4.5236016d0,pi2)
c
      tval = 0.5*nm
      q1= sin(tval)/cos(tval)
      p = atan(c(2)*q1)
      q = atan(c(3)*q1)
      v = p+q
      y = p-q
      wm = pi-2d0*atan(sin(q)*c(4)/sin(p))
      rom = hm-nm+y
      sin2 = sin(hm-pm)
      cos2 = cos(hm-pm)
      sin3 = sin(hm-hs)
      cos3 = cos(hm-hs)
      sin4 = 2d0*sin2*cos2
      cos4 = 2d0*cos2**2-1d0
      sin5 = 2d0*sin3*cos3
      cos5 = 2d0*cos3**2-1d0
      sin6 = sin5*cos2-cos5*sin2
      cos6 = cos5*cos2+sin5*sin2
      sin7 = sin5*cos2+cos5*sin2
      cos7 = cos5*cos2-sin5*sin2
      sin8 = sin5*sin1-cos5*cos1
      cos8 = cos5*cos1+sin5*sin1
      lm = rom + c(5)*sin2 +c(6)*sin4 +c(7)*sin6 +c(8)*sin5 +c(9)*sin7
     &         + c(10)*sin8*es + c(19)*sin1*es
      rm = 1d0 + c(11)*(c(12)*cos2 +c(13)*cos4 +c(14)*cos6 +c(15)*cos5
     &                 +c(16)*cos7 +c(17)*cos8*es)
      pm = acos(sin(lm)*sin(wm))
c
      tval = (0.5d0*pi+wm)*0.5d0
      y = sin(tval)/cos(tval)
      tval = (-0.5d0*pi+lm+pm)*0.5d0
      fim = 2d0*atan(y*sin(tval)/cos(tval))
      if(fim.lt.0d0)fim=fim+pi2
      rlonmo = mod(v+fim-gra,pi2)
      if(rlonmo.lt.0d0)rlonmo=rlonmo+pi2
      return
c
c set up constants
c
 200  in=1
      c1 = 1d0
      c2 = 2d0
      p5 = 0.5d0
      pi = 4d0*atan(c1)
      pi2 = 2d0*pi
      radian = pi/180d0
      degree = 180d0/pi
      ws = 23.452*radian
      wi =  5.145*radian
      em = 0.0549
      ratm = 0.074804
      ratm2 = ratm**2
      em2 = em**2
      tval = (0.5d0*pi+ws)*0.5d0
      c(1) = sin(tval)/cos(tval)
      c(2) = cos((wi-ws)*0.5d0)/cos((wi+ws)*0.5d0)
      c(3) = sin((wi-ws)*0.5d0)/sin((wi+ws)*0.5d0)
      tval = (wi-ws)*0.5d0
      c(4) = cos(tval)/sin(tval)
      c(5) = 2d0*em
      c(6) = 1.25d0*em2
      c(7) = (15d0/4d0 + 263d0*ratm/16d0)*ratm*em
      c(8) = (75d0*em2/(16d0*ratm) + 59d0*ratm/12d0 + 11d0/8d0)*ratm2
      c(9) = ratm2*em*17d0/8d0
      c(10) = 77d0*ratm2/16d0
      c(11) = 6d0/(6d0+ratm2)
      c(12) = em
      c(13) = em2
      c(14) = (329d0*ratm/64d0 + 15d0/8d0)*ratm*em
      c(15) = (15d0*em2/(4d0*ratm)+19d0*ratm/6d0 +1d0)*ratm2
      c(16) = 33d0*ratm2*em/16d0
      c(17) = 7d0*ratm2/2d0
      c(18) = sin(ws)
      c(19) = -3d0*ratm
      goto 100
c
      end

      subroutine skytime(hour,iday,month,iyear,plsgmt,hrgmt,centj)
      implicit real*8(a-h,o-z)
c
c23456789012345678901234567890123456789012345678901234567890123456789012
c
c  subroutine to calculate hrgmt the number of hours since 0000 on
c  1 Jan 1900 and the number of Julian centuries since 1200 on
c  31 Decemebr 1899.
c
c  Input:
c
c    hour   -  real*8  - hours since midnight
c    iday   -  integer - day of month (1=1st)
c    month  -  integer - number of month (1=January)
c    iyear  -  integer - year in range 1900 to 2000
c    plsgmt -  real*8  - time zone relative to Greenwich
c
c  Output
c
c    hrgmt  -  real*8 - hours since 0000 on 1st Jan 1900
c    centj  -  real*8 - Julian centuries since 1200 31st Dec 1899
c
      dimension mday(12)
      data  mday /0,31,59,90,120,151,181,212,243,273,304,334/
c
      if(iyear.lt.1900.or.iyear.ge.2100)then
        print *,' Program moonlight should only be used for'
        print *,' the period 1900 to 2100.'
        print *,' Programme stopping ...'
        stop
      endif
c
      nyear = iyear-1900
!      nleap = nyear/4
!      if(4*nleap.eq.nyear.and.month.le.2.and.nyear.ne.0)nleap=nleap-1
!      hrgmt = (nyear*365+nleap+mday(month)+iday-1)*24d0+hour-plsgmt
!      centj = (hrgmt+12d0)/(36525d0*24d0)

      call leapyear(iyear,ileap)
      call leapyears(iyear,nleapp)
      call leapyears(1900,n1900)
      nleapp = nleapp-460
!      print *,iyear,nyear,nleap,nleapp,ileap,n1900
      nleap = nleapp
      if(ileap.eq.1.and.month.gt.2)nleap = nleap + 1

      hrgmt = (nyear*365+nleap+mday(month)+iday-1)*24d0+hour-plsgmt
      centj = (hrgmt+12d0)/(36525d0*24d0)

      return
      end subroutine skytime
#endif

      subroutine newpag
      implicit real*8(a-h,o-z)
c
      character*1 latdrn,lngdrn,latsvc,lngsvc,ins*5,iew*4
      character*8 date
      data icode/1/
      save icode,  maxlin, lpage,  date
      save latsva, latsvb, lngsva, lngsvb, timsv, latsvc, lngsvc
      common latdeg,latmin,lngdeg,lngmin,timcor,lincnt,latdrn,lngdrn
c
      if(icode.gt.1)goto 2

   1  icode=2
      maxlin=40
      lpage=0
c     call date(date)
      date ='DD/MM/YY'

   3  lincnt = maxlin+1
      latsva=latdeg
      latsvb=latmin
      latsvc=latdrn
      lngsva=lngdeg
      lngsvb=lngmin
      lngsvc=lngdrn
      timsv=timcor

   2  if(latsva.ne.latdeg.or.latsvb.ne.latmin.or.latsvc.ne.latdrn
     &   .or.lngsva.ne.lngdeg.or.lngsvb.ne.lngmin
     &   .or.lngsvc.ne.lngdrn.or.timsv.ne.timcor)goto 3
      lincnt = lincnt+1
      if(lincnt.le.maxlin)goto 4

   5  lpage=lpage+1
      lincnt=0
      print 61
      print 6,lpage,date
  61  format(////)
   6  format(' Program- moonlight',10x,'Calculation of the relative',
     & ' amount of illumination each hour', 19x,'page-',i5,/,
     & ' Run on-',a8,32x,'on specified nights')
      print 7
   7  format(/,' Reference-',10x,'Austin, R.H., Phillips, B.F.,',
     & ' and Webb, D.J. (1976) A method for calculating',/,21x,
     & "moonlight illuminance at the Earth's surface. J. Appl.",
     & ' Ecol., 13(3), 741-748.',/)
c
      ins=latdrn
      if(ins.eq.'N')ins='North'
      if(ins.eq.'S')ins='South'
      print 8,latdeg,latmin,ins
   8  format(' Latitude -',i6,' Degrees  ',i2,' Minutes  ',a5)
      iew=lngdrn
      if(iew.eq.'E')iew='East'
      if(iew.eq.'W')iew='West'
      print 9,lngdeg,lngmin,iew
   9  format(' Longitude -',i5,' Degrees  ',i2,' Minutes  ',a4)
      tim=abs(timcor)
      if(timcor.ge.0d0)then
      print 12,tim
  12  format(' Number of hours that local time is in advance of',
     &       ' Greenwich- ',f8.4,' Hours.')
      else
        print 14,tim
  14  format(' Number of hours that local time is behind',
     &       ' Greenwich- ',f8.4,' Hours.')
      endif
      print 15
      print 16
      print 17
      print 18
      print 19
  15  format(//,
     &'    D a t e    Local     Angle in         Angle in       ',
     &'   Elongation     Distance to      percentage of illumina',
     &'nce')
  16  format(
     &'               time      degrees from     degrees from   ',
     &'   in degrees     moon in units    at full moon at zenith',
     &' at')
  17  format(
     &'               in 24     the Sun to       the Moon to    ',
     &'     (el)         of 1/(mean       distance equal to mean')
  18  format(
     &'               hour      the horizon.     the horizon.   ',
     &'                  equatorial       equatorial parallax.')
  19  format(
     &'               clock.        (as)             (am)         ',
     &'                parallax)',/,78x,'(rm)',/)

   4  return
      end

      subroutine addone(id,im,iy,jd,jm,jy)
      implicit real*8(a-h,o-z)
c
c  returns the date jd/jm/jy of the day following id/im/iy
c
      dimension km(12)
      data km/31,28,31,30,31,30,31,31,30,31,30,31/
c
      jd = id+1
      jm = im
      jy = iy
      if(im.lt.1.or.im.gt.12)then
        print *,' routine addone, month out of range'
        print *,' im = ',im
        print *,' stopping...'
        stop
      endif
c
      kd=km(im)
      if(im.ne.2)goto 3
      if((iy/4)*4.eq.iy)kd=kd+1
      if((iy/100)*100.eq.iy) kd=kd-1
      if((iy/400)*400.eq.iy) kd=kd+1
  3   if(id.lt.1.or.id.gt.kd)then
        print *,' routine addone, day out of range'
        print *,' id = ',id
        print *,' stopping...'
        stop
      endif
c
      if(id.ne.kd)return
      jd = 1
      jm = im+1
      if(im.lt.12)return
      jm=1
      jy=iy+1
      return
      end

      subroutine leapyear(iyear, ileap)
      implicit none
!
!  Subroutine to flag leap years.
!  If iyear is a leap year: ileap = 1
!                     else: ileap = 0
!
      integer, intent(in)  :: iyear
      integer, intent(out) :: ileap

      if(modulo(iyear,4).ne.0)then
        ileap = 0
      elseif(modulo(iyear,400).eq.0)then
        ileap = 1
      elseif(modulo(iyear,100).eq.0)then
        ileap = 0
      else
        ileap = 1
      endif
      return
      end subroutine leapyear

      subroutine leapyears(iyear, nleap)
      implicit none
!
!  Subroutine returns nleap the number of leapyears prior to
!  year iyear since year 1.
!
      integer, intent(in)  :: iyear
      integer, intent(out) :: nleap
      nleap = (iyear-1)/4 -(iyear-1)/100+(iyear-1)/400
      return
      end

      include "sunmoon.F"
