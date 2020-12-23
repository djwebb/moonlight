!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  File sunmoon.F
!
!  This file contains the routines used to determine the position
!  of the Sun and Moon.  These can be used independently of the
!  program moonlight.F
! ***
! ***  Version 2.4   (December 2020)
! ***  For unix systems running Fortran 90.
! ***
! ***  Copyright 1976, 2018, 2020 R.H. Austin, B.F. Phillips, D.J. Webb
! ***  Released under licence GPL-3.0-or-later
! ***
! ***  Reference:
! ***   Austin, R.H., Phillips, B.F. and Webb, D.J. (1976)  A method
! ***       for calculation moonlight illuminance at the Earth's
! ***       surface.  J. Appl. Ecol., 13(3), 741-748.
! ***
!=======================================================================
!
      subroutine esunmoon(as,am,el,rm,hour,iday,month,iyear,plsgmt,
     &                 pole, rlon)
      implicit real*8(a-h,o-z)
!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  subroutine to calculate elevation of the sun and moon above the horizon
!
!  Input:
!
!    hour   -  real*8  - hours since midnight
!    iday   -  integer - day of month (1=1st)
!    month  -  integer - number of month (1=January)
!    iyear  -  integer - year i.e. 2002, in range 1900 to 2100
!    plsgmt -  real*8  - time zone relative to Greenwich
!    pole   -  real*8  - polar angle of the place of interest in degrees
!                        (equals 90 - latitude north of equator)
!    rlon   -  real*8  - longitude of place of interest in degrees
!
!  Output:
!
!    as     -  real*8  - angle in degrees of sun above the horizon
!    am     -  real*8  - angle in degrees of moon above the horizon
!    el     -  real*8  - elongation in degrees
!    rm     -  real*8  - distance of the moon in units of 1/(mean
!                                   equatorial parallax)
!
!=======================================================================
!
      data pi/3.1415926532d0/
      data in /0/
      save in,pi,radian
!
      if(in.eq.0)then
        in = 1
        pi = 4d0*atan(1d0)
        radian = pi/180d0
      endif
!
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
!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  subroutine to calculate position of the moon and sun using
!  astronomical constants for epoch 1900
!
!  Input:
!
!    hrgmt  -  real*8  - hours since 0000 1st January 1900
!
!  Output:
!
!    gra    -  real*8  - right ascension of greenwich
!    ps     -  real*8  - polar angle of the sun
!    rlonsu -  real*8  - longitude of sun east of Greenwich
!    rs     -  real*8  - relative distance of the sun
!    pm     -  real*8  - polar angle of the moon
!    rlonmo -  real*8  - longitude of moon east of Greenwich
!    rm     -  real*8  - distance of the moon in units of 1/(mean
!                                   equatorial parallax)
!
!  Note:  1.  all angles in radians
!         2.  All internal calculations use real*16 for accuracy
!
!=======================================================================
!
      real*8  hrgmt,gra,ps,rlonsu,rs,pm,rlonmo,rm,t15
      real*16 c(19),ls,lm,nm,t
      data in/0/
      save in,c1,c2,p5,pi,pi2,radian,degree,ws,wi,em,ratm,ratm2,em2,
     &     tval,c
!
      if (in.eq.0)goto 200
!
 100  t = (hrgmt+12d0)/(36525d0*24d0)
!
!  calculate sun's position
!
      hs = mod((52d-7*t+628.3319509d0)*t+4.8816280d0,pi2)
      ps = mod((79d-7*t+0.0300053d0  )*t+4.9082295d0,pi2)
      es = (-126e-9*t-418e-7)*t+0.01675104
      cos1 = cos(hs-ps)
      sin1 = sin(hs-ps)
      ls = (2.5d0*es*cos1+2.)*es*sin1+hs
      rs = ((2d0*cos1**2-1d0)*es+cos1)*es+1d0
      ps = acos(sin(ls)*c(18))
!
      tval = (-0.5d0*pi+ls+ps)*0.5d0
      fis = 2d0*atan(c(1)*sin(tval)/cos(tval))
      t15 = mod(15d0*pi*hrgmt/180d0,pi2)
      gra = mod(t15-12d0*15d0*radian+hs,pi2)
      rlonsu = mod(fis-gra,pi2)
      if(t15.lt.0d0)t15 = t15+pi2
      if(gra.lt.0d0)gra = gra+pi2
      if(rlonsu.lt.0d0)rlonsu=rlonsu+pi2
!
!      print 11,hrgmt,t
!      print 12,hs*degree,ps*degree,es*degree
!      print 12,cos1,sin1,c(18)
!      print 12,es*degree,cos1,sin1,hs*degree
!      print 12,ls*degree,ps*degree
!      print 12,tval,fis,t15*degree
!      print 12,gra*degree,rlonsu*degree
!      print 12,90d0-ps*degree,rlonsu*degree
!  11  format(3f12.4)
!  12  format(4f12.6)
!
! calculate moon's position
!
      hm = mod((346d-7*t+8399.7092745d0)*t+4.7200089d0,pi2)
      pm = mod((-1801d-7*t+71.0180412d0)*t+5.8351526d0,pi2)
      nm = mod((363d-7*t-33.7571463d0)*t+4.5236016d0,pi2)
!
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
!
      tval = (0.5d0*pi+wm)*0.5d0
      y = sin(tval)/cos(tval)
      tval = (-0.5d0*pi+lm+pm)*0.5d0
      fim = 2d0*atan(y*sin(tval)/cos(tval))
      if(fim.lt.0d0)fim=fim+pi2
      rlonmo = mod(v+fim-gra,pi2)
      if(rlonmo.lt.0d0)rlonmo=rlonmo+pi2
      return
!
! set up constants
!
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
!
      end



      subroutine skytime(hour,iday,month,iyear,plsgmt,hrgmt,centj)
      implicit real*8(a-h,o-z)
!
!=======================================================================
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!  subroutine to calculate hrgmt the number of hours since 0000 on
!  1 Jan 1900 and the number of Julian centuries since 1200 on
!  31 Decemebr 1899.
!
!  Input:
!
!    hour   -  real*8  - hours since midnight
!    iday   -  integer - day of month (1=1st)
!    month  -  integer - number of month (1=January)
!    iyear  -  integer - year in range 1900 to 2000
!    plsgmt -  real*8  - time zone relative to Greenwich
!
!  Output
!
!    hrgmt  -  real*8 - hours since 0000 on 1st Jan 1900
!    centj  -  real*8 - Julian centuries since 1200 31st Dec 1899
!=======================================================================
!
      dimension mday(12)
      data  mday /0,31,59,90,120,151,181,212,243,273,304,334/
!
      if(iyear.lt.1000.or.iyear.ge.3000)then
        print *,' The errors in the program moonlight'
        print *,' may start being significat outside'
        print *,' the period 1000 to 3000 AD.'
        print *,' Programme stopping ...'
        stop
      endif
!
      nyear = iyear-1900
      nleap = nyear/4
      if(4*nleap.eq.nyear.and.month.le.2.and.nyear.ne.0)nleap=nleap-1
      hrgmt = (nyear*365+nleap+mday(month)+iday-1)*24d0+hour-plsgmt
      centj = (hrgmt+12d0)/(36525d0*24d0)
      return
      end

