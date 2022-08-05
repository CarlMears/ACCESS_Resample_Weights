    !Started from O:\ssmis\resampling\find_gains.f

    !reformatted/indented so I can read
    
    
    include 'o:\ssmis\routines\ssmis_module.f'
    include 'o:\ssmis\routines\ssmis_l2_module.f'
 
    program find_gains                                        
    use ssmis_module
    use ssmis_l2_module
    implicit none                                            

    integer(4) iorbit,iexist
    integer(4) ifreq,itarget,isat
    real(8) time_first
 
    write(*,*) 'enter ssmis number (16, 17, ...)'
!    read(*,*) isat
    isat=16
    if(isat.ne.16 .and. isat.ne.16) stop 'prog stopped, isat oob'
 
    call define_filenames(isat)
 
    !read in sample orbit
    open(6,file='ssmis_find_gains.lis')
    iorbit=14400
    call read_l1a_file(iorbit, iexist,time_first)
 
    do itarget=0,1
        do ifreq=1,3
            if(itarget.eq.1 .and. ifreq.ne.1) cycle  !only target is 19 ghz    pattern located at 37 ghz location
            call ssmis_gains(itarget,ifreq)
        enddo
    enddo
    
    stop 'norm end'
    end
                                 
 
                         
 
!     in this routine ifreq refers to antenna pattern and kfreq refers to location
    subroutine ssmis_gains(itarget,ifreq)
    use ssmis_module
    use ssmis_l2_module
    implicit none
 
    integer(4), parameter :: iscan0 =922
    integer(4), parameter:: nlat=1601
    integer(4), parameter:: nlon=2401
    real(4),  parameter :: xlat0= -8.
    real(4),  parameter :: xlon0=243.
    real(8), parameter :: coeff_a=2.d-6
    real(8), parameter :: coeff_b=4.d-4
    real(8), parameter :: thtfull(3)=(/2.00, 2.00, 1.22/)  !see memo2 in O:\ssmis\resampling
 
    character(100) filename
 
    real(8) xcel,xlook(3),boresight(3),costht,gain,range,cosbeta,sinbeta,bx,by,bz,xlatc,xlonc
    real(8) delta,cosdelta
    real(8) coeff_c,coeff_d,anglim,ththalf
 
    real(8) x_eci(3),y_eci(3),z_eci(3),l_eci(3)
    real(8) scpos_eci(3)
    real(8) observation_time,time_since_scan_start
    real(8) scpos0(3),scvel0(3)
    real(8) r,ru(3),sinlat,coslat,deltau,alpha,rearth
 
    real(8) cell0(3,nlon,nlat),cell(3,nlon,nlat)
    real(8) qsum,sumgain(3),gain0,deltasv,xmin,xmax
    real(4) gsum(nlon,nlat)

 
    integer(4) itarget,ifreq,jfreq,kfreq,igrp
    integer(4) iscan,iscan1,iscan2,icel,interval
    integer(4) ilat,ilon,istart,i
    integer(2) ilat2,ilon2
 
 
    data istart/1/
 
!     ========================================================================================================================
!     ============================================== begin starting block ====================================================
!     ========================================================================================================================
 
    if(istart.eq.1) then
    istart=0
 
    do ilat=1,nlat
    do ilon=1,nlon
    xlatc=0.01d0*(ilat-1) + xlat0
    xlonc=0.01d0*(ilon-1) + xlon0
    call fd_cell_vector(nlat,nlon,ilat,ilon,xlatc,xlonc, cell0,cell)
    enddo
    enddo
 
    do jfreq=1,3
    ththalf=0.5*thtfull(jfreq)
    anglim=6*ththalf
    coeff_d = -alog(0.5)/ththalf**2
    coeff_c =  sqrt(0.3*coeff_d)
 
    qsum=0; deltasv=-1
    do i=0,160000
    delta=0.0001d0*i
      if (delta.gt.anglim) exit
    gain = coeff_a + coeff_b*exp(-coeff_c*delta) +     exp(-coeff_d*delta*delta)
    if(i.eq.0) gain0=gain
    qsum=qsum + gain*sind(delta)
    if(abs(gain/gain0 -0.5).lt.0.0001) deltasv=delta
    enddo !i
 
    sumgain(jfreq)=qsum*6.283185d0*1.745329d-6      !two pi times integration step in radians (0.0001 deg)
 
    write(*,1001) jfreq, sumgain(jfreq), 2*deltasv,gain0,10*dlog10(gain)
    write(6,1001) jfreq, sumgain(jfreq), 2*deltasv,gain0,10*dlog10(gain)
 1001 format(i3,e15.5,f8.3,2f8.3)
    enddo !jfreq
 
    endif
 
!     ========================================================================================================================
!     ==============================================  end  starting block ====================================================
!     ========================================================================================================================
 
    if(itarget.eq.0) then
    kfreq=ifreq
    else
    kfreq=3  !for the targets always use location of 37 ghz
    endif

    if(kfreq.le.2) igrp=3
    if(kfreq.eq.3) igrp=4
 
    ththalf=0.5*thtfull(ifreq)
    anglim=6*ththalf
    coeff_d = -alog(0.5)/ththalf**2
    coeff_c =  sqrt(0.3*coeff_d)
 
    xmin=1.e30; xmax=-1.e30
 
    iscan1=iscan0 - 14
    iscan2=iscan0 + 14
    
    do 200 iscan=iscan1,iscan2

    if(iscan_flag(iscan).ne.0) stop 'pgm stopped,iscan_flag(iscan).ne.0'
 
    if(itarget.eq.1 .and. iscan.ne.iscan0) cycle

!     =================================================================================================================
!     ================== compute spacecraft (sc) y coordinate vector (y0) which varies very slowly with time ==========
!     =================================================================================================================

    scpos0=scpos(:,iscan)
    scvel0=scvel(:,iscan)
    call cross_norm(scvel0,scpos0, l_eci)
 
!     =================================================================================================================
!     ======================================== loop thru cells ========================================================
!     =================================================================================================================
 
    do 100 icel=1,maxcel_grp(igrp)
      gsum=0
 
    do 90 interval=0,9
 
    xcel= icel + (-0.45 + 0.1*interval)

      alpha=alpha_start(igrp) + (xcel-1)*(  144.d0/maxcel_grp(igrp))
    time_since_scan_start =   (xcel-1)*(0.7596d0/maxcel_grp(igrp))   !0.7596=(144/360)*1.899
     observation_time = scan_time(iscan) + time_since_scan_start
 
     scpos_eci = scpos0 + time_since_scan_start*scvel0
 
!     =================================================================================================================
!     ========================================== find spacecraft x, y, and z coordinate vectors =======================
!     =================================================================================================================

    r=sqrt(dot_product(scpos_eci,scpos_eci))
    ru=scpos_eci/r
      sinlat=ru(3)
      coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
      rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
    deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
     z_eci(1)=-ru(1)
    z_eci(2)=-ru(2)
    z_eci(3)=-ru(3)*(1+deltau)/(ffac+deltau)
    z_eci=z_eci/sqrt(dot_product(z_eci,z_eci))    

    call cross_norm(l_eci,z_eci, x_eci)
    call cross_norm(z_eci,x_eci, y_eci)

      cosbeta=cosd(beta(igrp))
      sinbeta=sind(beta(igrp))
      bx= sinbeta*cosd(alpha) !input alpha defined relative to the x axes,  alpha=-90 denotes bx=0, by=1
      by=-sinbeta*sind(alpha) !input alpha defined relative to the x axes , alpha=-90 denotes bx=0, by=1
      bz=cosbeta                
      boresight=bx*x_eci + by*y_eci + bz*z_eci

    qsum=0
 
    do ilat=1,nlat
    do ilon=1,nlon
 
    xlook= cell(:,ilon,ilat) - scpos_eci       !point from s/c down to cell
      range=sqrt(dot_product(xlook,xlook))
    xlook=xlook/range
 
      costht=-dot_product(xlook,cell0(:,ilon,ilat)) !near equator cell0 is essentially the same as geoid
    if(costht.le.0) stop 'pgm stopped,costht.le.0'
 
    cosdelta=dot_product(xlook,boresight)
    if(cosdelta.gt.1) cosdelta=1  !should never be near -1
    delta=acosd(cosdelta)
                                                                                     
      if (delta.le.anglim) then
    gain = coeff_a + coeff_b*exp(-coeff_c*delta) +     exp(-coeff_d*delta*delta)
    gain=gain*costht/(range*range)
    gsum(ilon,ilat)=gsum(ilon,ilat) + gain
    qsum           =qsum            + gain
    endif
 
    enddo !ilat
    enddo !ilon
 
    qsum=qsum*1110.d0**2    !area of cell is 1110 meters
    qsum=qsum/sumgain(ifreq)
    if(qsum.lt.xmin) xmin=qsum
    if(qsum.gt.xmax) xmax=qsum
 
   90 continue
 
    gsum=0.1*gsum !0.1 is dividing by 10 samples, i.e. interval=0,9    results in 10 samples per bin
 
    if(itarget.eq.0) then
      write(filename,9001) ifreq,iscan-iscan1+1,icel
 9001 format('c:\ssmis_gain_arrays\freq',i1,'\s',i2.2,'c',i3.3,'.dat')
      else
      write(filename,9002) ifreq,iscan-iscan1+1,icel
 9002 format('c:\ssmis_target_arrays\freq',i1,'\s',i2.2,'c',i3.3,'.dat')
      endif
 
      write(*,1003) filename
      write(6,1003) filename
 1003 format(1x,a100)
 
    call openbig(4,filename,'new')
 
    do ilat=1,nlat
    do ilon=1,nlon
    if(gsum(ilon,ilat).ne.0) then
    ilat2=ilat
    ilon2=ilon
    write(4) ilat2,ilon2,gsum(ilon,ilat)
    endif
    enddo
    enddo
 
    close(4)
 
  100 continue
  200 continue
 
      write(*,1004) itarget,ifreq,kfreq,xmin,xmax
      write(6,1004) itarget,ifreq,kfreq,xmin,xmax
 1004 format(' completed ',3i5,2f10.5)
      return
    end






    subroutine fd_cell_vector(nlat,nlon,ilat,ilon,xlat,xlon, cell0,cell)
 
    implicit none
 
c    earth dimensions match wgs84    
    real(8), parameter :: re=6378.137d3        !meters
      real(8), parameter :: rp=6356.752d3        !meters
 
    integer(4) nlat,nlon,ilat,ilon
    real(8) xlat,xlon,rearth
    real(8) cell0(3,nlon,nlat),cell(3,nlon,nlat)
    real(8) coslat,sinlat,coslon,sinlon
 
    coslat=cosd(xlat)
    sinlat=sind(xlat)
    coslon=cosd(xlon)
    sinlon=sind(xlon)
 
      cell0(1,ilon,ilat)=coslat*coslon
      cell0(2,ilon,ilat)=coslat*sinlon
      cell0(3,ilon,ilat)=sinlat
 
      rearth=re*rp/sqrt((rp*coslat)*(rp*coslat)+(re*sinlat)*(re*sinlat))
 
    cell(:,ilon,ilat)= rearth*cell0(:,ilon,ilat)
 
    return
    end
 

 
 
    include 'o:\ssmis\routines\filename_routines.f'
    include 'o:\ssmis\routines\define_filenames.f'
    include 'o:\ssmis\routines\read_l1a_file.f'
    include 'o:\ssmis\routines\ssmis_geolocation.f'
    include    'o:\amsr_l2\v05\routines\pre_nut_routines.f'
    include 'o:\ssmis\routines\math_routines.f'
     include 'o:\sun_moon\sun_vector.f'            
     include 'o:\sun_moon\moon_vector.f'            
    include 'x:\syspgm\fd_date_2000.f'
    include 'x:\syspgm\ck_file_time.f'
    include 'x:\syspgm\open_file_routines.f'
    include 'x:\syspgm\openbig.for'
 
 
    
