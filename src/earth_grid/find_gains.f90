!     earth rotation is turned off because it has a dif effect at dif points in orbit on the wts

    include '/mnt/oserver/o/amsr_l2/v04/routines/amsre_l2_module_v04.f'              
                                                                            
    program find_gains                                                                                                          
    use amsre_l2_module                                                      
                                                         
    implicit none                                               

    character(120) filename
    integer(4) iorbit,iexist
    integer(4) ifreq,ifreq1,ifreq2,itarget

    real(8) time_first
    integer :: ioerr,lu 
    character(len =80) :: iomsg


    write(*,*) 'enter 0 to do obs, enter 1 to do targets'
    read(*,*) itarget

    
    write(*,*) 'enter ifreq1,ifreq2'
    read(*,*) ifreq1,ifreq2
    
    write(filename,9001) itarget,ifreq1,ifreq2
 9001 format('find_gains_',i1,'_',i1,'_',i1,'.lis')                 
    open(6,file=filename)

    iorbit=3000

    call read_amsr_l0(iorbit, iexist,time_first)

    do ifreq=ifreq1,ifreq2  !ifreq=7 is 89a, ifreq=8 is 89b
    if(ifreq.eq.6) cycle
    if(itarget.eq.1 .and. ifreq.eq.4) cycle  !23 ghz is never a target
    call amsre_gains(itarget,ifreq)
    enddo                                                
                                                                    
    stop 'norm end'                                   
    end


!     in this routine ifreq refers to antenna pattern and kfreq refers to location
    subroutine amsre_gains(itarget,ifreq)
    use amsre_l2_module
    implicit none

!    these are values chosen to mimic nasda's geolocation calculations.
    integer(4),parameter :: max_obs_low =243
    integer(4),parameter :: max_obs_high=486
!    earth dimensions match wgs84          
    real(8), parameter :: re=6378.137d3
    real(8), parameter :: rp=6356.824d3
    real(8), parameter:: dif_epoch = -2556.5d0      !=jdofepoch  - 2451545.d0, jdofepoch = 2448988.5d0 !jan, 1 1993 ut=0 hours

    integer(4), parameter:: iscan0=1064   !I guess, the scan near the Equator for orbit 3000
    integer(4), parameter:: nlat=1601
    integer(4), parameter:: nlon=2101
    real(4),    parameter:: xlat0= -8
    real(4),    parameter:: xlon0=258.

!     the values of these amsr-e parameters come from e:/amsr/resampling_weights/pm1_amsr_mod.for
!         the values of th midori2   parameters come from e:/amsr/resampling_weights/adeosii_amsr_mod.for
!         for these antenna parameters, peter found one set for 89 ghz, which was the averaged of 89a and 89b
!         this value was put in the ifreq=8 slot
!         the ifreq=6 and 7 slots were set to 0
!         for midori2, the ifreq=6 and 7 slots contain the 50.3 and  52.8 patterns respectively.

!         here, i use slot 7 for 89a and slot 8 for 89b, so they have the same value
!         when i do midori-2, i will need to decide what to do about the 50 ghz channels (maybe put in average value into slot 6)
  
    real(8), parameter :: sin_anglim(8)=(/  0.105d0,  0.105d0,  0.0524d0, 0.0524d0, 0.0262d0, 0.000d0,  0.0262d0 ,  0.0262d0 /)
    real(8), parameter :: ant_approx_a(8)=(/4.343d-6, 2.096d-6,    1.890d-6, 1.623d-6, 7.248d-7, 0.000d0,  2.070d-6,   2.070d-6 /)
    real(8), parameter :: ant_approx_b(8)=(/6.892d-4, 4.059d-4,    3.727d-4, 7.251d-4, 3.051d-4, 0.000d0,  2.381d-4,   2.381d-4 /)
    real(8), parameter :: ant_approx_c(8)=(/0.503d0,  0.662d0,  1.391d0,  1.804d0,  1.964d0,  0.000d0,  4.593d0,    4.593d0  /)
    real(8), parameter :: ant_approx_d(8)=(/0.651d0,  1.345d0,     4.844d0,  4.721d0,  15.18d0,  0.000d0, 79.785d0,   79.785d0  /)

    character(80) filename

    real(8) high_res_step_time
    real(8) start_time(8),beta(8),azoffset(8)
    real(8) x_eci2000(3),y_eci2000(3),z_eci2000(3),x_eci(3),y_eci(3),z_eci(3),alpha
    real(8) scpos_eci2000(3),scpos_eci(3),scvel_eci2000(3),attdq_eci2000(4),days
    real(8) precession(3,3),nutation(3,3),np(3,3)
    real(8) observation_time, time_since_scan_start,xfac
    real(8) xcel,xlook(3),boresight(3),costht,gain,range,cosbeta,sinbeta,bx,by,bz,xlatc,xlonc
    real(8) delta,cosdelta
    real(8) coeff_a,coeff_b,coeff_c,coeff_d,anglim
    real(4) omega_avg,asum
    real(4) pitch_corr,rolloffset
    integer(4) itarget,ifreq,jfreq,kfreq
    integer(4) iscan,iscan1,iscan2,jcel,kscan,iasum,ibad,interval,max_obs
    integer(4) imode
    integer(4) ilat,ilon,istart,i
    integer(2) ilat2,ilon2

    real(8) cell0(3,nlon,nlat),cell(3,nlon,nlat)
    real(8) qsum,sumgain(8),gain0,deltasv,xmin,xmax
    real(4) gsum(nlon,nlat)

    integer :: ioerr,lu 
    character(len =80) :: iomsg


    data imode/1/  !1 means only due precession; otherwise do both precession and nutation
    data high_res_step_time/1.3d-3/
    data start_time/7*0.438425d0, 0.436150d0/

    data beta/47.67d0, 47.64d0, 4*47.57d0 , 47.57d0, 47.09d0/  
    data azoffset/-0.06d0, -0.17d0, 4*-0.17d0, -0.37d0, -0.10d0/  

    data rolloffset/0.09/
    data istart/1/

!   ========================================================================================================================
!   ============================================== begin starting block ====================================================
!   ========================================================================================================================


    !calculate vectors from denter of Earth to the locations
    !xlat0= -8, xlon0 = 258.
    if(istart.eq.1) then
        istart=0
        do ilat=1,nlat
            do ilon=1,nlon
                xlatc=0.01d0*(ilat-1) + xlat0
                xlonc=0.01d0*(ilon-1) + xlon0
                call fd_cell_vector(nlat,nlon,ilat,ilon,xlatc,xlonc, cell0,cell)
            enddo
        enddo

        !calculate anteann patterns in angle space
        do jfreq=1,8
            anglim=asind(sin_anglim(jfreq))
            coeff_a =  ant_approx_a(jfreq)
            coeff_b =  ant_approx_b(jfreq)
            coeff_c =  ant_approx_c(jfreq)
            coeff_d =  ant_approx_d(jfreq)

            qsum=0; deltasv=-1
            do i=0,70000
                delta=0.0001d0*i
                if (delta.le.anglim) then
                    gain = coeff_a + coeff_b*exp(-coeff_c*delta) +     exp(-coeff_d*delta*delta)
                    if(i.eq.0) gain0=gain
                else
                    gain=0
                endif
                qsum=qsum + gain*sind(delta)
                if(abs(gain/gain0 -0.5).lt.0.001) deltasv=delta
            enddo !i


            !print some checks
            sumgain(jfreq)=qsum*6.283185d0*1.745329d-6      !two pi times integration step in radians (0.0001 deg)

            write(*,1001) jfreq, sumgain(jfreq), 2*deltasv, &
                10*dlog10((coeff_a + coeff_b*exp(-coeff_c*anglim) + exp(-coeff_d*anglim*anglim))/gain0) 
            write(6,1001) jfreq, sumgain(jfreq), 2*deltasv, &
                 10*dlog10((coeff_a + coeff_b*exp(-coeff_c*anglim) + exp(-coeff_d*anglim*anglim))/gain0) 
                1001 format(i3,e15.5,f8.3,f8.2)
        enddo !jfreq
    endif

!   ========================================================================================================================
!   ==============================================  end  starting block ====================================================
!   ========================================================================================================================

    if(itarget.eq.0) then
        kfreq=ifreq
    else
        kfreq=5  !for the targets always use location of 37 ghz
    endif



    anglim=asind(sin_anglim(ifreq))
    coeff_a = ant_approx_a(ifreq)
    coeff_b = ant_approx_b(ifreq)
    coeff_c = ant_approx_c(ifreq)                           
    coeff_d = ant_approx_d(ifreq)

    if(kfreq.le.6) then
        max_obs=max_obs_low
        xfac=2
    else
        max_obs=max_obs_high
        xfac=1
    endif
    xmin=1.e30; xmax=-1.e30
    iscan1=iscan0 - 14
    iscan2=iscan0 + 14
     
    do 200 iscan=iscan1,iscan2 
        if(itarget.eq.1 .and. iscan.ne.iscan0) cycle
        ibad=0; iasum=0; asum=0
        do kscan=iscan-5,iscan+5
            if(kscan.lt.1 .or. kscan.gt.numscan) cycle
            if(iflag_l0(kscan).eq.10) cycle     !no science data
            if(omega(kscan).le.0) then         !bad omega
                ibad=ibad+1
            else
                iasum=iasum+1
                asum=asum+omega(kscan)
            endif
        enddo

        if(ibad.gt.1) iasum=0 !if there is more than one bad omega near scan, do not process
        if(iasum.eq.0 .and. iflag_l0(iscan).eq.0) iflag_l0(iscan)=8 
        if(iflag_l0(iscan).ne.0) stop 'i want all scans to be good, pgm stopped'
        omega_avg=asum/iasum
        pitch_corr=0


        do 100 jcel=1,max_obs
    
            gsum=0

            do 90 interval=0,9
                xcel=-0.45d0 +0.1d0*interval

                time_since_scan_start = (jcel-1+xcel)*xfac*high_res_step_time
                observation_time = ut1_1993(iscan) + time_since_scan_start
                days=observation_time/86400.d0 + dif_epoch 

                scpos_eci2000 = scpos(:,iscan) + time_since_scan_start*scvel(:,iscan)
                scvel_eci2000=  scvel(:,iscan)    !is only used for jcel=1, thus it does not need to change
                attdq_eci2000 = attdq(:,iscan) + time_since_scan_start*dattdq(:,iscan)

                call get_sc_axes(attdq_eci2000, x_eci2000,y_eci2000,z_eci2000)
                call roll_sc(rolloffset, x_eci2000,y_eci2000,z_eci2000)
                if(pitch_corr.ne.0) call pitch_sc(pitch_corr, x_eci2000,y_eci2000,z_eci2000)


                if(imode.eq.1) then
                    call fd_precession(days, np)
                else
                    call fd_precession(days, precession)
                    call fd_nutation(  days, nutation)
                    np=matmul(nutation,precession)
                endif

                scpos_eci=matmul(np,scpos_eci2000)
                x_eci    =matmul(np,x_eci2000)
                y_eci    =matmul(np,y_eci2000)
                z_eci    =matmul(np,z_eci2000)

                alpha=omega_avg*(time_since_scan_start + start_time(kfreq)) + azoffset(kfreq) - 180.

                cosbeta=cosd(beta(kfreq))                                                                      
                sinbeta=sind(beta(kfreq))                                                                      
                bx=sinbeta*cosd(-alpha) !input alpha defined relative to the -z axes                                       
                by=sinbeta*sind(-alpha) !input alpha defined relative to the -z axes                                       
                bz=cosbeta                 
                boresight(1)=bx*x_eci(1)+by*y_eci(1)+bz*z_eci(1)                            
                boresight(2)=bx*x_eci(2)+by*y_eci(2)+bz*z_eci(2)                            
                boresight(3)=bx*x_eci(3)+by*y_eci(3)+bz*z_eci(3)

                qsum=0

                do ilat=1,nlat
                    do ilon=1,nlon

                        xlook= cell(:,ilon,ilat) - scpos_eci       !point from s/c down to cell
                        range=sqrt(dot_product(xlook,xlook))
                        xlook=xlook/range

                        costht=-dot_product(xlook,cell0(:,ilon,ilat))
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

            90  continue
                gsum=0.1*gsum !0.1 is dividing by 10 samples, i.e. interval=0,9    results in 10 samples per bin
                if(itarget.eq.0) then
                write(filename,9001) ifreq,iscan-iscan1+1,jcel
            9001 format('c:/gain_arrays/freq',i1,'/s',i2.2,'c',i3.3,'.dat')
                else
                write(filename,9002) ifreq,iscan-iscan1+1,jcel
            9002 format('c:/target_arrays/freq',i1,'/s',i2.2,'c',i3.3,'.dat')
                endif
                write(*,1003) filename
                write(6,1003) filename
            1003 format(1x,a80)

            !call openbig(4,filename,'new')
            open(unit=4,file=filename, status='new',action='write',form='unformatted', access='stream',iostat=ioerr,iomsg=iomsg)

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

  100   continue
  200 continue

    write(*,1004) itarget,ifreq,kfreq,xmin,xmax
 1004 format(' completed ',3i5,2f10.5)
    return
    end subroutine amsre_gains

    subroutine fd_cell_vector(nlat,nlon,ilat,ilon,xlat,xlon, cell0,cell)

    implicit none

    real(8), parameter :: re=6378.137d3
    real(8), parameter :: rp=6356.824d3

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
                                                           
    include '/mnt/oserver/o/amsr_l0/filename_routines_v04.f'      
    include '/mnt/oserver/o/amsr_l0/read_binary_head.f'
    include '/mnt/oserver/o/amsr_l2/read_amsr_l0_v05.f'                             
    include '/mnt/oserver/o/amsr_l2/find_therm_v03.f'

    include '/mnt/oserver/o/amsr_l2/v04/routines/amsre_geolocation_v05.f'                 
    include '/mnt/oserver/o/amsr_l2/v04/routines/pre_nut_routines.f'
    include '/mnt/oserver/o/amsr_l2/gelocation_error/utc_ut1_routines.f'
    include '/mnt/oserver/o/sun_moon/sunloc1.f'            

    ! include 'x:/readgzip/gzip.f'
    ! include 'x:/syspgm/openbig.for'
    ! include 'x:/syspgm/openbig_try.f'
    ! include 'x:/syspgm/open_try_append.f'
    ! include 'x:/syspgm/systime.f'                  
    ! include 'x:/syspgm/fd_time87.f'
    ! include 'x:/syspgm/time87.for'              
    ! include 'x:/syspgm/xtime87.f'
