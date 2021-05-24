    !include 'o:\mwi\routines\mwi_sensor_module.f'   This info should be in simulation_param now....   

	program find_resampling_weights 
    !use mwi_sensor_module 
    use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
    use simulation_param, only: NBAND,FOVS_PER_SCAN,NFREQ,SCAN_PERIOD,ORBIT_PERIOD,SAMPLINGTIME
	implicit none

!   Frank's Comments:
!     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
!     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
    ! Carl: changed all references to mfreq to NBAND
    
    !integer(int32), parameter:: mfreq=6

	integer(int32), parameter:: nlat=1601
	integer(int32), parameter:: nlon=2401  !ssmis value, probably overkill for gmi conisered it low altitude 

	integer(int32), parameter:: npixel_max=100000  !max number of pixels saved for each observation footprint
	integer(int32), parameter:: nobs_max=1000      !max number of obs (source footprints) saved

    integer(int32), parameter :: kcel_max=58       !number of FOVS along scan in resampling area
      
!     for this program all gain resampling is for G11

    integer(int32), parameter :: icase_trg=1
    integer(int32), parameter :: icase_obs=1
     
    character(120) filename1,filename2

	integer(int16) ilat2(npixel_max,0:nobs_max),ilon2(npixel_max,0:nobs_max)   !list of ilat,ilon locations for footprint data
	integer(int32) npixel(0:nobs_max)                                          ! not sure yet
    real(real64) gobs(npixel_max,0:nobs_max)                                    ! list of gains at locations (this is how the footprints are stored)
    real(real64) gobs1(nlon,nlat)                                               ! gains in lon/lat space

	integer(int32) iscansv(nobs_max),icelsv(nobs_max)                          ! not sure yet
	integer(int32) iscan,icel,jcel,kcel,nobs,iobs,jobs,ilat,ilon,ipixel        ! various incideces
	integer(int32) ierror,iobssv
	integer(int32) ifreq,itrg
	integer(1) index_wt(-kcel_max:kcel_max,-14:14,maxcel)                  ! = 1 if a neighbor??

	real(real64) xsum
	real(real64), allocatable :: g(:,:),g_inv(:,:),testg(:,:),u(:),v(:),xwork(:),a(:) 
    real(real64) denom,xcoef,xmax,psum(0:4),smooth_fac
    
    !need to change this to allocatable

    real(real64),allocatable, weight(:,:,:) !maxcel,-kcel_max:kcel_max,-14:14,NBAND)  !index order is a bit weird, but this is what peter has been using
      
    integer(int32) iopt_cel,kscan,icel_trg
      
 	character(256) acmd(100)
	integer(int32) numarg
    character(6) awindow
    
    allocate(weight(maxcel,-kcel_max:kcel_max,-14:14,NBAND))

	call get_cmd_arg( numarg,acmd)
	if(numarg.ne.4) stop 'wrong number of command line arguments, pgm stopped'
	read(acmd(1),'(i6)') iopt_cel
	read(acmd(2),'(i6)') kscan
	read(acmd(3),'(i6)') itrg
	read(acmd(4),'(a6)') awindow

	write(filename1,9001) iopt_cel,kscan,itrg
	write(filename2,9002) iopt_cel,kscan,itrg
 9001 format('find_resampling_weights_',3i2.2 ,'.lis')
 9002 format('find_resampling_weights_',3i2.2 ,'.dat')

	open(6,   file=filename1,status='new')
	call openbig(4,filename2,       'new')
    
    
	call openbig(3,'O:\mwi\resampling\find_neighbors.dat','old')
	read(3) index_wt
	close(3)
	

	weight=0

	do ifreq=1,NBAND
	
        smooth_fac=1.0d-05
        
        !MWI specific!!
	    if(ifreq.eq.1 .and. itrg.ge.2) smooth_fac=1.0d-06  !11 ghz mapping onto 30,25,20 km,use smaller smoothing to presever 30 km.
 
	    do jcel=1,maxcel !this is the target 
            if(iopt_cel.eq.1) then
                icel_trg=2*jcel -1
            else
            if(jcel.eq.maxcel) cycle
                icel_trg=2*jcel
            endif
       
	        iobs=0
	        call read_gobs(1,npixel_max,nobs_max, itrg,icase_trg,kscan,icel_trg,iobs, npixel,ilat2,ilon2,gobs)	  !first argument = 1 denotes target gains
            ! the target gain is stored in the first obs position (0) in gobs -- remember gobs is (1:max_pixel,0:max_obs)
	
c     =============================================================================================================================
c     ========================================= compute iscan and icel that correspond to each iobs ===============================
c     =============================================================================================================================

	        iobs=0
	        do iscan=-14,14
	            do kcel=-kcel_max,kcel_max
	                if(index_wt(kcel,iscan,jcel).eq.0) cycle
	                icel=kcel + jcel
	                if(icel.lt.1 .or. icel.gt. maxcel) stop 'pgm stopped, icel oob'
	                iobs=iobs+1
	                if(iobs.gt.nobs_max) stop 'pgm stopped, iobs oob'
	                call read_gobs(0,npixel_max,nobs_max,ifreq,icase_obs,iscan,icel,iobs, npixel,ilat2,ilon2,gobs)  !first argument = 0 denotes obs gains
                    ! the source gains are stored in iobs positions 1..nobs in gobs

                    iscansv(iobs)=iscan
	                icelsv(iobs)=icel
	            enddo	!kcel
	        enddo	!iscan
	        nobs=iobs

            ! allocate the matrics for the linear algebra steps.
	        allocate(g(nobs,nobs),g_inv(nobs,nobs),testg(nobs,nobs), u(nobs), v(nobs), xwork(nobs), a(nobs))

	        do iobs=1,nobs
                gobs1=0
	            xsum=0
	            do ipixel=1,npixel(iobs)
	                ilat=ilat2(ipixel, iobs)
	                ilon=ilon2(ipixel, iobs)
	                gobs1(ilon,ilat)=gobs(ipixel,iobs)
	                xsum=xsum + gobs1(ilon,ilat)
	            enddo !ipixel
	
	            if(abs(xsum-1.d0).gt.1.e-12) then
	                write(*,*) xsum
	                stop 'error in normalization, pgm stopped'
	            endif

	            do jobs=0,nobs  !THIS loop included the target footprint jobs == 0

	                xsum=0
	                do ipixel=1,npixel(jobs)  !loop over the non-zero pixels in ob j
	                    ilat=ilat2(ipixel, jobs)
	                    ilon=ilon2(ipixel, jobs)
	                    xsum=xsum + gobs1(ilon,ilat)*gobs(ipixel,jobs)
	                enddo

	                if(jobs.eq.0) then
	                    v(iobs)=xsum
	                else
	                    g(iobs,jobs)=xsum
	                endif
	            enddo	 !jobs
	        enddo	 !iobs

            !not sure if i need this, g may already be symmetric
	        do iobs=1,nobs
	            do jobs=iobs,nobs
	                g(jobs,iobs)=g(iobs,jobs)
	            enddo
	        enddo


	        u=1
            ! add in the smoothing factor along the diagonal
	        do iobs=1,nobs
	            g(iobs,iobs)=g(iobs,iobs) + smooth_fac
	        enddo
 
            ! invert
	        call svdinv(g,nobs, g_inv,ierror)
	        if(ierror.ne.0) stop 'inversion did not converge, pgm stopped'

	        testg=matmul(g,g_inv)
	        do iobs=1,nobs
	            testg(iobs,iobs)=testg(iobs,iobs)-1
            enddo

            xwork=matmul(g_inv,u)
	        denom=dot_product(u,xwork)

            xwork=matmul(g_inv,v)
	        xcoef=dot_product(u,xwork)							    

	        xwork=v + (1-xcoef)*u/denom

	        a=matmul(g_inv,xwork)

	        xmax=-1.e30; psum=0
	        do iobs=1,nobs
	            iscan=iscansv(iobs)
	            icel=icelsv(iobs)
	            kcel=icel -	jcel
	            weight(jcel,kcel,iscan,ifreq)=a(iobs)
	
	            psum(0)=psum(0) + 1
 	            psum(1)=psum(1) + a(iobs)
 	            psum(2)=psum(2) + a(iobs)**2
 	            psum(3)=psum(3) + a(iobs)*iscan
 	            psum(4)=psum(4) + a(iobs)*icel
                if(a(iobs).gt.xmax) then
	                xmax=a(iobs)
	                iobssv=iobs
	            endif
	        enddo

	        psum(2)=sqrt(psum(2))
	        psum(3:4)=psum(3:4)/psum(1)

            write(*,6001) itrg,jcel,ifreq,nobs,maxval(abs(testg)),psum(1:2),xmax, &
                            iscansv(iobssv),icelsv(iobssv)-jcel,psum(3),psum(4)-jcel 
            write(6,6001) itrg,jcel,ifreq,nobs,maxval(abs(testg)),psum(1:2),xmax, &
                            iscansv(iobssv),icelsv(iobssv)-jcel,psum(3),psum(4)-jcel 
            6001 format(4i5,e15.5,3f13.8,2i5,2f8.3)

	        deallocate(g, g_inv, testg, u, v, xwork, a)                                                         

	    enddo	  !jcel
	enddo	  !ifreq

	write(4) weight
	close(4)
	close(6)
    deallocate(weight)
	stop 'norm end'								   
	end

    include 'O:\mwi\resampling\read_gobs.f'
    include 'O:\mwi\resampling\invert_matrix.f'
    include 'x:\syspgm\get_cmd_arg.f'
	include 'x:\syspgm\openbig.for'
