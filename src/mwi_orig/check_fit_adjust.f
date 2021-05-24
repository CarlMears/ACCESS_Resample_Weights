      include 'o:\mwi\routines\mwi_sensor_module.f'   

	program check_fit_adjust	 	  	 		   		   
	use mwi_sensor_module 	  	 		   		   
	implicit none

	integer(4), parameter:: nlat=1601
	integer(4), parameter:: nlon=2401  !ssmis value, probably overkill for gmi conisered it low altitude

	integer(4), parameter :: npixel_max=100000
	integer(4), parameter :: nobs_max=1000

      integer(4), parameter :: kcel_max=58  

c     for this program all gain resampling is for G11
      integer(4), parameter :: icase_trg=1
      integer(4), parameter :: icase_obs=1

      real(8), parameter :: pi=3.141593d0

c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6

      character(120) filename
      character(120) filename1,filename2

	integer(2) ilat2(npixel_max,0:nobs_max),ilon2(npixel_max,0:nobs_max)
	integer(4) npixel(0:nobs_max)
	real(8) gobs(npixel_max,0:nobs_max),gobs1(nlon,nlat),gobs_trg(nlon,nlat),gobs_sim(nlon,nlat)  

	integer(4) iscansv(nobs_max),icelsv(nobs_max)
	integer(4) iscan,icel,jcel,kcel,nobs,iobs,ilat,ilon,ipixel	   
	integer(4) ifreq,itrg
	integer(1) index_wt(-kcel_max:kcel_max,-14:14,maxcel)


      real(8) weight(maxcel,-kcel_max:kcel_max,-14:14,mfreq)  !index order is a bit weird, but this is what peter has been using
 	real(8) factor_noise(maxcel,mfreq)
 	real(8) factor_area (maxcel,mfreq)

	real(8) xsum
	real(8) trg_max,sim1_max

	integer(4) isum_trg,isum_sim
	real(8) area_trg,area_sim

      integer(4) iopt_cel,kscan,icel_trg
	
 	character(256) acmd(100)
	integer(4) numarg
	character(6) awindow

							  
	call get_cmd_arg( numarg,acmd)
	if(numarg.ne.4) stop 'wrong number of command line arguments, pgm stopped'
	read(acmd(1),'(i6)') iopt_cel
	read(acmd(2),'(i6)') kscan
	read(acmd(3),'(i6)') itrg
	read(acmd(5),'(a6)') awindow
	
	write(filename1,9001) iopt_cel,kscan,itrg
	write(filename2,9002) iopt_cel,kscan,itrg
 9001 format('check_fit_adjust_',3i2.2 ,'.lis')
 9002 format('check_fit_adjust_',3i2.2 ,'.dat')

	open(6,   file=filename1,status='new')
	call openbig(4,filename2,       'new')
	

	write(filename,9003) iopt_cel,kscan,itrg
 9003 format('O:\mwi\resampling\data\find_resampling_weights_adjust_',3i2.2 ,'.dat')
      write(*,*) filename
	call openbig(3,filename,'old')
	read(3) weight
	close(3)

      weight(:,:,:,5)=weight(:,:,:,4)  !set second 37 GHz to first 37 GHz
	
	call openbig(3,'O:\mwi\resampling\find_neighbors.dat','old')
	read(3) index_wt
	close(3)
	
	do jcel=1,maxcel      ! this is the target 

      if(iopt_cel.eq.1) then
      icel_trg=2*jcel -1
      else
      if(jcel.eq.maxcel) cycle
      icel_trg=2*jcel
      endif

	do ifreq=1,mfreq
	if(ifreq.eq.5) cycle  !second 37 ghz assumed to be same as first and there are not gain patterns for second 
	
	iobs=0
	call read_gobs_adjust(1,npixel_max,nobs_max,itrg,ifreq,icase_trg,kscan,icel_trg,iobs, npixel,ilat2,ilon2,gobs)	  !first argument = 1 denotes target gains

	iobs=0
	do iscan=-14,14
	do kcel=-kcel_max,kcel_max
	if(index_wt(kcel,iscan,jcel).eq.0) cycle
	icel=kcel + jcel
	if(icel.lt.1 .or. icel.gt. maxcel) stop 'pgm stopped, icel oob'
	iobs=iobs+1
	if(iobs.gt.nobs_max) stop 'pgm stopped, iobs oob'
	call read_gobs_adjust(0,npixel_max,nobs_max,itrg,ifreq,icase_obs,iscan,icel,iobs, npixel,ilat2,ilon2,gobs)  !first argument = 0 denotes obs gains
      iscansv(iobs)=iscan
	icelsv(iobs)=icel
	enddo	!iscan
	enddo	!kcel 
	nobs=iobs


	gobs_sim=0
	do iobs=0,nobs

	gobs1=0
	do ipixel=1,npixel(iobs)
	ilat=ilat2(ipixel, iobs)
	ilon=ilon2(ipixel, iobs)
	gobs1(ilon,ilat)=gobs(ipixel,iobs)
	enddo

	if(iobs.eq.0) then
	gobs_trg=gobs1
	else
      iscan=iscansv(iobs)
	icel=icelsv(iobs)
	kcel=icel-jcel
	gobs_sim=gobs_sim + weight(jcel,kcel,iscan,ifreq)*gobs1
	endif

	enddo	 !iobs
	
	trg_max=maxval(gobs_trg)
	sim1_max=maxval(gobs_sim)
 
 	xsum=0;
	isum_trg=0; isum_sim=0
	do ilat=1,nlat
	do ilon=1,nlon
	xsum=xsum + abs(gobs_sim(ilon,ilat)-gobs_trg(ilon,ilat))
	if(gobs_trg(ilon,ilat).gt.0.5*trg_max)  isum_trg=isum_trg + 1 
	if(gobs_sim(ilon,ilat).gt.0.5*sim1_max) isum_sim=isum_sim + 1
	enddo
	enddo

	area_trg=1.11*1.11*isum_trg
	area_sim=1.11*1.11*isum_sim

 	write(*,6001) itrg,jcel,ifreq,xsum,2*sqrt(area_trg/pi),2*sqrt(area_sim/pi)
	write(6,6001) itrg,jcel,ifreq,xsum,2*sqrt(area_trg/pi),2*sqrt(area_sim/pi)
 6001 format(3i5,15f10.5)

      factor_noise(jcel,ifreq)=xsum
      factor_area( jcel,ifreq)=2*sqrt(area_sim/pi)
      
	enddo	!ifreq
	enddo	!jcel

	write(4) factor_noise
	write(4) factor_area
	close(4)

	stop 'norm end'								   
	end





      include 'O:\mwi\resampling\read_gobs_adjust.f'
      include 'x:\syspgm\get_cmd_arg.f'
	include 'x:\syspgm\openbig.for'
