      include 'o:\mwi\routines\mwi_sensor_module.f'   

	program find_resampling_weights_leakage  
	use mwi_sensor_module 	  	 		   		   
	implicit none

c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6

	integer(4), parameter:: nlat=1601
	integer(4), parameter:: nlon=2401  !ssmis value, probably overkill for gmi conisered it low altitude

	integer(4), parameter:: npixel_max=100000
	integer(4), parameter:: nobs_max=1000
      integer(4), parameter :: kcel_max=58 
      
      real(8), parameter ::   smooth_fac=1.0d-06  
 
	integer(2) ilat2(npixel_max,0:nobs_max),ilon2(npixel_max,0:nobs_max) 
	integer(4) npixel(0:nobs_max)
	real(8) gobs(npixel_max,0:nobs_max),gobs1(nlon,nlat)

	integer(4) iscansv(nobs_max),icelsv(nobs_max)
	integer(4) iscan,icel,jcel,kcel,nobs,iobs,jobs,ilat,ilon,ipixel 
	integer(4) ierror,iobssv
	integer(4) icase,icase_trg,ifreq,icase_obs,icase_obs0
	integer(1) index_wt(-kcel_max:kcel_max,-14:14,maxcel)

	real(8) xsum
	real(8), allocatable :: g(:,:),g_inv(:,:),testg(:,:),u(:),v(:),xwork(:),a(:) 
	real(8) denom,xcoef,xmax,psum(0:4)
      real(8) weight( maxcel,-kcel_max:kcel_max,-14:14,mfreq)  !index order is a bit weird, but this is what peter has been using
      real(8) weight1(maxcel,-kcel_max:kcel_max,-14:14,mfreq)  
      real(8) weight2(maxcel,-kcel_max:kcel_max,-14:14,mfreq)  
      real(8) acoef(  maxcel,-kcel_max:kcel_max,-14:14,mfreq)  
      real(8) dcoef(  maxcel,-kcel_max:kcel_max,-14:14,mfreq)    

	real(8) beamwidth11(mfreq),sumgain11(mfreq)
	real(8) beamwidth22(mfreq),sumgain22(mfreq)
	real(8) xnorm_trg,xnorm_obs
	
	write(*,*) 'enter icase: 1 is 4stk, 2 is 3stk, 3 is 2stk'
	read(*,*) icase
      
	call openbig(3,'O:\mwi\antenna\data\theoretical_beamwidth_sumgain.dat','old')
	read(3) beamwidth11,sumgain11
	read(3) beamwidth22,sumgain22
	close(3)

	call openbig(3,'O:\mwi\resampling\find_neighbors.dat','old') 
	read(3) index_wt
	close(3)
	
	
      if(icase.eq.1) then
	open(6,   file='data\find_resampling_weights_leakage_4stk.lis',status='new')
      call openbig(4,'data\find_resampling_weights_leakage_4stk.dat','new')
      icase_trg= 1 !first stokes
      icase_obs0=3
      endif
      
      if(icase.eq.2) then
	open(6,   file='data\find_resampling_weights_leakage_3stk.lis',status='new')
      call openbig(4,'data\find_resampling_weights_leakage_3stk.dat','new')
      icase_trg= 2 !second stokes
      icase_obs0=5
      endif
      
      if(icase.eq.3) then
 	open(6,   file='data\find_resampling_weights_leakage_2stk.lis',status='new')
      call openbig(4,'data\find_resampling_weights_leakage_2stk.dat','new')
      icase_trg= 2 !second stokes
      icase_obs0=7
      endif

     
      do icase_obs=icase_obs0,icase_obs0+1
      
	weight=0
	
	do ifreq=1,mfreq
	if(ifreq.eq.5) cycle !skipping the second 37 ghz channel for now
	 
      call fd_xnorm(ifreq,icase_trg,  sumgain11,sumgain22, xnorm_trg)
      call fd_xnorm(ifreq,icase_obs,  sumgain11,sumgain22, xnorm_obs)

	do jcel=1,maxcel !this is the target 

c     in this program, the target gains is always the observation gain.
	iobs=0
	iscan=0
	call read_gobs(0,npixel_max,nobs_max,ifreq,icase_trg,iscan,jcel,iobs, npixel,ilat2,ilon2,gobs)  !first argument = 0 denotes obs gains which are used for targets in the pgm

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
      iscansv(iobs)=iscan
	icelsv(iobs)=icel
	enddo	!iscan
	enddo	!kcel 
	nobs=iobs

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


	do jobs=0,nobs

	xsum=0
	do ipixel=1,npixel(jobs)
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

c     not sure if i need this, g may already be symmetric
	do iobs=1,nobs
	do jobs=iobs,nobs
	g(jobs,iobs)=g(iobs,jobs)
	enddo
	enddo


	u=1

	do iobs=1,nobs
	g(iobs,iobs)=g(iobs,iobs) + smooth_fac
	enddo
 
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
	weight(jcel,kcel,iscan,ifreq)=(xnorm_trg/xnorm_obs)*a(iobs)
	
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

	write(*,6001) icase_obs,jcel,ifreq,nobs,maxval(abs(testg)),psum(1:2),xmax,
     & iscansv(iobssv),icelsv(iobssv)-jcel,psum(3),psum(4)-jcel 
	write(6,6001) icase_obs,jcel,ifreq,nobs,maxval(abs(testg)),psum(1:2),xmax,
     & iscansv(iobssv),icelsv(iobssv)-jcel,psum(3),psum(4)-jcel 
 6001 format(4i5,e15.5,3f13.8,2i5,2f8.3)

	deallocate(g, g_inv, testg, u, v, xwork, a)   
	
	enddo	  !jcel
	enddo	  !ifreq
	
	if(icase_obs.eq.icase_obs0)   weight1=weight
	if(icase_obs.eq.icase_obs0+1) weight2=weight
	enddo !icase_obs
	
	acoef=(weight1+weight2)/2
	dcoef=(weight2-weight1)/2  
	write(4) acoef
	write(4) dcoef
	close(4) 
	close(6)

	stop 'norm end'								   
	end


      subroutine fd_xnorm(ifreq,icase,sumgain11,sumgain22, xnorm)
      implicit none
      integer(4) ifreq,icase
 	real(8) sumgain11(6),sumgain22(6)
      real(8) xnorm
 	if(icase.eq.1 .or. icase.eq.3 .or. icase.eq.4) xnorm=sumgain11(ifreq)
 	if(icase.eq.2 .or. icase.eq.5 .or. icase.eq.6) xnorm=sumgain22(ifreq)
 	if(icase.eq.7) xnorm= 2*sumgain22(ifreq) - sumgain11(ifreq) 
 	if(icase.eq.8) xnorm=   sumgain11(ifreq)
 	return
 	end
 
 
      include 'O:\mwi\resampling\read_gobs.f'
      include 'O:\mwi\resampling\invert_matrix.f'
	include 'x:\syspgm\openbig.for'
