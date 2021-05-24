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
      integer(4), parameter :: kcel_max=58 
      
      real(8), parameter ::   smooth_fac=1.0d-06  
 
	integer(2) ilat2(npixel_max),ilon2(npixel_max) 
	integer(4) npixel
	real(8) gobs(npixel_max),gobs1(nlon,nlat)

	integer(4) iscan,jcel,ilat,ilon,ipixel 
	integer(4) icase,ifreq

	call openbig(4,'show_gains_on_earth.dat','replace')

	
      icase=4	 !m11, m22, m32, m41

	gobs1=0
	
c	do ifreq=1,mfreq
	do ifreq=1,1
	if(ifreq.eq.5) cycle !skipping the second 37 ghz channel for now
	 
	do jcel=1,maxcel,30
	write(*,*) jcel

c     in this program, the target gains is always the observation gain.
	iscan=0
	call read_gobs(npixel_max,ifreq,icase,iscan,jcel, npixel,ilat2,ilon2,gobs)  !first argument = 0 denotes obs gains which are used for targets in the pgm

	do ipixel=1,npixel
	ilat=ilat2(ipixel)
	ilon=ilon2(ipixel)
	gobs1(ilon,ilat)=gobs(ipixel)
	enddo !ipixel
	
	
	enddo	  !jcel
	enddo	  !ifreq
	write(4) gobs1
	close(4)
	stop 'norm end'								   
	end


 	
 	subroutine read_gobs(npixel_max,ifreq,icase,iscan,icel, npixel,ilat2,ilon2,gobs)
	implicit none
	character(80) filename
	integer(2) ilat2x,ilon2x
	integer(4) npixel_max,ifreq,icase,iscan,icel,i,ipixel

	integer(2) ilat2(npixel_max),ilon2(npixel_max)
	integer(4) npixel
	real(8) gobsx(4)  !m11, m22, m32, m41
	real(8) gobs(npixel_max)
	real(8) xsum
	

      write(filename,9001) ifreq,iscan+15,icel
 9001 format('O:\mwi\resampling\mwi_gain_arrays\freq',i1,'\s',i2.2,'c',i3.3,'.dat')


 	call openbig(3,filename,'old')
	ipixel=0; xsum=0
	do i=1,1234567890
	
	read(3,end=100) ilat2x,ilon2x,gobsx
	
	ipixel=ipixel+1
	if(ipixel.gt.npixel_max) stop 'pgm stopped, ipixel oob'
	ilat2(ipixel)=ilat2x
	ilon2(ipixel)=ilon2x
      gobs(ipixel)=gobsx(icase)
	npixel=ipixel
	xsum=xsum+gobs(ipixel)
	enddo
  100	continue
      close(3)
c      gobs(1:ipixel)=gobs(1:ipixel)/xsum
	return
	end
	


	include 'x:\syspgm\openbig.for'
