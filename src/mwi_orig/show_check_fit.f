      include 'o:\mwi\routines\mwi_sensor_module.f'   

	program show_check_fit	 	  	 		   		   
	use mwi_sensor_module 	  	 		   		   
	implicit none
	
      integer(4), parameter :: kcel_max=58    
      integer(4), parameter :: maxcel_out=2*maxcel-1   

c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6

	integer(4) itrg,ifreq,icel,kcel,iscan,jcel  
	integer(4) kscan,iopt_cel
	
	real(8) weight(maxcel,-kcel_max:kcel_max,-14:14,mfreq)  !index order is a bit weird, but this is what peter has been using
 	real(8) factor_noise(maxcel,mfreq)
 	real(8) factor_area( maxcel,mfreq)
	
	real(8) weight_out(maxcel_out,-kcel_max:kcel_max,-14:14,mfreq)  !index order is a bit weird, but this is what peter has been using
 	real(8) factor_noise_out(maxcel_out,mfreq,4)
 	real(8) factor_area_out( maxcel_out,mfreq,4)
	real(8) xnoise_fac_out(  maxcel_out,mfreq,4)
 	
	real(8) asum(0:2)
	
	
      call openbig(4,'show_check_fit.dat','replace')
      
 	do kscan=0,8
  	
	do itrg=1,4    !target resolution
  
c     ===================================================================
c     ===================== native cells ================================
c     ===================================================================
	iopt_cel=1
      call read_data(itrg,iopt_cel,kscan, weight,factor_noise,factor_area)
      	
	jcel=0
	do icel=1,maxcel_out,2
	jcel=jcel+1
	factor_noise_out(icel,:,itrg)=factor_noise(jcel,:)
	factor_area_out( icel,:,itrg)=factor_area( jcel,:)
	weight_out(icel,:,:,:)=weight(jcel,:,:,:) 
	enddo


c     ===================================================================
c     ===================== synthetic cells ================================
c     ===================================================================
	iopt_cel=2
      call read_data(itrg,iopt_cel,kscan, weight,factor_noise,factor_area)
	
	jcel=0
	do icel=2,maxcel_out-1,2
	jcel=jcel+1
	factor_noise_out(icel,:,itrg)=factor_noise(jcel,:)
	factor_area_out( icel,:,itrg)=factor_area( jcel,:)
	weight_out(icel,:,:,:)=weight(jcel,:,:,:) 
	enddo
	
	
	
	do ifreq=1,mfreq
	if(ifreq.eq.5) cycle  !second 37 ghz assumed to be same as first and there are not gain patterns for second
	

	do icel=1,maxcel_out
	
	if(weight_out(icel,0,0,ifreq).eq.0) stop 'pgm stopped, error1 in weights' 
	
	asum=0
	do kcel=-kcel_max,kcel_max
	do iscan=-14,14
	asum(0)=asum(0)+1
	asum(1)=asum(1) + weight_out(icel,kcel,iscan,ifreq)
	asum(2)=asum(2) + weight_out(icel,kcel,iscan,ifreq)**2
	enddo
	enddo
	xnoise_fac_out(icel,ifreq,itrg)=sqrt(asum(2))
	enddo  !icel

	enddo !ifreq
	enddo !itrg

	write(4) factor_noise_out
	write(4) factor_area_out
	write(4) xnoise_fac_out
	enddo !kscan

	close(4)

	stop 'norm end'								   
	end


      subroutine read_data(itrg,iopt_cel_in,kscan_in, weight,factor_noise,factor_area)
	use mwi_sensor_module 	  	 		   		   
      implicit none
      
      integer(4), parameter :: kcel_max=58  
c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6
        
      integer(4) itrg,iopt_cel_in,kscan_in
      integer(4) iopt_cel,kscan

	real(8) weight(maxcel,-kcel_max:kcel_max,-14:14,mfreq)  !index order is a bit weird, but this is what peter has been using
 	real(8) factor_noise(maxcel,mfreq)
 	real(8) factor_area( maxcel,mfreq)
      character(120) filename
      
      if(itrg.eq.1 .or. itrg.eq.3) then
      kscan=kscan_in
      iopt_cel=iopt_cel_in
      
      else
      iopt_cel=1 !elliptical targets have no synthetic cells
      if(kscan_in.ge.3 .and. kscan_in.le.5) then !near native sccan
      kscan=4  !native scan
      else
      kscan=5  !halfway synthetic scan
      endif
      
      endif
      
	write(filename,9001) iopt_cel,kscan,itrg 
 9001 format('data\check_fit_',3i2.2 ,'.dat')
	call openbig(3,filename,'old')
	read(3) factor_noise
	read(3) factor_area
	close(3)
	
	write(filename,9002) iopt_cel,kscan,itrg
 9002 format('data\find_resampling_weights_',3i2.2 ,'.dat')
	call openbig(3,filename,'old')
	read(3) weight
	close(3)
     
      return
      end


	include 'x:\syspgm\openbig.for'
