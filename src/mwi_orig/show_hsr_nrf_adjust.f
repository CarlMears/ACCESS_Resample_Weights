      include 'o:\mwi\routines\mwi_sensor_module.f'   

	program show_hsr_nrf 	  	 		   		   
	use mwi_sensor_module 	  	 		   		   
	implicit none
	
      integer(4), parameter :: kcel_max=58    
      integer(4), parameter :: maxcel_out=2*maxcel-1   
      integer(4), parameter :: ntrg=4   

c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6

	integer(4) itrg,ifreq,icel,jcel  
	integer(4) kscan,iopt_cel
	
 	real(8) factor_noise(maxcel,mfreq)
 	real(8) factor_area( maxcel,mfreq)
	
 	real(8) factor_area_out(maxcel_out,mfreq,ntrg)
 	real(8) factor_area_out_adjust(maxcel_out,mfreq,ntrg)
 	
	real(8) bsum(0:1),xmin,xmax

 	real(8) size_km(ntrg)
 	
 	real(8) d_beamwidth_out(maxcel_out,0:8,mfreq,ntrg)
 	
c 	data size_km/38.d0, 29.5d0, 24.5d0, 19.5d0/
 	data size_km/39.d0, 30.d0, 25.d0, 20.d0/
 	
	call openbig(3,'O:\mwi\antenna\show_adjust_beamwidths.dat','old')
	read(3) d_beamwidth_out
	close(3)


	
 	do kscan=0,8
  	
	do itrg=3,3 !1,ntrg
  
c     ===================================================================
c     ===================== native cells ================================
c     ===================================================================
	iopt_cel=1
      call read_data(1,itrg,iopt_cel,kscan, factor_noise,factor_area)
	jcel=0
	do icel=1,maxcel_out,2
	jcel=jcel+1
	factor_area_out( icel,:,itrg)=factor_area( jcel,:)
	enddo
c     ===================================================================
c     ===================== synthetic cells ================================
c     ===================================================================
	iopt_cel=2
      call read_data(1,itrg,iopt_cel,kscan, factor_noise,factor_area)
	jcel=0
	do icel=2,maxcel_out-1,2
	jcel=jcel+1
	factor_area_out( icel,:,itrg)=factor_area( jcel,:)
	enddo
	
  
c     ===================================================================
c     ===================== native cells ================================
c     ===================================================================
	iopt_cel=1
      call read_data(2,itrg,iopt_cel,kscan, factor_noise,factor_area)
	jcel=0
	do icel=1,maxcel_out,2
	jcel=jcel+1
	factor_area_out_adjust( icel,:,itrg)=factor_area( jcel,:)
	enddo
c     ===================================================================
c     ===================== synthetic cells ================================
c     ===================================================================
	iopt_cel=2
      call read_data(2,itrg,iopt_cel,kscan, factor_noise,factor_area)
	jcel=0
	do icel=2,maxcel_out-1,2
	jcel=jcel+1
	factor_area_out_adjust( icel,:,itrg)=factor_area( jcel,:)
	enddo	
	
	
	
	do ifreq=1,mfreq

	if(ifreq.eq.5) cycle  !second 37 ghz assumed to be same as first and there are not gain patterns for second
	if(ifreq.eq.1 .and. itrg.ge.3) cycle !footprint larger than target.  there is no requirement for this.

	bsum=0; xmin=1.e30; xmax=-1.e30
	do icel=1,maxcel_out

      if(iabs(icel-571).lt.270) then ! midscan
      else
      cycle 
      endif
c      if(icel.eq.maxcel_out) cycle



 	if(abs(d_beamwidth_out(icel,kscan,ifreq,itrg)).lt.0.02) cycle
 	if(abs(factor_area_out_adjust(icel,ifreq,itrg)-factor_area_out(icel,ifreq,itrg)).lt.0.20) cycle
 	
      bsum(0)=bsum(0) + 1	
      bsum(1)=bsum(1) - 
     & d_beamwidth_out(icel,kscan,ifreq,itrg)/(factor_area_out_adjust(icel,ifreq,itrg)-factor_area_out(icel,ifreq,itrg))
      cycle


      bsum(0)=bsum(0) + 1	
      bsum(1)=bsum(1) + factor_area_out_adjust(icel,ifreq,itrg)	



      if(factor_area_out_adjust(icel,ifreq,itrg).lt.xmin) xmin=factor_area_out_adjust(icel,ifreq,itrg)
      if(factor_area_out_adjust(icel,ifreq,itrg).gt.xmax) xmax=factor_area_out_adjust(icel,ifreq,itrg)
      
      if(factor_area_out_adjust(icel,ifreq,itrg).gt.size_km(itrg)) then
      write(*,1003) kscan,itrg,ifreq,icel,factor_area_out_adjust(icel,ifreq,itrg),size_km(itrg)
 1003 format(4i4,13f8.1) 
      endif
      
      enddo  !icel
      bsum(1)=bsum(1)/bsum(0)
      write(*,1005) kscan,itrg,ifreq,bsum(1)
 1005 format(3i6,f10.5)
      cycle
      
c     toggle write statements to simply output
      if(xmax.gt.size_km(itrg)) then
      write(*,1002) kscan,itrg,ifreq,bsum(1),xmin-size_km(itrg),xmax-size_km(itrg) !,xmax-xmin
 1002 format(3i4,13f8.1) 
      endif 	

	enddo !ifreq
	enddo !itrg

	enddo !kscan


	stop 'norm end'								   
	end


      subroutine read_data(icase,itrg,iopt_cel_in,kscan_in, factor_noise,factor_area)
	use mwi_sensor_module 	  	 		   		   
      implicit none
      
      integer(4), parameter :: kcel_max=58  
c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6
        
      integer(4) icase,itrg,iopt_cel_in,kscan_in
      integer(4) iopt_cel,kscan

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
      
      if(icase.eq.1) write(filename,9001) iopt_cel,kscan,itrg 
      if(icase.eq.2) write(filename,9002) iopt_cel,kscan,itrg 
 9001 format('data\check_fit_',3i2.2 ,'.dat')
 9002 format('data\check_fit_adjust_',3i2.2 ,'.dat')
	call openbig(3,filename,'old')
	read(3) factor_noise
	read(3) factor_area
	close(3)
     
      return
      end


	include 'x:\syspgm\openbig.for'
