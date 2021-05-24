      include 'o:\mwi\routines\mwi_sensor_module.f'   

	program make_processing_weights	 	  	 		   		   
	use mwi_sensor_module 	  	 		   		   
	implicit none

      integer(4), parameter :: kcel_max=58  

c     'o:\mwi\routines\mwi_sensor_module.f' has nfreq=5. new convention required nfreq=6
c     i did not want to change 'o:\mwi\routines\mwi_sensor_module.f', so i use mfreq
	integer(4), parameter:: mfreq=6
	
      character(120) filename,filename_ta_resample_wts

	integer(4) itrg,kscan,iopt_cel,icel,ifreq,index1,index2	   
      real(8) t

      real(8) weight(maxcel,-kcel_max:kcel_max,-14:14,mfreq)  
      
     	real(4), allocatable :: weight_38(:,:,:,:,:)
	real(4), allocatable :: weight_30(:,:,:,:,:)
	real(4), allocatable :: weight_25(:,:,:,:,:)
	real(4), allocatable :: weight_20(:,:,:,:,:)

      
	allocate (weight_38(18, maxcel, -kcel_max:kcel_max,-14:14,mfreq))
      allocate (weight_30( 2, maxcel, -kcel_max:kcel_max,-14:14,mfreq))  
	allocate (weight_25(18, maxcel, -kcel_max:kcel_max,-14:14,mfreq))
      allocate (weight_20( 2, maxcel, -kcel_max:kcel_max,-14:14,mfreq))   


      do itrg=1,4
      do kscan=0,8
      if(itrg.eq.2 .or. itrg.eq.4) then
      if(kscan.lt.4 .or. kscan.gt.5) cycle
      endif
      do iopt_cel=1,2
      if(itrg.eq.2 .or. itrg.eq.4) then
      if(iopt_cel.ne.1) cycle
      endif
      
	index1=kscan+1 + 9*(iopt_cel-1)
	index2=kscan-3

 
	write(filename,9003) iopt_cel,kscan,itrg
 9003 format('O:\mwi\resampling\data\find_resampling_weights_adjust_',3i2.2 ,'.dat')
      write(*,*) filename
	call openbig(3,filename,'old')
	read(3) weight
	close(3)
	
      weight(:,:,:,5)=weight(:,:,:,4)  !set second 37 GHz to first 37 GHz
	
c     check normalization
	do ifreq=1,mfreq
	do icel=1,maxcel
	t=sum(weight(icel,:,:,ifreq))
	if(iopt_cel.eq.2 .and. icel.eq.maxcel) then
	if(t.ne.0) stop 'error2'
	else
	if(abs(t-1).gt.1.e-10) stop 'error1'
	endif
	enddo
	enddo
	
	if(itrg.eq.1) then
	do ifreq=1,mfreq
	weight_38(index1,:,:,:,ifreq)=weight(:,:,:,ifreq) 
	enddo
	endif
	
	if(itrg.eq.2) then
	do ifreq=1,mfreq
	weight_30(index2,:,:,:,ifreq)=weight(:,:,:,ifreq) 
	enddo
	endif
	
	if(itrg.eq.3) then
	do ifreq=1,mfreq
	weight_25(index1,:,:,:,ifreq)=weight(:,:,:,ifreq)  
	enddo
	endif
	
	if(itrg.eq.4) then
	do ifreq=1,mfreq
	weight_20(index2,:,:,:,ifreq)=weight(:,:,:,ifreq) 
	enddo
	endif
	
      enddo  !iopt_cel
      enddo  !kscan
      enddo  !itrg

	filename_ta_resample_wts=  'O:\mwi\resampling\data\make_processing_weights_adjust.dat'
	call openbig(4,filename_ta_resample_wts,'new')
	write(4) weight_38
	write(4) weight_30
	write(4) weight_25
	write(4) weight_20
	close(4)


	stop 'norm end'								   
	end


	include 'x:\syspgm\openbig.for'
