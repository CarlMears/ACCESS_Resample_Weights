	program rsp_wt_draper 	  	 		   	      
	implicit none

      integer(4), parameter :: maxcel= 571 !extended scan
	integer(4), parameter :: kcel_max= 58
	integer(4), parameter :: kscan_max=14
	integer(4), parameter :: mfreq=6
	
	integer(4) ifreq
	real(8) arsp(maxcel,-kcel_max:kcel_max,-kscan_max:kscan_max,mfreq) !(A+B)/2
	real(8) drsp(maxcel,-kcel_max:kcel_max,-kscan_max:kscan_max,mfreq) !(B-A)/2 
      
	call openbig(3,'O:\mwi\resampling\data\find_resampling_weights_leakage_4stk.dat','old')
	read(3) arsp,drsp  
	close(3)
	
	call openbig(4,'O:\mwi\resampling\data\rsp_wt_draper_4stk.dat','new')

c     mfreq=6 is 11, 19, 24, 36, 37, 89  GHz
c     David only needs mfreq=1,2, and 4
c     Also David only need drsp
	do ifreq=1,4
	if(ifreq.eq.3) cycle !24 ghz
	write(4) drsp(:,:,:,ifreq)
	enddo
	close(4)
	
	call openbig(3,'O:\mwi\resampling\data\find_resampling_weights_leakage_3stk.dat','old')
	read(3) arsp,drsp  
	close(3)

	call openbig(4,'O:\mwi\resampling\data\rsp_wt_draper_3stk.dat','new')

c     mfreq=6 is 11, 19, 24, 36, 37, 89  GHz
c     David only needs mfreq=1,2, and 4
c     Also David only need drsp
	do ifreq=1,4
	if(ifreq.eq.3) cycle !24 ghz
	write(4) drsp(:,:,:,ifreq)
	enddo
	close(4)
	
      stop 'norm end'
      end
 


	include 'x:\syspgm\openbig.for'


