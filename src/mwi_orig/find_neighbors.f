c     May 27 2019 changed Nov 2 2020
c     maxcel changed from 553 to 571 to implement extended scan
c     new survey data files
	program find_neighbors	 	  	 		   		   
	implicit none											   

c     integer(4), parameter :: maxcel= 553 !regular scan
      integer(4), parameter :: maxcel= 571 !extended scan
      integer(4), parameter :: kcel_max=58  !derived by running this program trial and error for the 50 km neighbor distance
c50   real(4), parameter    :: diskm_max=50
      real(4), parameter    :: diskm_max=80
      
	integer(4) iscan,icel,jcel,kcel
	real(4) xlat(maxcel,-14:14),xlon(maxcel,-14:14)  
	real(4) diflat,diflon,coslat,dissq,diskm

	integer(1) index_wt(-kcel_max:kcel_max,-14:14,maxcel)
	integer(4) isum(maxcel)

	call openbig(3,'survey.dat','old')
	read(3) xlat,xlon
	close(3)
	
	open(6,file='find_neighbors.lis',status='new') 
	call openbig(4,'find_neighbors.dat','new')
	
	index_wt=0
	isum=0

	do jcel=1,maxcel  ! this is the target 
	 

	do iscan=-14,14
	do icel=1,maxcel
	diflat=xlat(icel,iscan)-xlat(jcel,0)  
	diflon=xlon(icel,iscan)-xlon(jcel,0)  !from survey i am sure there is no greenwich problem
	coslat=0.5*(cosd(xlat(icel,iscan))+cosd(xlat(jcel,0)))
	dissq=diflat**2 +(coslat*diflon)**2
	diskm=111*sqrt(dissq)

	kcel=icel-jcel
	
	if(diskm.gt.diskm_max) cycle

c50	if(iabs(kcel).gt.kcel_max) stop 'kcel oob, pgm stopped'
c50	if(iscan.eq.-14 .or. iscan.eq.14) stop 'iscan oob, pgm stopped' 
	if(iabs(kcel).gt.kcel_max) cycle

	index_wt(kcel,iscan,jcel)=1
	isum(               jcel)=isum(jcel)+1

	enddo	!icel
	enddo	!iscan

	enddo	!jcel  
	
	do jcel=1,maxcel  ! this is the target
	write(6,1001) jcel,isum(jcel)
 1001 format(9i6)
      enddo

	write(4) index_wt
	close(4)
	stop 'norm end'								   
	end


	include 'x:\syspgm\openbig.for'
