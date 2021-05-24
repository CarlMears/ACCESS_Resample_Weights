c     may 27 2019 changed nov 2 2020.  
c     1.  'o:\mwi\routines\mwi_sensor_module.f' was changed to implement extended scan
c     2.  reference iscan0 changed from 643 to 645
c     3.  altitude of orbit at equator changed from 823 to 833 km
c     4.  new 'o:\mwi\routines\read_l1a.f' with counter-clockwise scan and a few other little changes

      include 'o:\mwi\routines\mwi_sensor_module.f'   
      include 'o:\mwi\routines\mwi_simulator_module.f' 

	program survey     	 	 			 			
	use l2_module														   
	implicit none

	character(120) filename_l1a					    
	
	real(8) start_time

	real(8) secyr,secdy
	integer(4) lyear,idayjl,imon,idaymo
	integer(4) iorbit,iscan0,iwrite
	
	real(4) xlat0,xlon0,xlatsv(maxcel,29),xlonsv(maxcel,29)
	
	write(*,*) numscan0
	
	ksat=37  !mwi

	open(7,   file='survey.txt',status='new')
      open(6,   file='survey.lis',status='new')
	call openbig(4,'survey.dat',       'new')

	call define_filenames
      call readin_data_files
 
      iorbit=2	
	write(filename_l1a,9001) iorbit
 9001 format('O:\mwi\simulator\sc_orbits\600pm\r',i5.5)
	call openbig(2,filename_l1a,'old')

	call read_l1a(2,iorbit,filename_l1a, start_time)  !2 is ilu, 

      call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)

    	write(*,'(5i6,2f10.2)') iorbit,numscan,lyear,imon,idaymo,secdy/3600.
    	
    	if(abs(numscan-numscan0).gt.5) stop 'not a complete orbit, pgm stopped'

	call mwi_geolocation(iorbit)

	do iscan0=645-20,645+20
	call fd_survey(iscan0, xlat0,xlon0,xlatsv,xlonsv)
	
	iwrite=0
	if(iscan0.eq.645) iwrite=1 !645 determined from previous run of this program

	if(iwrite.eq.1) then
	if(abs(xlat0).gt.0.075) stop 'abs(xlat0).gt.0.075, pgm stopped'
      write(4) xlatsv,xlonsv
	write(7,7001) iorbit,iscan0,xlat0,xlon0
 7001 format(2i6,2f8.3)
      endif
      
	enddo  !iscan0


	stop 'norm end'
	end   
	
	
	
	subroutine fd_survey(iscan0, xlat0,xlon0,xlatsv,xlonsv)
	use l2_module														   
	implicit none
	
 	real(8) time_since_scan_start,observation_time,days,rotangle
	real(8) xmin1,xmax1,xmin2,xmax2
	real(8) ymin1,ymax1,ymin2,ymax2
      
      integer(4) iscan0
	real(4) xlat0,xlon0,xlatsv(maxcel,29),xlonsv(maxcel,29)
 
	integer(4) iscan1,iscan2,iscan_ref
	integer(4) iscan,icel
	real(4) xlat,xlon
	

      
	xmin1=1.e30; xmin2=1.e30; xmax1=-1.e30; xmax2=-1.e30
	ymin1=1.e30; ymin2=1.e30; ymax1=-1.e30; ymax2=-1.e30
 
c     11-37 ghz

	iscan1=iscan0 - 14
	iscan2=iscan0 + 14
	
	do iscan=iscan1,iscan2
	iscan_ref=iscan-iscan1+1
	observation_time = scan_time(iscan)  !scan time is ut1 time
	days=observation_time/86400.d0 - 0.5d0
 
	do icel=1,maxcel
      time_since_scan_start = sample_time*(icel-1)
	observation_time = scan_time(iscan) + time_since_scan_start 
	call get_gm_angle(observation_time,days, rotangle)
	
 	xlat=cellat(icel,iscan)
	xlon=cellon(icel,iscan)	 + rotangle  !convert from earth to inertial
	if(xlon.lt.  0) xlon=xlon+360	  
	if(xlon.ge.360) xlon=xlon-360    
	if(xlat.lt.xmin1) xmin1=xlat
	if(xlat.gt.xmax1) xmax1=xlat
	if(xlon.lt.xmin2) xmin2=xlon
	if(xlon.gt.xmax2) xmax2=xlon
	xlatsv(icel,iscan_ref)=xlat
	xlonsv(icel,iscan_ref)=xlon
	enddo  !icel
	enddo  !iscan
	
 
	xlat0=0.5*(xmin1+xmax1)
	xlon0=0.5*(xmin2+xmax2)
	
      write(6,6002) iscan0,xlat0,xlon0,xmin1,xmax1,xmin2,xmax2,xmax1-xmin1,xmax2-xmin2
 6002 format(i5,8f8.3)
 
	return
	end 
 
	  



	include 'o:\mwi\routines\define_filenames.f'    
	include 'o:\mwi\routines\readin_data_files.f' 
	include 'o:\mwi\routines\read_l1a.f' 
	
	include 'o:\mwi\routines\mwi_geolocation.f' 
	include 'o:\mwi\routines\fd_faraday_angle.f' 
      include 'o:\earth_magnetic_field\fd_earth_tec_field_allyears.f'
      include 'o:\earth_magnetic_field\fd_earth_tec_field_ubern_1hourly.f'
	include 'o:\amsr_l2\v05\routines\pre_nut_routines.f'
 	include 'o:\sun_moon\sun_vector.f'            
 	include 'o:\sun_moon\moon_vector.f'            

      include 'o:\earth_magnetic_field\find_earth_magnetic_field.f'
	

      include 'x:\mathlib\math_routines.f'

	include 'x:\syspgm\fd_date_2000.f'
	include 'x:\syspgm\ck_file_time.f'
	include 'x:\syspgm\openbig.for'
	include 'x:\syspgm\openbig_try.f'


	 