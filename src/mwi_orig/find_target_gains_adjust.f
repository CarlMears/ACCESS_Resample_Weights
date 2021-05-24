      include 'o:\mwi\routines\mwi_sensor_module.f'     
      include 'o:\mwi\routines\mwi_simulator_module.f' 

	program find_target_gains	  	 	  	 		   		   
	use l2_module														   
	implicit none
  
	character(120) filename,filename_l1a	 	 			    
	
	real(8) start_time

	real(8) secyr,secdy
	integer(4) lyear,idayjl,imon,idaymo 
	integer(4) iorbit
	
	integer(4) itrg,ifreq,iscan0
	real(4) xlat0,xlon0
	
      do itrg=1,4
      
      do ifreq=1,6
      if(ifreq.eq.5) cycle  !36 and 37 ghz currently the same. i dont bother doing ifreq=5 for now
      
	write(filename,9000) itrg,ifreq
 9000 format('O:\mwi\resampling\mwi_target_arrays_adjust\find_target_gains_',2i1,'.lis')
      open(6,file=filename, status='new')

	open(3,file='survey.txt')
	read(3,3001) iorbit,iscan0,xlat0,xlon0
 3001 format(2i6,2f8.3)
	close(3)	
	xlat0=0.5*nint(2*xlat0)
	xlon0=0.5*nint(2*xlon0)
	
	ksat=37  !mwi
	call define_filenames
      call readin_data_files

      iorbit=2	
	write(filename_l1a,9001) iorbit
 9001 format('O:\mwi\simulator\sc_orbits\600pm\r',i5.5)
	call openbig(2,filename_l1a,'old')

	call read_l1a(2,iorbit,filename_l1a, start_time)  !2 is ilu, 

      call fd_date_2000(start_time, secyr,lyear,idayjl,imon,idaymo,secdy)
    
    	write(*,'(i3,i6,i6,i5,2i3,f7.2)') ksat,iorbit,numscan,lyear,imon,idaymo,secdy/3600.
    	write(6,'(i3,i6,i6,i5,2i3,f7.2)') ksat,iorbit,numscan,lyear,imon,idaymo,secdy/3600.
	
	call mwi_geolocation(iorbit)
	
	call mwi_gains(itrg,ifreq,iscan0,xlat0,xlon0)
	close(6)
	
	enddo  !ifreq
	enddo  !itrg

	stop 'norm end'
	end   
	

 		
	

	subroutine mwi_gains(itrg,ifreq,iscan0,xlat0,xlon0)
	use l2_module														   
	implicit none
	
      integer(4), parameter :: maxcel_out=2*maxcel-1   
	
	integer(4), parameter:: nlat=1601
c	integer(4), parameter:: nlon=2101  !this was used for tmi, decided to increase to 2401, which is what ssmis used
	integer(4), parameter:: nlon=2401  !ssmis value, probably overkill for gmi conisered it low altitude

	real(8), parameter::  area=1110.d0**2	!area of cell is 1110 meters
	
	character(80) filename

	integer(4) itrg,ifreq,iscan0
	real(4) xlat0,xlon0
	real(4) xcel
	
	real(8) xlook(3),costht,domega,gain
	real(8) delta,cosdelta
	real(8) cell0(3,nlon,nlat),cell(3,nlon,nlat)
	real(8) qsum,xmin,xmax
	real(8) gain_domega(nlon,nlat)
	
	real(8) beamwidth(4),sumgain(4)
 	real(8) beamwidth_adjust(maxcel_out,0:8,6,4)

	integer(4) istart
	integer(4) iscan,icel
      integer(4) kscan,kscan1,kscan2,kkscan
	integer(4) ilat,ilon
	integer(2) ilat2,ilon2
	
	real(8) x_eci(3),y_eci(3),z_eci(3),alpha
	real(8) scpos_eci_cel(3),rotangle,days
	real(8) coslat,sinlat,cell_pos_eci0(3),cell_pos_eci(3) 
	
	real(8) sunvec(3)
	real(8) observation_time, time_since_scan_start
	real(8) x0(3),y0(3),z0(3)
	real(8) np(3,3)
      real(8) b(3),refl(3),xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt
	real(8) geolat,geolon
c     real(8) ha(3),va(3),xant,yant,elv,azm  !not need for these target gains

 	real(4) omega_avg
 	
	integer(4) ierror
	
	data istart/1/
	
c     ========================================================================================================================
c     ============================================== begin starting block ====================================================
c     ========================================================================================================================

      if(istart.eq.1) then
      istart=0
	call openbig(3,'O:\mwi\antenna\data\target_beamwidth_sumgain.dat','old')
	read(3) beamwidth
	read(3) sumgain
	close(3)
	
	call openbig(3,'O:\mwi\antenna\adjust_beamwidths.dat','old')
	read(3) beamwidth_adjust
	close(3)
	
	endif
	
	
	call fd_cell_vector(nlat,nlon,xlat0,xlon0, cell0,cell)

	iscan=iscan0

	xmin=1.e30; xmax=-1.e30
      
      if(iflag_scn(iscan).ne.0)  stop 'error1'
      
      omega_avg=6.*scan_rpm(iscan)  !deg/sec
	observation_time = scan_time(iscan)  !scan time is ut1 time
	days=observation_time/86400.d0 - 0.5d0
	call fd_precession(days, np)

	scpos_eci(:,iscan)=matmul(np,scpos_eci2000(:,iscan))
	scvel_eci(:,iscan)=matmul(np,scvel_eci2000(:,iscan))   

	z0= -scpos_eci(:,iscan)
	x0=  scvel_eci(:,iscan)
	call cross_norm(z0,x0, y0) 
	      
      if(itrg.eq.1 .or. itrg.eq.3) then
      kscan1=-4
      kscan2= 4
      else
      kscan1=0
      kscan2=1
      endif
      
	do kscan=kscan1,kscan2
	
      if(itrg.eq.1 .or. itrg.eq.3) then
      kkscan=kscan+4
      else
      if(kscan.eq.0)  kkscan=4
      if(kscan.eq.1)  kkscan=8
      endif
	
	do icel=1,maxcel_out
	xcel=1 + 0.5*(icel-1)     

      if(itrg.eq.1 .or. itrg.eq.3) then
      time_since_scan_start = sample_time*(xcel-1) + scan_period*kscan/9
      else
      time_since_scan_start = sample_time*(xcel-1) + scan_period*kscan/2.
      endif

	scpos_eci_cel = scpos_eci(:,iscan) + time_since_scan_start*scvel_eci(:,iscan)
	call get_sc_axes_from_rpy(scpos_eci_cel,y0,scrpy(1,iscan),scrpy(2,iscan),scrpy(3,iscan), x_eci,y_eci,z_eci)
   
c     drapers pattern already include pattern smearing from integration 
c     so i do not need the 	do interval=0,9 loop

	observation_time = scan_time(iscan) + time_since_scan_start 
	call get_gm_angle(observation_time,days, rotangle)

      alpha= alpha_start + omega_avg*sample_time*(xcel-1)

 	call fd_cell_parameters(3,re,rp,ffac,geosync_alt,
     &                        rotangle,scpos_eci_cel,x_eci,y_eci,z_eci,beta,alpha,sunvec,
     &                        b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)

	if(ierror.ne.0) stop 'ierror.ne.0, pgm stopped'

c     next two statements not needed for these gain targets
c     call cross_norm(b,z_eci, ha) 	 !ha=kx(-z)=(-b)x(-z)=bxz
c	call cross_norm(b,ha,    va)	 !va=haxk=bxha

c     ==============================================================================================================================	
c      compute the vector from the earth's center that passes through the current cell locations and goes to alttitude of satellite
c     ==============================================================================================================================	
      xlat=atand(ffac*tand(xlat))  !convert back to geocentric
	coslat=cosd(xlat)
	sinlat=sind(xlat)
      cell_pos_eci0(1)=coslat*cosd(xlon_rot+rotangle)
      cell_pos_eci0(2)=coslat*sind(xlon_rot+rotangle)
      cell_pos_eci0(3)=sinlat
	cell_pos_eci= sqrt(dot_product(scpos_eci_cel,scpos_eci_cel))*cell_pos_eci0  !earth center to sc altitude thru current cell

	gain_domega=0; qsum=0 
	 
	do ilat=1,nlat
	do ilon=1,nlon

	if(itrg.eq.1 .or. itrg.eq.3) then
	xlook= cell(:,ilon,ilat) - cell_pos_eci	   !point from s/c down to cell
	else
	xlook= cell(:,ilon,ilat) - scpos_eci_cel	   !point from s/c down to cell
	endif

      range=sqrt(dot_product(xlook,xlook))
	xlook=xlook/range
	
      costht=-dot_product(xlook,cell0(:,ilon,ilat)) 
	if(costht.le.0) stop 'pgm stopped,costht.le.0' 
	
	if(itrg.eq.1 .or. itrg.eq.3) then
	cosdelta=-dot_product(xlook,cell_pos_eci0)  !minus sign becuase cell_pos_eci0 is going up
	else
	cosdelta=dot_product(xlook,b)
	endif
	
	if(cosdelta.gt.1) cosdelta=1  !should never be near -1
	delta=acosd(cosdelta)
 
      if(delta.le.3*beamwidth_adjust(icel,kkscan,ifreq,itrg)) then  !see mathnotes
     	call fd_gain_targets(beamwidth_adjust(icel,kkscan,ifreq,itrg),delta, gain)
      domega=costht*area*sqrt(1-cell0(3,ilon,ilat)**2)/(range*range)  !coslat=sqrt(1-cell0(3,ilon,ilat)**2)
	gain=gain*domega
	gain_domega(ilon,ilat)=gain
	qsum=qsum + gain
	endif ! if(delta.le.3*beamwidth_adjust(icel,kkscan,ifreq,itrg)) 

	enddo !ilat
	enddo !ilon

	qsum=qsum/sumgain(itrg)
	if(qsum.lt.xmin) xmin=qsum
	if(qsum.gt.xmax) xmax=qsum

      write(filename,9001) itrg,ifreq,kscan+4,icel
 9001 format('O:\mwi\resampling\mwi_target_arrays_adjust\target',i1,'\freq',i1,'\s',i2.2,'c',i4.4,'.dat')

      write(*,1003) filename
 1003 format(1x,a80)

	call openbig(4,filename,'new')

	do ilat=1,nlat
	do ilon=1,nlon
	if(gain_domega(ilon,ilat).ne.0) then
	ilat2=ilat
	ilon2=ilon
	write(4) ilat2,ilon2,gain_domega(ilon,ilat)
	endif
	enddo  !ilon
	enddo  !ilat

	close(4)
	
	enddo  !icel
	enddo  !kscan

      write(*,1004) itrg,xmin,xmax
      write(6,1004) itrg,xmin,xmax
 1004 format(' completed ',i5,2f10.5)
      return
	end subroutine mwi_gains

 

	
	subroutine fd_cell_vector(nlat,nlon,xlat0,xlon0, cell0,cell)
	implicit none
 
c	earth dimensions match wgs84	
	real(8), parameter :: re=6378.137d3		!meters
      real(8), parameter :: rp=6356.752d3		!meters
 
	integer(4) nlat,nlon
	real(4) xlat0,xlon0
	real(8) cell0(3,nlon,nlat),cell(3,nlon,nlat)

	integer(4) ilat,ilon
	real(8) xlat,xlon,rearth
	real(8) coslat,sinlat,coslon,sinlon

	do ilat=1,nlat
	do ilon=1,nlon
	xlat=0.01d0*(ilat-1) + xlat0 - 8.0
	xlon=0.01d0*(ilon-1) + xlon0 -12.0
 
	coslat=cosd(xlat)
	sinlat=sind(xlat)
	coslon=cosd(xlon)
	sinlon=sind(xlon)
 
      cell0(1,ilon,ilat)=coslat*coslon
      cell0(2,ilon,ilat)=coslat*sinlon
      cell0(3,ilon,ilat)=sinlat
 
      rearth=re*rp/sqrt((rp*coslat)*(rp*coslat)+(re*sinlat)*(re*sinlat)) 
 
	cell(:,ilon,ilat)= rearth*cell0(:,ilon,ilat)
	
	enddo
	enddo
 
	return
	end


      include 'o:\mwi\antenna\fd_gain_targets.f'
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
	include 'x:\syspgm\openbig.for'


