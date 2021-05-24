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

	do iscan0=643-20,643+20
	call fd_survey(iscan0, xlat0,xlon0,xlatsv,xlonsv)
	
	iwrite=0
	if(iscan0.eq.643) iwrite=1 !924 determined from previous run of this program

    
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
	
 
	return
	end 
 
	  



	include 'o:\mwi\routines\define_filenames.f'    
	include 'o:\mwi\routines\readin_data_files.f' 
	include 'o:\mwi\routines\read_l1a.f' 
	
c	include 'o:\mwi\routines\mwi_geolocation.f' 
c     02/16/2019 changed 5/27/2019.  89 ghz now has same eia as other channels,
c     hence only one set of geoloc are required. 02/16/2019 version saved with date suffix


c     july 3 2018 changed feb 16 2019.
c     an additional call to get_sc_axes_from_rpy is put into the icel loop. 
c     this was an oversight.  x,y,z eci change during the course of a scan.
c     i check windsat, gmi, and amsr.  they all do it this correct way.  not sure why mwi was doing it wrong.
c     i also put a call to get_sc_axes_from_rpy when finding cold mirrror spillover.


c     jun 1 2018 changed july 3 2018.  variabiles frd_rad(icel,iscan,i89) and tec(icel,iscan,i89) are now computed

c     jan 3 2018 changed jun 1 2018.
c     1. ifreq is no longer an input
c     2. numcel is now the same for all freqs (it is maxcel)
c     3. sample time and alpha_start are now the same for all freqs
c     4. there are now just two values for beta (11-37 ghz, 89 ghz) rather than separate values for each freq
c     5. two sets of geolocation variable are found:  one for 11-37 ghz and the other for 89 ghz

	subroutine mwi_geolocation(iorbit)
	use l2_module
	implicit none
 
	real(8), parameter :: alphax=0  !dummy value for alpha
	real(8), parameter :: betax =0  !dummy value for beta

	integer(4) iorbit
	
	real(8) dot_product_unit8
	real(8) x_eci(3),y_eci(3),z_eci(3),alpha,r3,v3
	real(8) rotangle,days
	real(8) scpos_eci_cel(3),scpos_eci_cold(3)
	real(8) sunvec_eci2000(3),sunvec(3),sundis_km,moonvec_eci2000(3),moonvec(3),moondis_km
	real(8) observation_time, time_since_scan_start
	real(8) x0(3),y0(3),z0(3)
	real(8) np(3,3)
      real(8) b(3),refl(3),xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt
	real(8) sc_to_sun(3),dot1,dot2,cossunang
	real(8) geolat,geolon
 	real(4) omega_avg
      real(8) orbitsv
      real(8) ha(3),va(3)


	integer(4) iscan,icel
	integer(4) ierror
	
      integer(4) lyear,idayjl,imon,idaymo
      real(8) secyr,secdy
      	
      real(4) bmag(3) 
      real(8) cosr,sinr,bmag_eci(3)

      integer(4) isecdy !inserted 7/3/2018
 
	orbitsv=-1
      
	do iscan=1,numscan

      if(iflag_scn(iscan).ne.0)  cycle
      	
      omega_avg=6.*scan_rpm(iscan)  !deg/sec
	observation_time = scan_time(iscan)  !scan time is ut1 time
	days=observation_time/86400.d0 - 0.5d0
	call get_gm_angle(observation_time,days, rotangle)
	call fd_precession(days, np)

	scpos_eci(:,iscan)=matmul(np,scpos_eci2000(:,iscan))
	scvel_eci(:,iscan)=matmul(np,scvel_eci2000(:,iscan))   

	z0= -scpos_eci(:,iscan)
	x0=  scvel_eci(:,iscan)
	call cross_norm(z0,x0, y0) 
	call get_sc_axes_from_rpy(scpos_eci(:,iscan),y0,scrpy(1,iscan),scrpy(2,iscan),scrpy(3,iscan), x_eci,y_eci,z_eci)
	
	
c	call pitch_sc(pitch_corr, x_eci,y_eci,z_eci)
c	call  roll_sc(  -0.035d0, x_eci,y_eci,z_eci)

	
	

	r3=scpos_eci(3,iscan)/sqrt(dot_product(scpos_eci(:,iscan),scpos_eci(:,iscan)))
	v3=scvel_eci(3,iscan)/sqrt(dot_product(scvel_eci(:,iscan),scvel_eci(:,iscan)))
	zang(iscan)=datan2d(r3,v3) + 90.
	call fixang8(zang(iscan))
	
	orbit(iscan)=iorbit + zang(iscan)/360.d0
	if(iscan.lt.numscan0/2 .and. zang(iscan).gt.330.) orbit(iscan)=orbit(iscan)-1
	if(iscan.ge.numscan0/2 .and. zang(iscan).lt. 30.) orbit(iscan)=orbit(iscan)+1
	if(iscan.gt.1) then !probably need better check here
	if(orbit(iscan)-orbitsv.le.0) then
	write(*,1001) iscan,zang(iscan-1), zang(iscan)
 1001 format('!!orb err ', i6,10f8.3)
      stop 'pgm stopped'
      endif
      endif
	orbitsv=orbit(iscan)
	
c     ===============================================================================================================
c     =========================== compute sun and moon parameters ===================================================	
c     ===============================================================================================================

      call  sun_vector(days,  sunvec_eci2000, sundis_km)	 !returns vector in j2000 system, dis is km
	sunvec=matmul(np,sunvec_eci2000)
	sunvec_km(:,iscan)=sundis_km*sunvec

	sunlat(iscan)=asind(sunvec(3))
	sunlon(iscan)=atan2d(sunvec(2),sunvec(1)) - rotangle
	call fixang4( sunlon(iscan))
	sundis(iscan)=sundis_km/149.619d6
 
	sc_to_sun=sundis_km*sunvec - 1.d-3*scpos_eci(:,iscan)
	sc_to_sun=sc_to_sun/sqrt(dot_product(sc_to_sun,sc_to_sun))
 
	cossunang=-dot_product_unit8(z_eci,sc_to_sun)	!minus sign accounts for dif between windsat x,y,z and amsr x,y,z
      beta_sun(iscan)=acosd(cossunang)
 
	dot1= dot_product(x_eci,sc_to_sun)
	dot2=-dot_product(y_eci,sc_to_sun)      !minus sign accounts for dif between windsat x,y,z and amsr x,y,z
	alpha_sun(iscan)=atan2d(dot2,dot1)
	
      call moon_vector(days, moonvec_eci2000,moondis_km)	 !returns vector in j2000 system, dis is km
	moonvec=matmul(np,moonvec_eci2000)
	moonvec_km(:,iscan)=moondis_km*moonvec
	
c     ===============================================================================================================
c     ============================== compute spacecraft parameters ==================================================	
c     ===============================================================================================================
	
	
c     following call has boresight vector = sc_to_sun
c     if there is no earth intersection, ierror is set to 1 and this indicates sun if visible from spacecraft
 	call fd_cell_parameters(4,re,rp,ffac,geosync_alt,
     &                        rotangle,scpos_eci(:,iscan),x_eci,y_eci,z_eci,betax,alphax,sunvec,sc_to_sun,
     &                        refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)

      kflag_sun(iscan)=ierror
 


	call fd_cell_parameters(1,re,rp,ffac,geosync_alt,
     &                        rotangle,scpos_eci(:,iscan),x_eci,y_eci,z_eci,betax,alphax,sunvec,
     &                        b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)

	if(ierror.ne.0) iflag_scn(iscan)=ibset(iflag_scn(iscan),0)
	
	scloc(1,iscan)=xlat      !spacecraft geodetic latitude
	scloc(2,iscan)=xlon_rot  !spacecraft longitude east from greenwich
	scloc(3,iscan)=range     !spacecraft altitude in meters

	
cw      write(*,'(2i6,10f15.5)') iscan,ierror,zang(iscan),scloc(:,iscan),scrpy(:,iscan)

	
c     =========================================================================================================
c     =================================== compute earth magnetic field vector =================================
c     =========================================================================================================

      call fd_date_2000(scan_time(iscan), secyr,lyear,idayjl,imon,idaymo,secdy)
      call find_earth_magnetic_field(lyear,idayjl,scloc(:,iscan), bmag)
      bmag=1.d-4*bmag  !1.d-4 scales bmag to be near unity.  just for convenience

c	convert from earth to inertial
	cosr=cosd(rotangle)
	sinr=sind(rotangle)
	bmag_eci(1)=cosr*bmag(1) - sinr*bmag(2)
	bmag_eci(2)=cosr*bmag(2) + sinr*bmag(1)
	bmag_eci(3)=bmag(3)
	
c     bmag_sc is defined in a cooridinate systems with the z axis going away from earth
c     this is different than standard with  z axis pointing towards earth
	bmag_sc(1,iscan)= dot_product(bmag_eci,x_eci)
	bmag_sc(2,iscan)=-dot_product(bmag_eci,y_eci)
	bmag_sc(3,iscan)=-dot_product(bmag_eci,z_eci)





c     =========================================================================================================
c     =================================== compute locations of cells ========================================
c     =========================================================================================================

	do icel=1,maxcel

      time_since_scan_start = sample_time*(icel-1)
	observation_time = scan_time(iscan) + time_since_scan_start 
	call get_gm_angle(observation_time,days, rotangle)

	scpos_eci_cel = scpos_eci(:,iscan) + time_since_scan_start*scvel_eci(:,iscan)

 	call get_sc_axes_from_rpy(scpos_eci_cel,y0,scrpy(1,iscan),scrpy(2,iscan),scrpy(3,iscan), x_eci,y_eci,z_eci)
 	
 	


      alpha= alpha_start + omega_avg*time_since_scan_start
      
      if(icel.eq.277) then
      alpha=0  !for cell 277 alpha is actually-8.687329101576324E-003, here i set it to zero
      endif
      
 	call fd_cell_parameters(3,re,rp,ffac,geosync_alt,
     &                        rotangle,scpos_eci_cel,x_eci,y_eci,z_eci,beta,alpha,sunvec,
     &                        b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)
     
c     check signs to be consistent with draper
      call cross_norm(b,z_eci, ha) 	 !ha=kx(-z)=(-b)x(-z)=bxz
	call cross_norm(b,ha,    va)	 !va=haxk=bxha
c     when alpha=0, ha=-y_eci 
      if(icel.eq.277) then
      write(*,1011) zang(iscan),x_eci,y_eci,z_eci,scrpy(:,iscan),ha+y_eci,va
 1011 format(f9.5,5x,6(3f9.5,2x))
      endif

 
 
	if(ierror.ne.0) iflag_scn(iscan)=ibset(iflag_scn(iscan),2)

	cellat(icel,iscan)=xlat
	cellon(icel,iscan)=xlon_rot
      call fixang4( cellon(icel,iscan))  !needed because of real(8) to real(4) conversion
	celeia(icel,iscan)=thtinc  
	celazm(icel,iscan)=azim
	celsun(icel,iscan)=sunglt
	reflat(icel,iscan)=geolat
	reflon(icel,iscan)=geolon
      call fixang4( reflon(icel,iscan)) 
      
c     next 3 lines inserted 7/3/2018
      isecdy=int(secdy)
	call fd_faraday_angle(lyear,idayjl,isecdy,rotangle,scpos_eci(:,iscan),scloc(1,iscan),b, 
     &                      frd_rad(icel,iscan),tec(icel,iscan))
      
      enddo !icel
       
       

c     ===========================================================================================
c     ================================ compute location of cold spillover =======================
c     ===========================================================================================
      
	observation_time = scan_time(iscan) + time_cold_mirror
	call get_gm_angle(observation_time,days, rotangle)
	
	scpos_eci_cold = scpos_eci(:,iscan) + time_cold_mirror*scvel_eci(:,iscan)

 	call get_sc_axes_from_rpy(scpos_eci_cold,y0,scrpy(1,iscan),scrpy(2,iscan),scrpy(3,iscan), x_eci,y_eci,z_eci)

 	call fd_cell_parameters(2,re,rp,ffac,geosync_alt,
     &                        rotangle,scpos_eci(:,iscan),x_eci,y_eci,z_eci,beta,alpha_cold_mirror,sunvec,
     &                        b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)
     
	if(ierror.ne.0) iflag_scn(iscan)=ibset(iflag_scn(iscan),3)
 
 	cellat_cold(iscan)=xlat
	cellon_cold(iscan)=xlon_rot
      call fixang4( cellon_cold(iscan))  !needed because of real(8) to real(4) conversion
      
      enddo !iscan
      
 
	return
	end subroutine mwi_geolocation
 
 

 
 
 
 
 
 
	
 
c     x, y, z  s/c vectors
c     x vector points in general direction of s/c velocity
c     z vector points is s/c nadir and points towards earth for normal attitude
c     y vector is given by z cross x
 

c     icase=1 means compute satellite subpoint location
c     icase=2 means compute just cell locations given alpha and beta
c     icase=3 means compute all  cell parameters given alpha and beta
c     icase=4 means compute cell locations given boresight


	subroutine fd_cell_parameters(icase,re,rp,ffac,geosync_alt,
     &                              rotangle,scpos,x,y,z,beta,alpha,sunvec,
     &                              b,refl,xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon,ierror)
	implicit none

	real(8),    intent(in)    :: re,rp,ffac,geosync_alt,rotangle,scpos(3),x(3),y(3),z(3),beta,alpha,sunvec(3)
	integer(4), intent(in)    :: icase
	real(8),    intent(out)   :: refl(3),xlat,xlon_rot,range,rearth,thtinc,azim,pra,suninc,sunazm,sunglt,geolat,geolon
	integer(4), intent(out)   :: ierror
	real(8),    intent(inout) :: b(3)
 
	real(8) dot_product_unit8
 
	real(8) r,ru(3),sinlat,coslat,costht,xlon
	real(8) sinbeta,cosbeta
	real(8) bx,by,bz,cel1,cel2,cel3
	real(8) coslon,sinlon
	real(8) geoid(3),cosazim,cossunglt
	real(8) deltau
      real(8) ha(3),he(3),ve(3),dot1,dot2
      real(8) rcel(3),earth_to_geo(3),dist
 
	r=sqrt(dot_product(scpos,scpos))
	ru=scpos/r
 
	if(icase.eq.1) then	 !find pointing vector b(3) to geoid
      sinlat=ru(3)
      coslat=sqrt(ru(1)*ru(1)+ru(2)*ru(2))
      rearth=rp/sqrt(ffac*coslat**2+sinlat**2)  
	deltau=((r-rearth)/re)*sqrt(1 - (1-ffac)*sinlat**2/((ffac*coslat)**2+sinlat**2))
 	b(1)=-ru(1)
	b(2)=-ru(2)
	b(3)=-ru(3)*(1+deltau)/(ffac+deltau)
	b=b/sqrt(dot_product(b,b))
	
cx	b=z	

	endif
 
	if(icase.eq.2 .or. icase.eq.3) then	 !compute b(3) from alpha and beta
      cosbeta=cosd(beta)
      sinbeta=sind(beta)
      bx=sinbeta*cosd(alpha) 
      by=sinbeta*sind(alpha) 
      bz=cosbeta				
      b=bx*x + by*y + bz*z
	endif
 
c     for icase=4, b(3) is an input
	
      call loccel_range(r,ru(1),ru(2),ru(3),b(1),b(2),b(3),re,rp, cel1,cel2,cel3,rearth,range,ierror)
 
	if(ierror.ne.0) then
      refl=0
	xlat=0
	xlon_rot=0
	range=0
	rearth=0
	thtinc=0
	azim=0
	pra=0
	suninc=0
	sunazm=0
	sunglt=0
	return
	endif

	if(icase.eq.4) return  !case 4 is to see if s/c can see sun, so there is nothing else to computed
 
      coslat=sqrt(cel1*cel1+cel2*cel2)     !cosine of geocentric lat
      xlat=atand(cel3/(ffac*coslat))       !geodetic cell latitude
      xlon=atan2d(cel2,cel1)
      call fixang8( xlon)

	xlon_rot=xlon - rotangle	!convert from inertial to earth
      call fixang8( xlon_rot)

	if(icase.ne.3) return  !nothing else to do for case 1,2 and 4
 
 
c     =================================================================================================================
c     ======= compute earth inc. ang. (thtinc) and earth azimuth angle wrt. north (azim) ==============================
c     =================================================================================================================
 
      coslon=cosd(xlon)
      sinlon=sind(xlon)
      coslat=cosd(xlat)
      geoid(1)=coslat*coslon
      geoid(2)=coslat*sinlon
      geoid(3)=sind(xlat)
 
      costht=-dot_product_unit8(b,geoid)
      thtinc=acosd(costht)
 
	
	if(abs(costht).gt.0.9999999999d0) then	!too close to nadir to get azim
	azim=0
	else
      cosazim=(-sinlon*b(1)+coslon*b(2))/sqrt(1.-costht*costht)
      if(cosazim.lt.-1.) cosazim=-1.
      if(cosazim.gt. 1.) cosazim= 1.
      azim=acosd(cosazim)
      if(b(3)+costht*geoid(3).lt.0.) azim=-azim
      azim=90.-azim	 ! convert to relative to clockwise from north
      if(azim.lt.  0.) azim=azim+360.
      if(azim.gt.360.) azim=azim-360.
	endif
 
c     =================================================================================================================
c     ======= compute sun glint angle (sunglt), sun inc. ang. (suninc) and sun azimuth angle (sunazm) =================
c     =================================================================================================================
 
	refl=b + 2*costht*geoid
 
	cossunglt=dot_product_unit8(refl,sunvec)
      sunglt=acosd(cossunglt)
 
c     compute sun tht and azimuth angles
      costht=dot_product_unit8(sunvec,geoid)
      suninc=acosd(costht)
 
	if(abs(costht).gt.0.9999999999d0) then	!too close to nadir to get sunazm
	sunazm=0
	else
      cosazim=(-sinlon*sunvec(1)+coslon*sunvec(2))/sqrt(1.-costht*costht)
      if(cosazim.lt.-1.) cosazim=-1.
      if(cosazim.gt. 1.) cosazim= 1.
      sunazm=acosd(cosazim)
      if(sunvec(3)-costht*geoid(3).lt.0.) sunazm=-sunazm	 !minus sign used here becuase sun is pointing away from surface
      sunazm=90-sunazm	 ! convert to relative to clockwise from north
      call fixang8(sunazm)
	endif
 
 
 
c     =================================================================================================================
c     ============================= compute polarization rotation angle (pra) =========================================
c     =================================================================================================================
 
 
      call cross_norm(b,z, ha) 	 !ha=kx(-z)=(-b)x(-z)=bxz
 
	if(thtinc.gt.0.01) then
	call cross_norm(geoid,b, he) !he=kxg=(-b)xg=gxb
	else
	he=ha
	endif
 
	call cross_norm(b,he, ve)	 !ve=hexk=bxhe
 
	dot1=dot_product(ve,ha)
	dot2=dot_product(he,ha)
	pra=atan2d(dot1,dot2)
 

c     ========================================================================================================================
c     ==================================== 9/20/2013 update: geo lat/lon are computed ========================================
c     ========================================================================================================================
      
      rcel(1)=cel1
      rcel(2)=cel2
      rcel(3)=cel3
      
      dot1=dot_product(rcel,refl)
      
      dist=-rearth*dot1 + sqrt((rearth*dot1)**2 + (geosync_alt**2 - rearth**2))
      
      earth_to_geo=rearth*rcel + dist*refl
      
      dist=sqrt(earth_to_geo(1)**2 +earth_to_geo(2)**2 +earth_to_geo(3)**2)
      
      if(abs(dist-geosync_alt).gt.10.) stop 'earth_to_geo in error, pgm stopped' !exceeds 10 meter
      
      earth_to_geo=earth_to_geo/dist

	geolat=asind(earth_to_geo(3))
	geolon=atan2d(earth_to_geo(2),earth_to_geo(1)) - rotangle

      return
	end
 
 
 
 
 
 
      subroutine loccel_range(r,r1,r2,r3,s1,s2,s3,re,rp, cel1,cel2,cel3,rearth,range,ierror)
	implicit none
 
	real(8), intent(in)  :: r,r1,r2,r3,s1,s2,s3
	real(8), intent(out) :: cel1,cel2,cel3,rearth,range
	real(8) delta,a,cosx,b,c,celmag,arg
	real(8) re,rp
	integer(4) ierror
 
      delta=(re/rp)*(re/rp)-1.
      a=1.+delta*s3*s3
      cosx=-r1*s1-r2*s2-r3*s3
      b=cosx-delta*r3*s3
      c=1.-(re/r)*(re/r)+delta*r3*r3
 
	arg=b*b-a*c
	if(arg.lt.0) then !vector does not intercept earth
	ierror=1
	return
	else
	ierror=0
	endif
 
      range=(b-sqrt(arg))/a
      cel1=r1+range*s1
      cel2=r2+range*s2
      cel3=r3+range*s3
      celmag=sqrt(cel1*cel1+cel2*cel2+cel3*cel3)
      cel1=cel1/celmag
      cel2=cel2/celmag
      cel3=cel3/celmag
      rearth=r*celmag
      range=range*r
	if(range.le.0) ierror=1
      return
      end
 
 
 
 
				
c     this routine is modeled after the rotation matrix given in 'o:\tmi\level1\tmi_geoloc1.f'
c     this rotation is different than used for amsr.  amsr says it uses the '123' convention
c     the tmi/gmi handbooks say gmi uses 3-2-1 Euler
c     also i found that tmi/gmi uses a different x,y,z defintion.  
c     gmi and amsr x are the same (velocity), but the other two (y and z) are negative of each other.
c     tmi/gmi handbooks say: the Z-axis is toward the geocentric nadir
c     this is the same as my usual convention so it must be that for amsr jaxa uses opposite convention when defining rpy.

c     here the x,y,z are the amsr convention
c     so before doing the gmi rotation i apply a negative to y and z
c     after the rotation i apply a negative sign again to convert back to amsr convention
c     the rotation matrix is taken directly from 'o:\tmi\level1\tmi_geoloc1.f'
c     for tmi, i verified this works by testing a large rpy.
c     
 
	subroutine get_sc_axes_from_rpyxxx(scpos,y0,roll,pitch,yaw, x,y,z)
 
	implicit none
 
	!arguments
	real(8), intent(in):: scpos(3), y0(3), roll,pitch,yaw
	real(8), intent(out)::  x(3), y(3), z(3)
	
	!local variables
	real(8) x0(3), z0(3)
	real(8) att_mat(3,3),a0(3,3), a(3,3)
	real(8) cos_roll,sin_roll,cos_pitch,sin_pitch,cos_yaw,sin_yaw
 
  			
c     x0,y0,z0 are s/c axis for 0 roll,pitch,yaw
	z0=-scpos/sqrt(dot_product(scpos,scpos))
	call cross_norm(y0,z0, x0)

c     minus signs convert to tmi convention 
	a0(1,:)= x0
	a0(2,:)=-y0
	a0(3,:)=-z0
 
	cos_roll= cosd(roll)
	sin_roll= sind(roll)
	cos_pitch=cosd(pitch)
	sin_pitch=sind(pitch)
	cos_yaw=  cosd(yaw)
	sin_yaw=  sind(yaw)
 
 
c     goddard  3-2-1 convention fromm 'o:\tmi\level1\tmi_geoloc1.f'
 	att_mat(1,1)= cos_pitch*cos_yaw	      
	att_mat(1,2)=-cos_pitch*sin_yaw	      
	att_mat(1,3)= sin_pitch	      
	att_mat(2,1)= cos_roll*sin_yaw - sin_roll*sin_pitch*cos_yaw  	      
	att_mat(2,2)= cos_roll*cos_yaw + sin_roll*sin_pitch*sin_yaw        
	att_mat(2,3)= sin_roll*cos_pitch	      
	att_mat(3,1)=-sin_roll*sin_yaw - cos_roll*sin_pitch*cos_yaw 	      
	att_mat(3,2)=-sin_roll*cos_yaw + cos_roll*sin_pitch*sin_yaw 	      
	att_mat(3,3)= cos_roll*cos_pitch

	a=matmul(att_mat,a0)

c     minus signs convert back to amsr convention 
	x= a(1,:)
	y=-a(2,:)
	z=-a(3,:)
 
	return
	end
 
 
	subroutine get_sc_axes_from_rpy(scpos,y0,roll,pitch,yaw, x,y,z)
	implicit none

	!arguments
	real(8), intent(in):: scpos(3), y0(3), roll,pitch,yaw
	real(8), intent(out)::  x(3), y(3), z(3)
	
	!local variables
	real(8) x0(3), z0(3)
	real(8) att_mat(3,3),a0(3,3), a(3,3)
	real(8) cos_roll,sin_roll,cos_pitch,sin_pitch,cos_yaw,sin_yaw

  			  
c     x0,y0,z0 are s/c axis for 0 roll,pitch,yaw
	z0=-scpos/sqrt(dot_product(scpos,scpos))
	call cross_norm(y0,z0, x0)

	a0(1,:)=x0
	a0(2,:)=y0
	a0(3,:)=z0

	cos_roll= cosd(roll)
	sin_roll= sind(roll)  
	cos_pitch=cosd(pitch)
	sin_pitch=sind(pitch)
	cos_yaw=  cosd(yaw)
	sin_yaw=  sind(yaw)

c     following is the 'rpy' or  '123' convention

	att_mat(1,1)= cos_pitch*cos_yaw	      
	att_mat(1,2)= sin_roll*sin_pitch*cos_yaw + cos_roll*sin_yaw	      
	att_mat(1,3)=-cos_roll*sin_pitch*cos_yaw + sin_roll*sin_yaw	      
	att_mat(2,1)=-cos_pitch*sin_yaw	      
	att_mat(2,2)=-sin_roll*sin_pitch*sin_yaw + cos_roll*cos_yaw	      
	att_mat(2,3)= cos_roll*sin_pitch*sin_yaw + sin_roll*cos_yaw	      
	att_mat(3,1)= sin_pitch	      
	att_mat(3,2)=-sin_roll*cos_pitch	      
	att_mat(3,3)= cos_roll*cos_pitch

	a=matmul(att_mat,a0)

	x=a(1,:)
	y=a(2,:)
	z=a(3,:)

	return
	end

 
 
	subroutine pitch_sc(pitch, x,y,z)
 
	implicit none
	
	real(8) x(3), y(3), z(3)
	real(8) att_mat(3,3),a0(3,3), a(3,3)
	real(8) pitch,cos_pitch,sin_pitch
 
	if(pitch.eq.0) return
	
	a0(1,:)=x
	a0(2,:)=y
	a0(3,:)=z
 
	cos_pitch=cosd(pitch)
	sin_pitch=sind(pitch)
 
	att_mat(1,1)= cos_pitch	
	att_mat(1,2)= 0	
	att_mat(1,3)=-sin_pitch	
	att_mat(2,1)= 0	
	att_mat(2,2)= 1	
	att_mat(2,3)= 0	
	att_mat(3,1)= sin_pitch	
	att_mat(3,2)= 0	
	att_mat(3,3)= cos_pitch
 
	a=matmul(att_mat,a0)
 
	x=a(1,:)
	y=a(2,:)
	z=a(3,:)
 
	return
	end
 
 
	subroutine roll_sc(roll, x,y,z)
 
	implicit none
	
	real(8) x(3), y(3), z(3)
	real(8) att_mat(3,3),a0(3,3), a(3,3)
	real(8) roll,cos_roll,sin_roll
	
	if(roll.eq.0) return
 
	a0(1,:)=x
	a0(2,:)=y
	a0(3,:)=z
 
	cos_roll=cosd(roll)
	sin_roll=sind(roll)
 
	att_mat(1,1)= 1	
	att_mat(1,2)= 0	
	att_mat(1,3)= 0	
	att_mat(2,1)= 0
	att_mat(2,2)= cos_roll	
	att_mat(2,3)= sin_roll	
	att_mat(3,1)= 0
	att_mat(3,2)= -sin_roll	
	att_mat(3,3)= cos_roll
 
	a=matmul(att_mat,a0)
 
	x=a(1,:)
	y=a(2,:)
	z=a(3,:)
 
	return
	end
 





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


	 