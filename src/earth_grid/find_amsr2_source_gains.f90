program find_source_gains                                            
!use l2_module  
use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int16
use simulation_param, only: set_simulation_params,SimulationParameterData_f,NEXTRAFOV,NEXTRASCAN
!use simulation_param, only: NBAND,NPOLBAND,BAND_FREQS,POLBAND_FREQS,FOVS_PER_SCAN
use tdr_io, only: TypeTDR, tdr_alloc, tdr_free, tdr_read_nc
use math_routines, only: cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
use trig_degrees, only: sind,cosd,acosd,asind
use geolocation_core_routines, only: get_sc_axes_ecef,get_sc_axes,GeolocationScanlineData,prepare_virtual_fov,prepare_fov
use geolocation_core_routines, only: GeolocationCommonParameterData, GeolocationCommonInputData, &
                                    GeolocationSatelliteFOV, GeolocationGroundFOV, geolocate_onto_ground
use geolocation_types, only: geolocation_sensor_init
use pre_nut_routines, only: get_precession_matrix, get_gm_angle
use time_conversions, only: ccsds_utc_to_j2000_ut1,ccsds_utc_to_j2000_tt
use time_conversions, only: j2000_utc_to_j2000_ut1,j2000_utc_to_j2000_tt
use mwi_geolocation, only: set_gl_common_from_tdr,ellipsoid,init_dem_from_gl
use mwi_geolocation, only: geolocation_init_f, geolocation_runCompute_f
use geolocation_aux, only: geolocation_ancillary_load
use geolocation_extra_routines, only: find_subsat_point,convert_j2000_vec_to_ecef
use geolocation_core_routines, only: print_fov
use geolocation_types, only: geolocation_ancillary_init, geolocation_ancillary_free,  &
                            GeolocationAncillaryData_f,GeolocationParameterData_f
use dem, only: DigitalElevationMap
                                    
implicit none

type(TypeTDR) :: tdr_data
type(SimulationParameterData_f) :: sim_parameters
type(GeolocationScanlineData)   :: scanline
type(GeolocationSatelliteFOV)   :: satellite_fov
type(GeolocationGroundFOV)      :: ground_fov

type(GeolocationCommonParameterData) :: gl_parameters
type(GeolocationCommonInputData) :: gl_input

type(GeolocationParameterData_f) :: geoloc_parameters
type(GeolocationAncillaryData_f) :: geoloc_aux_data
type(DigitalElevationMap) :: ancillary_dem

character(120) filename,filename_tdr                    
integer(int32) iorbit
integer(int32) iscan0,iscan,center_fov,ifov,ilat,ilon,iscan1,iscan2
integer(int16) ilat2,ilon2
real(real32) xlat0,xlon0,xlat_center_last
real(real64) xlat,xlon
real(real64) time_this_scan_utc
real(real64) time_next_scan_utc
integer(int32), parameter:: nlat=1601
integer(int32), parameter:: nlon=2401 
real(real64), parameter :: SEC_PER_DAY = 86400.
real(real64) cell0(3,nlon,nlat),cell(3,nlon,nlat)
real(real64) :: xlook(3)
real(real64) :: range,delta,cosdelta,costht
logical :: found_ascending_crossing

integer(int32) :: footprint_size_int
integer(int32) :: iband,native_fov
character(10)  :: sat_name 

real(8), dimension(3) :: sc_pos_ecef,boresight_ecef

real(8),parameter :: Ident_3x3(3,3) = reshape((/1.0D0,0.0D0,0.0D0,  &
                                                0.0D0,1.0D0,0.0D0,  &
                                                0.0D0,0.0D0,1.0D0/),shape(Ident_3x3))


    !next block only for debugging
real(8) :: min_range
real(8) :: min_angle,min_lat,min_lon
integer(4) :: min_ilat,min_ilon
real(8) :: sc_lat, sc_lon, sc_alt,gm_angle
logical :: found
real(8),dimension(3) :: sc_x, sc_y, sc_z

!these describe the antenna patterns for AMSR2 (or maybe AMSRE !?)
real(8), parameter :: sin_anglim(8)=(/  0.105d0,  0.105d0,  0.0524d0, 0.0524d0, 0.0262d0, 0.000d0,  0.0262d0 ,  0.0262d0 /)
real(8), parameter :: ant_approx_a(8)=(/4.343d-6, 2.096d-6,    1.890d-6, 1.623d-6, 7.248d-7, 0.000d0,  2.070d-6,   2.070d-6 /)
real(8), parameter :: ant_approx_b(8)=(/6.892d-4, 4.059d-4,    3.727d-4, 7.251d-4, 3.051d-4, 0.000d0,  2.381d-4,   2.381d-4 /)
real(8), parameter :: ant_approx_c(8)=(/0.503d0,  0.662d0,  1.391d0,  1.804d0,  1.964d0,  0.000d0,  4.593d0,    4.593d0  /)
real(8), parameter :: ant_approx_d(8)=(/0.651d0,  1.345d0,     4.844d0,  4.721d0,  15.18d0,  0.000d0, 79.785d0,   79.785d0  /)
real(8) :: anglim

real(8) :: gain
real(real64) :: g(nlon,nlat)
real(real64) :: gsum,gmax

real(real64) :: max_gain
integer(int32) :: max_ilat,max_ilon,num_valid_pix,scan_num_for_file 

integer :: ioerr

sat_name = 'AMSR2'
footprint_size_int = 30

NEXTRAFOV = 1
NEXTRASCAN = 1

call geolocation_ancillary_init(geoloc_aux_data)
call geolocation_ancillary_load(geoloc_aux_data, "mwsdps-gl-anc-v3.nc", ioerr)  ! this loads the dem and land masks
if (ioerr /= 0) error stop 1

ioerr = geolocation_init_f(geoloc_parameters, geoloc_aux_data)  !transfers into rss verions of the sat parameters and geoloc data
if (ioerr /= 0) error stop 1

call init_dem_from_gl(geoloc_aux_data, ancillary_dem)

call set_simulation_params(trim(sat_name),sim_parameters)
ellipsoid%EarthEquatorialRadius = sim_parameters%re 
ellipsoid%EarthPolarRadius = sim_parameters%rp 

call tdr_alloc(tdr_data)
iorbit = 2

! read in the TDR file.  This file is constructed to already contain extra FOVs and Scans
! to produce the more finely spaced target locations.
! produced by simulated_tdr_extra_fovs.f90

write(filename_tdr,9001) trim(sat_name),iorbit,NEXTRAFOV,NEXTRASCAN
9001 format('/mnt/ops1p-ren/l/access/resampling/',a,'/tdr/r',i5.5,'.tdr.extra_fov_',i1,'_scan_',i1,'.h5')
!write(filename_tdr,9001) trim(sat_name),iorbit
!9001 format('/mnt/ops1p-ren/l/access/resampling/',a,'/tdr/r',i5.5,'.tdr.test.h5')
call tdr_read_nc(tdr_data, trim(filename_tdr), .false.)  !.false. means not fully polarimetric

! the TDR file contains all the information needed -- SC position, velocity, roll, pitch, yaw, etc.
! as well as the footprint locations compute by the TDR simulator.  
do iband = 1,5
    write(filename,9012) trim(sat_name),trim(sat_name),iband
    9012 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/source_gains/logs/find_source_gain_',a,'_band_',i2.2,'.log')
    open(7,file=filename, status='replace')

    !write(filename,9022) trim(sat_name),trim(sat_name),iband
    !9022 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/source_gains/logs/find_source_gain_',a,'_band_',i2.2,'.debug')
    !open(9,file=filename, status='replace')

    write(filename,9013) trim(sat_name),trim(sat_name),iband
    9013 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/source_gains/locs/find_source_gain_',a,'_band_',i2.2,'_locs.txt')
    open(8,file=filename, status='replace')


    !transfer data needed for the GL routines into the common data structures
    call set_gl_common_from_tdr(sim_parameters, gl_parameters, tdr_data, gl_input)
    
    
    !look for the ascending node crossing scan in the tdr file

    xlat_center_last = 85.5
    center_fov = (tdr_data%nfov + 1)/2
    do iscan = 1,tdr_data%nscan
        if (xlat_center_last < 0.0) then
            !if (tdr_data%sclat(iscan) >= 0.0) then
            if (tdr_data%latitude(center_fov,iscan) >= 0.0) then
                found_ascending_crossing = .true.
                iscan0 = iscan
                xlat0 = tdr_data%latitude(center_fov,iscan)
                xlon0 = tdr_data%longitude(center_fov,iscan)
                exit
            endif
        endif
        
        xlat_center_last = tdr_data%latitude(center_fov,iscan)
        !xlat_center_last = tdr_data%sclat(iscan)
    enddo

    ! these are determined from examining L/amsr2_sim_test/r00002.tdr.test  in panoply
    if (found_ascending_crossing) then
        write(6,44) xlat0,xlon0,iscan0
        write(7,44) xlat0,xlon0,iscan0
        44 format('Found Ascending Crossing, Lat = ',f7.2,' Lon = ',f7.2,' Scan Number = ',i5)
       
        write(6,55) tdr_data%sclat(iscan0),tdr_data%sclon(iscan0)
        55 format('SC Lat = ',f7.2,' Lon = ',f7.2)

        ! load cell0 amd cell arrays
        ! these are the unit vector and the vector pointing from earth center to lat/lon cell.

        scanline%time_j2000_ut1 = tdr_data%scan_time_utc(iscan0)
        call get_gm_angle(scanline%time_j2000_ut1, gm_angle)

        xlat0=0.5*nint(2*xlat0)
        xlon0=0.5*nint(2*(xlon0))

        do ilat = 1,nlat 
            do ilon = 1,nlon 
                xlat = 0.01D0*(ilat-1) + xlat0 - 8.0D0
                xlon = 0.01D0*(ilon-1) + xlon0 -12.0D0 
                call fd_cell_vector(nlat,nlon,ilat,ilon,xlat,xlon, cell0,cell)
            enddo
        enddo

        iscan1=iscan0 - 14*(sim_parameters%nextrascan + 1)
        iscan2=iscan0 + 14*(sim_parameters%nextrascan + 1)
         
        do iscan=iscan1,iscan2,(sim_parameters%nextrascan + 1)
        !do iscan=iscan0 ,iscan0+2,(sim_parameters%nextrascan + 1)

            ! prepare the scanline structure for this scan
            scanline%iscan = iscan
            scanline%iscan_native = iscan/2

            scan_num_for_file = 1 + (iscan-iscan1)/2

            scanline%time_j2000_tt = (j2000_utc_to_j2000_tt(tdr_data%scan_time_utc(iscan),  &
                                                        tdr_data%leap_seconds(iscan))) / SEC_PER_DAY
        
            call get_precession_matrix(scanline%time_j2000_tt, scanline%precession)
            
            ! The satellite position/velocity in ECI coordinates at the mean
            ! epoch of date at the start of the scan
           
            scanline%scpos = tdr_data%scpos(iscan,:)
            scanline%scvel = tdr_data%scvel(iscan,:)
            ! The solar position in ECI coordinates at the mean epoch of date
            scanline%sunvec = tdr_data%solar_vector_mics(iscan,:)

            ! The time to complete a full scan in seconds
            scanline%scan_period = 60.0/tdr_data%scan_rate(iscan)

            ! The attitude quaternions for the current scan and the next real scan (for interpolation)
            scanline%att_this_scan = tdr_data%scatt(iscan,:)
            scanline%att_next_scan = tdr_data%scatt(iscan+sim_parameters%nextrascan+1,:)

            call get_sc_axes(scanline%att_this_scan, scanline%precession,sc_x, sc_y, sc_z)

            ! The time in UTC for the current scan and the next real scan (for interpolation)
            time_this_scan_utc = tdr_data%scan_time_utc(iscan0)
            time_next_scan_utc = tdr_data%scan_time_utc(iscan0+sim_parameters%nextrascan+1)
         
            scanline%time_this_scan = time_this_scan_utc
            scanline%time_next_scan = time_next_scan_utc

            ! find the source gain for each location

            anglim=asind(sin_anglim(iband))


            do ifov = 1,tdr_data%nfov,(sim_parameters%NEXTRAFOV + 1)
                native_fov = (1+ifov)/(sim_parameters%NEXTRAFOV + 1)
                print *,'fov = ',ifov,' native FOV = ',native_fov
                ! calling with nextrascan and nextrafov = 0 because input data already includs these scans/fovs
                ! setting the ifov parameter to native fov so that the times work out right 
                satellite_fov = prepare_fov(scanline, gl_parameters, gl_input, native_fov, iband)
                !
                
                !call print_fov(satellite_fov,6)
            
                call find_subsat_point(ellipsoid, satellite_fov%scpos, satellite_fov%gm_angle, sc_lat, sc_lon, sc_alt, found)

                !print *,satellite_fov%gm_angle, sc_lat, sc_lon, sc_alt

               
                ground_fov = geolocate_onto_ground(satellite_fov, ellipsoid, ancillary_dem)
            
                 !rotate sc pos and boresight into ecef system

                call convert_j2000_vec_to_ecef(satellite_fov%scpos,Ident_3x3, satellite_fov%gm_angle,sc_pos_ecef)
                call convert_j2000_vec_to_ecef(satellite_fov%boresight,Ident_3x3, satellite_fov%gm_angle,boresight_ecef)

                min_range = 1.0d50
                min_angle = 360.0
                g = 0.0D0
                gsum = 0.0D0
                gmax = -1.0D0
                do ilat=1,nlat
                    do ilon=1,nlon

                        xlook=  cell(:,ilon,ilat) - sc_pos_ecef      !point from s/c down to cell, cell is in ecef

                        range=sqrt(dot_product(xlook,xlook))
                        if (range < min_range) min_range = range
                        xlook=xlook/range

                        costht=-dot_product(xlook,cell0(:,ilon,ilat))
                        if(costht.le.0) stop 'pgm stopped,costht.le.0'
                       
                        cosdelta=dot_product(xlook,boresight_ecef)
                        if(cosdelta.gt.1) then
                            cosdelta=1  !should never be near -1
                        endif
                        
                        delta=acosd(cosdelta)



                        if (delta.le.anglim) then
                            
                            gain = ant_approx_a(iband) + &
                                ant_approx_b(iband)*exp(-ant_approx_c(iband)*delta) + &
                                exp(-ant_approx_d(iband)*delta*delta)
                            gain = gain*costht/(range*range)
                            g(ilon,ilat) = gain
                            gsum = gsum + gain
                        endif

                        if (gain > gmax) then
                            gmax = gain
                            call find_subsat_point(ellipsoid, cell(:,ilon,ilat), 0.0D0, sc_lat, sc_lon, sc_alt, found)
                            min_lat = sc_lat
                            min_lon = sc_lon
                            min_ilat = ilat
                            min_ilon = ilon
                        endif
                enddo !ilat
            enddo !ilon

            write(8,222)iband,scan_num_for_file,native_fov,ground_fov%lat,ground_fov%lon,min_lat,min_lon, &
            min_lat -ground_fov%lat,min_lon-ground_fov%lon
            write(6,222)iband,scan_num_for_file,native_fov,ground_fov%lat,ground_fov%lon,min_lat,min_lon, &
                  min_lat -ground_fov%lat,min_lon-ground_fov%lon
            222 format(1x,i4.4,1x,i4.4,1x,i4.4,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5,1x,f12.5)

            if ((abs(ground_fov%lon - min_lon) > 0.025) .or. &
                (abs(ground_fov%lat - min_lat) > 0.025)) then
                print *,'Excess Lat/Lon Difference'
            endif

            write(filename,9002)iband,scan_num_for_file,native_fov
            9002 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/source_gains/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')

            open(unit=4,file=filename, status='replace',action='write',form='unformatted', access='stream')
            num_valid_pix = 0
            max_gain = 0.0

            do ilat=1,nlat
                do ilon=1,nlon
                    if(abs(g(ilon,ilat)) .gt. 1.0d-30) then
                        if (g(ilon,ilat) > max_gain) then
                            max_gain = g(ilon,ilat)
                            max_ilat = ilat
                            max_ilon = ilon
                        endif
                        ilat2=int(ilat,kind=2)
                        ilon2=int(ilon,kind=2)
                        write(4) ilat2,ilon2,g(ilon,ilat)
                        num_valid_pix = num_valid_pix + 1
                    endif
                enddo  !ilon
            enddo  !ilat
            close(4)

            write(6,10)num_valid_pix,trim(filename),max_ilat,max_ilon
            write(7,10)num_valid_pix,trim(filename),max_ilat,max_ilon
            10 format('Wrote ',i7,' valid pix to ',a,'  Max Gain at ',i4,2x,i4)
        enddo !ifov
    enddo !iscan

    close(6)
    close(7)
    close(8)


    else
        print *,'No Crossing Found'
    endif

enddo


    stop 'norm end'
    end   



    
    subroutine target_gains_simple(tdr_data,iscan0,kscan1,kscan2,nfov,xlat0,xlon0,beamwidth_km_int)
        
        use, intrinsic :: iso_fortran_env, only: real32, real64, int8, int32, ERROR_UNIT
        use tdr_io, only: TypeTDR
        use trig_degrees, only: cosd, sind         
        use wgs84, only: distance_between_lon_lats      
        use target_patterns, only: fd_gain_targets_km
        
        implicit none
    
        type(TypeTDR),intent(in) :: tdr_data
        
        integer(int32),intent(in)   :: iscan0,kscan1,kscan2,nfov
        real(real32),intent(in)      :: xlat0,xlon0
        integer(int32),intent(in)   :: beamwidth_km_int
        real(real32)                 :: beamwidth_km
        
        integer(int32), parameter:: nlat=1601
        integer(int32), parameter:: nlon=2401  !ssmis value, probably overkill for gmi considering its low altitude
    
        character(80) filename
    
        real(real64) gain 
    
        real(real64) qsum,xmin,xmax
        real(real64) gain_domega(nlon,nlat)
        real(real32) distance_km
        real(real32) :: xlat_grid,xlon_grid,xlon,xlat
        integer(int32) error,jscan,ifov
        integer(int32) kscan
        integer(int32) ilat,ilon,num_valid_pix
        integer(2) ilat2,ilon2
        integer :: ioerr
        integer(int32) :: max_ilon,max_ilat
        real(real64) :: max_gain
        
        beamwidth_km = real(beamwidth_km_int)

        ! in this version, the target locations are already contained in the extended TDR data file.
        ! for the idealized targets, we can just construct a circular gaussian centered at xlat, xlon
    
        gain_domega=0
        qsum=0 

        do kscan = kscan1,kscan2
            jscan = iscan0+kscan
            do ifov = 1,nfov
                xlat = tdr_data%latitude(ifov,jscan)
                xlon = tdr_data%longitude(ifov,jscan)
                print *,jscan,ifov,xlat,xlon
                max_gain = 0.0
                do ilat=1,nlat
                    do ilon=1,nlon
                        xlat_grid=0.01*(ilat-1) + xlat0 - 8.0
                        xlon_grid=0.01*(ilon-1) + xlon0 -12.0
                    
                         distance_km = distance_between_lon_lats(real(xlon),real(xlat),real(xlon_grid),real(xlat_grid),error)
    
                         if(distance_km .le. 3*beamwidth_km) then
                             gain = fd_gain_targets_km(beamwidth_km,distance_km)
                             gain_domega(ilon,ilat) = gain
                             qsum = qsum + gain
                             if (gain > max_gain) then
                                max_ilat = ilat
                                max_ilon = ilon
                                max_gain = gain
                             endif
                         else
                             gain_domega(ilon,ilat) = 0.0
                         endif
                    enddo !ilon
                enddo !ilat
               
    
                if(qsum.lt.xmin) xmin=qsum
                if(qsum.gt.xmax) xmax=qsum
    
                write(filename,9021)beamwidth_km_int,kscan+1,ifov
                9021 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/target_gains/circular_',i2.2,'km/s',i2.2,'c',i4.4,'.dat')

                open(unit=4,file=filename, status='new',action='write',form='unformatted', access='stream',iostat=ioerr)
                num_valid_pix = 0
                max_gain = 0.0

                do ilat=1,nlat
                     do ilon=1,nlon
                         if(abs(gain_domega(ilon,ilat)) .gt. 1.0d-30) then
                             ilat2=int(ilat,kind=2)
                             ilon2=int(ilon,kind=2)
                             write(4) ilat2,ilon2,gain_domega(ilon,ilat)
                             num_valid_pix = num_valid_pix + 1
                         endif
                     enddo  !ilon
                 enddo  !ilat
    
                 close(4)
                 write(6,10)num_valid_pix,filename,max_ilat,max_ilon
                 write(7,10)num_valid_pix,filename,max_ilat,max_ilon
                 10 format('Wrote ',i7,' valid pix to ',a,'  Max Gain at ',i4,2x,i4)
            enddo  !ifov
        enddo  !kscan
        return
        end subroutine target_gains_simple
    
    
    subroutine fd_cell_vector(nlat,nlon,ilat,ilon,xlat,xlon, cell0,cell)

        ! finds the vector points from the Earth center to the center of the cell at xlat,xlon
        ! vector is calculated in the ecef coordinate system
        ! cell contains the array of vectors
        ! cell0 contains the array of unit vectors

        use, intrinsic :: iso_fortran_env, only: real32, real64, int32
        use trig_degrees, only: sind,cosd

        implicit none

        real(real64), parameter :: re=6378.137d3
        real(real64), parameter :: rp=6356.824d3
    
        integer(int32) nlat,nlon,ilat,ilon
        real(real64) xlat,xlon,rearth
        real(real64) cell0(3,nlon,nlat),cell(3,nlon,nlat)
        real(real64) coslat,sinlat,coslon,sinlon
    
        coslat=cosd(xlat)
        sinlat=sind(xlat)
        coslon=cosd(xlon)
        sinlon=sind(xlon)
    
        cell0(1,ilon,ilat)=coslat*coslon                                
        cell0(2,ilon,ilat)=coslat*sinlon                                
        cell0(3,ilon,ilat)=sinlat
    
        rearth=re*rp/sqrt((rp*coslat)*(rp*coslat)+(re*sinlat)*(re*sinlat))                  
    
        cell(:,ilon,ilat)= rearth*cell0(:,ilon,ilat) 
    
        return
        end
    
                                                    
         
