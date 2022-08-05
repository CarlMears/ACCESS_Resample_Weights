    program find_target_gains                                            

    use simulation_param, only: set_simulation_params,SimulationParameterData_f,NEXTRAFOV,NEXTRASCAN
    use tdr_io, only: TypeTDR, tdr_alloc, tdr_free, tdr_read_nc
    use math_routines, only: cross_norm,dot_product_unit8,invert_3by3,fixang4,fixang8,minang4,minang8
                                         
    implicit none

    type(TypeTDR) :: tdr_data
    type(SimulationParameterData_f)    :: sim_parameters
    character(120) filename,filename_tdr                    
    integer(4) iorbit, extra_fov,extra_scan
    integer(4) iscan0,iscan,center_fov
    real(4) xlat0,xlon0,xlat_center_last

    logical :: found_ascending_crossing

    integer(4) :: footprint_size_int
    character(10)  :: sat_name 

    sat_name = 'AMSR2'
    footprint_size_int = 50

    NEXTRAFOV = 1
    NEXTRASCAN = 1
    
    call set_simulation_params(trim(sat_name),sim_parameters)
    call tdr_alloc(tdr_data)
    iorbit = 2

    ! read in the TDR file.  This file is constructed to already contain extra FOVs and Scans
    ! to produce the more finely spaced target locations.
    ! produced by simulated_tdr_extra_fovs.f90

    write(filename_tdr,9001) trim(sat_name),iorbit,NEXTRAFOV,NEXTRASCAN
    9001 format('/mnt/ops1p-ren/l/access/resampling/',a,'/tdr/r',i5.5,'.tdr.extra_fov_',i1,'_scan_',i1,'.h5')
    call tdr_read_nc(tdr_data, trim(filename_tdr), .false.)  !.false. means not fully polarimetric

    ! the TDR file contains all the information needed -- SC position, velocity, roll, pitch, yaw, etc.
    ! as well as the footprint locations compute by the TDR simulator.  
    write(filename,9002) trim(sat_name),footprint_size_int
    9002 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/target_gains/logs/find_target_gain_circular_',i2.2,'km.log')
    open(7,file=filename, status='replace')
    
    write(filename,9013) trim(sat_name),footprint_size_int
    9013 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/target_gains/locs/find_target_gain_circular_',i2.2,'km.loc')
    !open(8,file=filename, status='replace')
    !look for the ascending node crossing scan in the tdr file

    xlat_center_last = 85.5
    center_fov = (tdr_data%nfov + 1)/2
    do iscan = 1,tdr_data%nscan
        if (xlat_center_last < 0.0) then
            if (tdr_data%latitude(center_fov,iscan) >= 0.0) then
                found_ascending_crossing = .true.
                iscan0 = iscan
                xlat0 = tdr_data%latitude(center_fov,iscan)
                xlon0 = tdr_data%longitude(center_fov,iscan)
                exit
            endif
        endif
        xlat_center_last = tdr_data%latitude(center_fov,iscan)
    enddo

    ! these are determined from examining L/amsr2_sim_test/r00002.tdr.test  in panoply
    if (found_ascending_crossing) then
        write(6,44) xlat0,xlon0,iscan0
        write(7,44) xlat0,xlon0,iscan0
        44 format('Found Ascending Crossing, Lat = ',f7.2,' Lon = ',f7.2,' Scan Number = ',i5)
       

        xlat0=0.5*nint(2*xlat0)
        xlon0=0.5*nint(2*xlon0)
        

        call target_gains_simple(tdr_data,iscan0,0,1,tdr_data%nfov,xlat0,xlon0,footprint_size_int)
        close(6)
        close(7)
        close(8)

    else
        print *,'No Crossing Found'
    endif


    stop 'norm end'
    end   

    subroutine target_gains_simple(tdr_data,iscan0,kscan1,kscan2,nfov,xlat0,xlon0,beamwidth_km_int)
        
        use, intrinsic :: iso_fortran_env, only: real32, real64, int8, int32, ERROR_UNIT
        use tdr_io, only: TypeTDR
        use trig_degrees, only: cosd, sind         
        use wgs84, only: distance_between_lon_lats, distance_between_lon_lats_real64 
        use target_patterns, only: fd_gain_targets_km
        
        implicit none
    
        type(TypeTDR),intent(in) :: tdr_data
        
        integer(4),intent(in)   :: iscan0,kscan1,kscan2,nfov
        real(4),intent(in)      :: xlat0,xlon0
        integer(4),intent(in)   :: beamwidth_km_int
        real(4)                 :: beamwidth_km
        
        integer(4), parameter:: nlat=1601
        integer(4), parameter:: nlon=2401  !ssmis value, probably overkill for gmi considering its low altitude
    
        character(80) filename
    
        real(8) gain 
    
        real(8) qsum,xmin,xmax,distance_km_real64
        real(8) gain_domega(nlon,nlat)
        real(4) distance_km
        real(4) :: xlat_grid,xlon_grid,xlon,xlat
        integer(4) error,jscan,ifov
        integer(4) kscan
        integer(4) ilat,ilon,num_valid_pix
        integer(2) ilat2,ilon2
        integer :: ioerr
        integer(4) :: max_ilon,max_ilat
        real(8) :: max_gain
        
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
              
                write(8,222)kscan,ifov,xlat,xlon
                write(6,222)kscan,ifov,xlat,xlon
                222 format(1x,i4.4,1x,i4.4,1x,f12.5,1x,f12.5) 
                max_gain = 0.0
                do ilat=1,nlat
                    do ilon=1,nlon
                         xlat_grid=0.01*(ilat-1) + xlat0 - 8.0
                         xlon_grid=0.01*(ilon-1) + xlon0 -12.0
                    
                         distance_km = distance_between_lon_lats(xlon,xlat,xlon_grid,xlat_grid,error)
    
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
               
                write(filename,9001)beamwidth_km_int,kscan+1,ifov
                9001 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/target_gains/circular_',i2.2,'km/s',i2.2,'c',i4.4,'.dat')

                open(unit=4,file=filename, status='replace',action='write',form='unformatted', access='stream',iostat=ioerr)
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
    
    


