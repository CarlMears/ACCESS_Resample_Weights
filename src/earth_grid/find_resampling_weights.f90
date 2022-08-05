    
    program find_resampling_weights 

    use, intrinsic :: iso_fortran_env, only: int8, int16, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite

    !use simulation_param, only: NBAND,FOVS_PER_SCAN,NFREQ
    use simulation_param, only: set_simulation_params, SimulationParameterData_f
    use resampling_types, only: FootprintLocations,read_location_text_file
    use wgs84, only: distance_between_lon_lats
    use svd_invert, only: svdinv

    !gain pattern has two types -- CompactGrid is a form where only on-zero elements are
    !stored.  FullGrid is the entire 2D array with zero elements
    use gain_pattern, only: CompactGrid,FullGrid !,nlon,nlat,maxpts
    !some helper routines
    use gain_pattern, only: write_gain_full_grid,convert_to_full_Grid
    use gain_pattern, only: zero_full_grid,find_mean_grid_location
    !this calculates the overlap integrals
    use gain_pattern, only: compute_overlap_full_compact
    implicit none
    

    !integer(int32), parameter:: npixel_max=maxpts  !max number of pixels saved for each observation footprint
    integer(int32), parameter:: nobs_max=1500      !max number of obs (source footprints) saved
                                                   !probably will need to be bigger for larger fooprints
     
    character(120) filename1,target_location_file,source_location_file,weight_file, lisT_file!,filename2,
   
    !integer(int16) :: ilat2(npixel_max,0:nobs_max)
    !integer(int16) :: ilon2(npixel_max,0:nobs_max)   !list of ilat,ilon locations for footprint data
    !!integer(int32) npixel(0:nobs_max)                ! not sure yet
    
    integer(int32) iscansv(nobs_max),icelsv(nobs_max),iobs_numsv(nobs_max)                         ! not sure yet
    integer(int32) iscan,icel,nobs,iobs,jobs,ipixel        ! various incideces
    integer(int32) ierror,iobssv
    integer(int32) ifreq !,itrg

    integer(int16) ilat,ilon
  
    real(real64) xsum

    real(real64), allocatable :: g(:,:),g_inv(:,:),testg(:,:),u(:),v(:),xwork(:),a(:)

    real(real64), allocatable :: g_everywhere(:,:)  !array to store elements of g to prevent recalculation
    real(real64) denom,xcoef,xmax,psum(0:6),smooth_fac,delta_gain
      
    !integer(int32) kscan,icel_trg

    type(CompactGrid),dimension(0:nobs_max) :: gain_array

    type(FullGrid)     :: resampled_gain_full
    type(FullGrid)     :: source_gain_full

    character(len=10) sensor_name
    type(SimulationParameterData_f) :: sensor_data
    type(FootprintLocations)        :: target_locs
    type(FootprintLocations)        :: source_locs
    integer(int32) :: footprint_size_int

    integer(int32) :: target_scan_to_do,target_fov_to_do,i,distance_error
    real(real64)   :: target_lat,target_lon,source_lat,source_lon
    real(real32)   :: distance,distance_threshold
    logical        :: found_target
    integer(int32) :: num_valid_pix
    integer(int32) :: num_calculated,num_found 

    real(real64)   :: weights(243,29,485,2)  !fov_in, scan_in, fov_out, scan_out

    !summary stats for each target
    !integer(int32),dimension(485,2) :: n_obs_summary
    real(real64),dimension(485,2) :: tot_a_summary,tot_sqr_a_summary !,check_invert_summary
 
    

    !temp
    real(real64) :: max_diff

    sensor_name = 'AMSR2'
    footprint_size_int = 70
    distance_threshold = 3.0*footprint_size_int
    ifreq = 1
    
    call set_simulation_params(trim(sensor_name),sensor_data)

    write(list_file,9077) trim(sensor_name), footprint_size_int,ifreq
    9077 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a, &
                '/resample_weights_v3/stats/summary_stats_circular_',i2.2,'_band_',i2.2,'.txt')
    open(8,file=list_file, status='replace')

    write(target_location_file,9013) trim(sensor_name)!,footprint_size_int
    9013 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/target_gains_v2/locs/find_target_gain_circular_30km.loc')

    write(source_location_file,9023) trim(sensor_name),trim(sensor_name),ifreq
    9023 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/source_gains_v2/locs/find_source_gain_',a,'_band_',i2.2,'_locs.txt')
    
    call read_location_text_file(target_location_file,target_locs,2)
    call read_location_text_file(source_location_file,source_locs,3)



    allocate(g_everywhere(source_locs%num_footprints,source_locs%num_footprints))
    g_everywhere = ieee_value(0.0_real64, ieee_quiet_nan)

    do target_scan_to_do = 1,2
        do target_fov_to_do = 1,485
        found_target = .false.
        do i = 1,target_locs%num_footprints
            if ((target_locs%iscan(i) == target_scan_to_do-1) .and. (target_locs%ifov(i) == target_fov_to_do)) then
                found_target = .true.
                exit
            endif
        enddo

        if (.not. (found_target)) error stop
        
        target_lat = real(target_locs%lat(i),kind=real32)
        target_lon = target_locs%lon(i)
        print *,'Target Latitude: ',target_lat
        print *,'Target Longitude: ',target_lon

        !loop over source footprints, noting those that are close enough to the target footprint

        source_locs%use = 0
        nobs = 0
        do i = 1,source_locs%num_footprints
            source_lat = source_locs%lat(i)
            source_lon = source_locs%lon(i)
            
            distance = distance_between_lon_lats(real(source_lon,kind=real32), &
                                                real(source_lat,kind=real32), &
                                                real(target_lon,kind=real32), &
                                                real(target_lat,kind=real32),distance_error)

            if (distance < distance_threshold) then
                !write(6,333) source_locs%iscan(i),source_locs%ifov(i),source_lat,source_lon,distance
                !333 format(1x,i4,1x,i4,3(1x,f10.3))
                source_locs%use(i) = 1   ! setting this means that it will be used
                nobs = nobs + 1
            endif
        enddo

        write(6,334) nobs 
        334 format(i4,' observations will be included')
    

        smooth_fac=1.0d-05    
        !MWI specific!!
        
        !smooth_fac=1.0d-08
        !if(ifreq.eq.1 .and. itrg.ge.2) smooth_fac=1.0d-06  !11 ghz mapping onto 30,25,20 km,use smaller smoothing to presever 30 km.
    
        ! read in the target gain
        write(6,335)
        335 format('reading target gain pattern')

        ! eventually will loop here
        
        write(6,337)trim(sensor_name),target_scan_to_do, target_fov_to_do
        337 format('reading ',a,' target gain pattern, scan: ',i4,' fov: ',i4)

        select case(trim(sensor_name)) 
            case('AMSR2')
                        iobs=0
                        !call read_gobs_amsr2_v2(1,ifreq,target_scan_to_do,target_fov_to_do,gain_array(0))
                        call read_tar_amsr2_v3(footprint_size_int,ifreq,target_scan_to_do,target_fov_to_do,gain_array(0))
                    case default
                        print *,'Sensor ',sensor_name,' not implemented'
                        stop
        end select
                ! the target gain is stored in the first obs position (0) in gobs -- remember gobs is (1:max_pixel,0:max_obs)
        

    ! !     =============================================================================================================================
    ! !     ========================================= compute iscan and icel that correspond to each iobs ===============================
    ! !     =============================================================================================================================
        write(6,336)nobs,trim(sensor_name)
        336 format('reading ',i4,1x,a,' source gain patterns')
        iobs = 0
        do i = 1,source_locs%num_footprints         !loop over the entire source grid
            if (source_locs%use(i) == 0) cycle      !but skip the locations that are not used

            iscan = source_locs%iscan(i)            !indicec of the fovs that are used
            icel = source_locs%ifov(i)
         
            iobs=iobs+1
            if(iobs.gt.nobs_max) stop 'pgm stopped, iobs oob'
            select case(trim(sensor_name)) 
                case('AMSR2')
                    !read the source gain pattern 
                    call read_gobs_amsr2_v3(0,ifreq,iscan,icel,gain_array(iobs))
                case default
                    print *,'Sensor ',sensor_name,' not implemented'
                    stop
            end select

            iscansv(iobs)=iscan   
            icelsv(iobs)=icel
            iobs_numsv(iobs) = i
        enddo

        nobs=iobs

        ! allocate the matrics for the linear algebra steps.
        write(6,341)
        341 format('allocating arrays for svd calculations')
        allocate(g(nobs,nobs),g_inv(nobs,nobs),testg(nobs,nobs), u(nobs), v(nobs), xwork(nobs), a(nobs))
        
        write(6,743)
        743 format('calculating overlap integrals assuming symetry and avoiding recalculation')

        num_calculated = 0
        num_found = 0
        do iobs=1,nobs
            ! expand the first footprint to grid form
            
            call convert_to_full_grid(gain_array(iobs),source_gain_full,xsum)
            if(abs(xsum-1.d0).gt.1.e-12) then
                write(*,*) xsum
                stop 'error in normalization, pgm stopped'
            endif

            do jobs=0,iobs  !THIS loop includes the target footprint jobs == 0
                xsum=0
                if (jobs >= 1) then
                    if (ieee_is_finite(g_everywhere(iobs_numsv(iobs),iobs_numsv(jobs)))) then
                        num_found = num_found + 1
                        g(iobs,jobs) = g_everywhere(iobs_numsv(iobs),iobs_numsv(jobs))
                        cycle 
                    endif
                endif

                call compute_overlap_full_compact(source_gain_full,gain_array(jobs),xsum)

                if(jobs.eq.0) then
                    v(iobs)=xsum
                else
                    g(iobs,jobs)=xsum
                    g_everywhere(iobs_numsv(iobs),iobs_numsv(jobs)) = xsum
                endif
                num_calculated = num_calculated + 1
            enddo     !jobs
        enddo     !iobs

        write(6,166) num_calculated,num_found
        166 format('Calculated ',i6,' Elements, Reused ',i6,' Elements')

        !matrix us symmetric by definition, so reflect (we saved 50% of the calc by only doing one triangle)
        do iobs=1,nobs
            do jobs=1,iobs
                g(jobs,iobs)=g(iobs,jobs)
            enddo
        enddo

        u=1
        ! add in the smoothing factor along the diagonal
        do iobs=1,nobs
            g(iobs,iobs)=g(iobs,iobs) + smooth_fac
        enddo
    
        ! invert
        call svdinv(g,nobs, g_inv,ierror)
        if(ierror.ne.0) stop 'inversion did not converge, pgm stopped'

        testg=matmul(g,g_inv)
        do iobs=1,nobs
            testg(iobs,iobs)=testg(iobs,iobs)-1
        enddo

        max_diff = 0.0D0
        do iobs=1,nobs
            do jobs=1,iobs
                if (abs(testg(iobs,jobs)) > max_diff) max_diff = abs(testg(iobs,jobs))
            enddo
        enddo

        write(6,77) max_diff
        77 format('Max abs value of check matrix: ',f14.11)


        xwork=matmul(g_inv,u)
        denom=dot_product(u,xwork)

        xwork=matmul(g_inv,v)
        xcoef=dot_product(u,xwork)                                

        xwork=v + (1-xcoef)*u/denom

        a=matmul(g_inv,xwork)

        xmax=-1.e30; psum=0

        do iobs=1,nobs
            source_locs%weight(iobs_numsv(iobs)) = a(iobs)
        
            psum(0)=psum(0) + 1
            psum(1)=psum(1) + a(iobs)
            psum(2)=psum(2) + a(iobs)**2
            psum(3)=psum(3) + a(iobs)*iscansv(iobs)
            psum(4)=psum(4) + a(iobs)*icelsv(iobs)
            psum(5)=psum(5) + a(iobs)*source_locs%lon(iobs_numsv(iobs))
            psum(6)=psum(6) + a(iobs)*source_locs%lat(iobs_numsv(iobs))
            
            if(a(iobs).gt.xmax) then
                xmax=a(iobs)
                iobssv=iobs

                iobs_numsv(iobs) = i
            endif

            weights(icelsv(iobs),iscansv(iobs),target_fov_to_do,target_scan_to_do) = a(iobs)  !fov_in, scan_in, fov_out, scan_out
        enddo

        psum(2)=sqrt(psum(2))
        psum(3:6)=psum(3:6)/psum(1)

        write(6,801)
        801 format('Calculating resampled gain pattern')

        call zero_full_grid(resampled_gain_full)
        resampled_gain_full%xlat0 = 0.0
        resampled_gain_full%xlon0 = 0.0
        xsum = 0.0D0
        
        do iobs=1,nobs
            do ipixel=1,gain_array(iobs)%numpts
                ilat       = gain_array(iobs)%ilat(ipixel)
                ilon       = gain_array(iobs)%ilon(ipixel)
                delta_gain = gain_array(iobs)%gain(ipixel)*a(iobs)
                resampled_gain_full%gain(ilon,ilat)  = resampled_gain_full%gain(ilon,ilat) + delta_gain
                xsum=xsum + delta_gain
            enddo !ipixel
        enddo !iobs

        !write out the resampled gain pattern
        write(filename1,9001)footprint_size_int,ifreq,target_scan_to_do,target_fov_to_do
        9001 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/resample_gains_v3/circular_',i2.2,&
                 'km/band_',i2.2,'/s',i2.2,'c',i4.4,'.dat')

        write(6,802) filename1
        802 format('Writing resampled gain pattern to ',a)
        
        open(unit=4,file=filename1, status='replace',action='write',form='unformatted', access='stream')
        
        tot_a_summary(target_fov_to_do,target_scan_to_do) = psum(1)
        tot_sqr_a_summary(target_fov_to_do,target_scan_to_do) = psum(2)

        num_valid_pix = 0
        call write_gain_full_grid(resampled_gain_full,4,1.0d-8,num_valid_pix)

        write(6,6001) target_scan_to_do,target_fov_to_do,  &
                      nobs,maxval(abs(testg)), &
                      tot_a_summary(target_fov_to_do,target_scan_to_do), &
                      tot_sqr_a_summary(target_fov_to_do,target_scan_to_do), &
                      xmax
        write(8,6001) target_scan_to_do,target_fov_to_do,  &
                      nobs,maxval(abs(testg)), &
                      tot_a_summary(target_fov_to_do,target_scan_to_do), &
                      tot_sqr_a_summary(target_fov_to_do,target_scan_to_do), &
                      xmax
        6001 format(i2,i5,i5,e15.5,3f13.8,1x,f8.5,1x,f8.5,1x,f7.2)
        close(4)
        deallocate(g, g_inv, testg, u, v, xwork, a)   
      enddo                                                      
    enddo


    write(weight_file,555) trim(sensor_name),footprint_size_int,ifreq
    555 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/resample_weights_v3/circular_',i2.2,&
        'km/resample_weights_band_',i2.2,'.dat')

    open(unit=5,file=weight_file, status='replace',action='write',form='unformatted', access='stream')
    write(5) weights
    close(5)


    close(8)
    !deallocate(weight)
    ! need to zero and deallocated source locs
    stop 'norm end'                                   
    end