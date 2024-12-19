    
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
    use gain_pattern, only: allocate_full_grid,allocate_compact_grid
    use gain_pattern, only: write_gain_full_grid,convert_to_full_Grid
    use gain_pattern, only: zero_full_grid,find_mean_grid_location
    !this calculates the overlap integrals
    use gain_pattern, only: compute_overlap_full_compact

    !manages for source patterns to avoid excessive reads
    use gain_pattern, only: CompactGrid_Array_FAFO
    use gain_pattern, only: allocate_CompactGrid_Array_FAFO, &
                            insert_CompactGrid_Array_FAFO,   &
                            find_CompactGrid_Array_FAFO


    implicit none
    

    !integer(int32), parameter:: npixel_max=maxpts  !max number of pixels saved for each observation footprint
    integer(int32), parameter:: nobs_max=1500      !max number of obs (source footprints) saved
                                                   !probably will need to be bigger for larger fooprints
     
    character(120) filename1,target_location_file,source_location_file,weight_file, lisT_file,g_precompute_file!,filename2, 
    character(100) cmd_str
    integer(int16) iscansv(nobs_max),icelsv(nobs_max),iobs_numsv(nobs_max)                         ! not sure yet
    integer(int32) nobs,iobs,jobs,ipixel        ! various incideces
    integer(int16) iscan,icel
    integer(int32) ierror,iobssv
    integer(int32) ifreq !,itrg

    integer(int16) ilat,ilon
  
    real(real64) xsum

    real(real64), allocatable :: g(:,:),g_inv(:,:),testg(:,:),u(:),v(:),xwork(:),a(:)

    real(real64), allocatable :: g_everywhere(:,:)  !array to store elements of g to prevent recalculation
    real(real64) denom,xcoef,xmax,psum(0:6),smooth_fac,delta_gain

    type(CompactGrid) :: source_gain
    type(CompactGrid) :: target_gain
    type(CompactGrid_Array_FAFO)  :: gain_array
    integer(int32) :: max_num_pts_per_source_footprint

    type(FullGrid)     :: resampled_gain_full
    type(FullGrid)     :: target_gain_full

    character(len=10) sensor_name
    type(SimulationParameterData_f) :: sensor_data
    type(FootprintLocations)        :: target_locs
    type(FootprintLocations)        :: source_locs
    integer(int32) :: footprint_size_int
    integer(int32) :: max_pixels
    

    integer(int32) :: target_scan_to_do,target_fov_to_do,i,distance_error
    real(real64)   :: target_lat,target_lon,source_lat,source_lon
    real(real32)   :: distance,distance_threshold
    logical        :: found_target,found_source
    integer(int32) :: num_valid_pix
    integer(int32) :: num_calculated,num_non_zero

    integer(int32) :: source_index1,source_index2
    integer(int32) :: num_source_reused,num_source_loaded

    real(real64)   :: weights(486,57,485,2)  !fov_in, scan_in, fov_out, scan_out
    real(real64)   :: gtemp,temp
    integer(int32) :: num_passed,num_fixed
    !summary stats for each target
    !integer(int32),dimension(485,2) :: n_obs_summary
    real(real64),dimension(485,2) :: tot_a_summary,tot_sqr_a_summary !,check_invert_summary

    real(real64) :: max_a,min_a
    integer(int32) :: num_a_nonzero,num_a_positive,num_a_negative = 0
 
    !temp
    real(real64) :: max_diff

    write(*,*) command_argument_count()

    if (command_argument_count() .ne. 3) then
        write(*,*) 'please provide sensor name, footprint size and band number'
        stop
    endif

    call get_command_argument(1,cmd_str)
    sensor_name = trim(cmd_str(1:10))

    call get_command_argument(2,cmd_str)
    read(cmd_str,*) footprint_size_int

    call get_command_argument(3,cmd_str)
    read(cmd_str,*) ifreq

    write(*,*) 'Computing Resampling Weight for ',sensor_name
    write(*,*) 'Footprint Size = ',footprint_size_int,' km'
    write(*,*) 'Iband = ',ifreq

    if (ifreq .ne. 5) then
        write(*,*) 'Only band 5 is implemented'
        stop
    !endif

    select case(ifreq)
    case(0)
        max_num_pts_per_source_footprint = 32000
        !7 GHz footrpints are saved with a larger threshold = 0.0001 not 0.00001
    case(1)
        max_num_pts_per_source_footprint = 52000   
    case(2)
        max_num_pts_per_source_footprint = 12000
    case(3)
        max_num_pts_per_source_footprint = 13000
    case(4)
        max_num_pts_per_source_footprint = 5000
    case(5)
        max_num_pts_per_source_footprint = 5000
    end select

    call allocate_full_grid(resampled_gain_full)
    call allocate_full_grid(target_gain_full)
    call allocate_compact_grid(target_gain,300000)
    call allocate_compact_grid(source_gain,max_num_pts_per_source_footprint)

    ! do nobs=0,nobs_max
    !     call allocate_compact_grid(gain_array(nobs),max_num_pts_per_source_footprint)
    ! enddo

    call allocate_CompactGrid_Array_FAFO(gain_array,max_num_pts_per_source_footprint)

    distance_threshold = 2.0*footprint_size_int
    if (footprint_size_int .eq. 70) then
        print *,"Reducing distance threshold for 70 km footprints"
        distance_threshold = 1.7*footprint_size_int
    
    call set_simulation_params(trim(sensor_name),sensor_data)

    write(list_file,9077) trim(sensor_name), footprint_size_int,ifreq
    9077 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a, &
                '/resample_weights_v3/stats/summary_stats_circular_',i2.2,'_band_',i2.2,'.txt')
    open(8,file=list_file, status='replace')

    write(target_location_file,9013) trim(sensor_name)
    9013 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/target_gains_v2/locs/find_target_gain_circular_30km.loc')

    write(source_location_file,9023) trim(sensor_name),trim(sensor_name),ifreq
    9023 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/source_gains_v2/locs/find_source_gain_',a,'_band_',i2.2,'_locs.txt')
    
    call read_location_text_file(target_location_file,target_locs,2)
    call read_location_text_file(source_location_file,source_locs,3)

    allocate(g_everywhere(source_locs%num_footprints,source_locs%num_footprints))

    ! read in precomputed g everywhere
    write(g_precompute_file,9033) trim(sensor_name),trim(sensor_name),ifreq
    9033 format('/mnt/ops1p-ren/l/ACCESS/resampling/',a,'/source_gains_v3/g_everywhere/g_everywhere_',a,'_band_',i2.2,'.dat')
    open(unit=3, &
            file=g_precompute_file, &
            action='read',&
            form='unformatted', &
            access='stream')

    g_everywhere = 0.0D0
    num_non_zero = 0
    do i = 1,99999999
        read(3,end = 4444) source_index1,source_index2,gtemp
        if (g_everywhere(source_index1,source_index2) .ne. 0.0) then
            write(*,*) 'repeat',source_index1,source_index2
            exit
        endif
        g_everywhere(source_index1,source_index2) = gtemp
        num_non_zero = num_non_zero + 1
    enddo

    4444 continue
    write(*,*) 'Read in ',num_non_zero,' non zeros elements of G from'
    write(*,*) g_precompute_file
    num_passed = 0
    num_fixed = 0
    !matrix is symmetric by definition
    do iobs=1,source_locs%num_footprints
        do jobs=1,iobs
            if ((g_everywhere(jobs,iobs)-g_everywhere(iobs,jobs)) .gt. 0.00000001) then
                write(*,*) 'Symmetry Error in G'
                write(*,*) iobs,jobs
                write(*,*) g_everywhere(jobs,iobs),g_everywhere(iobs,jobs),g_everywhere(jobs,iobs)-g_everywhere(iobs,jobs)
                temp = 0.5*(g_everywhere(jobs,iobs)+g_everywhere(iobs,jobs))
                g_everywhere(jobs,iobs) = temp
                g_everywhere(iobs,jobs) = temp
                num_fixed = num_fixed + 1
            else
                num_passed=num_passed + 1
            endif
        enddo
    enddo

    write(*,*) 'num G passed symmetry check',num_passed
    write(*,*) 'num G failed symmetry check',num_fixed
    

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

        write(6,341)
        341 format('allocating arrays for svd calculations')
        allocate(g(nobs,nobs),g_inv(nobs,nobs),testg(nobs,nobs), u(nobs), v(nobs), xwork(nobs), a(nobs))
        


        smooth_fac=1.0d-05    
        if (footprint_size_int .eq. 70) then
            smooth_fac=1.0d-04  !try this to stabilize 70 km footprints
        endif
        
        ! read in the target gain
        write(6,335)
        335 format('reading target gain pattern')
        
        write(6,337)trim(sensor_name),target_scan_to_do, target_fov_to_do
        337 format('reading ',a,' target gain pattern, scan: ',i4,' fov: ',i4)

        select case(trim(sensor_name)) 
            case('AMSR2')
                        iobs=0
                        !call read_gobs_amsr2_v2(1,ifreq,target_scan_to_do,target_fov_to_do,gain_array(0))
                        call read_tar_amsr2_v3(footprint_size_int,ifreq,target_scan_to_do,target_fov_to_do,target_gain)
                    case default
                        print *,'Sensor ',sensor_name,' not implemented'
                        stop
        end select

        call convert_to_full_grid(target_gain,target_gain_full,xsum)
          
        write(6,336)nobs,trim(sensor_name)
        336 format('using ',i4,1x,a,' source gain patterns')
        iobs = 0
        num_source_loaded = 0
        num_source_reused = 0
        do i = 1,source_locs%num_footprints         !loop over the entire source grid
            if (source_locs%use(i) == 0) cycle      !but skip the locations that are not used

            iscan = source_locs%iscan(i)            !indicec of the fovs that are used
            icel = source_locs%ifov(i)

            ! check to see if footprint is already located in the array
            call find_CompactGrid_Array_FAFO(gain_array,iscan,icel,found_source,source_gain)

            if (.not. found_source) then
                select case(trim(sensor_name)) 
                case('AMSR2')
                    !read the source gain pattern 
                    call read_gobs_amsr2_v3(0,ifreq,iscan,icel,source_gain)
                case default
                    print *,'Sensor ',sensor_name,' not implemented'
                    stop
                end select

                call insert_CompactGrid_Array_FAFO(gain_array,source_gain,iscan,icel)
                
                num_source_loaded = num_source_loaded + 1
            else
                num_source_reused = num_source_reused + 1
            endif

            call compute_overlap_full_compact(target_gain_full,source_gain,xsum)
            iobs=iobs+1
            v(iobs)=xsum
            iscansv(iobs)=iscan   
            icelsv(iobs)=icel
            iobs_numsv(iobs) = i
        enddo

        write(*,*) 'Loaded ',num_source_loaded,' Reused ',num_source_reused
        nobs=iobs

        
        write(6,742)
        742 format('selecting relevant elements of g from g_everywhere')

        do iobs=1,nobs
            do jobs=1,nobs
                g(iobs,jobs) = g_everywhere(iobs_numsv(iobs),iobs_numsv(jobs))
            enddo
        enddo


        write(6,743)
        743 format('calculating target/source overlap integrals')

        num_calculated = 0

        !matrix is symmetric by definition
        do iobs=1,nobs
            do jobs=1,iobs
                if ((g(jobs,iobs)-g(iobs,jobs)) .gt. 0.000000001) then
                    write(*,*) 'Symmetry Error in G'
                    write(*,*) g(jobs,iobs),g(iobs,jobs),g(jobs,iobs)-g(iobs,jobs)
                    stop
                endif
            enddo
        enddo

        u=1
        ! add in the smoothing factor along the diagonal
        do iobs=1,nobs
            g(iobs,iobs)=g(iobs,iobs) + smooth_fac
        enddo
    
        ! invert
        call svdinv(g,nobs, g_inv,ierror)
        if(ierror.ne.0) then
            write(*,*) 'inversion did not converge, location skipped'
        endif

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

        !summary properties of a
        max_a = 0.0
        min_a = 0.0
        num_a_nonzero = 0
        num_a_positive = 0
        num_a_negative = 0

        do iobs=1,nobs
            if (a(iobs) .gt. max_a) then
                max_a = a(iobs)
            endif
            if (a(iobs) .lt. min_a) then
                min_a = a(iobs)
            endif
            if (a(iobs) .gt. 0.0) then
                num_a_nonzero = num_a_nonzero + 1
                num_a_positive = num_a_positive + 1
            endif
            if (a(iobs) .lt. 0.0) then
                num_a_nonzero = num_a_nonzero + 1
                num_a_negative = num_a_negative + 1
            endif
        enddo

        write(*,*) 'max weight= ',max_a
        write(*,*) 'min weight= ',min_a
        write(*,*) 'num non zero = ',num_a_nonzero
        write(*,*) 'num positive = ',num_a_positive
        write(*,*) 'num negative = ',num_a_negative

        psum(2)=sqrt(psum(2))
        psum(3:6)=psum(3:6)/psum(1)

        write(6,801)
        801 format('Calculating resampled gain pattern')

        call zero_full_grid(resampled_gain_full)
        resampled_gain_full%xlat0 = 0.0
        resampled_gain_full%xlon0 = 0.0
        xsum = 0.0D0
        
        do iobs=1,nobs
            iscan = iscansv(iobs)      !indicec of the fovs that are used
            icel = icelsv(iobs)
            call find_CompactGrid_Array_FAFO(gain_array,iscan,icel,found_source,source_gain)
            if (.not. found_source) then
                select case(trim(sensor_name)) 
                case('AMSR2')
                    !read the source gain pattern 
                    call read_gobs_amsr2_v3(0,ifreq,iscan,icel,source_gain)
                case default
                    print *,'Sensor ',sensor_name,' not implemented'
                    stop
                end select
            endif
            do ipixel=1,source_gain%numpts
                ilat       = source_gain%ilat(ipixel)
                ilon       = source_gain%ilon(ipixel)
                delta_gain = source_gain%gain(ipixel)*a(iobs)
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
        
        open(unit=17,file=filename1, status='replace',action='write',form='unformatted',access='stream')
        num_valid_pix = 0
        call write_gain_full_grid(resampled_gain_full,17,1.0d-8,num_valid_pix)
        close(17)


        tot_a_summary(target_fov_to_do,target_scan_to_do) = psum(1)
        tot_sqr_a_summary(target_fov_to_do,target_scan_to_do) = psum(2)

        num_valid_pix = 0

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