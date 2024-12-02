    
    program find_resampling_weights 

    use, intrinsic :: iso_fortran_env, only: int8, int16, int32, real32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite

    !use simulation_param, only: NBAND,FOVS_PER_SCAN,NFREQ
    !use simulation_param, only: set_simulation_params, SimulationParameterData_f
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

    use gain_pattern, only: CompactGrid_Array_FAFO
    use gain_pattern, only: allocate_CompactGrid_Array_FAFO, &
                            insert_CompactGrid_Array_FAFO,   &
                            find_CompactGrid_Array_FAFO

    use ssmis_gains, only: read_gain_obs_ssmis
    use ssmi_gains, only: read_gain_obs_ssmi
    
    implicit none
      
    character(120) source_location_file,output_file,output_file_text
                        
    integer(int16) iscan1,icel1,iscan2,icel2      ! various incideces
    integer(int32) num_to_do,num_computed,num_found,num_loaded,num_written,source1,source2
    integer(int32) output_unit
    integer(int32) ifreq
    integer(int32) max_num_pts_per_source_footprint
    integer(int32) iobs,jobs
    integer(int32) ksat
  
    real(real64) xsum

    real(real64), allocatable :: g_everywhere(:,:)  !array to store elements of g to prevent recalculation
    integer(int32), allocatable :: do_g_everywhere(:,:) 
    real(real64) :: g_mean,g_sum

    type(CompactGrid) ::  source1_gain
    type(CompactGrid) ::  source2_gain
    type(CompactGrid_Array_FAFO)  :: source_gain_array
    logical :: found_source

    type(FullGrid)     :: source1_gain_full
    type(FullGrid)     :: source2_gain_full
    type(FullGrid)     :: product_gain_full

    character(len=10) sensor_name
    !type(SimulationParameterData_f) :: sensor_data
    type(FootprintLocations)        :: source_locs

    integer(int32) :: distance_error
    real(real64)   :: source1_lat,source2_lat,source1_lon,source2_lon
    real(real32)   :: distance,distance_threshold,gain_threshold

    integer(4),dimension(6), parameter :: footprint_major_axis_AMSR2 = &
                    (/62, 42, 22, 26, 12, 5/)

    integer(4),dimension(6), parameter :: footprint_major_axis_SSMIS = &
                    (/0, 0, 73, 50, 41, 14/)

    integer(4),dimension(6), parameter :: footprint_major_axis_SSMI = &
                    (/0, 0, 69, 50, 37, 15/)

    real(real64) :: cpu_time_before,cpu_time_after
    real(real64) :: cpu_time_before2,cpu_time_after2
    real(real64) :: cpu_time_computing,cpu_time_loading


    max_num_pts_per_source_footprint = 80000
  
    call allocate_full_grid(source1_gain_full)
    call allocate_full_grid(source2_gain_full)
    call allocate_full_grid(product_gain_full)
    call allocate_compact_grid(source1_gain,max_num_pts_per_source_footprint)

    call allocate_CompactGrid_Array_FAFO(source_gain_array,max_num_pts_per_source_footprint)

    sensor_name = 'AMSR2'
    ksat = 13

    ifreq = 5
    if (sensor_name .eq. 'AMSR2') distance_threshold = 3.0*footprint_major_axis_AMSR2(ifreq+1)
    if (sensor_name .eq. 'SSMIS') distance_threshold = 3.0*footprint_major_axis_SSMIS(ifreq+1)
    if (sensor_name .eq. 'SSMI') distance_threshold = 3.0*footprint_major_axis_SSMI(ifreq+1)
    
    !call set_simulation_params(trim(sensor_name),sensor_data)

    if (sensor_name .eq. 'AMSR2') then
        write(source_location_file,9023) trim(sensor_name),trim(sensor_name),ifreq
        9023 format('/mnt/l/access/resampling/',a,'/source_gains_v3/locs/find_source_gain_',a,'_band_',i2.2,'_locs.txt')
        write(output_file,9033) trim(sensor_name),trim(sensor_name),ifreq
        9033 format('/mnt/l/access/resampling/',a,'/source_gains_v3/g_everywhere/g_everywhere_',a,'_band_',i2.2,'.v2.dat')
        write(output_file_text,9034) trim(sensor_name),trim(sensor_name),ifreq
        9034 format('/mnt/l/access/resampling/',a,'/source_gains_v3/g_everywhere/g_everywhere_',a,'_band_',i2.2,'.v2.txt')
    endif

    if (sensor_name .eq. 'SSMIS') then
        write(source_location_file,9123) trim(sensor_name),ksat,ksat,ifreq
        9123 format('/mnt/l/access/resampling/',a,'/f',i2.2,'/source_gains/locs/find_source_gain_f' &
                     ,i2.2,'_band_',i2.2,'_locs.txt')
        write(output_file,9133) trim(sensor_name),ksat,trim(sensor_name),ifreq
        9133 format('/mnt/l/access/resampling/',a,'/f',i2.2,'/source_gains/g_everywhere/g_everywhere_' &
                     ,a,'_band_',i2.2,'.dat')
        write(output_file_text,9134) trim(sensor_name),ksat,trim(sensor_name),ifreq
        9134 format('/mnt/l/access/resampling/',a,'/f',i2.2,'/source_gains/g_everywhere/g_everywhere_' &
        ,a,'_band_',i2.2,'.txt')
    endif

    if (sensor_name .eq. 'SSMI') then
        write(source_location_file,9112) trim(sensor_name),ksat,ksat,ifreq-2
        9112 format('/mnt/l/access/resampling/',a,'/f',i2.2,'/source_gains/locs/find_source_gain_f' &
                     ,i2.2,'_freq_',i2.2,'_locs.txt')
        write(output_file,9113) trim(sensor_name),ksat,trim(sensor_name),ifreq-2
        9113 format('/mnt/l/access/resampling/',a,'/f',i2.2,'/source_gains/g_everywhere/g_everywhere_' &
                     ,a,'_band_',i2.2,'.dat')
        write(output_file_text,9114) trim(sensor_name),ksat,trim(sensor_name),ifreq-2
        9114 format('/mnt/l/access/resampling/',a,'/f',i2.2,'/source_gains/g_everywhere/g_everywhere_' &
        ,a,'_band_',i2.2,'.txt')
    endif


    call read_location_text_file(source_location_file,source_locs,3)

    allocate(do_g_everywhere(source_locs%num_footprints,source_locs%num_footprints))
    do_g_everywhere = 0

    allocate(g_everywhere(source_locs%num_footprints,source_locs%num_footprints))
    g_everywhere = 0.0D0

    num_to_do = 0
    write(*,*) 'Finding source footprints closer than ',distance_threshold
    write(*,*) source_locs%num_footprints
    do source1 = 1,source_locs%num_footprints
        if (mod(source1,1000) .eq. 0) write(*,*) 'source1=',source1,' of ',source_locs%num_footprints
        source1_lat = source_locs%lat(source1)
        source1_lon = source_locs%lon(source1)

        do source2=source1,source_locs%num_footprints
            source2_lat = source_locs%lat(source2)
            source2_lon = source_locs%lon(source2)
            
            distance = distance_between_lon_lats(&
                real(source1_lon,kind=real32), &
                real(source1_lat,kind=real32), &
                real(source2_lon,kind=real32), &
                real(source2_lat,kind=real32),distance_error)
            
            if (distance < distance_threshold) then
                !compute G_ij
                do_g_everywhere(source1,source2) = 1
                num_to_do = num_to_do + 1
            endif
        enddo
    enddo

    print 899, num_to_do,source_locs%num_footprints*source_locs%num_footprints
    899 format('Need to calculate ',i9,' overlaps out of a total of ',i10,' array_elements')




    num_computed = 0
    do source1 = 1,source_locs%num_footprints

        source1_lat = source_locs%lat(source1)
        source1_lon = source_locs%lon(source1)

        iscan1 = source_locs%iscan(source1)
        icel1 = source_locs%ifov(source1)


        ! check to see if footprint is already located in the array
        call find_CompactGrid_Array_FAFO(source_gain_array,iscan1,icel1,found_source,source1_gain)

        if (.not. found_source) then
            select case(trim(sensor_name)) 
            case('AMSR2')
                !read the source gain pattern 
                call read_gobs_amsr2_v3(0,ifreq,iscan1,icel1,source1_gain)
            case('SSMIS')
                call read_gain_obs_ssmis(ksat,ifreq,iscan1,icel1,source1_gain)
            case('SSMI')
                call read_gain_obs_ssmi(ksat,ifreq-2,iscan1,icel1,source1_gain)
            case default
                print *,'Sensor ',sensor_name,' not implemented for read_gain_obs 1'
                stop
            end select
            call insert_CompactGrid_Array_FAFO(source_gain_array,source1_gain,iscan1,icel1)
        endif

        call convert_to_full_grid(source1_gain,source1_gain_full)

        !write(*,*) iscan1,icel1,'need to do',sum(do_g_everywhere(source1,:))
        call cpu_time(cpu_time_before)
        num_found=0
        num_loaded=0
        cpu_time_loading = 0.0
        cpu_time_computing = 0.0
        

        do source2=source1,source_locs%num_footprints
            if (do_g_everywhere(source1,source2) > 0) then
                iscan2 = source_locs%iscan(source2)
                icel2  = source_locs%ifov(source2)
                call cpu_time(cpu_time_before2)
                call find_CompactGrid_Array_FAFO(source_gain_array,iscan2,icel2,found_source,source2_gain)

                if (.not. found_source) then
                    select case(trim(sensor_name)) 
                    case('AMSR2')
                        !read the source gain pattern 
                        call read_gobs_amsr2_v3(0,ifreq,iscan2,icel2,source2_gain)
                    case('SSMIS')
                        call read_gain_obs_ssmis(ksat,ifreq,iscan2,icel2,source2_gain)
                    case('SSMI')
                        call read_gain_obs_ssmi(ksat,ifreq-2,iscan2,icel2,source2_gain)
                    case default
                        print *,'Sensor ',sensor_name,' not implemented for read gain obs 2'
                        stop
                    end select
                    call insert_CompactGrid_Array_FAFO(source_gain_array,source2_gain,iscan2,icel2)

                    num_loaded = num_loaded+1
                else
                    num_found = num_found+1
                endif
                call cpu_time(cpu_time_after2)
                cpu_time_loading = cpu_time_loading + cpu_time_after2-cpu_time_before2

                call cpu_time(cpu_time_before2)
                call compute_overlap_full_compact(source1_gain_full,source2_gain,xsum)
                call cpu_time(cpu_time_after2)
                cpu_time_computing = cpu_time_computing + cpu_time_after2-cpu_time_before2

                g_everywhere(source1,source2) = xsum
                g_everywhere(source2,source1) = xsum
                num_computed = num_computed + 1

                !print *,iscan1,icel1,iscan2,icel2,xsum
                !print *
            endif
        enddo
        call cpu_time(cpu_time_after)
        ! write(*,*) iscan1,icel1,'num found= ',num_found,' num_loaded = ',num_loaded, &
                ! cpu_time_after-cpu_time_before,cpu_time_loading,cpu_time_computing
        write(*,567)  iscan1,icel1,num_found,num_loaded,cpu_time_loading,cpu_time_computing
    567 format('iscan=',i3,' icel=',i4,' num found=',i4,' num_loaded=',i4,' cpu loading',f8.4,' cpu calc=',f8.4)
        !matrix is symmetric by definition
        do iobs=1,source_locs%num_footprints
            do jobs=1,iobs
                if ((g_everywhere(jobs,iobs)-g_everywhere(iobs,jobs)) .gt. 0.00000001) then
                    write(*,*) 'Symmetry Error in G'
                    write(*,*) iobs,jobs
                    write(*,*) g_everywhere(jobs,iobs),g_everywhere(iobs,jobs),g_everywhere(jobs,iobs)-g_everywhere(iobs,jobs)
                    stop
                endif
            enddo
        enddo
        !write(*,*) 'G passed symmetry check'
    enddo

    ! calculate the mean of the diagonal elements to set the citeria for writing the elements
    g_sum = 0.0D0

    do source1 = 1,source_locs%num_footprints
       g_sum = g_sum + g_everywhere(source1,source1)
    enddo
    
    g_mean = g_sum/source_locs%num_footprints  !mean of diagonal elements


    gain_threshold = 0.0001
    output_unit = 3

    open(unit=output_unit, &
         file=output_file, &
         !status='new',&
         action='write',&
         form='unformatted', &
         access='stream')

    open(unit=4, &
         file=output_file_text, &
         !status='new',&
         action='write',&
         form='formatted', &
         access='stream')
            
    num_written = 0
    do source1 = 1,source_locs%num_footprints
        do  source2 = 1,source_locs%num_footprints
            if (g_everywhere(source1,source2)/g_mean .gt. gain_threshold) then
                ! write out element
                write(output_unit)source1,source2,g_everywhere(source1,source2)
                num_written = num_written+1
                if (num_written .le. 1000) then
                    write(4,1111) source1,source2,g_everywhere(source1,source2)
                    1111 format(i4,",",i4,",",f11.8)
                endif
            endif
        enddo
    enddo

    close(output_unit)
    close(4)
    write(*,*) 'Wrote ',num_written,' elements to '
    write(*,*) output_file
    
    stop 'norm end'                                   
    end