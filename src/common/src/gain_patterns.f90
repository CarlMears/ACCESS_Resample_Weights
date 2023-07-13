! This module contains types and subroutines for antenna gain patterns on the Earth's surface.
!
! Gain patterns are stored in 2 ways:  
!   FullGrid (:,:) arrays of real64s
!   CompactGrid - lists of ilon,ilat, and gain for non-zero elements in the full grid.
!


module gain_pattern
  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int8, int16, ERROR_UNIT
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite

  implicit none 

  private
  public FullGrid, CompactGrid
  public allocate_full_grid,allocate_compact_grid

  public CompactGrid_Array_FAFO
  public allocate_CompactGrid_Array_FAFO,insert_CompactGrid_Array_FAFO,find_CompactGrid_Array_FAFO

  public find_mean_grid_location
  public convert_to_full_grid,convert_to_compact_grid
  public write_gain_compact_grid,write_gain_full_grid
  public read_gain_compact_grid
  public normalize_compact_grid, normalize_full_grid
  public zero_compact_grid, zero_full_grid
  public compute_overlap_full_compact
  public nlon,nlat


  ! These parameters define the gain grid -- in the future, could be variables.
  ! This definition is taken from Franks resampling calculations e.g. find_gains.f90
  ! The grid is large enough for all current satellites

  integer(int32),parameter  :: nlon = 2401
  integer(int32),parameter  :: nlat = 1601
  real(real64),parameter    :: dlon = 0.01D0
  real(real64),parameter    :: dlat = 0.01D0
  real(real64),parameter    :: lat_offset = ((nlat-1)/2)*dlat
  real(real64),parameter    :: lon_offset = ((nlon-1)/2)*dlon
  integer(int32),parameter  :: maxpts_default = 400000


  integer(int32),parameter :: num_grids_in_array = 3000

  type FullGrid
    real(real64), allocatable :: gain(:,:) !gain values for full grid
    real(real64)              :: xlat0    !latitude of the grid center
    real(real64)              :: xlon0    !longitude of the grid center
  end type FullGrid

  type CompactGrid
    integer(int32)                :: maxpts
    integer(int32)                :: numpts      !number of points with valid gains
    integer(int16),allocatable    :: ilon(:)     !longitude location in full grid
    integer(int16),allocatable    :: ilat(:)     !latitude location in full grid
    real(real64),allocatable      :: gain(:)     !gain at ilon,ilat
  end type CompactGrid

  type CompactGrid_Array_FAFO
    real(real64),dimension(num_grids_in_array)                  :: last_time_accessed
    integer(int16),dimension(num_grids_in_array)                :: iscan
    integer(int16),dimension(num_grids_in_array)                :: icel
    logical,dimension(num_grids_in_array)                       :: valid
    type(CompactGrid),dimension(num_grids_in_array)             :: grids
  endtype
                

contains

    subroutine allocate_full_grid(g_full)

        type(FullGrid),intent(inout)  :: g_full

        allocate(g_full%gain(nlon,nlat))
        g_full%gain = 0.0D0

    end subroutine allocate_full_grid

    subroutine deallocate_full_grid(g_full)

        type(FullGrid),intent(inout)  :: g_full

        deallocate(g_full%gain)

    end subroutine deallocate_full_grid

    subroutine allocate_compact_grid(g_compact,maxpts_in)

        type(CompactGrid),intent(inout) :: g_compact        !grid to be allocated
        integer(int32), intent(in), optional :: maxpts_in   !size of allocated arrays, optional

        !integer(int32) :: maxpts

        if (present(maxpts_in)) then
            g_compact%maxpts = maxpts_in
        else
            g_compact%maxpts = maxpts_default
        endif
            

        allocate(g_compact%ilon(g_compact%maxpts))
        allocate(g_compact%ilat(g_compact%maxpts))
        allocate(g_compact%gain(g_compact%maxpts))

        g_compact%ilon = -9999
        g_compact%ilat = -9999
        g_compact%gain = 0.0

    end subroutine allocate_compact_grid

    subroutine deallocate_compact_grid(g_compact)

        type(CompactGrid),intent(inout) :: g_compact        !grid to be allocated
        
            

        deallocate(g_compact%ilon)
        deallocate(g_compact%ilat)
        deallocate(g_compact%gain)


    end subroutine deallocate_compact_grid

    ! the next 3 routines maintain an array of compact grids.
    ! the array is least_recently_accessed first_out
    ! cpu_time() appears to return a value that is monotonic with increasing
    ! real world time, so it works here as a access time stamp

    subroutine allocate_CompactGrid_Array_FAFO(grid_array,maxpts_in)
        type(CompactGrid_Array_FAFO),intent(inout)  :: grid_array
        integer(int32),intent(in)                   :: maxpts_in

        integer(int32)      :: i
        real(real64)        :: cpu_time_now


        do i = 1,num_grids_in_array

            call allocate_compact_grid(grid_array%grids(i),maxpts_in)
            call cpu_time(cpu_time_now)
            grid_array%last_time_accessed(i) = cpu_time_now
            grid_array%iscan(i) = -99
            grid_array%icel(i) = -99
            grid_array%valid(i) = .false.
        enddo

    end subroutine allocate_CompactGrid_Array_FAFO

    subroutine insert_CompactGrid_Array_FAFO(grid_array,grid_to_insert,iscan,icel)
        type(CompactGrid_Array_FAFO),intent(inout)  :: grid_array
        type(CompactGrid),intent(in)                :: grid_to_insert
        integer(int16)                              :: iscan
        integer(int16)                              :: icel


        integer(int32)      :: i
        real(real64)        :: cpu_time_now
        !real(real64)        :: cpu_time_old                 

        !find the least recently access element
        i = minloc(grid_array%last_time_accessed,1)
        
        !cpu_time_old  = grid_array%last_time_accessed(i)

        grid_array%grids(i) = grid_to_insert
        grid_array%iscan(i) = iscan
        grid_array%icel(i)  = icel
        grid_array%valid(i) = .True.
        !set the time stamp
        call cpu_time(cpu_time_now)
        grid_array%last_time_accessed(i) = cpu_time_now


        !write(*,*) 'Inserted source gain for ',iscan,icel,'Last access=',cpu_time_old,'Now= ',cpu_time_now

    end subroutine insert_CompactGrid_Array_FAFO

    subroutine find_CompactGrid_Array_FAFO(grid_array,iscan,icel,found,grid_found)
        type(CompactGrid_Array_FAFO),intent(inout)  :: grid_array
        integer(int16),intent(in)                :: iscan
        integer(int16),intent(in)                :: icel
        logical,intent(out)                      :: found
        type(CompactGrid),intent(inout)          :: grid_found

        integer(int32)      :: i,j
        real(real64)        :: cpu_time_now


        found = .false.

        do i = 1,num_grids_in_array
            if (grid_array%valid(i)) then
                if ((grid_array%iscan(i) .eq. iscan) .and. &
                    (grid_array%icel(i) .eq. icel)) then
                        found = .true.
                        grid_found = grid_array%grids(i)
                        !set the time stamp to indicate access
                        call cpu_time(cpu_time_now)
                        grid_array%last_time_accessed(i) = cpu_time_now
                endif
            endif
        enddo

    end subroutine find_CompactGrid_Array_FAFO
 
    subroutine convert_to_full_grid(g_compact,g,sumgain)

        type(CompactGrid),intent(in)                  :: g_compact
        type(FullGrid),intent(inout)                  :: g 
        real(real64),intent(out),optional             :: sumgain

        integer(int32)                                :: n
        real(real64)                                  :: sum

        sum = 0.0D0
        g%gain = 0.0D0
        do n=1,g_compact%numpts
            g%gain(g_compact%ilon(n),g_compact%ilat(n)) = g_compact%gain(n)
            sum = sum + g_compact%gain(n)
        enddo
        if (present(sumgain)) then
            sumgain = sum
        endif

    end subroutine convert_to_full_grid

    subroutine convert_to_compact_grid(g,gain_threshold,g_compact)
        
        type(FullGrid),intent(in)                      :: g 
        real(real64),intent(in)                        :: gain_threshold
        type(CompactGrid),intent(inout)                :: g_compact

        integer(int32)                                :: n 
        integer(int16)                                :: ilon,ilat

        n = 0
        g_compact%ilat = -999
        g_compact%ilon = -999
        g_compact%gain = 0.0D0
        do ilat = 1,nlat 
            do ilon = 1,nlon
                if ((ieee_is_finite(g%gain(ilon,ilat))) .and. &
                    (abs(g%gain(ilon,ilat)) > gain_threshold)) then
                    n = n+1
                    g_compact%ilon(n)  = ilon
                    g_compact%ilat(n)  = ilat
                    g_compact%gain(n)  = g%gain(ilon,ilat)
                endif
            enddo
        enddo
        g_compact%numpts = n 
    
    end subroutine convert_to_compact_grid
        
    subroutine find_mean_grid_location(g,mean_lat,mean_lon,total_wt)

        type(FullGrid),intent(in)                     :: g      !gain pattern in lat/lon grid
        real(real64),intent(out)                      :: mean_lat,mean_lon,total_wt

        real(real64),dimension(0:2)                   :: psum
        integer(int32)                                :: ilat,ilon 

        real(real64)                                  :: xlat_grid,xlon_grid

        psum = 0.0D0

        do ilat=1,nlat
            do ilon=1,nlon
                xlat_grid = dlat*(ilat-1) + g%xlat0 - lat_offset
                xlon_grid = dlon*(ilon-1) + g%xlon0 - lon_offset
                psum(0) = psum(0) + g%gain(ilon,ilat)
                psum(1) = psum(1) + xlon_grid*g%gain(ilon,ilat)
                psum(2) = psum(2) + xlat_grid*g%gain(ilon,ilat)

            enddo  !ilon
        enddo  !ilat

        mean_lon = psum(1)/psum(0)
        mean_lat = psum(2)/psum(0)
        total_wt = psum(0)

    end subroutine find_mean_grid_location

    subroutine compute_overlap_full_compact(gain_full,gain_compact,xsum)

        type(FullGrid), intent(in)  :: gain_full
        type(CompactGrid), intent(in) :: gain_compact

        real(real64), intent(out)   :: xsum

        integer(int32) :: ipixel
        integer(int16) :: ilat,ilon

        xsum = 0.0D0

        do ipixel=1,gain_compact%numpts  !loop over the non-zero pixels in gain compact
            ilat = gain_compact%ilat(ipixel)
            ilon = gain_compact%ilon(ipixel)
            xsum = xsum + gain_full%gain(ilon,ilat)*gain_compact%gain(ipixel)
        enddo

    end subroutine compute_overlap_full_compact

    subroutine zero_full_grid(g)

        type(FullGrid),intent(inout)                     :: g      !gain pattern in lat/lon grid

        g%gain = 0.0D0

    end subroutine zero_full_grid

    subroutine normalize_full_grid(g)

        type(FullGrid),intent(inout)                     :: g      !gain pattern in lat/lon grid
        
        real(real64)                                  :: sum
        integer(int32)                                :: ilat,ilon 

        sum = 0.0D0

        do ilat=1,nlat
            do ilon=1,nlon
                sum = sum + g%gain(ilon,ilat)
            enddo  !ilon
        enddo  !ilat

        g%gain = g%gain/sum

    end subroutine normalize_full_grid

    subroutine zero_compact_grid(g_compact)

        type(CompactGrid),intent(inout)         :: g_compact      !gain pattern in compact grid

        g_compact%numpts = 0
        g_compact%ilat = -999
        g_compact%ilon = -999
        g_compact%gain = 0.0D0

    end subroutine zero_compact_grid

    subroutine normalize_compact_grid(g_compact)

        type(CompactGrid),intent(inout)         :: g_compact      !gain pattern in compact grid
        
        integer(int32)                          :: n
        real(real64)                            :: sum

        sum = 0.0D0

        do n = 1,g_compact%numpts
            sum = sum + g_compact%gain(n)
        enddo

        g_compact%gain = g_compact%gain/sum

    end subroutine normalize_compact_grid

    subroutine write_gain_full_grid(g,output_unit,gain_threshold,num_written)

        type(FullGrid),intent(in)      :: g
        integer,intent(in)             :: output_unit
        real(real64),intent(in)        :: gain_threshold
        integer(int32),intent(out)     :: num_written
        
        type(CompactGrid)              :: g_compact

        call allocate_compact_grid(g_compact)
        call convert_to_compact_grid(g,gain_threshold,g_compact)
        call write_gain_compact_grid(g_compact,output_unit,num_written)
        call deallocate_compact_grid(g_compact)

    end subroutine write_gain_full_grid

    subroutine write_gain_compact_grid(g_compact,output_unit,num_written)

        type(CompactGrid),intent(in)               :: g_compact
        integer,intent(in)                         :: output_unit
        integer(int32),intent(out)     :: num_written

        integer(int32)                             :: n

        do n = 1,g_compact%numpts
            write(output_unit)g_compact%ilat(n)-int(1,kind=int16) ,g_compact%ilon(n)-int(1,kind=int16), g_compact%gain(n)
            ! if (n .le. 20) then
            !     write(*,*)g_compact%ilat(n)-int(1,kind=int16) ,g_compact%ilon(n)-int(1,kind=int16), g_compact%gain(n)
            ! endif
        enddo
        num_written=g_compact%numpts

    end subroutine write_gain_compact_grid

    subroutine read_gain_compact_grid(input_unit,g_compact)

        integer,intent(in)                          :: input_unit
        type(CompactGrid),intent(inout)               :: g_compact

        integer(int32)                              :: n,i
        integer(int16)                              :: ilat,ilon
        real(real64)                                :: gain

        n = 0
        do i = 1, g_compact%maxpts
            read(input_unit,end = 99) ilat,ilon,gain
            n = n + 1
            g_compact%ilon(n) = ilon+int(1,kind=int16) ! +1 is to convert from python indexing
            g_compact%ilat(n) = ilat+int(1,kind=int16) 
            g_compact%gain(n) = gain
        enddo

        write(*,*) 'Error - end of file not reached!'
        write(*,*) 'i = ',i
        stop

        99 g_compact%numpts = n

    end subroutine read_gain_compact_grid
end module gain_pattern
