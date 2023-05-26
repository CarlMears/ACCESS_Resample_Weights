module resampling_types

  use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int8, ERROR_UNIT
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

  private
  public FootprintLocations, read_location_text_file
  type FootprintLocations
     integer(int32) :: num_footprints
     integer(int32), allocatable, dimension(:) :: iscan
     integer(int32), allocatable, dimension(:) :: ifov 
     real(real64), dimension(:), allocatable :: lat      
     real(real64), dimension(:), allocatable :: lon
     integer(int8), dimension(:), allocatable :: use
     real(real64), dimension(:), allocatable ::weight
  end type FootprintLocations

contains

  subroutine read_location_text_file(filename,locs,num_skip)

        character(len = *),intent(in)           :: filename
        type(FootprintLocations), intent(out) :: locs
        integer(int32),intent(in)               :: num_skip
        integer :: ioerr, lu, i,num_lines
        character(len=80) :: iomsg

        integer(int32)    :: kscan,kfov,iband
        real(real64)     :: lat, lon

        open(newunit=lu, file=filename, status='old', access='stream', &
            form='formatted', action='read', &
            iostat=ioerr, iomsg=iomsg)

        ! read thru once to count lines
        do i = 1,1000000
            if (num_skip == 3) then
              read(lu,'(1x,i4,1x,i4.4,1x,i4.4,1x,f12.5,1x,f12.5)', end = 99) iband,kscan, kfov, lat, lon
            else
              read(lu,'(1x,i4.4,1x,i4.4,1x,f12.5,1x,f12.5)', end = 99) kscan, kfov, lat, lon
            endif
        enddo

        99 num_lines = i-1

        allocate(locs%iscan(num_lines))
        allocate(locs%ifov(num_lines))
        allocate(locs%lat(num_lines))
        allocate(locs%lon(num_lines))
        allocate(locs%use(num_lines))
        allocate(locs%weight(num_lines))
        locs%num_footprints = num_lines

        locs%iscan = 0
        locs%ifov = 0
        locs%lat = -999.0D0
        locs%lon = -999.0D0
        locs%use = 0
        locs%weight = 0.0D0

        rewind(lu)

        do i = 1,num_lines
            if (num_skip == 3) then
              read(lu,'(1x,i4.4,1x,i4.4,1x,i4.4,1x,f12.5,1x,f12.5)', end = 99) iband, kscan, kfov, lat, lon
            else
              read(lu,'(1x,i4.4,1x,i4.4,1x,f12.5,1x,f12.5)', end = 99) kscan, kfov, lat, lon
            endif
            locs%iscan(i) = kscan
            locs%ifov(i) = kfov
            locs%lat(i) = lat 
            locs%lon(i) = lon 
        enddo

end subroutine read_location_text_file

end module resampling_types
