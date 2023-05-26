subroutine read_gobs_amsr2(icase,npixel_max,nobs_max,ifreq,iscan,icel,iobs, npixel,ilat2,ilon2,gobs)

    use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
    use gain_pattern, only: CompactGrid, read_gain_compact_grid
    implicit none

    integer(int32),intent(in) :: icase        ! ==1 if target case
    integer(int32),intent(in) :: npixel_max   ! max number of pixels >0 
    integer(int32),intent(in) :: nobs_max     ! max number obs to be included
    integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
    integer(int32),intent(in) :: iobs         ! index of the observation -- 0 means target
    integer(int32),intent(in) :: iscan        ! scan number
    integer(int32),intent(in) :: icel         ! fov number

    integer(int16), intent(out) :: ilat2(1:npixel_max)  !list of ilats
    integer(int16), intent(out) :: ilon2(1:npixel_max)  !list of ilons
    integer(int32), intent(out) :: npixel           !number of pixel for each gain pattern
    real(real64),intent(out)    :: gobs(1:npixel_max)   !gain for each >0 location
                                                                   !this holds the data for both the target and source footprints

    type(CompactGrid) :: g_compact
    character(80) filename
    integer(int16) ilat2x,ilon2x
    integer(int32) :: i,ipixel
    real(real64) :: gobsx
    real(real64) xsum

    integer :: ioerr
    character(len =80) :: iomsg

    integer :: beamwidth_km_int = 30

    if(icase.eq.1) then
        write(filename,9001)beamwidth_km_int,iscan,icel
        9001 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/target_gains/circular_',i2.2,'km/s',i2.2,'c',i4.4,'.dat')

    else
        write(filename,9002)ifreq,iscan,icel
        9002 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/source_gains/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
    endif

    open(unit=3,file=filename, status='old',action='read',form='unformatted', access='stream',iostat=ioerr,iomsg=iomsg)
                
    if (ioerr /= 0) then
        write(6,*) iomsg
        error stop 1
    endif

    call read_gain_compact_grid(3,g_compact)
     
    ipixel=0; xsum=0
    gobs = 0.0D0
    npixel = 0



    do i=1,1234567890
        read(3,end=100) ilat2x,ilon2x,gobsx
        ipixel=ipixel+1
        if(ipixel.gt.npixel_max) stop 'pgm stopped, ipixel oob'
        ilat2(ipixel) = ilat2x
        ilon2(ipixel) = ilon2x
        gobs( ipixel) = gobsx 
        npixel=ipixel
        xsum=xsum+gobs( ipixel)
    enddo
    100    continue
    close(3)

    ! enforce normailization
    gobs(1:ipixel)=gobs(1:ipixel)/xsum

    return
    end

    
    subroutine read_gobs_amsr2_v2(icase,ifreq,iscan,icel,g_compact)

        use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
        use gain_pattern, only: CompactGrid, read_gain_compact_grid, normalize_compact_grid,zero_compact_grid
        implicit none
    
        
        integer(int32),intent(in) :: icase        ! ==1 if target case
        integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
        integer(int32),intent(in) :: iscan        ! scan number
        integer(int32),intent(in) :: icel         ! fov number

        type(CompactGrid),intent(inout) :: g_compact

        character(80) filename
        integer :: ioerr
        character(len =80) :: iomsg
    
        if(icase.eq.1) then
            write(filename,9001)ifreq,iscan,icel
            9001 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/target_gains_v2/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
    
        else
            write(filename,9002)ifreq,iscan,icel
            9002 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/source_gains_v2/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
        endif
    
        open(unit=3,file=filename, status='old',action='read',form='unformatted', access='stream',iostat=ioerr,iomsg=iomsg)
                    
        if (ioerr /= 0) then
            write(6,*) iomsg
            error stop 1
        endif
    
        call zero_compact_grid(g_compact)
        call read_gain_compact_grid(3,g_compact)
        call normalize_compact_grid(g_compact)
        close(3)
    
        return
        end

        subroutine read_gobs_amsr2_v3(icase,ifreq,iscan,icel,g_compact)

            use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
            use gain_pattern, only: CompactGrid, read_gain_compact_grid, normalize_compact_grid,zero_compact_grid
            implicit none
        
            
            integer(int32),intent(in) :: icase        ! ==1 if target case
            integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
            integer(int16),intent(in) :: iscan        ! scan number
            integer(int16),intent(in) :: icel         ! fov number
    
            type(CompactGrid),intent(inout) :: g_compact
    
            character(80) filename
            
            integer :: ioerr
            character(len =80) :: iomsg
        
            if(icase.eq.1) then
                error stop 'v3 for targets makes no sense'
                write(filename,9001)ifreq,iscan,icel
                9001 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/target_gains_v2/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
            else
                write(filename,9002)ifreq,iscan,icel
                9002 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/source_gains_v3/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
            endif
        
            open(unit=3,file=filename, status='old',action='read',form='unformatted', access='stream',iostat=ioerr,iomsg=iomsg)
                        
            if (ioerr /= 0) then
                write(6,*) iomsg
                error stop 1
            endif
        
            call zero_compact_grid(g_compact)
            call read_gain_compact_grid(3,g_compact)
            call normalize_compact_grid(g_compact)
            close(3)
        
            return
            end

        subroutine read_tar_amsr2_v3(idiameter,ifreq,iscan,icel,g_compact)

                use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
                use gain_pattern, only: CompactGrid, read_gain_compact_grid, normalize_compact_grid,zero_compact_grid
                implicit none
            
                
                integer(int32),intent(in) :: idiameter        ! ==1 if target case
                integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
                integer(int32),intent(in) :: iscan        ! scan number
                integer(int32),intent(in) :: icel         ! fov number
        
                type(CompactGrid),intent(inout) :: g_compact
        
                character(100) filename
            
                integer :: ioerr
                character(len =80) :: iomsg
            
                write(filename,9002)idiameter,ifreq,iscan,icel
                    9002 format('/mnt/ops1p-ren/l/access/resampling/AMSR2/target_gains_v2/', &
                                i2.2,'km/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
            
                open(unit=3,file=filename, status='old',action='read',form='unformatted', access='stream',iostat=ioerr,iomsg=iomsg)
                            
                if (ioerr /= 0) then
                    write(6,*) iomsg
                    error stop 1
                endif
            
                call zero_compact_grid(g_compact)
                call read_gain_compact_grid(3,g_compact)
                call normalize_compact_grid(g_compact)
                close(3)
            
                return
                end