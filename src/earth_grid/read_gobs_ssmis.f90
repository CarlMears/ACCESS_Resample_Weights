module ssmis_gains

    implicit none
    private

    public read_gain_obs_ssmis,ead_gain_target_ssmis

contains
    subroutine read_gain_obs_ssmis(icase,ksat,ifreq,iscan,icel,g_compact)

        use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
        use gain_pattern, only: CompactGrid, read_gain_compact_grid, normalize_compact_grid,zero_compact_grid
        implicit none
    
        
        integer(int32),intent(in) :: icase        ! ==1 if target case
        integer(int32),intent(in) :: ksat
        integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
        integer(int16),intent(in) :: iscan        ! scan number
        integer(int16),intent(in) :: icel         ! fov number

        type(CompactGrid),intent(inout) :: g_compact

        character(80) filename
        
        integer :: ioerr
        character(len =80) :: iomsg
    
        if(icase.eq.1) then
            write(filename,9001)ksat,beamwidth_km_int,iscan,icel
            9001 format('/mnt/ops1p-ren/l/access/resampling/SSMIS/f',i2.2,'/target_gains/',i2.2,'km/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
        else
            write(filename,9002)ksat,ifreq,iscan,icel
            9002 format('/mnt/ops1p-ren/l/access/resampling/SSMIS/f',i2.2,'/source_gains/band_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
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
    end subroutine read_gain_obs_ssmis

    subroutine read_gain_target_ssmis(idiameter,ksat,ifreq,iscan,icel,g_compact)

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
    end subroutine read_gain_target_ssmis
end module ssmis_gains