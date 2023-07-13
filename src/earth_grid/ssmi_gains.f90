module ssmi_gains

    implicit none
    private

    public read_gain_obs_ssmi,read_gain_target_ssmi

contains
    subroutine read_gain_obs_ssmi(ksat,ifreq,iscan,icel,g_compact)

        use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
        use gain_pattern, only: CompactGrid, read_gain_compact_grid, normalize_compact_grid,zero_compact_grid

        implicit none
        
        integer(int32),intent(in) :: ksat
        integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
        integer(int16),intent(in) :: iscan        ! scan number
        integer(int16),intent(in) :: icel         ! fov number

        type(CompactGrid),intent(inout) :: g_compact

        character(80) filename
        
        integer :: ioerr
        character(len =80) :: iomsg
    
        write(filename,9002)ksat,ifreq,iscan,icel
        9002 format('/mnt/ops1p-ren/l/access/resampling/SSMI/f',i2.2,'/source_gains/freq_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
        !print *,filename
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
    end subroutine read_gain_obs_ssmi

    subroutine read_gain_target_ssmi(idiameter,ksat,ifreq,iscan,icel,g_compact)

        use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
        use gain_pattern, only: CompactGrid, read_gain_compact_grid, normalize_compact_grid,zero_compact_grid

        implicit none
    
        integer(int32),intent(in) :: idiameter    ! ==1 if target case
        integer(int32),intent(in) :: ksat         ! SSMIS satellite number 16,17,18   
        integer(int32),intent(in) :: ifreq        ! same as iband -- 1-8 for amsr2
        integer(int32),intent(in) :: iscan        ! scan number
        integer(int32),intent(in) :: icel         ! fov number

        type(CompactGrid),intent(inout) :: g_compact

        character(100) filename
    
        integer :: ioerr
        character(len =80) :: iomsg
    
        write(filename,9002)ksat,idiameter,ifreq,iscan,icel
        9002 format('/mnt/ops1p-ren/l/access/resampling/SSMI/f',i2.2,'/target_gains/', &
                        i2.2,'km/freq_',i2.2,'/s',i2.2,'c',i3.3,'.dat')
        !print *,filename
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
    end subroutine read_gain_target_ssmi
end module ssmi_gains