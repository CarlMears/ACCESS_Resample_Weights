subroutine read_gobs_amsre(icase,npixel_max,nobs_max,ifreq,iscan,icel,iobs, npixel,ilat2,ilon2,gobs)

    use, intrinsic :: iso_fortran_env, only: int16, int32, real32, real64
    implicit none

    
    integer(int32),intent(in) :: icase
    integer(int32),intent(in) :: npixel_max
    integer(int32),intent(in) :: nobs_max
    integer(int32),intent(in) :: ifreq
    integer(int32),intent(in) :: iobs
    integer(int32),intent(in) :: iscan
    integer(int32),intent(in) :: icel


    integer(int16), intent(inout) :: ilat2(npixel_max,0:nobs_max)
    integer(int16), intent(inout) :: ilon2(npixel_max,0:nobs_max)
    integer(int32), intent(inout) :: npixel(0:nobs_max)
    integer(real64),intent(inout) :: gobs(npixel_max,0:nobs_max)

    character(80) filename
    integer(int16) ilat2x,ilon2x
    integer(int32) :: i,ipixel
    real(real32) gobsx
    real(real64) xsum

    integer :: ioerr,lu 
    character(len =80) :: iomsg

    if(icase.eq.1) then
        write(filename,9001) ifreq,iscan+15,icel
        9001 format('c:\target_arrays\freq',i1,'\s',i2.2,'c',i3.3,'.dat')
    else
        write(filename,9002) ifreq,iscan+15,icel
        9002 format('/mnt/oserver/o/amsr_l2/v05/resampling/gain_array/freq',i1,'/s',i2.2,'c',i3.3,'.dat')
    endif

    !call openbig(3,filename,'old')
    open(unit=3,file=filename, status='old',action='read',form='unformatted', access='stream',iostat=ioerr,iomsg=iomsg)
    if (ioerr /= 0) then
        write(6,*) iomsg
        error stop 1
    endif
     
    ipixel=0; xsum=0
    
    do i=1,1234567890
        read(3,end=100) ilat2x,ilon2x,gobsx
        ipixel=ipixel+1
        if(ipixel.gt.npixel_max) stop 'pgm stopped, ipixel oob'
        ilat2(ipixel,iobs)=ilat2x
        ilon2(ipixel,iobs)=ilon2x
        gobs( ipixel,iobs)=gobsx 
        npixel(      iobs)=ipixel
        xsum=xsum+gobs( ipixel,iobs)
    enddo
    100    continue
    close(3)

    ! enforce normailization
    gobs(1:ipixel,iobs)=gobs(1:ipixel,iobs)/xsum

    return
    end