    subroutine get_l0_filename(iorbit, filename_l0,iexist) 

        implicit none

        integer(4),intent(in)       :: iorbit
        character(:),allocatable,intent(inout)  :: filename_l0
        integer(4),intent(out)      :: iexist
        !character*(*) filename_l0

        logical(4) lexist
        integer(4) iorbit1 
        integer(4) iorbit2

        iorbit1 = 1+5000*(floor((iorbit-1)/5000.0))
        iorbit2 = iorbit1 + 4999

        write(filename_l0,9001) iorbit1, iorbit2, iorbit
        9001 format("/mnt/ops1p/j/AMSR/AMSRE/L0/r",i5.5,"_",i5.5,"/r",i5.5,".gz")
        inquire(file=filename_l0,exist=lexist)
        if(lexist) then
            iexist=1
        else
            iexist=0
        endif
        return
    end subroutine get_l0_filename



!       subroutine get_l2b_filename(iorbit, filename_l2b,iexist) 
!       implicit none

!     character*(*) filename_l2b
!     character(80) adir
!     character(9) aorbit
!       logical(4) lexist
!     integer(4) iorbit,iexist

!       if(iorbit.lt.1 .or. iorbit.ge.25000) stop 'orbit oob in get_l2b_filename, pgm stopped'

!     inquire(file='c:\dell-2200b.txt',exist=lexist)

!     if(lexist) then
!     if(iorbit.ge.    1 .and. iorbit.le.13000) adir='h:\amsre_l2b_v04\r00001_13000\amsr_l2b_r'
!     if(iorbit.ge.13001 .and. iorbit.le.25000) adir='i:\amsre_l2b_v04\r13001_25000\amsr_l2b_r'
!     else
!     if(iorbit.ge.    1 .and. iorbit.le.13000) adir='\\dell-2200b\h\amsre_l2b_v04\r00001_13000\amsr_l2b_r'
!     if(iorbit.ge.13001 .and. iorbit.le.25000) adir='\\dell-2200b\i\amsre_l2b_v04\r13001_25000\amsr_l2b_r'
!     endif

!     write(aorbit,9001) iorbit
!  9001 format(i5.5,'.dat') 
!     filename_l2b=trim(adir)//aorbit

!     inquire(file=filename_l2b,exist=lexist)
!     if(lexist) then
!     iexist=1
!     else
!     iexist=0
!     endif

!     return
!     end
