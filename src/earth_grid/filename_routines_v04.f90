subroutine get_l0_filename(iorbit, filename_l0,iexist) 
    implicit none

    character*(*) filename_l0
    character(80) adir
    character(8) aorbit
    logical(4) lexist
    integer(4) iorbit,iexist

    inquire(file='c:\dell-2200b.txt',exist=lexist)

    if(iorbit.lt.1 .or. iorbit.ge.24500) stop 'orbit oob in get_l0_filename, pgm stopped'
    if(lexist) then
        if(iorbit.ge.    1 .and. iorbit.le. 3600) adir='e:\amsre_l0\r00001_03600\r'
        if(iorbit.ge. 3601 .and. iorbit.le.10500) adir='f:\amsre_l0\r03601_10500\r'
        if(iorbit.ge.10501 .and. iorbit.le.17500) adir='g:\amsre_l0\r10501_17500\r'
        if(iorbit.ge.17501 .and. iorbit.le.24500) adir='z:\amsre_l0\r17501_24500\r'
    else
        if(iorbit.ge.    1 .and. iorbit.le. 3600) adir='\\dell-2200b\e\amsre_l0\r00001_03600\r'
        if(iorbit.ge. 3601 .and. iorbit.le.10500) adir='\\dell-2200b\f\amsre_l0\r03601_10500\r'
        if(iorbit.ge.10501 .and. iorbit.le.17500) adir='\\dell-2200b\g\amsre_l0\r10501_17500\r'
        if(iorbit.ge.17501 .and. iorbit.le.24500) adir='\\dell-2200b\z\amsre_l0\r17501_24500\r'
    endif


    write(aorbit,9001) iorbit
 9001 format(i5.5,'.gz') 
    filename_l0=trim(adir)//aorbit
    inquire(file=filename_l0,exist=lexist)
    if(lexist) then
    iexist=1
    else
    iexist=0
    endif

    return
    end


      subroutine get_l2b_filename(iorbit, filename_l2b,iexist) 
      implicit none

    character*(*) filename_l2b
    character(80) adir
    character(9) aorbit
      logical(4) lexist
    integer(4) iorbit,iexist

      if(iorbit.lt.1 .or. iorbit.ge.25000) stop 'orbit oob in get_l2b_filename, pgm stopped'

    inquire(file='c:\dell-2200b.txt',exist=lexist)

    if(lexist) then
    if(iorbit.ge.    1 .and. iorbit.le.13000) adir='h:\amsre_l2b_v04\r00001_13000\amsr_l2b_r'
    if(iorbit.ge.13001 .and. iorbit.le.25000) adir='i:\amsre_l2b_v04\r13001_25000\amsr_l2b_r'
    else
    if(iorbit.ge.    1 .and. iorbit.le.13000) adir='\\dell-2200b\h\amsre_l2b_v04\r00001_13000\amsr_l2b_r'
    if(iorbit.ge.13001 .and. iorbit.le.25000) adir='\\dell-2200b\i\amsre_l2b_v04\r13001_25000\amsr_l2b_r'
    endif

    write(aorbit,9001) iorbit
 9001 format(i5.5,'.dat') 
    filename_l2b=trim(adir)//aorbit

    inquire(file=filename_l2b,exist=lexist)
    if(lexist) then
    iexist=1
    else
    iexist=0
    endif

    return
    end
