!file from mwi\routines\read_l1a

!this reads the l1a file from the simulator


!     jan 3 2018 change nov 2 2020
!     1.  scan_rpm(iscan) is set to exactly -26.5d0.  Apart from the minus sign, it was nearly this value anyway
!         the minus sign implements the counter-clockwise spin
!     2.  a very small scaling of 0.999923d0 is applied to scpos_eci and scvel_eci to make the altitude 
!         exactly equal 833 for iorbit=2, iscan=645

	subroutine read_l1a(ilu,iorbit,filename_l1a, start_time)
        use l2_module
        implicit none
 
        character(*) filename_l1a
        integer(4) ilu,iorbit

        real(8) time_last,start_time
        integer(4) iscan,iscan_last
        
        real(8) deltim
        
        call openbig(ilu,filename_l1a,'old')
        read(ilu)  iorbit_file,numscan,scan_time,scan_rpm,scpos_eci2000,scvel_eci2000,scrpy
        close(ilu)
      
!       =============== changes made on nov 2 2020 =======================
        scan_rpm=-26.5d0
!       scale to give exactly altitude of 833 for iorbit=2, iscan=645	 
        scpos_eci2000=0.999923d0*scpos_eci2000
        scvel_eci2000=0.999923d0*scvel_eci2000
!       ==================================================================
      
!     quality flags are currently set to zero 
        iflag_scn=0
	    iflag_cal=0

    	
	    if(iorbit.ne.iorbit_file) stop 'iorbit.ne.iorbit_file, pgm stopped'
	    if(numscan.gt.3500)  stop 'numscan oob, pgm stopped'
        start_time=scan_time(1)
	
	    iscan_last=-1
	    do iscan=1,numscan
	        if(iscan_last.ne.-1) then
	            deltim=scan_time(iscan)-time_last
	            if(abs((iscan-iscan_last)*scan_period-deltim).gt.0.1) then
	                write(*,6001) iorbit,iscan,iscan-iscan_last,scan_time(iscan),deltim
                    6001 format('time gap ',i8, 2i6, 2f18.4,2i4)
	            endif
	        endif
	        time_last=scan_time(iscan)
	        iscan_last=iscan
        enddo  !iscan
	return
	end subroutine




