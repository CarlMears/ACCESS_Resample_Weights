! forked from o:/ssmis/routines/fdl1b.f

!         6/18/2016 updated 6/20/2016.  more comments added
!         2/1/2010 changed 6/18/2016.  to fix f17 37v anomaly the argument iorbit is added to some routines
!         also call to set_iflag_chn is added

!         12/23/2009 version changed on 2/1/2010.  the following line was added:
!        call avg_cal_counts !this routine assumes future ssmis (18,19,20) will not have scan averaging aboard

!         not sure what happen between 12/16/2009 and 12/23/2009.  there was a lot of misc development during this period

!         12/8/2009 version changed on 12/16/2009. fd_ta_adjusted is now called after calling ssmis_geolocation.
!         for this current version of ssmis this should make no dif, but i think it is better practice to call geoloc first
!         for example, ssmi need celtht when doing ta adjustment

    subroutine fd_l1b(ioption,iorbit)  !ioption=1 means do 92 ghz
        use l2_module                                             
        implicit none

        integer(4), intent(in) :: ioption,iorbit

        real(4) diag(2,maxchn_max)
        integer(4) kgrp,nbad_counts(7)                                  

        call avg_cal_counts(iorbit) !this routine assumes future ssmis (18,19,20) will not have scan averaging aboard
        call ssmi_geolocation(0)    !0 denotes find s/c and sun/moon parameters
        call find_moon_angle        ! finds moon angles for all 6 location groups, yet-to-be specified groups are assumed to be igrp=4
        call ck_counts(iorbit)
        call ck_counts1(iorbit, nbad_counts)
        call set_iflag_chn(iorbit)  !always do the qc
        call fd_arm_temp

        if(ioption.eq.1) then
            kgrp=1                               !1 denotes kgrp=1, which is the img channels
            call prep_group(kgrp)                !1 denotes kgrp=1, which is the img channels
            call intrep_cold_counts(1,iorbit)    !1 denotes igrp=1 which is the first  part of the img channels            
            call intrep_cold_counts(2,iorbit)    !2 denotes igrp=2 which is the second part of the img channel
            call shift_cal_counts
            call ssmi_geolocation(2)             !2 denotes igrp=2, which is 91 ghz
            call fd_ta_adjusted(kgrp, diag)       

            !save kgrp=1 data so it wont be loss when get kgrp=2 data
            maxcel_hi=maxcel
            nch_hi=nch
            iflag_cal_hi=iflag_cal
            cellat_hi=cellat
            cellon_hi=cellon
            celtht_hi=celtht
            celphi_hi=celphi
            celsun_hi=celsun
            ta_hi=ta
        endif

        kgrp=2                               !2 denotes kgrp=2, which is the env channels
        call prep_group(kgrp)             
        call intrep_cold_counts(3,iorbit)    !3 denotes igrp=3 which is the first  part of the env channels            
        call intrep_cold_counts(4,iorbit)    !4 denotes igrp=4 which is the second part of the env channel
        call shift_cal_counts
        call ssmi_geolocation(4)             !4 denotes igrp=4, which is 37 ghz
        call fd_ta_adjusted(kgrp, diag)

        return
    end subroutine fd_l1b
