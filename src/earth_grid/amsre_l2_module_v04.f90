    module amsre_l2_module

    implicit none

    integer(4), parameter :: maxscan=4400
    integer(4), parameter :: maxcel=243
    integer(4), parameter :: maxscan_89=2*maxscan
    integer(4), parameter :: maxcel_89=2*maxcel

    character(100) land_filename,ice_filename,sst_filename,climate_filename,sss_filename,tbocean_filename

    real(8) timeofscan(maxscan),orbit(maxscan)
    real(4) omega(maxscan)
    real(8) scpos(3,maxscan),scvel(3,maxscan),attdq(4,maxscan),dattdq(4,maxscan)
    real(4) xscan(16,maxcel)

    real(4) scrpy(3,maxscan),scloc(3,maxscan)
    real(4) cellat(maxcel,maxscan),cellon(maxcel,maxscan),celtht(maxcel,maxscan),celphi(maxcel,maxscan),celsun(maxcel,maxscan)
    real(4) celrfi(maxcel,maxscan)
    real(4) cellat_89(maxcel_89,maxscan_89),cellon_89(maxcel_89,maxscan_89),cellat_cold(maxscan),cellon_cold(maxscan)
    real(4) sunlat(maxscan),sunlon(maxscan),sundis(maxscan)
    real(4) therm(8,maxscan),therm_sps(16,maxscan),avg_count(2,16,maxscan),therm_avg(8,maxscan)
    real(4) ta(maxcel,maxscan,12),ta_89(maxcel_89,maxscan_89,2),tar(maxcel,maxscan,22)
    real(4) ta_cold(16,maxscan),ta_hot(16,maxscan),beta_term(maxcel,maxscan,2),ta_linear(maxcel,maxscan,2)
    real(4) sstcl(maxcel,maxscan),wincl(maxcel,maxscan),vapcl(maxcel,maxscan),phir(maxcel,maxscan),ep(0:9,maxcel,maxscan)
    real(4) wgcm(maxcel,maxscan), sst_reg(maxcel,maxscan)
    integer(4) numscan,num_count(2,16,maxscan)

    integer(2) leap_sec(maxscan)
    integer(2) earthcounts(16,maxcel_89,maxscan),csmcounts(16,32,maxscan),htscounts(16,32,maxscan),rx_offset_gain(16, 2,maxscan)
    integer(2) itbsss(10,maxcel,maxscan)

    integer(1) iflag_l0(maxscan),iflag_rx(16,maxscan),iflag_cal(16,maxscan),isurcel(5,maxcel,maxscan),ice_flag1(maxcel,maxscan)
    integer(1) ice_flag2(maxcel,maxscan),ioob_flag(4,maxcel,maxscan),iret_flag(4,maxcel,maxscan),irainadj(4,maxcel,maxscan)

    integer(1) ice_map(1440,720,2)

    end module amsre_l2_module

