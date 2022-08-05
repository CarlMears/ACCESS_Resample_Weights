!     5/27/2019 changed 10/31/2020  sign for alpha_cold_mirror changes to implement counter clockwise scan
!     not sure if alpha_cold_mirror is a placeholder or a real value

!     7/3/2018 changed 5/27/2019.  89 ghz now has same eia as other channels,
!     hence only one set of geoloc are required. 7/3/2018 version saved with date suffix
!     see comment cc1 for the arrays that are affected by this


!     6/1/2018 change 7/3/2018. frd_rad(maxcel,maxscan,2) and tec(maxcel,maxscan,2) added 


!     4/19/2018 changed 6/1/2018,  geolocation variables are now given just for 11-37 ghz and 89 ghz rather than for every freq
!     nch is now 10, not 12

    module l2_module
    use sensor_module
    implicit none

!    earth dimensions match wgs84    
    real(8),    parameter :: re=6378.137d3        !meters
    real(8),    parameter :: rp=6356.752d3        !meters
    real(8),    parameter :: ffac=(rp/re)*(rp/re)
    real(8),    parameter :: geosync_alt   =42164.0e3      !satellite orbits, montenbruck and gill, page 294
      
    real(8),    parameter :: orbit_period=6082.58d0  !orbital period
      
    integer(4), parameter :: numscan0=nint(orbit_period/scan_period)  !typical number of scans per orbit

    real(8),    parameter :: time_cold_mirror   =1.      !time for cold mirror observation relative to begin of scan
!     real(8),    parameter :: alpha_cold_mirror  =100.      !azimth angle for cold mirror observation relative to begin of scan
    real(8),    parameter :: alpha_cold_mirror  =-100.      !azimth angle for cold mirror observation relative to begin of scan
     
    real(4), parameter :: ta_earth_min(nch)=(/ 55., 55., 55., 55., 55., 55., 55., 55., 55., 55./)
    real(4), parameter :: ta_earth_max(nch)=(/330.,330.,330.,330.,330.,330.,330.,330.,330.,330./)

    character(120) ice_filename,sss_filename,filename_percent_land
    character(120) filename_apc,filename_rsp

    !     ====================================================================================================
    !     ========================= tables read in by readin_data_files ======================================
    !     ====================================================================================================
    real(4) freq(nfreq)

    character(1) percent_land_tab(4320,2160,5)  !land table
    character(1) amask_ice(45,180,12)           !monthly climatology sea ice map
    integer(1) isal_corr(10,360,180,12)
    real(4) delta_apc(2,nfreq),chi_apc(2,nfreq)
    real(4) resample_wt(maxcel,-32:32,-14:14,9,2)  

    integer(4) ksat
    
    !     ====================================================================================================
    !     ========================================== mwi l1a variables  ======================================
    !     ====================================================================================================

    integer(4) iorbit_file,numscan
    integer(4) iflag_goddard(10,maxscan)
    real(8) scan_time(maxscan),scan_rpm(maxscan)
    real(8) scpos_eci2000(3,maxscan),scvel_eci2000(3,maxscan),scrpy(3,maxscan) 
    integer(2) iflag_scn(maxscan)
    integer(2) iflag_cal(nch,maxscan)


    !     ====================================================================================================
    !     ====================================== geolocation variables  ======================================
    !     ====================================================================================================
    real(8) orbit(maxscan),zang(maxscan),scpos_eci(3,maxscan),scvel_eci(3,maxscan)
    integer(1) kflag_sun(maxscan)

    real(4) scloc(3,maxscan)
    real(4) sunlat(maxscan),sunlon(maxscan),sundis(maxscan),sunvec_km(3,maxscan),moonvec_km(3,maxscan)
    real(4) alpha_sun(maxscan),beta_sun(maxscan)
    real(4) bmag_sc(3,maxscan)    
    
    !  for the following arrays, the last subscript 2 is removed 
    real(4)  cellat(maxcel,maxscan),cellon(maxcel,maxscan),celeia(maxcel,maxscan)
    real(4)  celazm(maxcel,maxscan),celsun(maxcel,maxscan)
    real(4)  reflat(maxcel,maxscan),reflon(maxcel,maxscan)     !geostationary lat and lon of reflected ray
    real(4) frd_rad(maxcel,maxscan),   tec(maxcel,maxscan)     !added july 3 2018
      
    real(4) cellat_cold(maxscan),cellon_cold(maxscan)
    end module l2_module
    
    

