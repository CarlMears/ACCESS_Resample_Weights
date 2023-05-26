! Routines to set the simulation parameters
!
! The parameters here are combined from geolocation parameters, and common_sim_data in the MWI simulator

module simulation_param
  use, intrinsic :: iso_fortran_env, only: int32, real32, real64
  !use geolocation_types, only: GeolocationParameterData_f
  !use cfgio_mod, only: cfg_t, parse_cfg
  implicit none
  private
  public :: set_simulation_params,SimulationParameterData_f
  public :: RE,RP,FFAC,FRA_REF_FREQ
  public :: MAXSCAN,FOVS_PER_SCAN,NBAND,NPOLBAND,NCHAN,NCHAN_SIM,NEXTRAFOV,NEXTRASCAN
  public :: BAND_FREQS,POLBAND_FREQS,NFREQ
  public :: SCAN_PERIOD,ORBIT_PERIOD,SAMPLINGTIME

  real(real64) :: RE    !6378.137d3 ! meters
  real(real64) :: RP    !6356.752d3 ! meters
  real(real64) :: FFAC  !=(RP/RE)**2

  real(real32)   :: FRA_REF_FREQ

  integer(int32) :: FOVS_PER_SCAN 
  integer(int32) :: MAXSCAN
  integer(int32) :: NFREQ
  integer(int32) :: NBAND
  integer(int32) :: NPOLBAND
  integer(int32) :: NCHAN
  integer(int32) :: NCHAN_SIM 
  integer(int32) :: NEXTRAFOV
  integer(int32) :: NEXTRASCAN

  real(real64)   :: SCAN_PERIOD
  real(real64)   :: ORBIT_PERIOD
  real(real64)   :: SAMPLINGTIME

  real(real32), dimension(:), pointer :: BAND_FREQS != [10.85, 18.85, 23.8, 36.75, 37.3, 89.0]
  real(real32), dimension(:), pointer :: POLBAND_FREQS != [10.85, 18.85, 36.75]

    ! Simulation parameter data
  type SimulationParameterData_f
     ! earth dimensions match wgs84
     real(real64) :: RE    !6378.137d3 ! meters
     real(real64) :: RP    !6356.752d3 ! meters
     real(real64) :: FFAC  !=(RP/RE)**2

     ! Faraday rotation angle reference frequency in GHz
     real(real32) :: FRA_REF_FREQ !=1.413

     ! -------------------------------
     ! sensor-specific constants
     integer(int32) :: FOVS_PER_SCAN != 553
     integer(int32) :: MAXSCAN != 4400

     !----------------------------------
     ! variables for adding extra FOVS -- these are set from the command line
     integer(int32) :: NEXTRAFOV
     integer(int32) :: NEXTRASCAN 

     ! Number of frequencies: 10.85, 18.85, 23.80, 36.75, 37.3, 89 GHz
     integer(int32) :: NFREQ != 6

     ! Number of frequency bands
     integer(int32) :: NBAND != 6
 
     ! Number of polarimetric bands (10.85, 18.85, 36.75 GHz)
     integer(int32) :: NPOLBAND != 3
 
     ! Number of channels, or frequency/polarization combinations
     integer(int32) :: NCHAN != 17
     integer(int32) :: NCHAN_SIM != 18

     real(real32), dimension(:), pointer :: BAND_FREQS != [10.85, 18.85, 23.8, 36.75, 37.3, 89.0]
     real(real32), dimension(:), pointer :: POLBAND_FREQS != [10.85, 18.85, 36.75]
   
     integer, dimension(:), pointer :: CHAN_SIM_TO_TB_FREQ
     integer, dimension(:), pointer :: CHAN_SIM_TO_TB_POL

     real(real64) :: SCAN_PERIOD != 60.0d0 / 26.5d0
     real(real64) :: ORBIT_PERIOD !=6082.58d0
     integer(int32) :: NUMSCAN0  !=nint(ORBIT_PERIOD/SCAN_PERIOD)

    ! The sample time (or integration time or tau) for each measurement is 1.8 ms
     real(real64) :: SAMPLE_TIME != 1.8d-3
     real(real32) :: SamplingTime
     real(real32), dimension(:), pointer :: OffNadirAngle
     real(real32), dimension(:), pointer :: FirstSampleScanAngle
     real(real32), dimension(:), pointer :: FirstSampleDelay
     real(real32), dimension(:), pointer :: PolarizationTilt
     real(real32) :: EarthEquatorialRadius
     real(real32) :: EarthPolarRadius
     real(real32), dimension(:), pointer :: ColdSkyScanAngle
     real(real32) :: SolarRadioFlux

     integer(int32) :: NAMBIG
  end type SimulationParameterData_f

contains

  ! Set the constants
  subroutine set_simulation_params(sensor_name,data)
    
    character(len=*), intent(in) :: sensor_name
    type(SimulationParameterData_f), intent(inout) :: data

      ! earth dimensions match wgs84
    RE=6378.137d3 ! meters
    data%RE=RE ! meters

    RP=6356.752d3 ! meters
    data%RP=RP ! meters

    FFAC=(data%RP/data%RE)**2
    data%FFAC=FFAC

    ! Faraday rotation angle rreference frequency in GHz
    FRA_REF_FREQ=1.413
    data%FRA_REF_FREQ=FRA_REF_FREQ

    data%NAMBIG = 4

    select case(trim(sensor_name))
      case('MWI')
          ! -------------------------------
        ! MWI-specific constants
        FOVS_PER_SCAN = 553
        data%FOVS_PER_SCAN = FOVS_PER_SCAN

        MAXSCAN = 4400
        data%MAXSCAN = MAXSCAN
  
        ! Number of frequencies: 10.85, 18.85, 23.80, 36.75, 37.3, 89 GHz
        NFREQ = 6
        data%NFREQ = NFREQ
  
        ! Number of frequency bands
        NBAND = 6
        data%NBAND = NBAND
  
        ! Number of polarimetric bands (10.85, 18.85, 36.75 GHz)
        NPOLBAND = 3
        data%NPOLBAND = NPOLBAND
  
        ! Number of channels, or frequency/polarization combinations
        NCHAN = 17
        data%NCHAN = NCHAN
  
        ! Number of channels, or frequency/polarization combinations (for
        ! simulation purposes, an extra channel at 23 GHz h-pol is added)
        NCHAN_SIM = 18
        data%NCHAN_SIM = NCHAN_SIM

        allocate(BAND_FREQS(NBAND), &
                 POLBAND_FREQS(NPOLBAND))
        allocate(data%BAND_FREQS(data%NBAND), &
                 data%POLBAND_FREQS(data%NPOLBAND))

        BAND_FREQS = [10.85, 18.85, 23.8, 36.75, 37.3, 89.0]
        POLBAND_FREQS = [10.85, 18.85, 36.75]
        data%BAND_FREQS = BAND_FREQS
        data%POLBAND_FREQS = POLBAND_FREQS
   
        ! ----------------------------------
        ! Mapping of MWI channels to frequency band/polarization indices. In
        ! order: 10.85 v/h/3/4, 18.85 v/h/3/4, 23.8 v/h, 36.75 v/h/3/4, 37.3
        ! v/h, 89 v/h. (NOTE that the 23.8 h-pol channel is included here
        ! for simulation purposes.)
        allocate(data%CHAN_SIM_TO_TB_FREQ(data%NCHAN_SIM),  &
                 data%CHAN_SIM_TO_TB_POL(data%NCHAN_SIM))
  
      
  
        data%CHAN_SIM_TO_TB_FREQ = [ &
         1, 1, 1, 1, & ! 10.85 GHz
         2, 2, 2, 2, & ! 18.85 GHz
         3, 3, & ! 23.8 GHz
         4, 4, 4, 4, & ! 36.75 GHz
         5, 5, & ! 37.3 GHz
         6, 6] ! 89 GHz
        data%CHAN_SIM_TO_TB_POL = [ &
         1, 2, 3, 4, & ! 10.85 GHz
         1, 2, 3, 4, & ! 18.85 GHz
         1, 2, & ! 23.8 GHz
         1, 2, 3, 4, & ! 36.75 GHz
         1, 2, & ! 37.3 GHz
         1, 2] ! 89 GHz
  
        ! 26.5 RPM is 2.26 s per scan
        data%SCAN_PERIOD = 60.0d0 / 26.5d0
        data%ORBIT_PERIOD=6082.58d0
  
        ! typical number of scans per orbit
        data%NUMSCAN0=nint(data%ORBIT_PERIOD/data%SCAN_PERIOD)
 
        ! Maximum number of wind ambiguities
    
        ! The sample time (or integration time or tau) for each
        ! measurement is 1.8 ms. The units here are in milliseconds.

        data%SamplingTime = 1.8
      
        allocate(data%OffNadirAngle(data%NBAND), &
                 data%FirstSampleScanAngle(data%NBAND), &
                 data%FirstSampleDelay(data%NBAND), &
                 data%PolarizationTilt(data%NPOLBAND))
        ! Off-nadir angle of the boresight in degrees for each band (same for all bands for MWI)
        data%OffNadirAngle(:) = 44.3
        data%FirstSampleScanAngle(:) = -78.9912
        data%FirstSampleDelay(:) = 0.0

        ! These values are unknown so set to zero
        data%PolarizationTilt(:) = 0.0

        ! This is based on data from David Draper
        allocate(data%ColdSkyScanAngle(data%NBAND))
        data%ColdSkyScanAngle(:) = [158.145, 198.402, 198.389, 181.177, 181.177, 181.111]

     
        data%SolarRadioFlux = 110.

    case('AMSR2')

        ! -------------------------------
      ! AMSR2-specific constants
      FOVS_PER_SCAN = 243      
      data%FOVS_PER_SCAN = FOVS_PER_SCAN 

      MAXSCAN = 4400
      data%MAXSCAN = MAXSCAN

      ! these are set using the command line for now.
      data%NEXTRAFOV = NEXTRAFOV
      data%NEXTRASCAN = NEXTRASCAN

      ! Number of frequencies: 10.85, 18.85, 23.80, 36.75, 37.3, 89 GHz
      NFREQ = 7
      data%NFREQ = NFREQ

      ! Number of frequency bands  - included 89 a + b
      NBAND = 8
      data%NBAND = NBAND

      ! Number of polarimetric bands (10.85, 18.85, 36.75 GHz)
      NPOLBAND = 0
      data%NPOLBAND = NPOLBAND

      ! Number of channels, or frequency/polarization combinations
      NCHAN = 16
      data%NCHAN = NCHAN

      ! Number of channels, or frequency/polarization combinations (for
      ! simulation purposes, an extra channel at 23 GHz h-pol is added)
      NCHAN_SIM = 16
      data%NCHAN_SIM = NCHAN_SIM

      allocate(BAND_FREQS(NBAND))
      allocate(data%BAND_FREQS(data%NBAND))
              
      BAND_FREQS = [6.925,7.3,10.65,18.7,23.8,36.5,89.0,89.0]
      data%BAND_FREQS = BAND_FREQS
      

      ! ----------------------------------
      ! Mapping of MWI channels to frequency band/polarization indices. In
      ! order: 10.85 v/h/3/4, 18.85 v/h/3/4, 23.8 v/h, 36.75 v/h/3/4, 37.3
      ! v/h, 89 v/h. (NOTE that the 23.8 h-pol channel is included here
      ! for simulation purposes.)
      allocate(data%CHAN_SIM_TO_TB_FREQ(data%NCHAN_SIM),  &
               data%CHAN_SIM_TO_TB_POL(data%NCHAN_SIM))

    

      data%CHAN_SIM_TO_TB_FREQ = [ &
       1, 1, & ! 6.925 GHz
       2, 2, & ! 7.3 GHz
       3, 3, & ! 10.65 GHz
       4, 4, & ! 18.7 GHz
       5, 5, & ! 23.8 GHz
       6, 6, & ! 36.5 GHz
       7, 7, & ! 89.0 Ghz (a)
       8, 8] ! 89.0 GHz (b)
      data%CHAN_SIM_TO_TB_POL = [ &
       1, 2, & ! 6.925 GHz
       1, 2, & ! 7.3 GHz
       1, 2, & ! 10.65 GHz
       1, 2, & ! 18.7 GHz
       1, 2, & ! 23.8 GHz
       1, 2, & ! 36.5 Ghz
       1, 2, & ! 89 GHz
       1, 2]  ! 89 GHz

      ! 40.0 RPM is 1.5 s per scan
      SCAN_PERIOD = 60.0d0 / 40.0d0
      data%SCAN_PERIOD = SCAN_PERIOD

      ORBIT_PERIOD = 5929.564d0 ! from jan 2020 TLE
      data%ORBIT_PERIOD = ORBIT_PERIOD

      ! typical number of scans per orbit
      data%NUMSCAN0=nint(data%ORBIT_PERIOD/data%SCAN_PERIOD)

      ! Maximum number of wind ambiguities
  
      ! The sample time (or integration time or tau) for each
      ! measurement is 1.8 ms. The units here are in milliseconds.

      SAMPLINGTIME = 2.6
      data%SAMPLINGTIME = 2.6
    
      allocate(data%OffNadirAngle(data%NBAND), &
               data%FirstSampleScanAngle(data%NBAND), &
               data%FirstSampleDelay(data%NBAND), &
               data%PolarizationTilt(data%NCHAN_SIM))
      ! Off-nadir angle of the boresight in degrees for each band (same for all bands for MWI)
      data%OffNadirAngle(1:7) = 47.505  !nominal from user manual
      data%OffNadirAngle(8) = 47.115 !needs to be something different for the in between scan

      data%FirstSampleScanAngle(1:6) = -75.083 !360.0/data%scan_period * 0.4375 -0.083 - 180.0
      data%FirstSampleScanAngle(7) = -75.395   !360.0/data%scan_period * 0.4375 -0.395 - 180.0
      data%FirstSampleScanAngle(8) = -75.725   !360.0/data%scan_period * -0.4375 -0.725 - 180.0


      data%FirstSampleDelay(:) = 0.0

      ! These values are unknown so set to zero
      data%PolarizationTilt(:) = 0.0

      ! This is based on data from David Draper
      allocate(data%ColdSkyScanAngle(data%NBAND))
      data%ColdSkyScanAngle(:) = [158.145, 198.402, 198.389, 181.177, 181.177, 181.111]

   
      data%SolarRadioFlux = 110.
      
    case('SSMIS')

      ! -------------------------------
      ! AMSR2-specific constants
      FOVS_PER_SCAN = 243      
      data%FOVS_PER_SCAN = FOVS_PER_SCAN 

      MAXSCAN = 4400
      data%MAXSCAN = MAXSCAN

      ! these are set using the command line for now.
      data%NEXTRAFOV = NEXTRAFOV
      data%NEXTRASCAN = NEXTRASCAN

      NFREQ = 7
      data%NFREQ = NFREQ

      ! Number of frequency bands  - included 89 a + b
      NBAND = 8
      data%NBAND = NBAND

      ! Number of polarimetric bands
      NPOLBAND = 0
      data%NPOLBAND = NPOLBAND

      !  Number of channels, or frequency/polarization combinations
      NCHAN = 16
      data%NCHAN = NCHAN

      ! Number of channels, or frequency/polarization combinations (for
      ! simulation purposes, an extra channel at 23 GHz h-pol is added)
      NCHAN_SIM = 16
      data%NCHAN_SIM = NCHAN_SIM

      allocate(BAND_FREQS(NBAND))
      allocate(data%BAND_FREQS(data%NBAND))
            
      BAND_FREQS = [6.925,7.3,10.65,18.7,23.8,36.5,89.0,89.0]
      data%BAND_FREQS = BAND_FREQS
    

      ! ----------------------------------
      ! Mapping of MWI channels to frequency band/polarization indices. In
      ! order: 10.85 v/h/3/4, 18.85 v/h/3/4, 23.8 v/h, 36.75 v/h/3/4, 37.3
      ! v/h, 89 v/h. (NOTE that the 23.8 h-pol channel is included here
      ! for simulation purposes.)
      allocate(data%CHAN_SIM_TO_TB_FREQ(data%NCHAN_SIM),  &
               data%CHAN_SIM_TO_TB_POL(data%NCHAN_SIM))

  

    data%CHAN_SIM_TO_TB_FREQ = [ &
     1, 1, & ! 6.925 GHz
     2, 2, & ! 7.3 GHz
     3, 3, & ! 10.65 GHz
     4, 4, & ! 18.7 GHz
     5, 5, & ! 23.8 GHz
     6, 6, & ! 36.5 GHz
     7, 7, & ! 89.0 Ghz (a)
     8, 8] ! 89.0 GHz (b)
    data%CHAN_SIM_TO_TB_POL = [ &
     1, 2, & ! 6.925 GHz
     1, 2, & ! 7.3 GHz
     1, 2, & ! 10.65 GHz
     1, 2, & ! 18.7 GHz
     1, 2, & ! 23.8 GHz
     1, 2, & ! 36.5 Ghz
     1, 2, & ! 89 GHz
     1, 2]  ! 89 GHz

    ! 40.0 RPM is 1.5 s per scan
    SCAN_PERIOD = 60.0d0 / 40.0d0
    data%SCAN_PERIOD = SCAN_PERIOD

    ORBIT_PERIOD = 5929.564d0 ! from jan 2020 TLE
    data%ORBIT_PERIOD = ORBIT_PERIOD

    ! typical number of scans per orbit
    data%NUMSCAN0=nint(data%ORBIT_PERIOD/data%SCAN_PERIOD)

    ! Maximum number of wind ambiguities

    ! The sample time (or integration time or tau) for each
    ! measurement is 1.8 ms. The units here are in milliseconds.

    SAMPLINGTIME = 2.6
    data%SAMPLINGTIME = 2.6
  
    allocate(data%OffNadirAngle(data%NBAND), &
             data%FirstSampleScanAngle(data%NBAND), &
             data%FirstSampleDelay(data%NBAND), &
             data%PolarizationTilt(data%NCHAN_SIM))
    ! Off-nadir angle of the boresight in degrees for each band (same for all bands for MWI)
    data%OffNadirAngle(1:7) = 47.505  !nominal from user manual
    data%OffNadirAngle(8) = 47.115 !needs to be something different for the in between scan

    data%FirstSampleScanAngle(1:6) = -75.083 !360.0/data%scan_period * 0.4375 -0.083 - 180.0
    data%FirstSampleScanAngle(7) = -75.395   !360.0/data%scan_period * 0.4375 -0.395 - 180.0
    data%FirstSampleScanAngle(8) = -75.725   !360.0/data%scan_period * -0.4375 -0.725 - 180.0


    data%FirstSampleDelay(:) = 0.0

    ! These values are unknown so set to zero
    data%PolarizationTilt(:) = 0.0

    ! This is based on data from David Draper
    allocate(data%ColdSkyScanAngle(data%NBAND))
    data%ColdSkyScanAngle(:) = [158.145, 198.402, 198.389, 181.177, 181.177, 181.111]

 
    data%SolarRadioFlux = 110.
    case default
      print *,'Sensor ',sensor_name,' not implemented'
      stop
    end select

  end subroutine set_simulation_params

end module simulation_param
