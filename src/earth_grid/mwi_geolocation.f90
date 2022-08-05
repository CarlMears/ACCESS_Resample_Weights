! ---------------------------------------------------------------------------
!
! Ball Aerospace & Technologies Corp.
! ---------------------------------------------------------------------------
! WSF-M PROJECT: MWSDPS - MicroWave Sensor Data Processing Software
! ---------------------------------------------------------------------------
! This is an unpublished work, 2020 Ball Aerospace & Technologies
! Corp. and Remote Sensing Systems Inc. All rights are reserved,
! subject to the U.S. Government license described herein. This work
! was developed in whole or in part with U.S. Government sponsorship
! and was delivered to the Government pursuant to DFARS
! 252.227-13(b)(1) with unlimited rights, including the rights to
! reproduce, prepare derivative works, distribute copies to the
! public, and perform publicly and display publicly, by or on behalf
! of the Government.
! ---------------------------------------------------------------------------
! Government Purpose Rights Notice Contract Number: FA8810-18-C-0002
! Contractor Name: Ball Aerospace & Technologies Corp. Contractor
! Address: 1600 Commerce Street, Boulder, Colorado 80301
! ---------------------------------------------------------------------------
! The Government's rights to use, modify, reproduce, release, perform,
! display, or disclose these technical data are restricted by
! paragraph (b)(2) of the Rights in Technical Data Non-commercial
! Items clause contained in the above identified contract. Any
! reproduction of technical data or portions thereof marked with this
! legend must also reproduce the markings. DISTRIBUTION STATEMENT D:
! Distribution authorized to the Department of Defense and U.S. DoD
! contractors only. Reason: Critical Technology, Export Controlled.
! Date of determination: 15 Nov 2017. Other requests shall be referred
! to SMC/RSK, El Segundo, CA, 90245-2808. WARNING - This document
! contains technical data whose export is restricted by the Arms
! Export Control Act (Title 22, U.S.C., Sec 2751, et seq.) or the
! Export Administration Act of 1979 (Title 50, U.S.C., App. 2401 et
! seq.), as amended. Violations of these export laws are subject to
! severe criminal penalties. Disseminate in accordance with provisions
! of DoD Directive 5230.25.
!
! ---------------------------------------------------------------------------
!
! WSF-M MWI Geolocation algorithm
module mwi_geolocation
  use, intrinsic :: iso_fortran_env, only: real32, real64, int8, int32, ERROR_UNIT
  use, intrinsic :: iso_c_binding, only: c_float, c_double, c_char, c_int, c_ptr, c_f_pointer, C_NULL_CHAR
  use dem, only: DigitalElevationMap
  use faraday_angle, only: compute_ionosphere_pierce_point, compute_faraday_angle, find_tec_ntcmgl
  use geolocation_core_routines, only: &
       EarthFigure, GeolocationCommonParameterData, GeolocationCommonInputData, &
       GeolocationScanlineData, GeolocationSatelliteFOV, GeolocationGroundFOV, &
       prepare_scanline, prepare_fov, geolocate_onto_ground, get_sc_axes, attitude_outside_range
  use geolocation_extra_routines, only: find_subsat_point, &
       intersects_earth, find_polrot_angle, convert_sc_to_ecef, find_cold_view_outputs
  use geomag, only: find_earth_magnetic_field
  use land_fraction, only: LandFractionMap
  use mwsdps_constants, only: nBand, nStokesBand, &
       nLonDEM, nLatDEM, nResLF, nLatLF, nLonLF, BAND_FREQS, ERRNO_OK, &
       GL_QUALITY_EXCLUSION_BIT, GL_QUALITY_DEGRADATION_BIT, GL_QUALITY_VALIDATION_BIT, &
       GL_EXCLUSION_NO_INTERSECTION_BIT, &
       GL_DEGRADATION_GEOMAG_EXPIRED_BIT, GL_DEGRADATION_INCIDENCE_OOB_BIT, GL_DEGRADATION_ALTITUDE_OOB_BIT, &
       GL_DEGRADATION_SPIN_OOB_BIT, GL_DEGRADATION_ATTITUDE_OOB_BIT, &
       GL_VALIDATION_LAND_BIT, GL_VALIDATION_COASTLINE_BIT, &
       CHECK_INCIDENCE_MIN, CHECK_INCIDENCE_MAX, &
       CHECK_ALTITUDE_MIN, CHECK_ALTITUDE_MAX, &
       CHECK_SPIN_MIN, CHECK_SPIN_MAX
        
  use geolocation_types, only: GeolocationParameterData, &
       GeolocationParameterData_f, &
       GeolocationAncillaryData, &
       GeolocationAncillaryData_f, &
       GeolocationSensorData, &
       GeolocationSensorData_f, &
       GeolocationOutputData, &
       GeolocationOutputData_f, &
       geolocation_sensor_convert, geolocation_output_convert
  use simulation_param, only: SimulationParameterData_f
  use mwsdps_logger, only: log_info, log_warning
  use pre_nut_routines, only: get_gm_angle

  use tdr_io, only: TypeTDR

  implicit none
  private
  public :: geolocation_init_f, geolocation_runCompute_f, set_gl_common_from_tdr,ellipsoid

  ! Private internal state
  type(GeolocationParameterData_f) :: parameters
  type(GeolocationAncillaryData_f) :: ancillary
  type(DigitalElevationMap) :: ancillary_dem
  type(LandFractionMap) :: ancillary_lf
  type(EarthFigure) :: ellipsoid

contains

  ! -----------------------------------------------------------
  function geolocation_init(parameter_data, ancillary_data) bind(C)
    integer(c_int) :: geolocation_init
    type(GeolocationParameterData), intent(in) :: parameter_data
    type(GeolocationAncillaryData), intent(in) :: ancillary_data

    real(real32), dimension(:), pointer :: one_dim
    real(real32), dimension(:, :), pointer :: two_dims
    real(real32), dimension(:, :, :), pointer :: three_dims

    ! Copy parameters to module-private data
    parameters%SamplingTime = parameter_data%SamplingTime
    parameters%EarthEquatorialRadius = parameter_data%EarthEquatorialRadius
    parameters%EarthPolarRadius = parameter_data%EarthPolarRadius
    parameters%SolarRadioFlux = parameter_data%SolarRadioFlux
    call c_f_pointer(parameter_data%OffNadirAngle, one_dim, [nBand])
    parameters%OffNadirAngle => one_dim(:)
    call c_f_pointer(parameter_data%PolarizationTilt, one_dim, [nStokesBand])
    parameters%PolarizationTilt => one_dim(:)
    call c_f_pointer(parameter_data%ColdSkyScanAngle, one_dim, [nBand])
    parameters%ColdSkyScanAngle => one_dim(:)

    ! Copy ancillary data to module-private data
    call c_f_pointer(ancillary_data%LatitudeDEM, one_dim, [nLatDEM])
    ancillary%LatitudeDEM => one_dim(:)
    call c_f_pointer(ancillary_data%LongitudeDEM, one_dim, [nLonDEM])
    ancillary%LongitudeDEM => one_dim(:)
    call c_f_pointer(ancillary_data%LatitudeLF, one_dim, [nLatLF])
    ancillary%LatitudeLf => one_dim(:)
    call c_f_pointer(ancillary_data%LongitudeLF, one_dim, [nLonLF])
    ancillary%LongitudeLF => one_dim(:)
    call c_f_pointer(ancillary_data%ElevationMean, two_dims, [nLonDEM, nLatDEM])
    ancillary%ElevationMean => two_dims(:, :)
    call c_f_pointer(ancillary_data%ElevationMin, two_dims, [nLonDEM, nLatDEM])
    ancillary%ElevationMin => two_dims(:, :)
    call c_f_pointer(ancillary_data%ElevationMax, two_dims, [nLonDEM, nLatDEM])
    ancillary%ElevationMax => two_dims(:, :)
    call c_f_pointer(ancillary_data%LandFraction, three_dims, [nLonLF, nLatLF, nResLF])
    ancillary%LandFraction => three_dims(:, :, :)

    ! Copy from the MWSDPS ancillary data to the internal DEM and land fraction data
    call init_dem_from_gl(ancillary, ancillary_dem)
    call init_lf_from_gl(ancillary, ancillary_lf)

    ! Store Earth figure as meters
    ellipsoid%EarthEquatorialRadius = real(parameters%EarthEquatorialRadius, real64) * 1d3
    ellipsoid%EarthPolarRadius = real(parameters%EarthPolarRadius, real64) * 1d3

    geolocation_init = ERRNO_OK
  end function geolocation_init

  ! -----------------------------------------------------------
  function geolocation_init_f(parameter_data, ancillary_data)
    integer(c_int) :: geolocation_init_f
    type(GeolocationParameterData_f), intent(in) :: parameter_data
    type(GeolocationAncillaryData_f), intent(in) :: ancillary_data

    ! Copy parameters and ancillary data to module-private data
    parameters = parameter_data
    ancillary = ancillary_data

    ! Copy from the MWSDPS ancillary data to the internal DEM and land fraction data
    call init_dem_from_gl(ancillary, ancillary_dem)
    call init_lf_from_gl(ancillary, ancillary_lf)

    ! Store Earth figure as meters
    ellipsoid%EarthEquatorialRadius = real(parameters%EarthEquatorialRadius, real64) 
    ellipsoid%EarthPolarRadius = real(parameters%EarthPolarRadius, real64) 

    geolocation_init_f = ERRNO_OK
  end function geolocation_init_f

  ! -----------------------------------------------------------
  function geolocation_runCompute(input_data, output_data) bind(C)
    integer(c_int) :: geolocation_runCompute
    type(GeolocationSensorData), intent(in) :: input_data
    type(GeolocationOutputData), intent(inout) :: output_data

    type(GeolocationSensorData_f) :: input_data_f
    type(GeolocationOutputData_f) :: output_data_f

    type(SimulationParameterData_f) :: parameter_data

    call geolocation_sensor_convert(input_data, input_data_f)
    call geolocation_output_convert(output_data, output_data_f)
    !this routine is not currently functional because imulationParameterData_f does not exist on CPP side of code
   
    geolocation_runCompute = geolocation_runCompute_f(parameter_data,input_data_f, output_data_f)
  end function geolocation_runCompute

  ! -----------------------------------------------------------
  function geolocation_runCompute_f(parameter_data,input_data, output_data)
    integer(c_int) :: geolocation_runCompute_f
    type(SimulationParameterData_f), intent(in) :: parameter_data
    type(GeolocationSensorData_f),   intent(in) :: input_data
    type(GeolocationOutputData_f),   intent(inout) :: output_data

    integer :: iscan, iband, ifov, ipolband, ires
    real(real64) :: gm_angle
    real(real64), dimension(3, 3) :: mics
    real(real64), dimension(3) :: sc_x, sc_y, sc_z
    real(real64), dimension(3) :: moonvec, sc_to_moon
    real(real64), dimension(3) :: sunvec, sc_to_sun
    real(real64) :: moon_distance, sun_distance
    logical :: found, mag_invalid, mag_expired
    real(real64) :: sc_lat, sc_lon, sc_alt
    real(real64) :: lat, lon, pra
    real(real64), dimension(3) :: ipp_bmag, ipp_loc, sc_mag
    real(real64) :: tec, faraday_rotation_angle

    type(GeolocationCommonParameterData) :: gl_parameters
    type(GeolocationCommonInputData) :: gl_input
    type(GeolocationScanlineData) :: scanline
    type(GeolocationSatelliteFOV) :: satellite_fov
    type(GeolocationGroundFOV) :: ground_fov

    character(len=200) :: log_message

    ! If the land fraction is greater than 50%, the IsLand validation flag is set
    real(real32), parameter :: LAND_THRESHOLD = 0.5
    ! If the land fraction is between 0.25% and 99.75%, the Coastline validation flag is set
    real(real32), parameter :: COASTLINE_THRESHOLD = 0.0025

    output_data%nScan = input_data%nScan
    output_data%nFOV = input_data%nFOV
    output_data%nBand = input_data%nBand
    output_data%nPolBand = input_data%nPolBand

    ! Prepare the two common types
    call set_gl_common(parameters, gl_parameters, input_data, gl_input)

    scans: do iscan=1, input_data%nScan
       write (log_message, '(A, " ", i0, "/", i0)') "Geolocating scan", iscan, input_data%nScan
       call log_info(trim(log_message))

       scanline = prepare_scanline(gl_input, iscan)

       ! Compute spacecraft parameters
       call get_gm_angle(scanline%time_j2000_ut1, gm_angle)
       call find_subsat_point(ellipsoid, scanline%scpos, gm_angle, sc_lat, sc_lon, sc_alt, found)
       if (found) then
          ! spacecraft geodetic latitude, longitude east of prime
          ! meridian, and altitude in meters
          output_data%SpacecraftLatitude(iscan) = real(sc_lat, real32)
          output_data%SpacecraftLongitude(iscan) = real(sc_lon, real32)
          output_data%SpacecraftAltitude(iscan) = real(sc_alt, real32)

          ! Check that the altitude is within range
          if (sc_alt < CHECK_ALTITUDE_MIN .or. sc_alt > CHECK_ALTITUDE_MAX) then
             write (log_message, '(A, " ", f0.2)') "Altitude out of range:", sc_alt
             call log_warning(trim(log_message))
             call set_degradation_flag(GL_DEGRADATION_ALTITUDE_OOB_BIT, &
                  output_data%DegradationConditionType(:, iscan, :), &
                  output_data%QualityFlag(:, iscan, :))
          end if
       else
          ! Spacecraft doesn't point to Earth! For all bands/FOVs for
          ! this scan, set the NoIntersection exclusion flag.
          call log_warning("No Earth intersection found")
          call set_exclusion_flag(GL_EXCLUSION_NO_INTERSECTION_BIT, &
               output_data%ExclusionConditionType(:, iscan, :), &
               output_data%QualityFlag(:, iscan, :))
       end if

       ! Find the SACS-S coordinate frame at the mean epoch of date
       ! (nominally oriented pointing toward geodetic nadir)
       call get_sc_axes(real(input_data%SpacecraftAttitude(iscan, :), real64), scanline%precession, &
            sc_x, sc_y, sc_z)

       ! This matrix is used to convert vectors from ECI to SACS-S
       ! (aka MWI instrument coordinate system, MICS)
       mics(1, :) = sc_x
       mics(2, :) = sc_y
       mics(3, :) = sc_z

       ! Warn if the attitude is too far from nominal
       if (attitude_outside_range(ellipsoid, scanline%scpos, scanline%scvel, sc_x, sc_z)) then
          ! The flag is set for every band and sample for this scan.
          ! Another check below is performed in case the attitude is
          ! out of range for a single band/sample.
          call log_warning("Attitude out of range")
          call set_degradation_flag(GL_DEGRADATION_ATTITUDE_OOB_BIT, &
               output_data%DegradationConditionType(:, iscan, :), &
               output_data%QualityFlag(:, iscan, :))
       end if

       ! Query the geomagnetic field model at the spacecraft location
       ! and transform it to the SACS-S frame
       call find_earth_magnetic_field(parameters, scanline%cdate, [sc_lat, sc_lon, sc_alt], &
            sc_mag, mag_invalid, mag_expired)
       if (mag_invalid) then
          call log_warning("Geomagnetic model inputs are invalid")
       end if
       if (mag_expired) then
          ! The magnetic field model is expired, so set the
          ! GeoMagModelExpired degradation condition
          call log_warning("Geomagnetic model expired")
          call set_degradation_flag(GL_DEGRADATION_GEOMAG_EXPIRED_BIT, &
               output_data%DegradationConditionType(:output_data%nFOV, iscan, :nBand), &
               output_data%QualityFlag(:output_data%nFOV, iscan, :nBand))
       end if
       output_data%ModelMagneticFieldAtSensor(iscan, :) = real(matmul(mics, sc_mag), real32)

       ! Convert lunar angle from ECI J2000 to mean epoch of date, and
       ! normalized in the SACS-S frame. Since the lunar vector is in
       ! meters and is large, to avoid precision loss, precession is
       ! applied only after it is normalized.
       moon_distance = norm2(input_data%LunarPositionECI(iscan, :))
       moonvec = input_data%LunarPositionECI(iscan, :) / moon_distance
       moonvec = matmul(scanline%precession, moonvec) * moon_distance
       sc_to_moon = moonvec - scanline%scpos
       sc_to_moon = sc_to_moon / norm2(sc_to_moon)
       output_data%LunarVector(iscan, :) = real(matmul(mics, sc_to_moon), real32)

       ! Convert solar vector to mean epoch of date, and normalized in
       ! the SACS-S frame. Since the solar vector is in meters and is
       ! quite large, to avoid precision loss, precession is applied
       ! only after it is normalized.
       sun_distance = norm2(input_data%SolarPositionECI(iscan, :))
       sunvec = input_data%SolarPositionECI(iscan, :) / sun_distance
       sunvec = matmul(scanline%precession, sunvec) * sun_distance
       sc_to_sun = sunvec - scanline%scpos
       sc_to_sun = sc_to_sun / norm2(sc_to_sun)
       output_data%SolarVector(iscan, :) = real(matmul(mics, sc_to_sun), real32)

       ! The solar eclipse state is detected when the direction from
       ! the spacecraft to the sun intersects the Earth ellipsoid
       if (intersects_earth(ellipsoid, scanline%scpos, sc_to_sun)) then
          output_data%SolarEclipseState(iscan) = int(1, int8)
       else
          output_data%SolarEclipseState(iscan) = int(0, int8)
       end if

       ! Check that the scan rate is within range
       if (input_data%ScanRate(iscan) < CHECK_SPIN_MIN .or. input_data%ScanRate(iscan) > CHECK_SPIN_MAX) then
          write (log_message, '(A, " ", f0.2)') "Scan rate out of range:", input_data%ScanRate(iscan)
          call log_warning(trim(log_message))
          call set_degradation_flag(GL_DEGRADATION_SPIN_OOB_BIT, &
               output_data%DegradationConditionType(:, iscan, :), &
               output_data%QualityFlag(:, iscan, :))
       end if

       fovs: do ifov = 1, output_data%nFOV
          bands: do iband = 1, input_data%nBand
             satellite_fov = prepare_fov(scanline, gl_parameters, gl_input, ifov, iband)

             ! Warn if the attitude is too far from nominal
             if (attitude_outside_range(ellipsoid, satellite_fov%scpos, scanline%scvel, &
                  satellite_fov%sc_x, satellite_fov%sc_z)) then
                call log_warning("Attitude out of range")
                call set_degradation_flag(GL_DEGRADATION_ATTITUDE_OOB_BIT, &
                     output_data%DegradationConditionType(ifov, iscan, iband), &
                     output_data%QualityFlag(ifov, iscan, iband))
             end if

             ground_fov = geolocate_onto_ground(satellite_fov, ellipsoid, ancillary_dem)
             if (.not. ground_fov%found) then
                ! No intersection, so set the NoIntersection exclusion
                ! flag and move to the next FOV
                call log_warning("No Earth intersection found")
                call set_exclusion_flag(GL_EXCLUSION_NO_INTERSECTION_BIT, &
                     output_data%ExclusionConditionType(ifov, iscan, iband), &
                     output_data%QualityFlag(ifov, iscan, iband))
                cycle
             end if

             output_data%Latitude(ifov, iscan, iband)=real(ground_fov%lat, real32)
             output_data%Longitude(ifov, iscan, iband)=real(ground_fov%lon, real32)
             ! Only the elevation is stored for the first band
             if (iband == 1) then
                output_data%Elevation(ifov, iscan) = real(ground_fov%alt, real32)
             end if

             output_data%EarthIncidenceAngle(ifov, iscan, iband)=real(ground_fov%inc, real32)
             output_data%EarthAzimuthAngle(ifov, iscan, iband)=real(ground_fov%azi, real32)

             ! Check that the incidence angle is within range
             if (ground_fov%inc < CHECK_INCIDENCE_MIN .or. ground_fov%inc > CHECK_INCIDENCE_MAX) then
                write (log_message, '(A, " ", f0.2)') "Incidence angle out of range:", ground_fov%inc
                call log_warning(trim(log_message))
                call set_degradation_flag(GL_DEGRADATION_INCIDENCE_OOB_BIT, &
                     output_data%DegradationConditionType(ifov, iscan, iband), &
                     output_data%QualityFlag(ifov, iscan, iband))
             end if

             ! These values are only computed/stored for the 10.85 GHz
             ! band
             if (iband == 1) then
                output_data%SolarAzimuthAngle(ifov, iscan)=real(ground_fov%sunazi, real32)
                output_data%SolarZenithAngle(ifov, iscan)=real(ground_fov%suninc, real32)
                output_data%SunGlintAngle(ifov, iscan)=real(ground_fov%sunglt, real32)

                ! Find the ionosphere pierce point and the magnetic field at
                ! that location
                call compute_ionosphere_pierce_point(parameters, scanline%cdate, satellite_fov%gm_angle, &
                     satellite_fov%scpos, sc_lat, satellite_fov%boresight, ipp_loc, ipp_bmag, found, mag_expired)
                if (.not. found) then
                   ! No intersection with the ionosphere, so set the
                   ! NoIntersection exclusion flag and move to the
                   ! next FOV. However, if we have an Earth
                   ! intersection above, then it should be impossible
                   ! to not hit the ionosphere...so this seems
                   ! unlikely to happen.
                   call set_exclusion_flag(GL_EXCLUSION_NO_INTERSECTION_BIT, &
                        output_data%ExclusionConditionType(ifov, iscan, :nBand), &
                        output_data%QualityFlag(ifov, iscan, :nBand))
                   cycle
                end if
                if (mag_expired) then
                   ! The magnetic field model is expired, so set the
                   ! GeoMagModelExpired degradation condition
                   call log_warning("Geomagnetic model expired")
                   call set_degradation_flag(GL_DEGRADATION_GEOMAG_EXPIRED_BIT, &
                        output_data%DegradationConditionType(ifov, iscan, :nBand), &
                        output_data%QualityFlag(ifov, iscan, :nBand))
                end if
                output_data%IonospherePiercePoint(ifov, iscan, :)=real(ipp_loc, real32)
                output_data%IonosphereMagneticField(ifov, iscan, :)=real(ipp_bmag, real32)

                ! TEC is in units of 1e16 electrons/m^2, aka a TECU
                call find_tec_ntcmgl(parameters, scanline%cdate, ipp_loc(1), ipp_loc(2), tec)

                ! Scale the vertically integrated TEC to account for
                ! the altitude of MWI being in the middle of the TEC column
                tec = tec * 0.75
                output_data%TotalElectronContent(ifov, iscan) = real(tec, real32)
             end if

             ! Land fraction map
             select case (iband)
             case (1)
                ! 10.85 GHz => 30 km resolution
                ires = 4
             case (2, 3)
                ! 18.85 GHz and 23.8 GHz => 20 km resolution
                ires = 2
             case default
                ! 36.75, 37.3, and 89 GHz => 15 km resolution
                ires = 1
             end select
             output_data%LandFraction(ifov, iscan, iband) = ancillary_lf%query(ires, ground_fov%lat, ground_fov%lon)

             ! Validation flag: IsLand is set if land fraction >= 50%
             if (output_data%LandFraction(ifov, iscan, iband) >= LAND_THRESHOLD) then
                call set_validation_flag(GL_VALIDATION_LAND_BIT, &
                     output_data%ValidationConditionType(ifov, iscan, iband), &
                     output_data%QualityFlag(ifov, iscan, iband))
             end if

             ! Validation flag: Coastline is set if land fraction between 0.25% and 99.75%
             if (output_data%LandFraction(ifov, iscan, iband) >= COASTLINE_THRESHOLD .and. &
                  output_data%LandFraction(ifov, iscan, iband) < 1. - COASTLINE_THRESHOLD) then
                call set_validation_flag(GL_VALIDATION_COASTLINE_BIT, &
                     output_data%ValidationConditionType(ifov, iscan, iband), &
                     output_data%QualityFlag(ifov, iscan, iband))
             end if

             if (output_data%nPolBand > 0) then
               ! only do this for polarimetric sensors?
               select case (iband)
                  case (1)
                     ! 10.85 GHz
                     ipolband = 1
                  case (2)
                     ! 18.85 GHz
                     ipolband = 2
                  case (4)
                     ! 36.75 GHz
                     ipolband = 3
                  case default
                     ! 23.8, 37.3, and 89 GHz
                     ipolband = 0
                  end select
               else
                  ipolband = 0
               endif
               

             ! These are only computed for fully polarimetric bands
             if (ipolband > 0) then
                call find_polrot_angle(satellite_fov%boresight, -satellite_fov%sc_z, &
                     ground_fov%surf_normal, ground_fov%inc < 0.01, pra)
                output_data%PolarizationRotationAngle(ifov, iscan, ipolband)=real(pra, real32) &
                     + parameters%PolarizationTilt(ipolband)

                call compute_faraday_angle(parameters, real(BAND_FREQS(iband), real64), &
                     satellite_fov%scpos, sc_lat, tec, ipp_bmag, satellite_fov%boresight, &
                     faraday_rotation_angle, found)
                if (found) then
                   output_data%FaradayRotationAngle(ifov, iscan, ipolband) = real(faraday_rotation_angle, real32)
                end if
             end if
          end do bands
       end do fovs

       ! For diagnostic purposes, output the satellite position,
       ! velocity, and attitude in ECEF coordinates at the mean epoch
       ! of date. (The precession matrix is already computed for this
       ! scan, but the GMST needs to be recomputed.)
       block
         real(real32), dimension(3) :: scpos_ecef, scvel_ecef
         real(real32), dimension(4) :: scatt_ecef
         real(real32), dimension(3) :: scpos_j2000, scvel_j2000
         real(real32), dimension(4) :: scatt_j2000

         ! These small arrays are needed to avoid runtime temporaries
         scpos_j2000(:) = input_data%SpacecraftPosition(iscan, :)
         scvel_j2000(:) = input_data%SpacecraftVelocity(iscan, :)
         scatt_j2000(:) = input_data%SpacecraftAttitude(iscan, :)

         call get_gm_angle(scanline%time_j2000_ut1, gm_angle)
         call convert_sc_to_ecef(scpos_j2000, scvel_j2000, scatt_j2000, &
              scanline%precession, gm_angle, &
              scpos_ecef, scvel_ecef, scatt_ecef)
         output_data%SpacecraftPositionECEF(iscan, :) = scpos_ecef(:)
         output_data%SpacecraftVelocityECEF(iscan, :) = scvel_ecef(:)
         output_data%SpacecraftAttitudeECEF(iscan, :) = scatt_ecef(:)
       end block

       do iband = 1,input_data%nBand
          call find_cold_view_outputs(scanline, parameters, ellipsoid, &
               input_data, iscan, iband, lat, lon)
          output_data%LatitudeColdView(iscan, iband) = real(lat, real32)
          output_data%LongitudeColdView(iscan, iband) = real(lon, real32)

          ! Look up the land fraction at the cold view location
          select case (iband)
          case (1)
             ! 10.85 GHz => 30 km resolution
             ires = 4
          case (2, 3)
             ! 18.85 GHz and 23.8 GHz => 20 km resolution
             ires = 2
          case default
             ! 36.75, 37.3, and 89 GHz => 15 km resolution
             ires = 1
          end select
          output_data%LandFractionColdView(iscan, iband) = ancillary_lf%query(ires, lat, lon)
       end do
    end do scans

    geolocation_runCompute_f = ERRNO_OK
  end function geolocation_runCompute_f

  ! -----------------------------------------------------------
  function geolocation_free() bind(C)
    integer(c_int) :: geolocation_free

    ! Deallocate module-private memory
    associate(a => ancillary_dem)
      deallocate(a%ElevationMean, a%ElevationMin, a%ElevationMax)
    end associate
    associate(a => ancillary_lf)
      deallocate(a%LandFraction, a%LatitudeLF, a%LongitudeLF)
    end associate

    geolocation_free = ERRNO_OK
  end function geolocation_free

  ! -----------------------------------------------------------
  ! Construct a DEM type by copying the geolocation ancillary data.
  subroutine init_dem_from_gl(anc_data, dem_data)
    type(GeolocationAncillaryData_f), intent(in) :: anc_data
    type(DigitalElevationMap), intent(inout) :: dem_data

    allocate(dem_data%ElevationMean, source=anc_data%ElevationMean)
    allocate(dem_data%ElevationMax, source=anc_data%ElevationMax)
    allocate(dem_data%ElevationMin, source=anc_data%ElevationMin)
  end subroutine init_dem_from_gl

  ! -----------------------------------------------------------
  ! Construct a LandFractionMap type by copying the geolocation ancillary data.
  subroutine init_lf_from_gl(anc_data, lf_data)
    type(GeolocationAncillaryData_f), intent(in) :: anc_data
    type(LandFractionMap), intent(inout) :: lf_data

    allocate(lf_data%LandFraction, source=anc_data%LandFraction)
    allocate(lf_data%LatitudeLF, source=anc_data%LatitudeLF)
    allocate(lf_data%LongitudeLF, source=anc_data%LongitudeLF)
  end subroutine init_lf_from_gl

  ! -----------------------------------------------------------
  ! Prepare the two geolocation common types. These contain a subset
  ! of the full parameters/input data for the geolocation algorithm.
  subroutine set_gl_common(params_full, params_sub, input_full, input_sub)
    type(GeolocationParameterData_f), intent(in), target :: params_full
    type(GeolocationCommonParameterData), intent(out) :: params_sub
    type(GeolocationSensorData_f), intent(in), target :: input_full
    type(GeolocationCommonInputData), intent(out) :: input_sub

    params_sub%SamplingTime = params_full%SamplingTime
    params_sub%OffNadirAngle => params_full%OffNadirAngle

    input_sub%nScan = input_full%nScan
    input_sub%nFOV = input_full%nFOV
    input_sub%nBand = input_full%nBand
    input_sub%SpacecraftPosition => input_full%SpacecraftPosition
    input_sub%SpacecraftVelocity => input_full%SpacecraftVelocity
    input_sub%SpacecraftAttitude => input_full%SpacecraftAttitude
    input_sub%SolarPositionECI => input_full%SolarPositionECI
    input_sub%ScanRate => input_full%ScanRate
    input_sub%ScanStartTimeUTC => input_full%ScanStartTimeUTC
    input_sub%LeapSeconds => input_full%LeapSeconds
    input_sub%DeltaUT1UTC => input_full%DeltaUT1UTC
    input_sub%FirstSampleDelay => input_full%FirstSampleDelay
    input_sub%FirstSampleScanAngle => input_full%FirstSampleScanAngle
  end subroutine set_gl_common

  subroutine set_gl_common_from_tdr(sim_params, params_sub, tdr_data, input_sub)

   

   type(SimulationParameterData_f), intent(in), target :: sim_params
   type(GeolocationCommonParameterData), intent(out) :: params_sub
   type(TypeTDR), intent(in), target :: tdr_data
   type(GeolocationCommonInputData), intent(out) :: input_sub

   params_sub%SamplingTime = sim_params%SamplingTime
   params_sub%OffNadirAngle => sim_params%OffNadirAngle

   input_sub%nScan = tdr_data%nScan
   input_sub%nFOV = tdr_data%nFOV
   input_sub%nBand = tdr_data%nBand
   input_sub%SpacecraftPosition => tdr_data%scpos
   input_sub%SpacecraftVelocity => tdr_data%scvel
   input_sub%SpacecraftAttitude => tdr_data%scatt
   input_sub%SolarPositionECI => tdr_data%Solar_Position
   input_sub%ScanRate => tdr_data%Scan_Rate
   input_sub%ScanStartTimeUTC => tdr_data%Scan_Time_UTC
   input_sub%LeapSeconds => tdr_data%Leap_Seconds
   input_sub%DeltaUT1UTC => tdr_data%Delta_UT1UTC
   input_sub%FirstSampleDelay => tdr_data%First_Sample_Delay
   input_sub%FirstSampleScanAngle => tdr_data%First_Sample_Scan_Angle
end subroutine set_gl_common_from_tdr



  ! -----------------------------------------------------------
  ! Set an exclusion condition flag.
  !
  ! exclusion_bit: the bit number (0 to 7) to set in exclusion_type
  ! exclusion_type: the exclusion type flag is updated
  ! quality_flag: the general quality flag is updated to indicate an exclusion condition
  !
  ! The valid exclusion_bit values are the GL_EXCLUSION_*_BIT constants in this module.
  elemental subroutine set_exclusion_flag(exclusion_bit, exclusion_type, quality_flag)
    integer(int8), intent(in) :: exclusion_bit
    integer(int8), intent(inout) :: quality_flag, exclusion_type

    quality_flag = ibset(quality_flag, GL_QUALITY_EXCLUSION_BIT)
    exclusion_type = ibset(exclusion_type, exclusion_bit)
  end subroutine set_exclusion_flag

  ! -----------------------------------------------------------
  ! Set a degradation condition flag.
  !
  ! degradation_bit: the bit number (0 to 7) to set in degradation_type
  ! degradation_type: the degradation type flag is updated
  ! quality_flag: the general quality flag is updated to indicate a degradation condition
  !
  ! The valid degradation_bit values are the GL_DEGRADATION_*_BIT constants in this module.
  elemental subroutine set_degradation_flag(degradation_bit, degradation_type, quality_flag)
    integer(int8), intent(in) :: degradation_bit
    integer(int8), intent(inout) :: quality_flag, degradation_type

    quality_flag = ibset(quality_flag, GL_QUALITY_DEGRADATION_BIT)
    degradation_type = ibset(degradation_type, degradation_bit)
  end subroutine set_degradation_flag

  ! -----------------------------------------------------------
  ! Set a validation condition flag.
  !
  ! validation_bit: the bit number (0 to 7) to set in validation_type
  ! validation_type: the validation type flag is updated
  ! quality_flag: the general quality flag is updated to indicate a validation condition
  !
  ! The valid validation_bit values are the GL_VALIDATION_*_BIT constants in this module.
  elemental subroutine set_validation_flag(validation_bit, validation_type, quality_flag)
    integer(int8), intent(in) :: validation_bit
    integer(int8), intent(inout) :: quality_flag, validation_type

    quality_flag = ibset(quality_flag, GL_QUALITY_VALIDATION_BIT)
    validation_type = ibset(validation_type, validation_bit)
  end subroutine set_validation_flag

end module mwi_geolocation
