module wgs84

    use, intrinsic :: iso_fortran_env, only: real64
    implicit none

    !all constants public by default
    public

!    earth dimensions match wgs84    
    real(real64),    parameter :: re=6378.137d3        !meters
    real(real64),    parameter :: rp=6356.752d3        !meters
    real(real64),    parameter :: ffac=(rp/re)*(rp/re)
    real(real64),    parameter :: geosync_alt   =42164.0d3      !satellite orbits, montenbruck and gill, page 294

contains


    real(4) function distance_between_lon_lats(lon1,lat1,lon2,lat2,error)
        use, intrinsic :: iso_fortran_env, only: real32, real64, int32, ERROR_UNIT   
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_quiet_nan,ieee_value
        use trig_degrees, only: sind, cosd

            real(4),intent(IN)            :: lon1
            real(4),intent(IN)            :: lat1
            real(4),intent(IN)            :: lon2
            real(4),intent(IN)            :: lat2
            integer(4),intent(INOUT)    :: error

            !real(4)                        :: EARTH_RADIUS = 6345.7
            !real(4)                        :: x1,y1,z1,x2,y2,z2
            real(4)                        :: nan_f32
            real(8)                        :: distance_64
            !real(4)                        :: distance_between_lon_lats3

            real(4)                        :: lon1_local,lon2_local

            nan_f32 = ieee_value(0.0_real32, ieee_quiet_nan)
            lon1_local = lon1
            lon2_local = lon2

            if (lon1_local .lt. 0) lon1_local = lon1_local + 360.0
            if (lon2_local .lt. 0) lon2_local = lon2_local + 360.0

            if((lon1_local .lt. 0.0) .or. (lon1_local .gt. 360.0) .or. &
               (lon2_local .lt. 0.0) .or. (lon2.gt. 360.0) .or. &
               (lat1 .lt. -90.0) .or. (lat1 .gt. 90.0) .or. &
               (lat2 .lt. -90.0) .or. (lat2 .gt. 90.0) .or. &
               (.not. (ieee_is_finite(lon1_local))) .or. &
               (.not. (ieee_is_finite(lon2_local))) .or. &
               (.not. (ieee_is_finite(lat1))) .or. &
               (.not. (ieee_is_finite(lat2)))) then

                error = -1
                distance_between_lon_lats = nan_f32
                return
            endif

            distance_64 = distance_between_lon_lats_real64(real(lon1_local,kind=real64), &
                                                                         real(lat1,kind=real64), &
                                                                         real(lon2_local,kind=real64), &
                                                                         real(lat2,kind=real64), &
                                                                         error)

            distance_between_lon_lats = real(distance_64,kind=real32)

        end function distance_between_lon_lats

        real(real64) function distance_between_lon_lats_real64(lon1,lat1,lon2,lat2,error)
            use, intrinsic :: iso_fortran_env, only: real32, real64, int32, ERROR_UNIT   
            use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_quiet_nan,ieee_value
            use trig_degrees, only: sind, cosd

            real(real64),intent(IN)            :: lon1
            real(real64),intent(IN)            :: lat1
            real(real64),intent(IN)            :: lon2
            real(real64),intent(IN)            :: lat2
            integer(4),intent(INOUT)    :: error

            real(real64)                        :: EARTH_RADIUS = 6345.7
            real(real64)                        :: x1,y1,z1,x2,y2,z2
            real(real64)                        :: cos_a
            real(real64)                        :: nan_f64
            real(real64)                        :: lon1_local,lon2_local


            nan_f64 = ieee_value(0.0_real64, ieee_quiet_nan)
            lon1_local = lon1
            lon2_local = lon2
            if (lon1_local .lt. 0) lon1_local = lon1_local + 360.0
            if (lon2_local .lt. 0) lon2_local = lon2_local + 360.0

           
            if((lon1_local .lt. 0.0) .or. (lon1_local .gt. 360.0) .or. &
            (lon2_local .lt. 0.0) .or. (lon2_local .gt. 360.0) .or. &
            (lat1 .lt. -90.0) .or. (lat1 .gt. 90.0) .or. &
            (lat2 .lt. -90.0) .or. (lat2 .gt. 90.0)) then

                error = -1
                distance_between_lon_lats_real64 = nan_f64
                return
            endif
        
            x1 = cosd(lat1)*cosd(lon1_local)
            y1 = cosd(lat1)*sind(lon1_local)
            z1 = sind(lat1)

            x2 = cosd(lat2)*cosd(lon2_local)
            y2 = cosd(lat2)*sind(lon2_local)
            z2 = sind(lat2)

            cos_a = x1*x2 + y1*y2 + z1*z2

            if (ieee_is_finite(cos_a)) then
                if (cos_a.gt.1.00) then
                    cos_a = 1.00      !fixes rounding errors since not using Real*8 
                endif
                distance_between_lon_lats_real64 = EARTH_RADIUS * acos(cos_a)
            else
                distance_between_lon_lats_real64 = nan_f64
                error = -1
            endif

        return
    end function distance_between_lon_lats_real64

    !
end module wgs84