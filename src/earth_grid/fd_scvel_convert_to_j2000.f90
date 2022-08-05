! Convert ECEF position/time samples to position/velocity/attitude in ECI J2000
  !
  ! The scan_time is in seconds since 2000-01-01 00:00:00 UTC, which
  ! is 12 hours earlier than the J2000 epoch.
  subroutine fd_scvel_convert_to_j2000(numscan, scan_time, scpos, &
    rpm, scpos_eci2000, scvel_eci2000, attitude)
    
    integer(int32), intent(in) :: numscan
    real(real64), dimension(:), intent(in) :: scan_time
    real(real64), dimension(:, :), intent(in) :: scpos
    real(real64), dimension(:), intent(out) :: rpm
    real(real64), dimension(:, :), intent(out) :: scpos_eci2000, scvel_eci2000, attitude

    real(real64), parameter :: d_time=2.d0 !arbitary, could have be 1 sec

    integer :: memstat
    integer(int32) :: iscan

    real(real64) :: rotangle, rotangle_step
    real(real64), dimension(3, 3) :: np_inv, rot_mat
    real(real64), dimension(3) :: scpos_step, scpos_step_eci2000
    real(real64), dimension(:, :), allocatable :: vel
    real(real64), dimension(3) :: x_source, y_source, z_source
    real(real64), dimension(3) :: x_dest, y_dest, z_dest
    real(real64) :: scan_time_ut1_s, scan_time_tt_day

    if (numscan < 2) error stop 'Too few scans'

    allocate(vel(3, size(scan_time)), stat=memstat)
    if (memstat /= 0) error stop 'Cannot allocate space for vel'

    ! ===============================================================================
    ! == compute earth-referenced velocity vector using postion vectors ==============
    ! The time is in seconds since an epoch without including leap
    ! seconds. So it's essentially TAI time with some offset.
    ! ===============================================================================
    rpm(1) = 60.d0/(scan_time(2)-scan_time(1))
    vel(:,1)=(scpos(:,2)-scpos(:,1))/(scan_time(2)-scan_time(1))
    do iscan=2,numscan-1
        rpm(iscan) = 60.d0/(scan_time(iscan+1)-scan_time(iscan))
        vel(:,iscan)=(scpos(:,iscan+1)-scpos(:,iscan-1))/(scan_time(iscan+1)-scan_time(iscan-1))
    end do
    rpm(numscan) = 60.d0/(scan_time(numscan)-scan_time(numscan-1))
    vel(:,numscan)=(scpos(:,numscan)-scpos(:,numscan-1))/(scan_time(numscan)-scan_time(numscan-1))

    scans: do iscan=1,numscan
        ! ======================================
        ! ===  convert time ====================
        ! ======================================

        ! The input times are SI seconds since 2000-01-01 00:00:00 UTC.
        ! We need to convert that to seconds since 2000-01-01 12:00:00
        ! UT1 for the GMST and to days since 2000-01-01 12:00:00 TT for
        ! the precession matrix.
        scan_time_ut1_s = time_2000_to_ut1(scan_time(iscan))
        scan_time_tt_day = time_2000_to_tt(scan_time(iscan))

        ! ======================================
        ! ===  convert to ECI J2000 system =====
        ! ======================================
        call get_gm_angle(scan_time_ut1_s, rotangle)
        ! Because the precession matrix is orthonormal, the inverse
        ! precession (converting from some epoch of mean date to J2000)
        ! is simply the matrix transpose
        call get_precession_matrix(scan_time_tt_day, np_inv)
        np_inv = transpose(np_inv)

        call get_gm_angle(scan_time_ut1_s + d_time, rotangle_step)

        scpos_step=scpos(:,iscan) + d_time*vel(:,iscan)

        call convert_earthfix_to_2000(rotangle,np_inv,scpos(:,iscan), scpos_eci2000(:,iscan))
        call convert_earthfix_to_2000(rotangle_step,np_inv,scpos_step, scpos_step_eci2000)

        scvel_eci2000(:,iscan)=(scpos_step_eci2000-scpos_eci2000(:,iscan))/d_time

        ! ======================================
        ! ===  determine attitude quaternion ===
        ! ======================================

        ! The "source" coordinate system is such that the X, Y, and Z
        ! axes are parallel to (but offset from) the ECI axes. The
        ! "destination" coordinate system has -Z pointing to geodetic
        ! nadir, X points along the flight direction, and Y is parallel
        ! to the orbit normal vector.
        !
        ! What we want is to determine the quaternion that rotates from
        ! "source" to "destination". All coordinates are in ECI J2000
        ! coordinates.
        call orient_frame_parallel_to_eci(x_source, y_source, z_source)
        call orient_frame_to_nadir(scpos_eci2000(:, iscan), scvel_eci2000(:, iscan), np_inv, &
            x_dest, y_dest, z_dest)

        ! The attitude quaternion is just the rotation between the two.
        ! So construct the rotation matrix and then extract the
        ! quaternion.
        call rotation_matrix_from_coords(x_source, y_source, z_source, &
            x_dest, y_dest, z_dest, rot_mat)
        attitude(:, iscan) = rotation_matrix_to_quaternion(rot_mat)
    end do scans
    deallocate(vel)
end subroutine fd_scvel_convert_to_j2000





