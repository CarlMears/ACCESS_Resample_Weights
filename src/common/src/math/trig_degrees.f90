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
! Wrappers to call trigonometric functions in terms of degrees instead
! of radians
module trig_degrees
  use, intrinsic :: iso_fortran_env, only: real32, real64
  implicit none
  private
  public :: sind, cosd, tand, dsind, dcosd, dtand
  public :: asind, acosd, atand, dasind, dacosd, datand
  public :: atan2d, datan2d
  public :: RAD2DEG_F32, RAD2DEG_F64
  public :: DEG2RAD_F32, DEG2RAD_F64

  real(real32), parameter :: PI_F32 = 4.0_real32 * atan(1.0_real32)
  real(real64), parameter :: PI_F64 = 4.0_real64 * atan(1.0_real64)

  real(real32), parameter :: DEG2RAD_F32 = PI_F32 / 180.0_real32
  real(real64), parameter :: DEG2RAD_F64 = PI_F64 / 180.0_real64

  real(real32), parameter :: RAD2DEG_F32 = 180.0_real32 / PI_F32
  real(real64), parameter :: RAD2DEG_F64 = 180.0_real64 / PI_F64

  interface sind
    module procedure rsind, dsind
  end interface sind

  interface cosd
    module procedure rcosd, dcosd
  end interface cosd

  interface tand
    module procedure rtand, dtand
  end interface tand

  interface asind
    module procedure rasind, dasind
  end interface asind

  interface acosd
    module procedure racosd, dacosd
  end interface acosd

  interface atand
    module procedure ratand, datand
  end interface atand

  interface atan2d
    module procedure ratan2d, datan2d
  end interface atan2d

contains

  ! ------------------------------------------

  pure elemental function rsind(x)
    real(real32), intent(in) :: x
    real(real32) :: rsind
    rsind = sin(x * DEG2RAD_F32)
  end function rsind

  pure elemental function rcosd(x)
    real(real32), intent(in) :: x
    real(real32) :: rcosd
    rcosd = cos(x * DEG2RAD_F32)
  end function rcosd

  pure elemental function rtand(x)
    real(real32), intent(in) :: x
    real(real32) :: rtand
    rtand = tan(x * DEG2RAD_F32)
  end function rtand

  ! ------------------------------------------

  pure elemental function dsind(x)
    real(real64), intent(in) :: x
    real(real64) :: dsind
    dsind = sin(x * DEG2RAD_F64)
  end function dsind

  pure elemental function dcosd(x)
    real(real64), intent(in) :: x
    real(real64) :: dcosd
    dcosd = cos(x * DEG2RAD_F64)
  end function dcosd

  pure elemental function dtand(x)
    real(real64), intent(in) :: x
    real(real64) :: dtand
    dtand = tan(x * DEG2RAD_F64)
  end function dtand

  ! ------------------------------------------

  pure elemental function rasind(x)
    real(real32), intent(in) :: x
    real(real32) :: rasind
    rasind = RAD2DEG_F32 * asin(x)
  end function rasind

  pure elemental function racosd(x)
    real(real32), intent(in) :: x
    real(real32) :: racosd
    racosd = RAD2DEG_F32 * acos(x)
  end function racosd

  pure elemental function ratand(x)
    real(real32), intent(in) :: x
    real(real32) :: ratand
    ratand = RAD2DEG_F32 * atan(x)
  end function ratand

  ! ------------------------------------------

  pure elemental function dasind(x)
    real(real64), intent(in) :: x
    real(real64) :: dasind
    dasind = RAD2DEG_F64 * asin(x)
  end function dasind

  pure elemental function dacosd(x)
    real(real64), intent(in) :: x
    real(real64) :: dacosd
    dacosd = RAD2DEG_F64 * acos(x)
  end function dacosd

  pure elemental function datand(x)
    real(real64), intent(in) :: x
    real(real64) :: datand
    datand = RAD2DEG_F64 * atan(x)
  end function datand

  ! ------------------------------------------

  pure elemental function ratan2d(x, y)
    real(real32), intent(in) :: x, y
    real(real32) :: ratan2d
    ratan2d = RAD2DEG_F32 * atan2(x, y)
  end function ratan2d

  pure elemental function datan2d(x, y)
    real(real64), intent(in) :: x, y
    real(real64) :: datan2d
    datan2d = RAD2DEG_F64 * atan2(x, y)
  end function datan2d

  ! ------------------------------------------

end module trig_degrees
