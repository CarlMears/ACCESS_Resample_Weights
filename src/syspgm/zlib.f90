! Fortran wrapper to zlib

module zlib_wrapper
  use, intrinsic :: iso_c_binding, only: c_char, c_int, C_NULL_CHAR, c_size_t, c_sizeof
  use, intrinsic :: iso_fortran_env, only: int32, int8
  implicit none
  private
  public :: fastunzip

  ! From gzip_inflate.c, a small zlib wrapper
  interface
     ! int32_t gzip_inflate(const char *fname, int8_t *buff, size_t buff_len, size_t *nbytes);
     function gzip_inflate(fname, buff, buff_len, nbytes) bind(c)
       use, intrinsic ::  iso_c_binding
       use, intrinsic ::  iso_fortran_env
       implicit none
       character(c_char), dimension(*), intent(in) :: fname
       ! integer(c_int8_t), dimension(*), intent(out) :: buff
       character(c_char), dimension(*), intent(out) :: buff
       integer(c_size_t), value, intent(in) :: buff_len
       integer(c_size_t), intent(out) :: nbytes
       integer(c_int32_t) :: gzip_inflate
     end function gzip_inflate
  end interface

contains

  ! Decompress the gzipped file in-memory
  !
  ! filename: path to gzipped file
  ! isize: length of uncompressed file (in bytes)
  ! abuf: decompressed data
  ! ierr: return code
  !
  ! Note this is not a very general-purpose function since it's tuned
  ! for the the L1A buffer size/shape.
  subroutine fastunzip(filename, isize, abuf, ierr)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: isize
    integer(int32), intent(out) :: ierr
    integer, parameter :: reclen = 39816
    character(len=reclen), dimension(:), intent(inout) :: abuf

    integer(c_size_t) :: nbytes, abuf_bytes
    abuf_bytes = size(abuf, 1, c_size_t) * storage_size(abuf) / 8

    ierr = gzip_inflate(trim(filename) // C_NULL_CHAR, abuf, abuf_bytes, nbytes)
    isize = int(nbytes, int32)
  end subroutine fastunzip
end module zlib_wrapper
