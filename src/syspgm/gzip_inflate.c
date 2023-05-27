// Small zlib wrapper to decompress a gzipped file in-memory
//
// This is meant to make things easier from the Fortran side, so
// Fortran doesn't have to worry about the gzFile struct and pointers
// and such.

#include <stdio.h>
#include <stdint.h>
#include <zlib.h>

// Public function prototypes
int32_t gzip_inflate(const char *fname, int8_t *buff, size_t buff_len, size_t *nbytes);

// Private function prototypes
static int uncompressed_size(const char *fname, size_t *nbytes);

// Open the file and decompress it into the given buffer
//
// fname: the gzipped file to open
// buff: a buffer of bytes that will be written
// buff_len: the length of the given buffer
// nbytes: the number of uncompressed bytes
//
// Return 0 if successful, 1 if failed (e.g., the buffer was too
// small for the file)
int32_t gzip_inflate(const char *fname, int8_t *buff, size_t buff_len, size_t *nbytes) {

  int err;
  gzFile file;

  // Read the uncompressed size from the file footer (assuming we have
  // a valid gzip file)
  err = uncompressed_size(fname, nbytes);
  if (err != 0 || buff_len < *nbytes)
	return 1;

  // Open the file and use a large (128 kiB) internal buffer
  file = gzopen(fname, "rb");
  if (file == NULL)
	return 1;
  err = gzbuffer(file, 1024 * 128);
  if (err != 0)
	return 1;

  // Decompress it into the provided buffer
  err = gzread(file, buff, (unsigned int) *nbytes);
  if (err < 0 || (size_t) err < *nbytes)
	return 1;

  // All done
  if (gzclose_r(file) != Z_OK)
	return 1;
  else
	return 0;
}


// A gzip-format file stores the length of the uncompressed data
// (modulo 2**32) at the end of the file
int uncompressed_size(const char *fname, size_t *nbytes) {
  FILE *file;

  uint32_t trail_size;

  file = fopen(fname, "rb");
  if (file == NULL)
	return 1;
  if (fseek(file, -4, SEEK_END) == -1)
	return 1;
  if (fread(&trail_size, sizeof(uint32_t), 1, file) != 1)
	return 1;
  if (fclose(file) != 0)
	return 1;

  *nbytes = trail_size;
  return 0;
}
