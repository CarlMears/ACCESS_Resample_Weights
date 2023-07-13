// Simple implementation of the MWSDPS logger interface
#include <stdio.h>

// Prototypes copied from MWSDPS
void clog_info(const char *message);
void clog_warning(const char *message);
void clog_error(const char *message);
void clog_test(const char *message);

// Implementations below
void clog_info(const char *message) {
  printf("INFO: %s\n", message);
}

void clog_warning(const char *message) {
  printf("WARNING: %s\n", message);
}

void clog_error(const char *message) {
  printf("ERROR: %s\n", message);
}

void clog_test(const char *message) {
  printf("TEST: %s\n", message);
}
