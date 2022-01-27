#ifndef FILESYSTEM_HPP
#define FILESYSTEM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "type.hpp"
#include <fcntl.h>

#ifdef _WIN32
#include "unistdWindows.hpp"
#else
#include <unistd.h>
#endif


inline bool file_exists(std::string filename) {
  struct stat st;
  return stat(filename.c_str(), &st)==0;
}

inline long file_size(std::string filename) {
  struct stat st;
  assert(stat(filename.c_str(), &st)==0);
  return st.st_size;
}

#endif