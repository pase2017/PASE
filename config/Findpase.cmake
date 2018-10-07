# Findpase.cmake
#
# Input:
#    pase_root: installation path
#    use
#
# Output:
#    pase_FOUND       : TRUE if found, FALSE otherwise
#    pase_INCLUDE_DIRS: include   paths
#    pase_LIBRARIES   : libraries paths
#

# Not provide a default installation path yet
if(NOT pase_root)
#  set(pase_root ~/Software)
endif()

# find path of pase_config.h
find_path(pase_INCLUDE_DIRS 
          NAMES pase_config.h 
          PATHS ${pase_root}/include 
          NO_DEFAULT_PATH)

# find path of libpase.a
if(use_hypre)
  list(APPEND pase_lib_all "hypre")
endif()
list(APPEND pase_lib_all "kernel")

foreach(lib ${pase_lib_all})
  find_library(pase_lib_${lib}
               NAMES libpase_${lib}.a
               PATHS ${pase_root}/lib
               NO_DEFAULT_PATH)
  list(APPEND pase_LIBRARIES ${pase_lib_${lib}})
endforeach()

# set pase_FOUND
if(pase_INCLUDE_DIRS AND pase_LIBRARIES)
  set(pase_FOUND TRUE)
else()
  set(pase_FOUND FALSE)
  message(FATAL_ERROR "CANNOT find a valid pase path.")
endif()
