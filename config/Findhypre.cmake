# Findhypre.cmake
#
# Input:
#    hypre_root: installation path
#
# Output:
#    hypre_FOUND       : TRUE if found, FALSE otherwise
#    hypre_INCLUDE_DIRS: include   paths
#    hypre_LIBRARIES   : libraries paths
#

# Not provide a default installation path yet
if(NOT hypre_root)
#  set(hypre_root ~/Software)
endif()

# find path of HYPRE.h and HYPRE_config.h
find_path(hypre_INCLUDE_DIRS 
          NAMES HYPRE.h HYPRE_config.h 
          PATHS ${hypre_root}/include 
          NO_DEFAULT_PATH)

# find path of libHYPRE.a
find_library(hypre_LIBRARIES
             NAMES libHYPRE.a
             PATHS ${hypre_root}/lib
             NO_DEFAULT_PATH)

# set hypre_FOUND
if(hypre_INCLUDE_DIRS AND hypre_LIBRARIES)
  set(hypre_FOUND TRUE)
else()
  set(hypre_FOUND FALSE)
  message(FATAL_ERROR "CANNOT find a valid hypre path.")
endif()
