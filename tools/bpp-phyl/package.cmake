# CMake package file for Bio++ PhylLib
# Authors:
#   Francois Gindraud (2017)
# Created: 08/03/2017

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was package.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

if (NOT bpp-phyl_FOUND)
  # Deps
  find_package (bpp-core 4.1.0 REQUIRED)
  find_package (bpp-seq 12.0.0 REQUIRED)
  # Add targets
  include ("${CMAKE_CURRENT_LIST_DIR}/bpp-phyl-targets.cmake")
  # Append targets to convenient lists
  set (BPP_LIBS_STATIC "${BPP_LIBS_STATIC}" bpp-phyl-static)
  set (BPP_LIBS_SHARED "${BPP_LIBS_SHARED}" bpp-phyl-shared)
  # Print some path info for targets
  get_property (static-location TARGET bpp-phyl-static PROPERTY LOCATION)
  get_property (shared-location TARGET bpp-phyl-shared PROPERTY LOCATION)
  get_property (header-location TARGET bpp-phyl-static PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  message (STATUS "bpp-phyl 12.0.0 found:")
  message (STATUS "  static lib: ${static-location}")
  message (STATUS "  shared lib: ${shared-location}")
  message (STATUS "  includes: ${header-location}")
  unset (static-location)
  unset (shared-location)
  unset (header-location)
endif (NOT bpp-phyl_FOUND)
