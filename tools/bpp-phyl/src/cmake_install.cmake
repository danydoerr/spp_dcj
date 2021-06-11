# Install script for directory: /home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/libbpp-phyl.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbpp-phyl.so.12.0.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbpp-phyl.so.12"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbpp-phyl.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/libbpp-phyl.so.12.0.0"
    "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/libbpp-phyl.so.12"
    "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/libbpp-phyl.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbpp-phyl.so.12.0.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbpp-phyl.so.12"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbpp-phyl.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/lib:"
           NEW_RPATH "")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE DIRECTORY FILES "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/Bpp" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/bpp-phyl/bpp-phyl-targets.cmake")
    file(DIFFERENT EXPORT_FILE_CHANGED FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/bpp-phyl/bpp-phyl-targets.cmake"
         "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/CMakeFiles/Export/lib/cmake/bpp-phyl/bpp-phyl-targets.cmake")
    if(EXPORT_FILE_CHANGED)
      file(GLOB OLD_CONFIG_FILES "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/bpp-phyl/bpp-phyl-targets-*.cmake")
      if(OLD_CONFIG_FILES)
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/bpp-phyl/bpp-phyl-targets.cmake\" will be replaced.  Removing files [${OLD_CONFIG_FILES}].")
        file(REMOVE ${OLD_CONFIG_FILES})
      endif()
    endif()
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/bpp-phyl" TYPE FILE FILES "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/CMakeFiles/Export/lib/cmake/bpp-phyl/bpp-phyl-targets.cmake")
  if("${CMAKE_INSTALL_CONFIG_NAME}" MATCHES "^([Rr][Ee][Ll][Ww][Ii][Tt][Hh][Dd][Ee][Bb][Ii][Nn][Ff][Oo])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/bpp-phyl" TYPE FILE FILES "/home/cchauve/ownCloud/DCJ-SPP/tools/Bio++/bpp-phyl/src/CMakeFiles/Export/lib/cmake/bpp-phyl/bpp-phyl-targets-relwithdebinfo.cmake")
  endif()
endif()

