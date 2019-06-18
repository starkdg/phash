cmake_minimum_required (VERSION 3.13)

include(ExternalProject)

project(pHash)

set(CMAKE_PROJECT_VERSION 0.9.7)

set(CMAKE_BUILD_TYPE Release CACHE STRING "build type" FORCE)
set(USE_IMAGE_HASH ON CACHE BOOLEAN "compile image hashes")
set(USE_AUDIO_HASH ON CACHE BOOLEAN "compile audio hashes")
set(USE_VIDEO_HASH ON CACHE BOOLEAN "compile video hash")
set(USE_TEXT_HASH ON CACHE BOOLEAN "compile text hash")
set(USE_PTHREAD ON CACHE BOOLEAN "use pthread library")
set(USE_MPG123 ON CACHE BOOLEAN "use mpg123 library")

configure_file(src/pHash.h.cmake src/pHash.h)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)

add_subdirectory(bindings)

if (USE_IMAGE_HASH)
  find_library(libpng png)
  if (libpng-NOTFOUND)
	message(SEND_ERROR "no libpng found")
  endif()
  find_library(libjpeg jpeg)
  if (libjpeg-NOTFOUND)
	message(SEND_ERROR "no libjpeg found")
  endif()
  find_library(libtiff tiff)
  if (libtiff-NOTFOUND)
	message(SEND_ERROR "no libtiff found")
  endif()
  find_library(libpthread pthread)
  if (libpthread-NOTFOUND)
	message(FATAL_ERROR "no pthread library found")
  endif()
  find_library(librt rt)
  if (librt-NOTFOUND)
	message(FATAL_ERROR "no rt library found")
  endif()

  list(APPEND LIB_SRCS src/pHash.cpp)
  list(APPEND LIB_DEPS png jpeg tiff pthread rt)
endif()

if (USE_AUDIO_HASH)
  find_library(libsnd sndfile)
  if (libsnd-NOTFOUND)
	message(FATAL_ERROR "no libsndfile found")
  endif()
  find_library(libsamplerate samplerate)
  if (libsamplerate-NOTFOUND)
	message(FATAL_ERROR "no libsamplerate found")
  endif()
  find_library(libmpg123 mpg123)
  if (libmpg123-NOTFOUND)
	message(FATAL_ERROR "no libmpg123 found")
  endif()
  find_library(libvorbis vorbis)
  if (libvorbis-NOTFOUND)
	message(FATAL_ERROR "no libvorbis found")
  endif()
  find_library(libogg ogg)
  if (libogg-NOTFOUND)
	message(FATAL_ERROR "no libogg found")
  endif()
  
  list(APPEND LIB_SRCS "src/audiophash.cpp")
  list(APPEND LIB_DEPS ${libsnd} ${libsamplerate} ${libmpg123} ${libvorbis} ${libogg})
endif()


add_library(pHash SHARED ${LIB_SRCS})
add_library(pHash-static STATIC ${LIB_SRCS})

target_link_libraries(pHash ${LIB_DEPS})
target_link_libraries(pHash	${LIBS_DEPS})


install(TARGETS pHash DESTINATION lib)
install(TARGETS pHash-static DESTINATION lib)
install(FILES src/pHash.h DESTINATION include)


if(HAVE_IMAGE_HASH)
add_executable(TestDCT "examples/test_imagephash.cpp")
target_link_libraries (TestDCT pHash)


add_executable(TestBMB "examples/test_bmbimagehash.cpp")
target_link_libraries (TestBMB pHash)


add_executable(TestMH "examples/test_mhimagehash.cpp")
target_link_libraries (TestMH pHash)
endif()



if(HAVE_AUDIO_HASH)
add_executable(TestAudio "examples/test_audiophash.cpp")
target_link_libraries (TestAudio pHash)
endif()

if(HAVE_VIDEO_HASH)
add_executable(TestVideoHash "examples/test_dctvideohash.cpp")
target_link_libraries (TestVideoHash pHash)
endif()


include (InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set (CPACK_PACKAGE_VERSION_MAJOR "${pHash_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${pHash_VERSION_MINOR}")
include (CPack)
