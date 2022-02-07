#---------------------------------------------------
# Molmodel 
#
# Creates SimTK library, base name=SimTKmolmodel.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d) and
# provision is made for an optional "namespace" (ns)
# and version number (vn).
#
# Windows:
#   [ns_]SimTKmolmodel[_vn][_d].dll
#   [ns_]SimTKmolmodel[_vn][_d].lib
#   [ns_]SimTKmolmodel[_vn]_static[_d].lib
# Unix:
#   lib[ns_]SimTKmolmodel[_vn][_d].so
#   lib[ns_]SimTKmolmodel[_vn]_static[_d].a
#
# All libraries are installed in 
#   %ProgramFiles%\SimTK\lib    (Windows)
#   /usr/local/SimTK/lib[64]    (Linux, Mac)
#
# Also creates an OpenMM plugin DLL that is used at
# runtime to determine whether OpenMM is available.
# That DLL is named
#   OpenMMPlugin[_d].dll
#   libOpenMMPlugin[_d].so
#   libOpenMMPlugin[_d].dylib
# And there is no static version.
#----------------------------------------------------

#cmake_minimum_required(VERSION 3.9.4)
cmake_minimum_required(VERSION 3.5.1)

project(Molmodel C CXX)

set(MOLMODEL_MAJOR_VERSION 3)
set(MOLMODEL_MINOR_VERSION 0)
set(MOLMODEL_PATCH_VERSION 0)

set(MOLMODEL_COPYRIGHT_YEARS "2006-12")

# underbar separated list of dotted authors, no spaces or commas
set(MOLMODEL_AUTHORS "Christopher.Bruns_Michael.Sherman")

# Report the version number to the CMake UI. Don't include the 
# build version if it is zero.
set(PATCH_VERSION_STRING)
IF(MOLMODEL_PATCH_VERSION)
    SET(PATCH_VERSION_STRING ".${MOLMODEL_PATCH_VERSION}")
ENDIF()

set(MOLMODEL_VERSION 
    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}${PATCH_VERSION_STRING}"
    CACHE STRING 
    "This is the version that will be built (can't be changed here)." 
    FORCE)

# This is the suffix if we're building and depending on versioned libraries.
set(VN "_${MOLMODEL_VERSION}")

set(BUILD_BINARY_DIR ${CMAKE_BINARY_DIR}
    CACHE PATH 
    "The Molmodel build (not the install) puts all the libraries and executables together here (with /Release, etc. appended on some platforms).")

# Make everything go in the same binary directory. (These are CMake-defined
# variables.)
set(EXECUTABLE_OUTPUT_PATH ${BUILD_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${BUILD_BINARY_DIR})

# Static libraries, tests, and examples won't be built unless this 
# is set.
SET(BUILD_STATIC_LIBRARIES FALSE CACHE BOOL
"Build '_static' versions of libraries in addition to dynamic libraries?")

# Use this to generate a private set of libraries whose names
# won't conflict with installed versions.
set(BUILD_USING_NAMESPACE "" CACHE STRING
	"All library names will be prefixed with 'xxx_' if this is set to xxx.")

set(BUILD_UNVERSIONED_LIBRARIES TRUE CACHE BOOL
 "Build library names, and assume dependency names, with no version numbers?")

set(BUILD_VERSIONED_LIBRARIES FALSE CACHE BOOL
 "Build library names, and assume dependency names, with version numbers?")

set(NS)
if(BUILD_USING_NAMESPACE)
    set(NS "${BUILD_USING_NAMESPACE}_")
endif()

set(MOLMODEL_LIBRARY_NAME ${NS}SimTKmolmodel CACHE STRING
"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
FORCE)


# Permit use of custom FindOpenMM and FindSimbody modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

set(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Control volume of build output" )

set(OpenMM_ON OFF)
find_package(OpenMM)
if(OpenMM_FOUND)
    set(OpenMM_ON ON)
endif()
set(MOLMODEL_USE_OpenMM ${OpenMM_ON} CACHE BOOL
    "Control whether Molmodel builds the OpenMM Plugin (requires that OpenMM is installed on the build machine)." )

# Caution: this variable is automatically created by the CMake
# ENABLE_TESTING() command, but we'll take it over here for
# our own purposes too.
set(BUILD_TESTING ON CACHE BOOL
	"Control building of Molmodel test programs." )

set(BUILD_EXAMPLES ON CACHE BOOL
	"Control building of Molmodel example programs." )

# Turning this off reduces the build time (and space) substantially,
# but you may miss the occasional odd bug. Also currently on Windows it
# is easier to debug the static tests than the DLL-liked ones.
SET( BUILD_TESTING_STATIC ON CACHE BOOL
    "If building static libraries, build static test and example programs too?" )

SET( BUILD_TESTING_SHARED ON CACHE BOOL
    "If building test or example programs, include dynamically-linked ones?" )

#
# Create a platform name useful for some platform-specific stuff.
IF(WIN32)
    SET(PLATFORM_NAME Windows)
    SET(NATIVE_COPY_CMD copy)
ELSEIF(APPLE)
    SET(PLATFORM_NAME Mac)
    SET(NATIVE_COPY_CMD cp)
ELSE()
    SET(PLATFORM_NAME Linux)
    SET(NATIVE_COPY_CMD cp)
ENDIF()

# In addition to the platform name we need to know the Application Binary
# Interface (ABI) we're building for. Currently that is either x86, meaning
# 32 bit Intel instruction set, or x64 for 64 bit Intel instruction set.

IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    SET(PLATFORM_ABI x64)
ELSE()
    SET(PLATFORM_ABI x86)
ENDIF()

SET(BUILD_PLATFORM "${PLATFORM_NAME}:${PLATFORM_ABI}" CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)

# 
# Make sure "sandbox" input & output directories are set. During
# SimTK Core build, the build system will set them. Otherwise, we'll
# set them to sensible local values.
# If SimTK_INSTALL_PREFIX is set then
# it is a sandbox installation area, otherwise we want to install
# in /usr/local/SimTK or %ProgramFiles%\SimTK. 
#
# As of Ubuntu 12, 64 bit libraries are in lib, not lib64, so we are 
# not setting this on any platform now. If you are building on an
# older system you might have to set it manually.
#
SET(LIB64) # suffix will be set only for 64 bit gcc builds
IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang") # Linux, Mac, Cygwin
#IF(${CMAKE_C_COMPILER} MATCHES "gcc") # Linux, Mac, Cygwin

    IF(NOT SimTK_INSTALL_PREFIX)
        SET(SimTK_INSTALL_PREFIX "/usr/local" 
            CACHE PATH "Install directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

    #IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
        #SET( LIB64 64 )
    #ENDIF()

ELSE() # Windows
    # On Win64, there are two environment variables, ProgramFiles(x86) for
    # 32 bit binaries, and ProgramW6432 for 64 bit binaries. The variable
    # ProgramFiles returns one or the other depending on whether the requesting
    # executable is 32 or 64 bits. That's CMake in this case, but we want
    # the right directory for the target. CMAKE_INSTALL_DIR after 2.8.3 gets
    # this right but still isn't what we want.
    IF( ${CMAKE_SIZEOF_VOID_P} EQUAL 8 )
	# 64 bit target on Win64
	set(PROGFILEDIR "$ENV{ProgramW6432}")
    ELSE() # Target is 32 bit
	set(PROGFILEDIR "$ENV{ProgramFiles(x86)}") # present if 64bit Windows
	if (NOT PROGFILEDIR)
	    set(PROGFILEDIR "$ENV{ProgramFiles}") # on 32bit Windows
	endif()
    ENDIF()

    IF(NOT SimTK_INSTALL_PREFIX)
	#HOREA
	MESSAGE( WARNING "I should take in consideration another default parameter for : SimTK_INSTALL_PREFIX" )
        SET(SimTK_INSTALL_PREFIX "${PROGFILEDIR}/SimTK" 
            CACHE PATH "Install directory")
    ENDIF(NOT SimTK_INSTALL_PREFIX)

ENDIF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")

# CMake installs go into CMAKE_INSTALL_PREFIX, which is always
# set to something incorrect by CMake. We'll use SimTK_INSTALL_PREFIX
# instead which is passed in by the build system or set above.
SET(CMAKE_INSTALL_PREFIX ${SimTK_INSTALL_PREFIX} CACHE STRING 
    "Install path prefix." FORCE)
# User should set SimTK_INSTALL_PREFIX, so suppress cmake's
mark_as_advanced(CMAKE_INSTALL_PREFIX)

# On Mac, build universal binaries that work on OS X 10.5 (Leopard) or later.
# TODO: we had problems using the OSX lapack for the 64-bit (x86_64) 
# architecture so we're only building the 32-bit (i386) architecture by
# default for now.
IF(APPLE)
    IF(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
	#SET(CMAKE_OSX_ARCHITECTURES "i386;x86_64" CACHE STRING "The processor architectures to build for" FORCE)
        SET(CMAKE_OSX_ARCHITECTURES "i386" CACHE STRING "The processor architectures to build for" FORCE)
        #10.5 sdk doesn't exist anymore on Mountain Lion. We let the system decide which is best.
        #SET(CMAKE_OSX_DEPLOYMENT_TARGET "10.5" CACHE STRING "The minimum version of OS X to support" FORCE)
        #SET(CMAKE_OSX_SYSROOT "/Developer/SDKs/MacOSX10.5.sdk" CACHE STRING "The SDK to build against" FORCE)
    ENDIF(NOT CMAKE_OSX_DEPLOYMENT_TARGET)
    SET(LIB64 )
ENDIF(APPLE)

IF(UNIX AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE Release CACHE STRING "Debug or Release build" FORCE)
ENDIF (UNIX AND NOT CMAKE_BUILD_TYPE)


# This will be initialized to the environment variable of the same name
# if it is set, otherwise it will be empty.
SET(ENV{SimTK_INSTALL_DIR} ${SimTK_INSTALL_PREFIX})
SET(SimTK_INSTALL_DIR $ENV{SimTK_INSTALL_DIR} CACHE
  PATH "Where to look for the Simbody libraries; you can change this.")

# This should find Simbody in SimTK_SDK if it is set, otherwise hunt for it.
find_package(Simbody REQUIRED) # CONFIG
include_directories(${Simbody_INCLUDE_DIR})
link_directories(${Simbody_LIB_DIR})

## Choose the maximum level of x86 instruction set that the compiler is 
## allowed to use. SSE2 is ubiquitous enough now that we don't mind
## abandoning machines that can't handle those instructions. SSE3 migh
## also be reasonable by now (April 2009) so this default should be
## revisited soon. This can be set to a different value by the person
## running CMake.
set(BUILD_INST_SET "sse2"  # use SSE2 instruction set by default
    CACHE STRING "CPU instruction level compiler is permitted to use.")
MARK_AS_ADVANCED( BUILD_INST_SET )

## When building in any of the Release modes, tell gcc to use full optimization and
## to generate SSE2 floating point instructions. Here we are specifying *all* of the
## Release flags, overriding CMake's defaults.
## Watch out for optimizer bugs in particular gcc versions!

IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")
    # # string(TOLOWER ${BUILD_INST_SET} GCC_INST_SET)
    # SET(GCC_INST_SET "-march=native -mtune=native") # use the native architecture, not some generic sse2 optimizations

    # # Get the gcc version number in major.minor.build format
    # execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
    #                 OUTPUT_VARIABLE GCC_VERSION)

    # # Unrolling fixed-count loops is a useful optimization for Simmatrix.
    # SET(GCC_OPT_ENABLE "-funroll-loops ${ROBO_EXTRA_FLAGS}")

    # # If you know of optimization bugs that affect SimTK in particular
    # # gcc versions, this is the place to turn off those optimizations.
    # SET(GCC_OPT_DISABLE)
    
    # # # We know Gcc 4.4.3 on Ubuntu 10 is buggy and that Snow Leopard's
    # # # 4.2.1 is not. To be safe for now we'll assume anything over 4.3
    # # # should have these disabled.
    # # if (GCC_VERSION VERSION_GREATER 4.3 OR GCC_VERSION VERSION_EQUAL 4.3)
    # #     SET(GCC_OPT_DISABLE 
	# # "-fno-tree-vrp -fno-strict-aliasing -fno-guess-branch-probability")
    # # endif()
    
    # # C++
    # SET(BUILD_CXX_FLAGS_DEBUG          "-g ${GCC_INST_SET}" 
    #   CACHE STRING "g++ Debug build compile flags")
    # SET(BUILD_CXX_FLAGS_RELEASE        
    #   "-DNDEBUG -O3 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}" 
    #   CACHE STRING "g++ Release build compile flags")
    # SET(BUILD_CXX_FLAGS_RELWITHDEBINFO 
    #   "-DNDEBUG -O3 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}" 
    #   CACHE STRING "g++ RelWithDebInfo build compile flags")
    # SET(BUILD_CXX_FLAGS_MINSIZEREL     "-DNDEBUG -Os ${GCC_INST_SET}" 
    #   CACHE STRING "g++ MinSizeRel build compile flags")

    # # C
    # SET(BUILD_C_FLAGS_DEBUG            "-g ${GCC_INST_SET}" 
    #   CACHE STRING "gcc Debug build compile flags")
    # SET(BUILD_C_FLAGS_RELEASE          
    #   "-DNDEBUG -O3 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}" 
    #   CACHE STRING "gcc Release build compile flags")
    # SET(BUILD_C_FLAGS_RELWITHDEBINFO   
    #   "-DNDEBUG -O3 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}" 
    #   CACHE STRING "gcc RelWithDebInfo build compile flags")
    # SET(BUILD_C_FLAGS_MINSIZEREL       "-DNDEBUG -Os ${GCC_INST_SET}" 
    #   CACHE STRING "gcc MinSizeRel build compile flags")

    # # C++
    # SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    # SET(CMAKE_CXX_FLAGS_DEBUG ${BUILD_CXX_FLAGS_DEBUG}
    #     CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    # SET(CMAKE_CXX_FLAGS_RELEASE ${BUILD_CXX_FLAGS_RELEASE}
    #     CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    # SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${BUILD_CXX_FLAGS_RELWITHDEBINFO}
    #     CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    # SET(CMAKE_CXX_FLAGS_MINSIZEREL ${BUILD_CXX_FLAGS_MINSIZEREL}
    #     CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    # # C
    # SET(CMAKE_C_FLAGS_DEBUG ${BUILD_C_FLAGS_DEBUG}
    #     CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    # SET(CMAKE_C_FLAGS_RELEASE ${BUILD_C_FLAGS_RELEASE}         
    #     CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    # SET(CMAKE_C_FLAGS_RELWITHDEBINFO ${BUILD_C_FLAGS_RELWITHDEBINFO}
    #     CACHE STRING "Can't change here -- see BUILD_C..." FORCE)
    # SET(CMAKE_C_FLAGS_MINSIZEREL ${BUILD_C_FLAGS_MINSIZEREL}
    #     CACHE STRING "Can't change here -- see BUILD_C..." FORCE)

ENDIF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")

## When building in any of the Release modes, tell VC++ cl compiler to use intrinsics 
## (i.e. sqrt instruction rather than sqrt subroutine) with flag /Oi.

IF(WIN32 AND ${CMAKE_C_COMPILER} MATCHES "cl")
    STRING(TOUPPER ${BUILD_INST_SET} CL_INST_SET)

    ## C++
    SET(BUILD_CXX_FLAGS_DEBUG        
	"/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Debug build compile flags")
    SET(BUILD_CXX_FLAGS_RELEASE        
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ Release build compile flags")
    SET(BUILD_CXX_FLAGS_RELWITHDEBINFO 
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ RelWithDebInfo build compile flags")
    SET(BUILD_CXX_FLAGS_MINSIZEREL 
	"/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC++ MinSizeRel build compile flags")

    ## C
    SET(BUILD_C_FLAGS_DEBUG        
	"/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC Debug build compile flags")
    SET(BUILD_C_FLAGS_RELEASE        
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC Release build compile flags")
    SET(BUILD_C_FLAGS_RELWITHDEBINFO 
	"/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC RelWithDebInfo build compile flags")
    SET(BUILD_C_FLAGS_MINSIZEREL 
	"/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- /arch:${CL_INST_SET}" 
        CACHE STRING "VC MinSizeRel build compile flags")

    ## C++
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    SET(CMAKE_CXX_FLAGS_DEBUG ${BUILD_CXX_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE ${BUILD_CXX_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${BUILD_CXX_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL ${BUILD_CXX_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)

    ## C
    SET(CMAKE_C_FLAGS_DEBUG ${BUILD_C_FLAGS_DEBUG}
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    SET(CMAKE_C_FLAGS_RELEASE ${BUILD_C_FLAGS_RELEASE}
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    SET(CMAKE_C_FLAGS_RELWITHDEBINFO ${BUILD_C_FLAGS_RELWITHDEBINFO}
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)
    SET(CMAKE_C_FLAGS_MINSIZEREL ${BUILD_C_FLAGS_MINSIZEREL}
        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)

ENDIF(WIN32 AND ${CMAKE_C_COMPILER} MATCHES "cl")

MARK_AS_ADVANCED(
    BUILD_CXX_FLAGS_DEBUG
    BUILD_CXX_FLAGS_RELEASE
    BUILD_CXX_FLAGS_RELWITHDEBINFO
    BUILD_CXX_FLAGS_MINSIZEREL
    BUILD_C_FLAGS_DEBUG
    BUILD_C_FLAGS_RELEASE
    BUILD_C_FLAGS_RELWITHDEBINFO
    BUILD_C_FLAGS_MINSIZEREL
)

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
set(MOLMODEL_SOURCE_SUBDIRS .)


# Collect up information about the version of the molmodel library we're building
# and make it available to the code so it can be built into the binaries.

# Get the subversion revision number if we can
# It's possible that WIN32 installs use svnversion through cygwin
# so we'll try for both svnversion.exe and svnversion. Note that
# this will result in warnings if all you have is Tortoise without
# Cygwin, and your "about" string will say "unknown" rather than
# providing the SVN version of the source.

# If there is an actual reason for the SVNVERSION vs. SVNVERSION_EXE
# nonsense that was here before, come talk to me.  --Chris Bruns
# 1) CMake already understands that Windows executables have ".exe" at the end.
# 2) Providing a PATHS argument to find_program is *NOT* hard coding a single path.
# Please read the CMake documentation for more details.
find_program(SVNVERSION_EXECUTABLE svnversion PATHS "C:/cygwin/bin")
if(SVNVERSION_EXECUTABLE)
    exec_program(${SVNVERSION_EXECUTABLE}
              # Works better on Windows to set working directory rather than
              # passing argument to svnversion
              ${CMAKE_SOURCE_DIR} # cwd for run
              OUTPUT_VARIABLE OUT)
    set(MOLMODEL_SVN_REVISION "${OUT}" CACHE STRING "Molmodel svn revision number" FORCE)
else(SVNVERSION_EXECUTABLE)
    message(STATUS 
         "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)"
        )
    set(MOLMODEL_SVN_REVISION "unknown" CACHE STRING "Molmodel svn revision number")
endif(SVNVERSION_EXECUTABLE)
mark_as_advanced(MOLMODEL_SVN_REVISION)

# Remove colon from build version, for easier placement in directory names
STRING(REPLACE ":" "_" MOLMODEL_SVN_REVISION ${MOLMODEL_SVN_REVISION})

ADD_DEFINITIONS(-DSimTK_MOLMODEL_LIBRARY_NAME=${MOLMODEL_LIBRARY_NAME}
                -DSimTK_MOLMODEL_MAJOR_VERSION=${MOLMODEL_MAJOR_VERSION}
                -DSimTK_MOLMODEL_MINOR_VERSION=${MOLMODEL_MINOR_VERSION}
		-DSimTK_MOLMODEL_PATCH_VERSION=${MOLMODEL_PATCH_VERSION})

# CMake quotes automatically when building Visual Studio projects but we need
# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes
# in Visual Studio which end up in the binary.

IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES FALSE)
ELSE (${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES TRUE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

##TODO: doesn't work without quotes in nightly build
set(NEED_QUOTES TRUE)

IF(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_MOLMODEL_SVN_REVISION="${MOLMODEL_SVN_REVISION}"
                   -DSimTK_MOLMODEL_COPYRIGHT_YEARS="${MOLMODEL_COPYRIGHT_YEARS}"
                   -DSimTK_MOLMODEL_AUTHORS="${MOLMODEL_AUTHORS}")
ELSE(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_MOLMODEL_SVN_REVISION=${MOLMODEL_SVN_REVISION}
                   -DSimTK_MOLMODEL_COPYRIGHT_YEARS=${MOLMODEL_COPYRIGHT_YEARS}
                   -DSimTK_MOLMODEL_AUTHORS=${MOLMODEL_AUTHORS})
ENDIF(NEED_QUOTES)

# -DSimTK_MOLMODEL_TYPE has to be defined in the target subdirectories.
# -Dmolmodel_EXPORTS defined automatically when Windows DLL build is being done.

set(SHARED_TARGET ${MOLMODEL_LIBRARY_NAME})
set(STATIC_TARGET ${MOLMODEL_LIBRARY_NAME}_static)
set(SHARED_TARGET_VN ${MOLMODEL_LIBRARY_NAME}${VN})
set(STATIC_TARGET_VN ${MOLMODEL_LIBRARY_NAME}${VN}_static)

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
    SET(STATIC_TARGET ${STATIC_TARGET}_d)
    SET(SHARED_TARGET_VN ${SHARED_TARGET_VN}_d)
    SET(STATIC_TARGET_VN ${STATIC_TARGET_VN}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

## Test against the unversioned libraries if they are being built;
## otherwise against the versioned libraries.
IF(BUILD_UNVERSIONED_LIBRARIES)
	SET(TEST_SHARED_TARGET ${SHARED_TARGET})
	SET(TEST_STATIC_TARGET ${STATIC_TARGET})
ELSE(BUILD_UNVERSIONED_LIBRARIES)
	SET(TEST_SHARED_TARGET ${SHARED_TARGET_VN})
	SET(TEST_STATIC_TARGET ${STATIC_TARGET_VN})
ENDIF(BUILD_UNVERSIONED_LIBRARIES)

IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")
    SET(ADDITIONAL_LINK_LIBRARIES)
ELSE(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")
    ## Assume Microsoft Visual Studio
	## ws2_32 is for socket calls in VMD connection
    SET(ADDITIONAL_LINK_LIBRARIES ws2_32)
ENDIF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")

# These are all the places to search for header files which are
# to be part of the API.
set(API_INCLUDE_DIRS) # start empty
FOREACH(subdir ${MOLMODEL_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${subdir}/include 
                         ${subdir}/include/molmodel 
                         ${subdir}/include/molmodel/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
set(API_REL_INCLUDE_FILES)   # start these out empty
set(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
set(SOURCE_FILES) # empty
set(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${MOLMODEL_SOURCE_SUBDIRS})
    FILE(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*.c
                         ${subdir}/src/*/*.cpp ${subdir}/src/*/*.c)
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h )
    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in
    ## SimTK/include if Molmodel was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/src)

#
# Allow automated build and dashboard.
#
INCLUDE (Dart)

IF (BUILD_TESTING)
    #IF (UNIX AND NOT CYGWIN AND NOT APPLE)
    #  IF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
    #    ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
    #    LINK_LIBRARIES(gcov)
    #  ENDIF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
    #ENDIF (UNIX AND NOT CYGWIN AND NOT APPLE)

    #
    # Testing
    #
    ENABLE_TESTING()

ENDIF (BUILD_TESTING)

INCLUDE(ApiDoxygen.cmake)

# libraries are installed from their subdirectories; headers here

# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
#HOREA this hack is to provide to gmodmodel access to some internal headers... practically... a hack
FILE(GLOB CORE_HEADERS_APP	   src/*h                       */src/*h )
FILE(GLOB TOP_HEADERS      include/molmodel/*.h          */include/molmodel/*.h)
FILE(GLOB INTERNAL_HEADERS include/molmodel/internal/*.h */include/molmodel/internal/*.h)
INSTALL_FILES(/include/                 FILES ${CORE_HEADERS})
#HOREA
INSTALL_FILES(/include/			FILES ${CORE_HEADERS_APP})
INSTALL_FILES(/include/molmodel/         FILES ${TOP_HEADERS})
INSTALL_FILES(/include/molmodel/internal FILES ${INTERNAL_HEADERS})

# Install documents.
FILE(GLOB TOPLEVEL_DOCS doc/*.pdf doc/*.txt)
INSTALL(FILES ${TOPLEVEL_DOCS} DESTINATION doc)

# These are at the end because we want them processed after
# all the various variables have been set above.

IF(BUILD_STATIC_LIBRARIES)
        ADD_SUBDIRECTORY( staticTarget )
ENDIF()
ADD_SUBDIRECTORY( sharedTarget )

IF (MOLMODEL_USE_OpenMM)
    ADD_SUBDIRECTORY( OpenMMPlugin )
ENDIF (MOLMODEL_USE_OpenMM)

IF( BUILD_EXAMPLES )
  ADD_SUBDIRECTORY( examples )
ENDIF()

IF( BUILD_TESTING )
  ADD_SUBDIRECTORY( tests )
ENDIF()


