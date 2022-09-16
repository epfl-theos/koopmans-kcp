cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
# CMake < 3.13 cannot install external targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2152

# CMake < 3.14 generates incorrect dependency graphs with
# alias targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2521

# CMake < 3.20 does not support NVHPC compiler id
# https://cmake.org/cmake/help/latest/release/3.20.html#compilers

# Policy CMP0048: The project() command manages VERSION variables
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

project(koopmans-kcp
    VERSION 0.1
    DESCRIPTION "Implementation of Koopmans functionals with full orbital optimization in Quantum ESPRESSO v4.1"
    LANGUAGES Fortran C)

if(${qe_BINARY_DIR} STREQUAL ${qe_SOURCE_DIR})
    message(FATAL_ERROR "QE source folder cannot be safely used as a build folder!")
endif()

# CMake < v3.18 cannot discover the ARM Performance Library
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
    if(CMAKE_VERSION VERSION_LESS 3.18.0)
        message("-- CMake versions less than 3.18 cannot automatically discover the ARM Performance Library!")
    endif()
endif()

##########################################################
# Define the paths for static libraries and executables
##########################################################
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${qe_BINARY_DIR}/lib 
    CACHE
    PATH "Single output directory for building all libraries.")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${qe_BINARY_DIR}/bin 
    CACHE
    PATH "Single output directory for building all executables.")

###########################################################
# Build helpers
###########################################################
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include(cmake/kcpHelpers.cmake)

###########################################################
# Build Type
# Ensure that a specific, default build type is set when
# none has been explicitly set by the user
###########################################################
qe_ensure_build_type("Release")

###########################################################
# Modules
###########################################################
include(CheckFunctionExists)
# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

###########################################################
# Build Options
###########################################################
include(CMakeDependentOption)
option(QE_ENABLE_MPI
    "enable distributed execution support via MPI" ON)
option(QE_ENABLE_TRACE
    "enable execution tracing output" OFF)
option(QE_ENABLE_MPI_INPLACE
    "enable inplace MPI calls (ignored when QE_ENABLE_MPI=OFF)" OFF)
option(QE_ENABLE_MPI_MODULE
    "use MPI via Fortran module instead of mpif.h header inclusion" OFF)
option(QE_ENABLE_BARRIER
    "enable global synchronization between execution units" OFF)
option(QE_LAPACK_INTERNAL
    "enable internal reference LAPACK" OFF)
option(QE_ENABLE_SCALAPACK
    "enable SCALAPACK execution units" OFF)
cmake_dependent_option(QE_ENABLE_SCALAPACK_QRCP
    "enable SCALAPACK QRCP in pw2wannier90 (requires SCALAPACK>=2.1.0 or Intel MKL>=2020)"
    OFF "QE_ENABLE_SCALAPACK" OFF)
option(QE_ENABLE_STATIC_BUILD
    "enable fully static build of executables" OFF)
set(QE_FFTW_VENDOR "AUTO" CACHE
    STRING "select a specific FFTW library [Intel_DFTI, Intel_FFTW3, ArmPL, IBMESSL, FFTW3, Internal]")
set(QE_ENABLE_SANITIZER "none" CACHE STRING "none,asan,ubsan,tsan,msan")

# TODO change all ifdefs throughout code base to match
# cmake options
# TODO symbols beginning with '__' followed by a capital
# character are reserved for standard library use (at
# least in C, not sure about Fortran), change all feature
# macros to avoid weird behaviors

# Disable all configuration headers used to be generated
# by configure (see <qe>/include/)
qe_add_global_compile_definitions(QE_NO_CONFIG_H)

if(QE_ENABLE_TRACE)
    qe_add_global_compile_definitions(__TRACE)
endif()
if(QE_ENABLE_MPI_INPLACE)
    qe_add_global_compile_definitions(__USE_INPLACE_MPI)
endif()
if(QE_ENABLE_MPI_MODULE)
    qe_add_global_compile_definitions(__MPI_MODULE)
endif()
if(QE_ENABLE_BARRIER)
    qe_add_global_compile_definitions(__USE_BARRIER)
endif()
if(QE_ENABLE_MPI)
    # OMPI_SKIP_MPICXX: skip CXX APIs on openmpi, cause trouble to C APIs
    qe_add_global_compile_definitions(__MPI OMPI_SKIP_MPICXX)
    if(QE_ENABLE_MPI_GPU_AWARE)
        qe_add_global_compile_definitions(__GPU_MPI)
    endif()
endif()
if(QE_ENABLE_SCALAPACK)
    qe_add_global_compile_definitions(__SCALAPACK)
endif()

# Feature checks
check_function_exists(mallinfo HAVE_MALLINFO)
if(HAVE_MALLINFO)
    qe_add_global_compile_definitions(HAVE_MALLINFO)
endif()

# Check options consistency
if(QE_ENABLE_STATIC_BUILD AND BUILD_SHARED_LIBS)
    message(FATAL_ERROR "Full static build of QE executables requires static QE internal libraries. QE_ENABLE_STATIC_BUILD and BUILD_SHARED_LIBS cannot be both ON")
endif()
if(QE_ENABLE_SCALAPACK AND NOT QE_ENABLE_MPI)
    message(FATAL_ERROR "SCALAPACK requires MPI support, enable it with '-DQE_ENABLE_MPI=ON' or disable SCALAPACK with '-DQE_ENABLE_SCALAPACK=OFF'")
endif()
if(QE_ENABLE_MPI_GPU_AWARE AND NOT (QE_ENABLE_CUDA AND QE_ENABLE_MPI))
    message(FATAL_ERROR "GPU aware MPI requires both MPI and CUDA features enabled")
endif()
# if(QE_ENABLE_HDF5 AND NOT QE_ENABLE_MPI)
#    message(FATAL_ERROR "HDF5 requires MPI support, enable it with '-DQE_ENABLE_MPI=ON' or disable HDF5 with '-DQE_ENABLE_HDF5=OFF'")
# endif()

# Add optional sanitizers ASAN, UBSAN, MSAN
set(VALID_SANITIZERS "none" "asan" "ubsan" "tsan" "msan")
# Perform sanitizer option check, only works in debug mode
if(NOT QE_ENABLE_SANITIZER IN_LIST VALID_SANITIZERS)
  message(FATAL_ERROR "Invalid -DQE_ENABLE_SANITIZER=${QE_ENABLE_SANITIZER}, value must be one of ${VALID_SANITIZERS}")
else()
  message(STATUS "Enable sanitizer QE_ENABLE_SANITIZER=${QE_ENABLE_SANITIZER}")
endif()
# only GNU works right now
if(NOT QE_ENABLE_SANITIZER STREQUAL "none" AND NOT CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  message(FATAL_ERROR "-DQE_ENABLE_SANITIZER=${QE_ENABLE_SANITIZER} only works with the GNU compiler")
endif()

############################################################
# C preprocessor
# Note: reply on the compiler preprocessor whenever possible
############################################################
if(DEFINED ENV{CPP})
    set(QE_CPP_DEFAULT $ENV{CPP})
else()
    set(QE_CPP_DEFAULT cpp)
endif()
# QE_CPP_DEFAULT is only effective when cached QE_CPP doesn't exist.
set(QE_CPP ${QE_CPP_DEFAULT} CACHE
    STRING "C preprocessor for qe_preprocess_source in qeHelpers.cmake")
find_program(QE_CPP_FULL_PATH NAMES ${QE_CPP} DOC "C preprocessor full path")
if(QE_CPP_FULL_PATH)
    message(STATUS "C preprocessor used by qe_preprocess_source in qeHelpers.cmake: ${QE_CPP_FULL_PATH}")
else()
    set(QE_CPP_SAVED ${QE_CPP})
    unset(QE_CPP CACHE)
    message(FATAL_ERROR "C preprocessor ${QE_CPP_SAVED} not found. Pass a working one to CMake via QE_CPP!")
endif()

###########################################################
# language standard requirements
###########################################################
# TODO need to require all compilers using the same one
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    set(CMAKE_C_STANDARD 11)
    set(CMAKE_C_STANDARD_REQUIRED ON)
    set(CMAKE_C_EXTENSIONS OFF)
endif()

###########################################################
# check Fortran compiler -isystem option support
###########################################################
include(CheckFortranCompilerFlag)
check_fortran_compiler_flag("-isystem ." Fortran_ISYSTEM_SUPPORTED)
if(NOT Fortran_ISYSTEM_SUPPORTED AND NOT DEFINED CMAKE_NO_SYSTEM_FROM_IMPORTED)
  set(CMAKE_NO_SYSTEM_FROM_IMPORTED ON)
endif()

############################################################
# Compiler vendor specific options
############################################################
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    include(GNUFortranCompiler)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
    include(NVFortranCompiler)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
    include(IBMFortranCompiler)
endif()

if(QE_ENABLE_STATIC_BUILD)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
endif()

###########################################################
# OpenMP
# The following targets will be defined:
add_library(qe_openmp_fortran INTERFACE)
add_library(qe_openmp_c INTERFACE)
qe_install_targets(qe_openmp_fortran qe_openmp_c)
###########################################################
if(QE_ENABLE_OPENMP)
    find_package(OpenMP REQUIRED Fortran C)
    target_link_libraries(qe_openmp_fortran INTERFACE OpenMP::OpenMP_Fortran)
    target_link_libraries(qe_openmp_c INTERFACE OpenMP::OpenMP_C)
endif(QE_ENABLE_OPENMP)

###########################################################
# MPI
# The following targets will be defined:
add_library(qe_mpi_fortran INTERFACE)
qe_install_targets(qe_mpi_fortran)
###########################################################
if(QE_ENABLE_MPI)
    find_package(MPI REQUIRED Fortran)
    target_link_libraries(qe_mpi_fortran
        INTERFACE MPI::MPI_Fortran)
    message(STATUS "MPI settings used by CTest")
    message("     MPIEXEC_EXECUTABLE : ${MPIEXEC_EXECUTABLE}")
    message("     MPIEXEC_NUMPROC_FLAG : ${MPIEXEC_NUMPROC_FLAG}")
    message("     MPIEXEC_PREFLAGS : ${MPIEXEC_PREFLAGS}")
    string(REPLACE ";" " " MPIEXEC_PREFLAGS_PRINT "${MPIEXEC_PREFLAGS}")
    message("   Tests run as : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} <NUM_PROCS> ${MPIEXEC_PREFLAGS_PRINT} <EXECUTABLE>")
endif(QE_ENABLE_MPI)

###########################################################
# Git
###########################################################
find_package(Git 2.13 REQUIRED)
if(EXISTS ${qe_SOURCE_DIR}/.git)
  message(STATUS "Source files are cloned from a git repository.")
  set(IS_GIT_PROJECT 1)
  include(GitInfo)
else()
  message(STATUS "Source files are not cloned from a git repository.")
endif()

###########################################################
# Lapack
# The following targets will be defined:
add_library(qe_lapack INTERFACE)
qe_install_targets(qe_lapack)
###########################################################
if(NOT QE_LAPACK_INTERNAL)
    if(NOT BLA_VENDOR)
        if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64.*")
            message(STATUS "Trying to find LAPACK from Intel MKL")
            if(QE_ENABLE_OPENMP)
                SET(BLA_VENDOR Intel10_64lp)
            else()
                SET(BLA_VENDOR Intel10_64lp_seq)
            endif()
            find_package(LAPACK)
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "i686.*|i386.*|x86.*")
            message(STATUS "Trying to find LAPACK from Intel MKL - 32bit")
            SET(BLA_VENDOR Intel10_32)
            find_package(LAPACK)
        elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)")
            message(STATUS "Trying to find LAPACK from ARM Performance Library")
            if(QE_ENABLE_OPENMP)
                SET(BLA_VENDOR Arm_mp)
            else()
                SET(BLA_VENDOR Arm)
            endif()
            find_package(LAPACK)
        endif()
        if(NOT LAPACK_FOUND)
            message(STATUS "Trying to find alternative LAPACK libraries")
            SET(BLA_VENDOR All)
            if(QE_ENABLE_OPENMP)
                set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_Fortran_FLAGS})
            endif()
            find_package(LAPACK)
            unset(CMAKE_REQUIRED_LINK_OPTIONS)
        endif()
    else()
        if(QE_ENABLE_OPENMP)
            set(CMAKE_REQUIRED_LINK_OPTIONS ${OpenMP_Fortran_FLAGS})
        endif()
        find_package(LAPACK)
        unset(CMAKE_REQUIRED_LINK_OPTIONS)
    endif()
    if(LAPACK_FOUND)
        list(APPEND _lapack_libs
            ${BLAS_LIBRARIES}
            ${BLAS_LINKER_FLAGS}
            ${LAPACK_LIBRARIES}
            ${LAPACK_LINKER_FLAGS})
        if(QE_ENABLE_OPENMP)
            list(APPEND _lapack_libs ${OpenMP_Fortran_LIBRARIES})
        endif()
        list(REMOVE_DUPLICATES "${_lapack_libs}")
        message(STATUS "Found LAPACK: ${_lapack_libs}")
        target_link_libraries(qe_lapack INTERFACE ${_lapack_libs})
        set(CMAKE_REQUIRED_LIBRARIES ${_lapack_libs})
        check_fortran_function_exists(zhpev ZHPEV_FOUND)
        unset(CMAKE_REQUIRED_LIBRARIES)
        if(NOT ZHPEV_FOUND)
          unset(ZHPEV_FOUND CACHE)
          message(FATAL_ERROR "Incomplete LAPACK! function zhpev not found!")
        endif()
    else()
        message(FATAL_ERROR "Failed to find a complete set of external BLAS/LAPACK library by FindLAPACK. "
                            "Variables controlling FindLAPACK can be found at CMake online documentation. "
                            "Alternatively, '-DQE_LAPACK_INTERNAL=ON' may be used to enable reference LAPACK "
                            "at a performance loss compared to optimized libraries.")
    endif()
else()
    message(WARNING "Internal reference LAPACK is enabled! It is less performant than vendor optimized libraries.")
    if(CMAKE_Fortran_COMPILER_ID MATCHES "XL")
        message(FATAL_ERROR "IBM XL compilers cannot build internal LAPACK with QE "
                            "due to the conflict in flags for free vs fixed format. "
                            "Please use an optimized LAPACK or build internal reference LAPACK separately.")
    endif()
    message(STATUS "Installing LAPACK via submodule")
    qe_git_submodule_update(external/lapack)
    add_subdirectory(external/lapack)
    target_link_libraries(qe_lapack INTERFACE lapack)
    # make lapack ready for other external libraries like mbd
    set(LAPACK_LIBRARIES lapack)
endif()

###########################################################
# SCALAPACK
# The following targets will be defined:
add_library(qe_scalapack INTERFACE)
qe_install_targets(qe_scalapack)
###########################################################
if(QE_ENABLE_SCALAPACK)
    find_package(SCALAPACK REQUIRED QUIET)
    message(STATUS "Found SCALAPACK: ${SCALAPACK_LIBRARIES};${SCALAPACK_LINKER_FLAGS}")
    target_link_libraries(qe_scalapack
        INTERFACE
            ${SCALAPACK_LIBRARIES}
            ${SCALAPACK_LINKER_FLAGS})
    if(QE_ENABLE_SCALAPACK_QRCP)
        include(CheckFortranFunctionExists)
        set(CMAKE_REQUIRED_LIBRARIES "${SCALAPACK_LIBRARIES}")
        check_fortran_function_exists("pzgeqpf" SCALAPACK_PZGEQPF_WORKS)
        unset(CMAKE_REQUIRED_LIBRARIES)
        if(SCALAPACK_PZGEQPF_WORKS)
            message(STATUS "Found pzgeqpf, add ScaLAPACK pzgeqpf macro")
            qe_add_global_compile_definitions(__SCALAPACK_QRCP)
        else()
            message(FATAL_ERROR "QE_ENABLE_SCALAPACK_QRCP requested but the current ScaLAPACK installation doesn't contain pzgeqpf!")
        endif()
    endif()
endif(QE_ENABLE_SCALAPACK)

###########################################################
# Components
###########################################################
add_subdirectory(Modules)
add_subdirectory(clib)
add_subdirectory(flib)
add_subdirectory(iotk)
add_subdirectory(AFC90)
add_subdirectory(CPV)

###########################################################
# Pkg-config
###########################################################
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/koopmans-kcp.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/koopmans-kcp.pc
    @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/koopmans-kcp.pc
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)

###########################################################
# Exports
###########################################################
install(EXPORT kcpTargets
        FILE kcpTargets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/koopmans-kcp)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    kcpConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

configure_file(cmake/kcpConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/kcpConfig.cmake @ONLY)

install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/kcpConfigVersion.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/kcpConfig.cmake
        DESTINATION
            ${CMAKE_INSTALL_LIBDIR}/cmake/koopmans-kcp)