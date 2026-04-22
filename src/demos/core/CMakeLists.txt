#--------------------------------------------------------------
# Add executables

set(DEMOS
  demo_CH_coords
  demo_CH_linalg
  demo_CH_matrix_ref
  demo_CH_math
  demo_CH_archive
  demo_CH_buildsystem
  demo_CH_powertrain
  demo_CH_functions
  demo_CH_solver
  demo_CH_EulerAngles
  demo_CH_filesystem
  demo_CH_baseline_contact
  demo_CH_sdf_point_contact
)

# OpenVDB probe demo is built separately below if OpenVDB is found
list(APPEND DEMOS_WITH_OPENVDB demo_CH_openvdb_probe demo_CH_sdf_point_contact_openvdb demo_CH_sdf_patch_contact_openvdb demo_CH_sdf_patch_force_openvdb)

#--------------------------------------------------------------

# Set include paths, compiler & linker flags, and libraries

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

list(APPEND LIBS Chrono_core)

if(CH_ENABLE_MODULE_PARDISO_MKL)
  list(APPEND LIBS "Chrono_pardisomkl")
endif()

if(CH_ENABLE_MODULE_MUMPS)
    list(APPEND LIBS "Chrono_mumps")
endif()

#--------------------------------------------------------------

message(STATUS "Add demo programs for CORE module")

foreach(PROGRAM ${DEMOS})
    add_executable(${PROGRAM}  "${PROGRAM}.cpp")
    source_group(""  FILES "${PROGRAM}.cpp")

    set_target_properties(${PROGRAM} PROPERTIES FOLDER demos)
    set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
    target_link_libraries(${PROGRAM} ${LIBS})
    if(MSVC)
      set_target_properties(${PROGRAM} PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
    endif()

    install(TARGETS ${PROGRAM} CONFIGURATIONS Release DESTINATION ${CH_INSTALL_DEMO})
endforeach(PROGRAM)

#--------------------------------------------------------------
# OpenVDB probe demo (separate target with OpenVDB + TBB + Imath deps)
#--------------------------------------------------------------

find_package(OpenVDB QUIET)
if(OpenVDB_FOUND)
  message(STATUS "OpenVDB found: ${OpenVDB_VERSION}")
  message(STATUS "  Include dir: ${OpenVDB_INCLUDE_DIRS}")
  message(STATUS "  Libraries:   ${OpenVDB_LIBRARIES}")

  # Locate TBB libraries for linking
  find_path(TBB_LIB_DIR NAMES tbb12.lib tbb12_debug.lib
    PATHS "C:/vcpkg/installed/x64-windows/lib" "C:/vcpkg/installed/x64-windows/debug/lib"
    NO_DEFAULT_PATH
  )
  if(TBB_LIB_DIR)
    message(STATUS "TBB lib dir: ${TBB_LIB_DIR}")
  endif()

  # Use release OpenVDB lib for Release builds
  set(OPENVDB_RELEASE_LIB "C:/vcpkg/installed/x64-windows/lib/openvdb.lib")
  set(OPENVDB_LIBS Imath-3_2)
  if(EXISTS "${OPENVDB_RELEASE_LIB}")
    list(APPEND OPENVDB_LIBS "${OPENVDB_RELEASE_LIB}")
    message(STATUS "Using release OpenVDB lib")
  else()
    list(APPEND OPENVDB_LIBS ${OpenVDB_LIBRARIES})
  endif()

  foreach(PROGRAM ${DEMOS_WITH_OPENVDB})
    add_executable(${PROGRAM} "${PROGRAM}.cpp")
    source_group("" FILES "${PROGRAM}.cpp")
    set_target_properties(${PROGRAM} PROPERTIES FOLDER demos)
    set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")

    # Link libraries: Chrono libs, then Imath, then OpenVDB (order matters for MSVC)
    target_link_libraries(${PROGRAM} ${LIBS} ${OPENVDB_LIBS})
    target_include_directories(${PROGRAM} PRIVATE ${OpenVDB_INCLUDE_DIRS})

    # Define IMATH_DLL so that Imath symbols are imported from DLL
    target_compile_definitions(${PROGRAM} PRIVATE IMATH_DLL)

    # Add TBB library search directory
    if(TBB_LIB_DIR)
      target_link_directories(${PROGRAM} PRIVATE "${TBB_LIB_DIR}")
    endif()

    if(MSVC)
      set_target_properties(${PROGRAM} PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
    endif()

    install(TARGETS ${PROGRAM} CONFIGURATIONS Release DESTINATION ${CH_INSTALL_DEMO})
  endforeach()
else()
  message(STATUS "OpenVDB not found - skipping OpenVDB demo programs")
endif()
