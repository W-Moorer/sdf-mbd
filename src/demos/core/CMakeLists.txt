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
  demo_CH_field_contact_regression
  demo_CH_field_contact_feature_switching
  demo_CH_field_contact_multislip_interlock
  demo_CH_chrono_native_contact_baseline
  demo_CH_friction_validation
)

# OpenVDB probe demo is built separately below if OpenVDB is found
list(APPEND DEMOS_WITH_OPENVDB demo_CH_openvdb_probe demo_CH_sdf_point_contact_openvdb demo_CH_sdf_patch_contact_openvdb demo_CH_sdf_patch_force_openvdb demo_CH_sdf_patch_persistence_openvdb demo_CH_sdf_patch_constitutive_openvdb demo_CH_sdf_patch_constitutive_cases_openvdb demo_CH_sdf_patch_constitutive_tangential_openvdb demo_CH_sdf_patch_constitutive_stickslip_openvdb demo_CH_sdf_patch_constitutive_frame_openvdb demo_CH_sdf_patch_constitutive_transport_limited_openvdb demo_CH_sdf_patch_transport_ablation_openvdb demo_CH_sdf_patch_adaptive_carryover_openvdb demo_CH_sdf_patch_tangential_history_v1_openvdb demo_CH_field_contact_primitives_openvdb demo_CH_field_contact_split_merge_openvdb demo_CH_field_contact_mesh_sdf_nonconvex_openvdb demo_CH_field_contact_matched_mesh_native_openvdb demo_CH_field_contact_free_dynamics_openvdb demo_CH_field_contact_cad_mesh_case_study_openvdb demo_CH_field_contact_chrono_original_mesh_baseline_openvdb demo_CH_simple_gear_openvdb demo_CH_rev_joint_clearance_openvdb)

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

if(CH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS)
  add_test(NAME demo_CH_field_contact_regression
           COMMAND "$<TARGET_FILE:demo_CH_field_contact_regression>")
  set_tests_properties(demo_CH_field_contact_regression PROPERTIES
                       LABELS "collision;field_contact")

  add_test(NAME demo_CH_friction_validation
           COMMAND "$<TARGET_FILE:demo_CH_friction_validation>"
                   "${CMAKE_SOURCE_DIR}/out/friction_validation")
  set_tests_properties(demo_CH_friction_validation PROPERTIES
                       LABELS "collision;field_contact;friction"
                       TIMEOUT 60)

  find_package(Python3 COMPONENTS Interpreter QUIET)
  if(Python3_Interpreter_FOUND)
    add_test(NAME demo_CH_field_contact_acceptance_regression
             COMMAND "${Python3_EXECUTABLE}"
                     "${CMAKE_SOURCE_DIR}/scripts/check_field_contact_acceptance.py"
                     --project-root "${CMAKE_SOURCE_DIR}"
                     --feature-demo "$<TARGET_FILE:demo_CH_field_contact_feature_switching>"
                     --multislip-demo "$<TARGET_FILE:demo_CH_field_contact_multislip_interlock>"
                     --native-demo "$<TARGET_FILE:demo_CH_chrono_native_contact_baseline>"
                     --compare-script "${CMAKE_SOURCE_DIR}/scripts/compare_chrono_native_field.py")
    set_tests_properties(demo_CH_field_contact_acceptance_regression PROPERTIES
                         LABELS "collision;field_contact;acceptance"
                         TIMEOUT 300)
  else()
    message(STATUS "Python3 interpreter not found - skipping field-contact acceptance regression test")
  endif()
endif()

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

  add_executable(demo_CH_paper_examples_openvdb "${CMAKE_SOURCE_DIR}/paper_example/paper_examples_openvdb.cpp")
  source_group("" FILES "${CMAKE_SOURCE_DIR}/paper_example/paper_examples_openvdb.cpp")
  set_target_properties(demo_CH_paper_examples_openvdb PROPERTIES FOLDER demos)
  set_property(TARGET demo_CH_paper_examples_openvdb PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:demo_CH_paper_examples_openvdb>")
  target_link_libraries(demo_CH_paper_examples_openvdb ${LIBS} ${OPENVDB_LIBS})
  target_include_directories(demo_CH_paper_examples_openvdb PRIVATE ${OpenVDB_INCLUDE_DIRS})
  target_compile_definitions(demo_CH_paper_examples_openvdb PRIVATE IMATH_DLL)
  if(TBB_LIB_DIR)
    target_link_directories(demo_CH_paper_examples_openvdb PRIVATE "${TBB_LIB_DIR}")
  endif()
  if(MSVC)
    set_target_properties(demo_CH_paper_examples_openvdb PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
  endif()
  install(TARGETS demo_CH_paper_examples_openvdb CONFIGURATIONS Release DESTINATION ${CH_INSTALL_DEMO})

  if(CH_ENABLE_FIELD_CONTACT_REGRESSION_TESTS)
    find_package(Python3 COMPONENTS Interpreter QUIET)
    if(Python3_Interpreter_FOUND)
      add_test(NAME demo_CH_field_contact_matched_mesh_native_baseline
               COMMAND "$<TARGET_FILE:demo_CH_field_contact_matched_mesh_native_openvdb>")
      set_tests_properties(demo_CH_field_contact_matched_mesh_native_baseline PROPERTIES
                           LABELS "collision;field_contact;matched_mesh"
                           FIXTURES_SETUP matched_mesh_baseline
                           TIMEOUT 300)

      add_test(NAME demo_CH_field_contact_matched_mesh_regression
               COMMAND "${Python3_EXECUTABLE}"
                       "${CMAKE_SOURCE_DIR}/scripts/check_matched_mesh_baseline.py"
                       --project-root "${CMAKE_SOURCE_DIR}")
      set_tests_properties(demo_CH_field_contact_matched_mesh_regression PROPERTIES
                           LABELS "collision;field_contact;matched_mesh;acceptance"
                           FIXTURES_REQUIRED matched_mesh_baseline
                           TIMEOUT 60)

      add_test(NAME demo_CH_field_contact_free_dynamics_response
               COMMAND "$<TARGET_FILE:demo_CH_field_contact_free_dynamics_openvdb>")
      set_tests_properties(demo_CH_field_contact_free_dynamics_response PROPERTIES
                           LABELS "collision;field_contact;free_dynamics"
                           FIXTURES_SETUP free_dynamics_response
                           TIMEOUT 300)

      add_test(NAME demo_CH_field_contact_free_dynamics_regression
               COMMAND "${Python3_EXECUTABLE}"
                       "${CMAKE_SOURCE_DIR}/scripts/check_free_dynamics_response.py"
                       --project-root "${CMAKE_SOURCE_DIR}")
      set_tests_properties(demo_CH_field_contact_free_dynamics_regression PROPERTIES
                           LABELS "collision;field_contact;free_dynamics;acceptance"
                           FIXTURES_REQUIRED free_dynamics_response
                           TIMEOUT 60)

      add_test(NAME demo_CH_field_contact_cad_mesh_case_study
               COMMAND "$<TARGET_FILE:demo_CH_field_contact_cad_mesh_case_study_openvdb>")
      set_tests_properties(demo_CH_field_contact_cad_mesh_case_study PROPERTIES
                           LABELS "collision;field_contact;cad_mesh"
                           FIXTURES_SETUP cad_mesh_case_study
                           TIMEOUT 300)

      add_test(NAME demo_CH_field_contact_cad_mesh_case_study_regression
               COMMAND "${Python3_EXECUTABLE}"
                       "${CMAKE_SOURCE_DIR}/scripts/check_cad_mesh_case_study.py"
                       --project-root "${CMAKE_SOURCE_DIR}")
      set_tests_properties(demo_CH_field_contact_cad_mesh_case_study_regression PROPERTIES
                           LABELS "collision;field_contact;cad_mesh;acceptance"
                           FIXTURES_REQUIRED cad_mesh_case_study
                           TIMEOUT 60)

      add_test(NAME demo_CH_field_contact_chrono_original_mesh_baseline
               COMMAND "$<TARGET_FILE:demo_CH_field_contact_chrono_original_mesh_baseline_openvdb>")
      set_tests_properties(demo_CH_field_contact_chrono_original_mesh_baseline PROPERTIES
                           LABELS "collision;field_contact;chrono_original_mesh"
                           FIXTURES_SETUP chrono_original_mesh_baseline
                           TIMEOUT 300)

      add_test(NAME demo_CH_field_contact_chrono_original_mesh_regression
               COMMAND "${Python3_EXECUTABLE}"
                       "${CMAKE_SOURCE_DIR}/scripts/check_chrono_original_mesh_baseline.py"
                       --project-root "${CMAKE_SOURCE_DIR}")
      set_tests_properties(demo_CH_field_contact_chrono_original_mesh_regression PROPERTIES
                           LABELS "collision;field_contact;chrono_original_mesh;acceptance"
                           FIXTURES_REQUIRED chrono_original_mesh_baseline
                           TIMEOUT 60)

      add_test(NAME demo_CH_paper_examples_dynamic
               COMMAND "$<TARGET_FILE:demo_CH_paper_examples_openvdb>")
      set_tests_properties(demo_CH_paper_examples_dynamic PROPERTIES
                           LABELS "collision;field_contact;paper_examples"
                           FIXTURES_SETUP paper_examples_dynamic
                           TIMEOUT 420)

      add_test(NAME demo_CH_paper_examples_regression
               COMMAND "${Python3_EXECUTABLE}"
                       "${CMAKE_SOURCE_DIR}/paper_example/check_paper_examples.py"
                       --project-root "${CMAKE_SOURCE_DIR}")
      set_tests_properties(demo_CH_paper_examples_regression PROPERTIES
                           LABELS "collision;field_contact;paper_examples;acceptance"
                           FIXTURES_REQUIRED paper_examples_figures
                           TIMEOUT 60)

      add_test(NAME demo_CH_paper_examples_figures
               COMMAND "${Python3_EXECUTABLE}"
                       "${CMAKE_SOURCE_DIR}/paper_example/plot_paper_examples.py"
                       --project-root "${CMAKE_SOURCE_DIR}")
      set_tests_properties(demo_CH_paper_examples_figures PROPERTIES
                           LABELS "collision;field_contact;paper_examples;figures"
                           FIXTURES_REQUIRED paper_examples_dynamic
                           FIXTURES_SETUP paper_examples_figures
                           TIMEOUT 60)
    else()
      message(STATUS "Python3 interpreter not found - skipping OpenVDB field-contact regression tests")
    endif()
  endif()
else()
  message(STATUS "OpenVDB not found - skipping OpenVDB demo programs")
endif()
