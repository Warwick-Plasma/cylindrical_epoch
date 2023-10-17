if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_definitions(-DHASH_DEBUG)
endif()
set(VALGRIND OFF CACHE BOOL "Compile with valgrind support")
if(VALGRIND)
    add_definitions(-DVALGRIND)
endif()

set(SOURCES
   sdf_control.c sdf_derived.c sdf_extension_util.c sdf_helper.c
   sdf_input.c sdf_input_cartesian.c sdf_input_point.c sdf_input_station.c
   sdf_modify.c sdf_output.c sdf_util.c stack_allocator.c
   )

add_library(sdfc STATIC ${SOURCES} commit_info.h)
if(PARALLEL)
   set_source_files_properties(
       ${SOURCES}
       PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif()
target_include_directories(
    sdfc
    INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/commit_info.h
   COMMAND sh gen_commit_string.sh ${CMAKE_CURRENT_BINARY_DIR}
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

install(
    TARGETS sdfc
    EXPORT SDFCConfig
    COMPONENT develop
    LIBRARY DESTINATION lib
    INCLUDES DESTINATION include
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/SDFCConfigVersion.cmake
    VERSION ${SDFC_VERSION}
    COMPATIBILITY SameMajorVersion
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/SDFCConfigVersion.cmake
    DESTINATION lib/SDFC
    COMPONENT develop
)

export(TARGETS sdfc FILE ${CMAKE_CURRENT_BINARY_DIR}/SDFCConfig.cmake)
install(EXPORT SDFCConfig DESTINATION lib/SDFC COMPONENT develop)
