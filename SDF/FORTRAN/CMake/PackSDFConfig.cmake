find_package(PythonInterp)
find_file(PACK_PY pack.py
   PATHS ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/../src
   DOC "Path to pack.py")
find_file(PACK_CMAKE_IN pack.cmake.in build_scripts
   PATHS ${CMAKE_CURRENT_LIST_DIR}
   DOC "Path to pack.cmake.in")
find_file(CONFIGURE_PACK_TARGET ConfigurePackTarget.cmake
   PATHS ${CMAKE_CURRENT_LIST_DIR}
   DOC "Path to ConfigurePackTarget.cmake")
set(PACK_FILES ${PACK_PY} ${PACK_CMAKE_IN} ${CONFIGURE_PACK_TARGET}
   CACHE STRING "Python and CMake scripts for packing Fortran source")
