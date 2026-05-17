cmake_minimum_required(VERSION 3.16)

project(ligra_openmp LANGUAGES CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib") # output directory of static libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin") # output directory of shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin") # output directory of executables

find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

set(LIGRA_APPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/apps")
set(LIGRA_TARGETS
    encoder
    decoder
    BFS
    BC
    BellmanFord
    Components
    Components-Shortcut
    Radii
    PageRank
    PageRankDelta
    BFSCC
    BFS-Bitvector
    KCore
    MIS
    Triangle
    CF
)

set(LIGRA_TARGET_PATHS)
foreach(LIGRA_TARGET IN LISTS LIGRA_TARGETS)
    list(APPEND LIGRA_TARGET_PATHS "${LIGRA_APPS_DIR}/${LIGRA_TARGET}")
endforeach()

add_custom_target(ligra_openmp ALL
    COMMAND "${CMAKE_COMMAND}" -E env --unset=CILK --unset=MKLROOT OPENMP=1
            "${MAKE_EXECUTABLE}" -C "${LIGRA_APPS_DIR}" clean
    COMMAND "${CMAKE_COMMAND}" -E env --unset=CILK --unset=MKLROOT OPENMP=1
            "${MAKE_EXECUTABLE}" -C "${LIGRA_APPS_DIR}" -j
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different ${LIGRA_TARGET_PATHS} "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
    COMMENT "Build Ligra apps with OpenMP"
)
