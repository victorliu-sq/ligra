cmake_minimum_required(VERSION 3.16)

project(ligra_openmp LANGUAGES CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib") # output directory of static libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin") # output directory of shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin") # output directory of executables

find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

set(LIGRA_APPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/apps")

function(add_ligra_openmp_app TARGET_NAME APP_NAME)
    add_custom_target(${TARGET_NAME}
        COMMAND "${CMAKE_COMMAND}" -E env --unset=CILK --unset=MKLROOT OPENMP=1
                "${MAKE_EXECUTABLE}" -C "${LIGRA_APPS_DIR}" "${APP_NAME}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${LIGRA_APPS_DIR}/${APP_NAME}"
                "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${APP_NAME}"
        COMMENT "Build Ligra ${APP_NAME} with OpenMP"
    )
endfunction()

add_ligra_openmp_app(ligra_bfs BFS)
add_ligra_openmp_app(ligra_bc BC)
add_ligra_openmp_app(ligra_cc Components)
add_ligra_openmp_app(ligra_sssp SSSP)
add_ligra_openmp_app(ligra_multisssp MultiSSSP)
add_ligra_openmp_app(ligra_multisssp_throughput MultiSSSP_Throughput)
add_ligra_openmp_app(ligra_bfs_throughput BFS_Throughput)
add_ligra_openmp_app(ligra_bc_throughput BC_Throughput)
add_ligra_openmp_app(ligra_sssp_throughput SSSP_Throughput)

add_custom_target(ligra_openmp ALL
    DEPENDS ligra_bfs ligra_bc ligra_sssp ligra_multisssp
    COMMENT "Build selected Ligra OpenMP apps: BFS, BC, SSSP, MultiSSSP"
)
