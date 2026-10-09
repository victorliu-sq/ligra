cmake_minimum_required(VERSION 3.16)

project(ligra_cilk LANGUAGES CXX)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/lib") # output directory of static libraries
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin") # output directory of shared libraries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin") # output directory of executables

find_program(MAKE_EXECUTABLE NAMES make gmake REQUIRED)

set(LIGRA_APPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/apps")

# The Ligra apps are built with OpenCilk (clang -fopencilk, Cilk's
# work-stealing runtime, which Ligra was designed for) as <App>_cilk.
# OpenCilk is installed by deps/targets/scripts/ligra/install_opencilk.sh.
set(LIGRA_OPENCILK_HOME "${CMAKE_SOURCE_DIR}/../deps/opencilk" CACHE PATH "OpenCilk install prefix")
if(NOT EXISTS "${LIGRA_OPENCILK_HOME}/bin/clang++")
    message(FATAL_ERROR "OpenCilk not found at ${LIGRA_OPENCILK_HOME}; run "
                        "deps/targets/scripts/ligra/install_opencilk.sh (or set LIGRA_OPENCILK_HOME)")
endif()

function(add_ligra_cilk_app TARGET_NAME APP_NAME)
    add_custom_target(${TARGET_NAME}
        COMMAND "${MAKE_EXECUTABLE}" -C "${LIGRA_APPS_DIR}" "${APP_NAME}_cilk" "OPENCILK_HOME=${LIGRA_OPENCILK_HOME}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${LIGRA_APPS_DIR}/${APP_NAME}_cilk"
                "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${APP_NAME}_cilk"
        COMMENT "Build Ligra ${APP_NAME} with OpenCilk"
    )
endfunction()

add_ligra_cilk_app(ligra_bfs_cilk BFS)
add_ligra_cilk_app(ligra_sssp_cilk SSSP)
add_ligra_cilk_app(ligra_multisssp_cilk MultiSSSP)
add_ligra_cilk_app(ligra_bfs_throughput_cilk BFS_Throughput)
add_ligra_cilk_app(ligra_sssp_throughput_cilk SSSP_Throughput)
add_ligra_cilk_app(ligra_multisssp_throughput_cilk MultiSSSP_Throughput)
# X-Blossom-Pro maximum matching written as a Ligra app (apps/XBlossomPro.C).
add_ligra_cilk_app(ligra_xblossom_pro_cilk XBlossomPro)

add_custom_target(ligra_cilk ALL
    DEPENDS ligra_bfs_cilk ligra_sssp_cilk ligra_multisssp_cilk
            ligra_bfs_throughput_cilk ligra_sssp_throughput_cilk ligra_multisssp_throughput_cilk
            ligra_xblossom_pro_cilk
    COMMENT "Build the Ligra apps with OpenCilk: BFS, SSSP, MultiSSSP (+ _Throughput), XBlossomPro"
)
