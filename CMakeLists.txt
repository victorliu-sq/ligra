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

# <App>_cilkscale: the same app with Cilkscale, the work and span of each Compute().
function(add_ligra_cilkscale_app TARGET_NAME APP_NAME)
    add_custom_target(${TARGET_NAME}
        COMMAND "${MAKE_EXECUTABLE}" -C "${LIGRA_APPS_DIR}" "${APP_NAME}_cilkscale" "OPENCILK_HOME=${LIGRA_OPENCILK_HOME}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${LIGRA_APPS_DIR}/${APP_NAME}_cilkscale"
                "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${APP_NAME}_cilkscale"
        COMMENT "Build Ligra ${APP_NAME} with OpenCilk and Cilkscale"
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
# Its reuse-measurement variants: XB-Ligra (no tree reuse), XB-Pro-Ligra with counters.
add_ligra_cilk_app(ligra_xblossom_noreuse_cilk XBlossomPro_noreuse)
add_ligra_cilk_app(ligra_xblossom_pro_counters_cilk XBlossomPro_counters)
# XB-Pro's own design (path mutex, per-edge allocations, serial copies), for measure/xb-cpu-opt.
add_ligra_cilk_app(ligra_xblossom_pro_xbpro_cilk XBlossomPro_xbpro)
# XB-Pro's path table (a vector per node) instead of the default triplet records, for the design comparisons.
add_ligra_cilk_app(ligra_xblossom_pro_vector_cilk XBlossomPro_vector)

# Iterative algorithms (SEP-Graph's classification): PageRank and Triangle Counting.
add_ligra_cilk_app(ligra_pagerank_cilk PageRank)
add_ligra_cilk_app(ligra_triangle_cilk Triangle)

add_ligra_cilkscale_app(ligra_bfs_cilkscale BFS)
add_ligra_cilkscale_app(ligra_sssp_cilkscale SSSP)
add_ligra_cilkscale_app(ligra_multisssp_cilkscale MultiSSSP)
add_ligra_cilkscale_app(ligra_xblossom_pro_cilkscale XBlossomPro)
add_ligra_cilkscale_app(ligra_pagerank_cilkscale PageRank)
add_ligra_cilkscale_app(ligra_triangle_cilkscale Triangle)

add_custom_target(ligra_cilk ALL
    DEPENDS ligra_bfs_cilk ligra_sssp_cilk ligra_multisssp_cilk
            ligra_bfs_throughput_cilk ligra_sssp_throughput_cilk ligra_multisssp_throughput_cilk
            ligra_xblossom_pro_cilk ligra_xblossom_noreuse_cilk ligra_xblossom_pro_counters_cilk
            ligra_xblossom_pro_xbpro_cilk ligra_xblossom_pro_vector_cilk
            ligra_bfs_cilkscale ligra_sssp_cilkscale ligra_multisssp_cilkscale ligra_xblossom_pro_cilkscale
            ligra_pagerank_cilk ligra_triangle_cilk ligra_pagerank_cilkscale ligra_triangle_cilkscale
    COMMENT "Build the Ligra apps with OpenCilk: BFS, SSSP, MultiSSSP (+ _Throughput), XBlossomPro, PageRank, Triangle"
)
