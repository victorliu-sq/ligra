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
# X-Blossom-Pro maximum matching written as a Ligra app (apps/XBlossomPro.C).
add_ligra_openmp_app(ligra_xblossom_pro XBlossomPro)
# The same with dynamic OpenMP scheduling (load-balance experiment).
add_ligra_openmp_app(ligra_xblossom_pro_dyn XBlossomPro_dyn)

# The same apps built with OpenCilk (work stealing) as <App>_cilk, when
# OpenCilk is installed (deps/targets/scripts/ligra/install_opencilk.sh).
set(LIGRA_OPENCILK_HOME "${CMAKE_SOURCE_DIR}/../deps/opencilk" CACHE PATH "OpenCilk install prefix")
function(add_ligra_opencilk_app TARGET_NAME APP_NAME)
    add_custom_target(${TARGET_NAME}
        COMMAND "${MAKE_EXECUTABLE}" -C "${LIGRA_APPS_DIR}" "${APP_NAME}_cilk" "OPENCILK_HOME=${LIGRA_OPENCILK_HOME}"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}"
        COMMAND "${CMAKE_COMMAND}" -E copy_if_different
                "${LIGRA_APPS_DIR}/${APP_NAME}_cilk"
                "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${APP_NAME}_cilk"
        COMMENT "Build Ligra ${APP_NAME} with OpenCilk"
    )
endfunction()
if(EXISTS "${LIGRA_OPENCILK_HOME}/bin/clang++")
    add_ligra_opencilk_app(ligra_bfs_cilk BFS)
    add_ligra_opencilk_app(ligra_bc_cilk BC)
    add_ligra_opencilk_app(ligra_sssp_cilk SSSP)
    add_ligra_opencilk_app(ligra_multisssp_cilk MultiSSSP)
    add_ligra_opencilk_app(ligra_xblossom_pro_cilk XBlossomPro)
    add_custom_target(ligra_cilk DEPENDS ligra_bfs_cilk ligra_bc_cilk ligra_sssp_cilk ligra_multisssp_cilk
                                         ligra_xblossom_pro_cilk)
else()
    add_custom_target(ligra_cilk COMMAND "${CMAKE_COMMAND}" -E echo
        "ligra_cilk: OpenCilk not found at ${LIGRA_OPENCILK_HOME}; skipped")
endif()

add_custom_target(ligra_openmp ALL
    DEPENDS ligra_bfs ligra_bc ligra_sssp ligra_multisssp
    COMMENT "Build selected Ligra OpenMP apps: BFS, BC, SSSP, MultiSSSP"
)
