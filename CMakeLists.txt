cmake_minimum_required(VERSION 4.1...4.3)

project(heat_bvp LANGUAGES CXX)

include(FetchContent)

set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
    tomlplusplus
    GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
    GIT_TAG        v3.4.0
    GIT_SHALLOW    TRUE
)

FetchContent_Declare(
    exprtk
    GIT_REPOSITORY https://github.com/ArashPartow/exprtk.git
    GIT_TAG        master
    GIT_SHALLOW    TRUE
)

FetchContent_MakeAvailable(tomlplusplus exprtk)

add_library(exprtk INTERFACE)
target_include_directories(exprtk INTERFACE
    ${exprtk_SOURCE_DIR}
)
add_library(exprtk::exprtk ALIAS exprtk)

add_library(heat_core STATIC
    src/io/config_loader.cpp
    src/io/writer.cpp
    src/math/bvp_solver.cpp
    src/app/app.cpp
)
target_include_directories(heat_core PUBLIC include)
target_link_libraries(heat_core PUBLIC exprtk::exprtk tomlplusplus::tomlplusplus)
target_compile_features(heat_core PUBLIC cxx_std_23)

add_executable(heat_app
    src/main.cpp
)
target_link_libraries(heat_app PRIVATE heat_core)

