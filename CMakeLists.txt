cmake_minimum_required(VERSION 3.27)
project(DCE_cpp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


include_directories("/home/polya/cppLocal/eigen-3.4.0")

include_directories("/home/polya/cppLocal/msgpack-c-cpp_master/include")

set(Boost_NO_SYSTEM_PATHS YES)
if (Boost_NO_SYSTEM_PATHS)
    set(Boost_ROOT "/home/polya/cppLocal/boostBuild")
    set(Boost_INCLUDE_DIR "/home/polya/cppLocal/boostBuild/include")
    set(Boost_LIBRARY_DIR "/home/polya/cppLocal/boostBuild/lib")
endif (Boost_NO_SYSTEM_PATHS)

add_executable(DCE main.cpp timeEvolution.cpp
        )
find_package(Boost  REQUIRED filesystem)



if(Boost_FOUND)
    message(STATUS "boost root is: ${Boost_ROOT}")
    message(STATUS "boost include path is : ${Boost_INCLUDE_DIR}")
    message(STATUS "boost library path is : ${Boost_LIBRARY_DIR}")
    message(STATUS "boost libraries is : ${Boost_LIBRARIES}")
    include_directories(${Boost_INCLUDE_DIR})
    link_directories(${Boost_LIBRARY_DIR})
else()
    message(WARNING "boost not found.")
endif()

set_source_files_properties(
        timeEvolution.cpp
        PROPERTIES
        COMPILE_OPTIONS "-Ofast"
)

file(GLOB LIBS_A "${Boost_LIBRARY_DIR}/*filesystem.a")
file(GLOB LIBS_SO "${Boost_LIBRARY_DIR}/*filesystem.so")

target_link_libraries(DCE ${LIBS_A} ${LIBS_SO})
