cmake_minimum_required(VERSION 3.7...3.18)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(CHyPSS VERSION 0.1
               DESCRIPTION "MFEM Learning Code"
               LANGUAGES CXX)

find_package(MPI)
message("${MPI_CXX_LIBRARIES}")


set(CHYPSS_SOURCE "./chypss")


add_executable(scale_test)
set_property(TARGET scale_test PROPERTY CXX_STANDARD 11)
target_include_directories(scale_test PRIVATE /usr/local/include)
target_include_directories(scale_test PRIVATE /Users/rmc298/lib_builds/hypre/include)
target_include_directories(scale_test PRIVATE "${MPI_INCLUDE_PATH}")
target_link_libraries(scale_test PRIVATE "${MPI_CXX_LIBRARIES}")
target_link_libraries(scale_test PRIVATE /usr/local/lib/libmfem.a)
target_link_libraries(scale_test PRIVATE /Users/rmc298/lib_builds/hypre/lib/libHYPRE.a)
target_link_libraries(scale_test PRIVATE /Users/rmc298/lib_builds/hypre/lib/libHYPRE.a)
target_link_libraries(scale_test PRIVATE /usr/local/lib/libmetis.dylib)

add_subdirectory("${CHYPSS_SOURCE}")
