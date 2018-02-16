cmake_minimum_required(VERSION 3.9)
project(RunningStatistics)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package( Threads REQUIRED )
MESSAGE(STATUS "Threads Lib: " ${CMAKE_THREAD_LIBS_INIT})
MESSAGE(STATUS "Using PTHreads?: "  ${CMAKE_USE_PTHREADS_INIT})

add_executable(rs  rs.cc)

add_executable(rsmt rs.cc)
target_compile_options(rsmt PRIVATE "-D_RUNNING_STATISTICS_MT_")
target_link_libraries(rsmt ${CMAKE_THREAD_LIBS_INIT})