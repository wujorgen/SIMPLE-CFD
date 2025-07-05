cmake_minimum_required(VERSION 3.15...4.0)

project(
	SIMPLECavity
	VERSION 0.1
	LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

add_executable(
	SIMPLE
	main.cpp 
	MomentumEquations.cpp 
	Corrections.cpp 
	Boundary.cpp
	Utilities.cpp
	InputReader.cpp
	SIMPLE.cpp
)

find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
	target_link_libraries(SIMPLE PRIVATE Eigen3::Eigen)
else()
	if (APPLE)
		target_include_directories(SIMPLE PRIVATE /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/)
	elseif(UNIX)
		target_include_directories(SIMPLE PRIVATE /usr/local/include/eigen-3.4.0)
	endif()
endif()

target_compile_options(SIMPLE PRIVATE -march=native)


