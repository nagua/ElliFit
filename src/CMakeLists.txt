project(elliFit)


set(CMAKE_EXPORT_COMPILE_COMMANDS 1)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED 14)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
add_definitions(-D_USE_MATH_DEFINES)

find_package(Eigen3 REQUIRED CONFIG)

set (SRC_LIST
	main.cpp
  )

set (HDR_LIST
  ellipse.h
  ellifit.h
  )

add_executable (${PROJECT_NAME} ${SRC_LIST} ${HDR_LIST})
target_link_libraries (${PROJECT_NAME} ${OpenCV_LIBS} Eigen3::Eigen)
