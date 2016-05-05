cmake_minimum_required(VERSION 3.4)

enable_language(Fortran)

project(oceanbot)

# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_Fortran_FLAGS "-O -traceback")

include_directories("/projects/v45/fortran/modules/")
link_directories("/projects/v45/fortran/modules/")

# Create an executable file 
add_executable("oceanbot" "oceanbot.f90")

target_link_libraries(oceanbot modules netcdff)
