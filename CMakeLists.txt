cmake_minimum_required(VERSION 3.16)

project(Object_ChatGPT_MININEC LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(Object_ChatGPT_MININEC main.cpp
    geometry.hpp
    geometry.cpp
    kernel.hpp
    kernel_mininec.cpp
    impedance_assembler.hpp
    impedance_assembler.cpp
    mom_solver.hpp
    mom_solver.cpp
    kernel_mininec.hpp
    impedance_matrix.hpp)

include(GNUInstallDirs)
install(TARGETS Object_ChatGPT_MININEC
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
