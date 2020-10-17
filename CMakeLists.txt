cmake_minimum_required(VERSION 3.16)

project(JacobiEllepticCurve LANGUAGES C)
find_package(OpenSSL REQUIRED)
add_executable(JacobiEllepticCurve main.c curve.c curve.h)
target_link_libraries(JacobiEllepticCurve OpenSSL::Crypto)
