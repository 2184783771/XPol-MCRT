cmake_minimum_required(VERSION 3.16)
project(XPolMCRT LANGUAGES CXX)

add_executable(xpol_mcrt
  corona-lamppost.cpp
  corona-slab.cpp
  corona-toroidal.cpp
  corona.cpp
  detector.cpp
  disk-reflection.cpp
  main.cpp
  run-simulation.cpp
  run-report.cpp
  model.cpp
  model-parse-basic.cpp
  model-parse-xrb.cpp
  observe.cpp
  photon.cpp
  specmath.cpp
)

set_property(DIRECTORY PROPERTY VS_STARTUP_PROJECT xpol_mcrt)

set_target_properties(xpol_mcrt PROPERTIES
  CXX_STANDARD 11
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  OUTPUT_NAME "XPol-MCRT"
  VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

set(VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}" CACHE PATH
    "Working directory for Visual Studio debugging" FORCE)
target_compile_definitions(xpol_mcrt PRIVATE
  XPOL_SOURCE_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

if(MSVC)
  target_compile_definitions(xpol_mcrt PRIVATE _CRT_SECURE_NO_WARNINGS)
  target_link_options(xpol_mcrt PRIVATE /STACK:67108864)
endif()
