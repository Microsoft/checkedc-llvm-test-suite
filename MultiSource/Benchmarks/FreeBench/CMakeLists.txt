add_subdirectory(distray)
add_subdirectory(fourinarow)
add_subdirectory(mason)
add_subdirectory(neural)
add_subdirectory(pcompress2)
add_subdirectory(pifft)
if(NOT "${ARCH}" STREQUAL "Sparc")
  add_subdirectory(analyzer)
endif()
