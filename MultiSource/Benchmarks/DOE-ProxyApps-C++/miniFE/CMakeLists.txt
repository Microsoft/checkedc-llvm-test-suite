list(APPEND CXXFLAGS -I. -DMINIFE_SCALAR=double  -DMINIFE_LOCAL_ORDINAL=int -DMINIFE_GLOBAL_ORDINAL=int -DMINIFE_CSR_MATRIX -DREDSTORM) # -DREDSTORM for mkdir
list(APPEND LDFLAGS -lm)
set(RUN_OPTIONS -nx 64 -ny 64 -nz 64)
llvm_multisource(miniFE)
target_include_directories(miniFE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
