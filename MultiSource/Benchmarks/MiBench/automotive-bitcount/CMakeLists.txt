list(APPEND LDFLAGS -lm)
set(RUN_OPTIONS 1125000)
llvm_multisource(automotive-bitcount)
