list(APPEND LDFLAGS -lm)
list(APPEND CFLAGS -std=gnu99)
set(RUN_OPTIONS 12500 5)
llvm_multisource(IndirectAddressing-flt)
