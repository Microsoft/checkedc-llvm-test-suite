list(APPEND LDFLAGS -lm)
list(APPEND CFLAGS -std=gnu99)
set(RUN_OPTIONS 12500 14)
llvm_multisource(IndirectAddressing-dbl)
