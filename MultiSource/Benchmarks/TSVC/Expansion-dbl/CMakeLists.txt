list(APPEND LDFLAGS -lm)
list(APPEND CFLAGS -std=gnu99)
set(RUN_OPTIONS 4160 14)
llvm_multisource(Expansion-dbl)
