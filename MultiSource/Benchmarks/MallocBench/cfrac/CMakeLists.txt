list(APPEND CPPFLAGS -DNOMEMOPT)
list(APPEND LDFLAGS -lm)
if(LARGE_PROBLEM_SIZE)
  set(RUN_OPTIONS 376529507064485971626234556089714053987)
else()
  set(RUN_OPTIONS 41757646344123832613190542166099121)
endif()
llvm_multisource(cfrac
  cfrac.c
  pops.c
  pconst.c
  pio.c
  pabs.c
  pneg.c
  pcmp.c
  podd.c
  phalf.c
  padd.c
  psub.c
  pmul.c
  pdivmod.c
  psqrt.c
  ppowmod.c
  atop.c
  ptoa.c
  itop.c
  utop.c
  ptou.c
  errorp.c
  pfloat.c
  pidiv.c
  pimod.c
  picmp.c
  primes.c
  pcfrac.c
  pgcd.c
)
