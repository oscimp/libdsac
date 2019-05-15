# Include header
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include/dsac")

# Create library dynamic
add_library(dsac SHARED
  crossspectrum.c
  fir.c
  fft.c
  cordic.c
  fileHandler.c
  gnuplot_i.c
  nco.c
  piosbr_mod.c
  reg_lin.c
  tf.c
  t_pnm.c
  linspace.c
  logspace.c
  logalizit.c
  conversions.c
  fake_signals.c
)

set_target_properties(dsac
  PROPERTIES SOVERSION 1.0
)
