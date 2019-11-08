# Include header
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../include/dsac")

# Create library dynamic
add_library(dsac SHARED
  conversions.c
  cordic.c
  crossspectrum.c
  fake_signals.c
  fft.c
  file_handler.c
  fir.c
  gnuplot_i.c
  linspace.c
  logalizit.c
  logspace.c
  nco.c
  piosbr_mod.c
  # piosbr.c
  reg_lin.c
  t_pnm.c
  tf.c
)

set_target_properties(dsac
  PROPERTIES SOVERSION 1.0
)
