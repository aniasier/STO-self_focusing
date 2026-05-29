#!/bin/bash
gfortran  -Ofast  -m64  -ffixed-line-length-none -fdefault-real-8  -fdefault-double-8  -std=legacy \
          -Wfatal-errors -ffpe-summary=none -fopenmp -ffree-form \
          -I${MKLROOT}/include  \
          -fpic STO_self_focusing_simple_model_1D.f90 \
	  -Wl,--start-group \
	  ${MKLROOT}/lib/intel64/libmkl_gf_lp64.so \
          ${MKLROOT}/lib/intel64/libmkl_gnu_thread.so \
          ${MKLROOT}/lib/intel64/libmkl_core.so \
          -Wl,--end-group \
          -lgomp -lpthread -lm -ldl
