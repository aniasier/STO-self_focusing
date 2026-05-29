# --- Variables ---
FC = gfortran
TARGET = STO.exe

# Your specific flags
FFLAGS = -Ofast -m64 -ffixed-line-length-none -fdefault-real-8 \
         -fdefault-double-8 -std=legacy -Wfatal-errors \
         -ffpe-summary=none -fopenmp -ffree-form -fpic

INCLUDES = -I$(MKLROOT)/include

# MKL Linking
MKL_LIBS = -Wl,--start-group \
           $(MKLROOT)/lib/intel64/libmkl_gf_lp64.so \
           $(MKLROOT)/lib/intel64/libmkl_gnu_thread.so \
           $(MKLROOT)/lib/intel64/libmkl_core.so \
           -Wl,--end-group \
           -lgomp -lpthread -lm -ldl

# Object files
OBJ = dielectric.o poisson_solver.o main.o

# --- Rules ---

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ) $(MKL_LIBS)

dielectric.o: dielectric.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c dielectric.f90

poisson_solver.o: poisson_solver.f90
	$(FC) $(FFLAGS) $(INCLUDES) -c poisson_solver.f90

main.o: main.f90 dielectric.o poisson_solver.o
	$(FC) $(FFLAGS) $(INCLUDES) -c main.f90

clean:
	rm -f *.o *.mod $(TARGET)