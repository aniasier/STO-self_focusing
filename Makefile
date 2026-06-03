# --- Variables ---
FC = gfortran
TARGET = STO.exe

SRC_DIR = src
OBJ_DIR = obj
MOD_DIR = mod

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

SRC = $(SRC_DIR)/constants.f90 \
      $(SRC_DIR)/utils.f90 \
      $(SRC_DIR)/dielectric.f90 \
      $(SRC_DIR)/poisson_solver.f90 \
      $(SRC_DIR)/main.f90

# Object files
OBJ = $(OBJ_DIR)/constants.o \
      $(OBJ_DIR)/utils.o \
      $(OBJ_DIR)/dielectric.o \
      $(OBJ_DIR)/poisson_solver.o \
      $(OBJ_DIR)/main.o

# --- Rules ---

all: directories $(TARGET)

directories:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(MOD_DIR)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) -o $(TARGET) $(OBJ) $(MKL_LIBS)

$(OBJ_DIR)/constants.o: $(SRC_DIR)/constants.f90
	$(FC) $(FFLAGS) $(INCLUDES) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/utils.o: $(SRC_DIR)/utils.f90 $(OBJ_DIR)/constants.o
	$(FC) $(FFLAGS) $(INCLUDES) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/dielectric.o: $(SRC_DIR)/dielectric.f90 $(OBJ_DIR)/constants.o
	$(FC) $(FFLAGS) $(INCLUDES) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/poisson_solver.o: $(SRC_DIR)/poisson_solver.f90 $(OBJ_DIR)/constants.o $(OBJ_DIR)/dielectric.o
	$(FC) $(FFLAGS) $(INCLUDES) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

$(OBJ_DIR)/main.o: $(SRC_DIR)/main.f90 $(OBJ_DIR)/constants.o $(OBJ_DIR)/dielectric.o $(OBJ_DIR)/poisson_solver.o $(OBJ_DIR)/utils.o
	$(FC) $(FFLAGS) $(INCLUDES) -J$(MOD_DIR) -I$(MOD_DIR) -c $< -o $@

DEBUG_FLAGS = -O0 -g -cpp -DDEBUG \
              -m64 -ffixed-line-length-none \
              -fdefault-real-8 -fdefault-double-8 \
              -std=legacy -Wfatal-errors \
              -ffpe-summary=none -fopenmp -ffree-form -fpic

debug:
	$(MAKE) clean
	$(MAKE) FFLAGS="$(DEBUG_FLAGS)"
clean:
	rm -rf $(OBJ_DIR) $(MOD_DIR) $(TARGET)

