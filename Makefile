SERVER_LIB_ROOT := /home/zxli_1/Lib_90_new
LOCAL_LIB_ROOT  := /home/yyk/Lib_90_new

DETECTED_LIB_ROOT := $(firstword $(wildcard $(SERVER_LIB_ROOT) $(LOCAL_LIB_ROOT)))
LIB_ROOT ?= $(if $(DETECTED_LIB_ROOT),$(DETECTED_LIB_ROOT),$(SERVER_LIB_ROOT))

MPIIFORT := $(shell command -v mpiifort 2>/dev/null)
MPIIFX   := $(shell command -v mpiifx   2>/dev/null)
MPIFORT  := $(shell command -v mpifort  2>/dev/null)
MPIF90   := $(shell command -v mpif90   2>/dev/null)
GFORTRAN := $(shell command -v gfortran 2>/dev/null)
IFORT    := $(shell command -v ifort    2>/dev/null)
IFX      := $(shell command -v ifx      2>/dev/null)

ifneq ($(filter command% environment,$(origin FC)),)
  FC_CMD := $(FC)
else ifneq ($(MPIIFORT),)
  ifneq ($(IFORT),)
    FC_CMD := $(MPIIFORT) -fc=ifort
  else ifneq ($(IFX),)
    FC_CMD := $(MPIIFORT) -fc=ifx
  else
    FC_CMD := $(MPIIFORT)
  endif
else ifneq ($(MPIIFX),)
  ifneq ($(IFX),)
    FC_CMD := $(MPIIFX) -fc=ifx
  else
    FC_CMD := $(MPIIFX)
  endif
else ifneq ($(MPIFORT),)
  FC_CMD := $(MPIFORT)
else ifneq ($(MPIF90),)
  FC_CMD := $(MPIF90)
else ifneq ($(GFORTRAN),)
  FC_CMD := $(GFORTRAN)
endif

ifneq ($(filter clean,$(MAKECMDGOALS)),clean)
  ifneq ($(filter help,$(MAKECMDGOALS)),help)
    ifeq ($(strip $(FC_CMD)),)
      $(error No Fortran compiler found. Set FC explicitly, for example: make FC=mpiifort)
    endif
  endif
endif

SRC_DIR := src
BUILD_DIR := build
OBJ_DIR := $(BUILD_DIR)/obj
MOD_DIR := $(BUILD_DIR)/mod
TARGET := $(BUILD_DIR)/bosonDQMC.out

FFLAGS ?= -O3 -I$(LIB_ROOT)/Modules -I$(MOD_DIR)
SUFFIX ?=
LDFLAGS ?= -mkl
PYTHON ?= python3
MPI_NP ?= 1
RUN_DIR ?= runs/examples/triangle_3x2
BENCHMARK_RUN_DIR ?= benchmarks/fixtures/triangle_3x2_mc_outputs
BENCHMARK_REFERENCE ?= benchmarks/references/triangle_3x2.json

ifneq (,$(findstring gfortran,$(FC_CMD)))
  MOD_FLAG ?= -J$(MOD_DIR)
else
  MOD_FLAG ?= -module $(MOD_DIR)
endif

LDLIBS := $(LIB_ROOT)/Modules/modules_90.a \
          $(LIB_ROOT)/MyEis/libeis.a \
          $(LIB_ROOT)/MyNag/libnag.a \
          $(LIB_ROOT)/MyLin/liblin.a \
          $(LIB_ROOT)/Ran/libran.a

SOURCES := \
  calc_basic.f90 \
  lattice.f90 \
  fields.f90 \
  non_interact.f90 \
  operator_Hubbard.f90 \
  model.f90 \
  process_matrix.f90 \
  multiply.f90 \
  obser_equal.f90 \
  obser_tau.f90 \
  stabilization.f90 \
  localU.f90 \
  dynamics.f90 \
  local_sweep.f90 \
  fourier_trans.f90 \
  main.f90

OBJECTS := $(addprefix $(OBJ_DIR)/,$(SOURCES:.f90=.o))

.PHONY: all build clean print-config run-example benchmark benchmark-ed help

all: build

build: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC_CMD) $(LDFLAGS) -o $@ $(OBJECTS) $(LDLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR) $(MOD_DIR)
	$(FC_CMD) $(SUFFIX) -c $(FFLAGS) $(MOD_FLAG) -o $@ $<

$(OBJ_DIR) $(MOD_DIR):
	mkdir -p $@

print-config:
	@echo "FC=$(FC_CMD)"
	@echo "LIB_ROOT=$(LIB_ROOT)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"
	@echo "MOD_FLAG=$(MOD_FLAG)"
	@echo "TARGET=$(TARGET)"
	@echo "RUN_DIR=$(RUN_DIR)"
	@echo "BENCHMARK_RUN_DIR=$(BENCHMARK_RUN_DIR)"
	@echo "BENCHMARK_REFERENCE=$(BENCHMARK_REFERENCE)"

run-example: build
	bash scripts/run_local.sh $(RUN_DIR) $(MPI_NP)

benchmark:
	$(PYTHON) benchmarks/compare.py --reference $(BENCHMARK_REFERENCE) --run-dir $(BENCHMARK_RUN_DIR)

benchmark-ed:
	cd benchmarks/ed && cp params_triangle_3x2.txt params.txt && $(PYTHON) EDtriangle_symm_NEblock.py

clean:
	rm -rf $(BUILD_DIR)
	rm -f *.mod *.o *.lst *.opt-report bosonDQMC.out

help:
	@echo "Targets:"
	@echo "  make build          Build build/bosonDQMC.out"
	@echo "  make run-example    Run the default example with mpirun"
	@echo "  make benchmark      Compare fixture or run output to ED reference"
	@echo "  make benchmark-ed   Optionally recompute ED reference in benchmarks/ed"
	@echo "  make clean          Remove generated build artifacts"
	@echo "  make print-config   Print compiler and path configuration"
