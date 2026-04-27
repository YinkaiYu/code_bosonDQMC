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
  ifeq ($(strip $(FC_CMD)),)
    $(error No Fortran compiler found. Set FC explicitly, for example: make FC=mpiifort)
  endif
endif

FFLAGS ?= -c -O3 -I$(LIB_ROOT)/Modules
SUFFIX ?=
LDFLAGS ?= -mkl

LDLIBS := $(LIB_ROOT)/Modules/modules_90.a \
          $(LIB_ROOT)/MyEis/libeis.a \
          $(LIB_ROOT)/MyNag/libnag.a \
          $(LIB_ROOT)/MyLin/liblin.a \
          $(LIB_ROOT)/Ran/libran.a
#         $(LIB_ROOT)/LaPack/lapack.a \
#         $(LIB_ROOT)/Blas/libblas.a

.PHONY: all clean print-config

all:
	$(MAKE) -f Compile FC="$(FC_CMD)" LF="$(LDFLAGS)" FLAGS="$(FFLAGS)" LIBS="$(LDLIBS)" SUFFIX="$(SUFFIX)"

print-config:
	@echo "FC=$(FC_CMD)"
	@echo "LIB_ROOT=$(LIB_ROOT)"
	@echo "FFLAGS=$(FFLAGS)"
	@echo "LDFLAGS=$(LDFLAGS)"

clean:
	$(MAKE) -f Compile clean
	@rm -f *.mod *.lst *.opt-report
