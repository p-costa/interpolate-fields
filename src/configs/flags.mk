ifeq ($(strip $(FFLAGS_DEBUG)),1)

ifeq "$(FCOMP)" "GNU"
FFLAGS := -O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
endif

ifeq "$(FCOMP)" "INTEL"
FFLAGS := -O0 -g -traceback -fpe0 -stand f18
endif

ifeq "$(FCOMP)" "NVIDIA"
FFLAGS := -O0 -g -traceback -Mstandard -Minform=inform -Mbackslash -Mbounds -Mchkptr -Mchkstk
endif
  
endif
ifeq ($(strip $(FFLAGS_DEBUG_MAX)),1)

ifeq "$(FCOMP)" "GNU"
FFLAGS := -O0 -Wall -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 -std=f2018
endif

endif

ifeq ($(strip $(FFLAGS_OPT)),1)

ifeq "$(FCOMP)" "GNU"
FFLAGS := -O3
endif

ifeq "$(FCOMP)" "INTEL"
FFLAGS := -O3
endif

ifeq "$(FCOMP)" "NVIDIA"
FFLAGS := -O3
endif
  
endif

ifeq ($(strip $(FFLAGS_OPT_MAX)),1)

ifeq "$(FCOMP)" "GNU"
FFLAGS := -Ofast -march=native
endif

ifeq "$(FCOMP)" "INTEL"
FFLAGS := -fast -xHost
endif

ifeq "$(FCOMP)" "NVIDIA"
FFLAGS := -fast
endif
  
endif

ifeq ($(strip $(CPP_DEBUG)),1)
DEFINES += -D_DEBUG
endif

ifeq ($(strip $(NON_UNIFORM_Z)),1)
DEFINES += -D_NON_UNIFORM_Z
endif
