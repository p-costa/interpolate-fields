ifeq ($(strip $(DEBUG)),1)

ifeq "$(ARCH)" "GNU"
FFLAGS := -O0 -g -fbacktrace -Wall -Wextra -pedantic -fcheck=all -finit-real=snan -ffpe-trap=invalid -std=f2018
endif

ifeq "$(ARCH)" "INTEL"
FFLAGS := -O0 -g -traceback -fpe0 -stand f18
endif

ifeq "$(ARCH)" "NVIDIA"
FFLAGS := -O0 -g -traceback -Mstandard -Minform=inform -Mbackslash -Mbounds -Mchkptr -Mchkstk
endif
  
endif
ifeq ($(strip $(DEBUG_MAX)),1)

ifeq "$(ARCH)" "GNU"
FFLAGS := -O0 -Wall -Wextra -Wimplicit-interface -Wno-unused-function -fPIC -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -finit-real=snan -finit-integer=-99999999 -std=f2018
endif

endif

ifeq ($(strip $(OPT)),1)

ifeq "$(ARCH)" "GNU"
FFLAGS := -O3
endif

ifeq "$(ARCH)" "INTEL"
FFLAGS := -O3
endif

ifeq "$(ARCH)" "NVIDIA"
FFLAGS := -O3
endif
  
endif

ifeq ($(strip $(OPT_MAX)),1)

ifeq "$(ARCH)" "GNU"
FFLAGS := -Ofast -march=native
endif

ifeq "$(ARCH)" "INTEL"
FFLAGS := -fast -xHost
endif

ifeq "$(ARCH)" "NVIDIA"
FFLAGS := -fast
endif
  
endif
