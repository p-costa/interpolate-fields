ifeq "$(FCOMP)" "GNU"
FC = mpifort
endif
ifeq "$(FCOMP)" "INTEL"
FC = mpiifort
endif
ifeq "$(FCOMP)" "NVIDIA"
FC = mpifort
endif
