ifeq "$(ARCH)" "GNU"
FC = mpifort
endif
ifeq "$(ARCH)" "INTEL"
FC = mpiifort
endif
ifeq "$(ARCH)" "NVIDIA"
FC = mpifort
endif
