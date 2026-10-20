#
# HDF5 (`USE_HDF5=1` in build.conf): parallel HDF5 with the Fortran interface; the defaults are those of the
# OpenMPI build of Debian/Ubuntu (package libhdf5-openmpi-dev); elsewhere, set HDF5_INCS and HDF5_LIBS, e.g.
# `make HDF5_INCS=-I$HDF5_DIR/include HDF5_LIBS="-L$HDF5_DIR/lib -lhdf5_fortran -lhdf5"`
#
ifeq ($(strip $(USE_HDF5)),1)
DEFINES += -D_USE_HDF5
HDF5_INCS ?= -I/usr/include/hdf5/openmpi
HDF5_LIBS ?= -lhdf5_openmpi_fortran -lhdf5_openmpi
INCS += $(HDF5_INCS)
LIBS += $(HDF5_LIBS)
endif
