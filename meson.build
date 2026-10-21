project('interpit', 'fortran', meson_version: '>=0.49')

lib_sources = files('src/bound.f90',   \
                    'src/common.f90',  \
                    'src/decomp.f90',  \
                    'src/grid.f90',    \
                    'src/interp.f90',  \
                    'src/io.f90',      \
                    'src/io_hdf5.f90', \
                    'src/param.f90')

# the sources are preprocessed (optional HDF5 support)
fc = meson.get_compiler('fortran')
cpp_flags = {'gcc': '-cpp', 'intel': '-fpp', 'intel-llvm': '-fpp', 'nvidia_hpc': '-Mpreprocess'}
add_project_arguments(cpp_flags.get(fc.get_id(), '-cpp'), language: 'fortran')

deps = [dependency('mpi', language: 'fortran', required: true)]
if get_option('hdf5')
  deps += dependency('hdf5', language: 'fortran', required: true)
  add_project_arguments('-D_USE_HDF5', language: 'fortran')
endif

lib = static_library('interpit', lib_sources, dependencies: deps)
executable('interpit', 'app/interpolate-fields.f90', link_with: lib, dependencies: deps)
executable('task-checkpoints', 'app/task-checkpoints.f90', link_with: lib, dependencies: deps)
