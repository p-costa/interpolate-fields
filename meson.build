project('interpit', 'fortran', meson_version: '>=0.49')

sources = files('src/bound.f90',  \
                'src/common.f90', \
                'src/io.f90',     \
                'app/interpolate-fields.f90')

mpif    = dependency('mpi', language: 'fortran', required: true)

executable('interpit', sources, dependencies: [mpif])
