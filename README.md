# interpolate-fields
This program interpolates DNS data to a new grid using linear interpolation. It can run and store the interpolated data in a massively-parallel setting using MPI I/O. It handles the grids and file formats of [CaNS](https://github.com/CaNS-World/CaNS) and [CaNS-EIGEN](https://github.com/CaNS-World/CaNS-EIGEN): each of the input and output grids can be uniform, non-uniform along z (CaNS), or non-uniform along all directions (CaNS-EIGEN). Information about the boundary conditions of the system must be prescribed too. It is a cleaned up version of the routines used to create the initial conditions for particle-laden turbulent channel transport with 10^10 Eulerian grid points.

## Building
Build with `make` (compiler and flags in `src/configs/build.conf`); the executables are `exe/interpit` (the interpolation) and `exe/task-checkpoints` (merging and splitting checkpoints saved with one file per task, see below). The program can also be built with [fpm](https://fpm.fortran-lang.org) (`FPM_FC=mpifort fpm build`) or [meson](https://mesonbuild.com) (`meson setup build && ninja -C build`). There are no compile-time options for the precision or the grids: the precision of the input is detected from the files, and that of the output is set at runtime.

HDF5 input and output (the format of CaNS and CaNS-EIGEN with `io_backend = 'hdf5'`) is optional, and needs parallel HDF5 with the Fortran interface: set `USE_HDF5=1` in `src/configs/build.conf` (or `make USE_HDF5=1`); the include and library flags default to those of the OpenMPI build of Debian/Ubuntu (package `libhdf5-openmpi-dev`), and can be set with `HDF5_INCS` and `HDF5_LIBS` (see `src/configs/libs.mk`). With meson, `meson setup build -Dhdf5=true`; with fpm, e.g. `FPM_FC=mpifort fpm build --flag "-D_USE_HDF5 -I/usr/include/hdf5/openmpi" --link-flag "-Wl,--no-as-needed -lhdf5_openmpi_fortran -lhdf5_openmpi"`. Without HDF5 support, the program stops with a message if an HDF5 file is requested.

## Usage
```
mpirun -n <ntasks> ./exe/interpit [namelist file]
```
The case is set in a namelist file (default: `interp.nml` in the working directory), described below. The program:

* reads, for each variable, the input checkpoint, in single or double precision (detected from the file size, or from the HDF5 dataset), or several checkpoints tiled periodically into a larger domain (see `&tiling`);
* imposes the boundary conditions on the ghost points, as CaNS does;
* optionally, subtracts the volume-weighted mean or a constant offset (see `&operations`);
* interpolates the field onto the output grid by trilinear interpolation, in double precision, with weights from the actual coordinates of the variable's location; output points that coincide with input points reproduce them exactly;
* optionally, subtracts the residual mean of the interpolated field, and adds a uniform velocity to u, v and w (change of reference frame);
* writes the output in the current CaNS checkpoint format, one file per variable, `<output_prefix>_<var>.bin` (the field followed by the time and the time step) or `<output_prefix>_<var>.h5`, in the chosen precision.

Each MPI task owns a block of the output grid and reads exactly the block of the input field that it needs (plus, when the output grid is coarser, the input points between the blocks, for the statistics), wrapping around periodic directions, so that any process grid can be used, with z decomposed or not, and without divisibility constraints. The variables are processed one at a time: the memory footprint per task is that of one input block, plus one output block in double precision (and, for single-precision output, a conversion buffer of at most 64 MiB). The program prints, per variable, the minimum, maximum and volume-weighted mean of the input and output fields, and, at the end, the time spent reading, interpolating and writing, and the memory of the field arrays.

## Input namelist

### `&interp`
| parameter | default | meaning |
|---|---|---|
| `dims(1:3)` | `1, 0, 0` | MPI process grid; zeros are filled by `MPI_DIMS_CREATE` (the default gives x-aligned pencils, as the files are stored with x fastest) |
| `variables(:)` | `'u', 'v', 'w', 'p'` | variables to process, any subset or order: `u`, `v`, `w`, `p`, scalars `s_001`, `s_002`, ... as in CaNS, or other variables described in a `&variable` group |
| `input_format` | `'per-variable'` | `'per-variable'`: current CaNS checkpoints, `<input_prefix>_<var>.bin`; `'hdf5'`: HDF5 checkpoints, `<input_prefix>_<var>.h5`; `'legacy'`: single-file checkpoint `input_file` with u, v, w, p and any scalars, one after the other, then the time and time step |
| `input_prefix` | | prefix of the per-variable input files, e.g. `'data/fld'` for `data/fld_u.bin`, ... |
| `input_file` | | legacy single-file checkpoint |
| `output_format` | `'per-variable'` | `'per-variable'` (`<output_prefix>_<var>.bin`) or `'hdf5'` (`<output_prefix>_<var>.h5`), independently of the input |
| `output_prefix` | | prefix of the output files (which must differ from the input files) |
| `output_precision` | `'input'` | `'single'`, `'double'`, or `'input'` (that of each input file) |
| `is_compress_output` | `F` | HDF5 output: chunked and gzip-compressed datasets, as in CaNS (`is_use_compression`) |
| `is_reset_time` | `F` | set the time and time step of the output to zero (e.g., for an initial condition); otherwise, they are those of the input |

### `&grid_in` and `&grid_out`
The input and output grids, set independently, per direction:

| parameter | meaning |
|---|---|
| `ng(1:3)` | number of grid cells |
| `l(1:3)` | domain length; needed for uniform directions only (if given with a grid file, it is checked against it) |
| `grid_files(1:3)` | grid file per direction; `''` (default) for a uniform grid |

The grid files are those written by CaNS (`grid.bin`, along z) and CaNS-EIGEN (`grid_x.bin`, `grid_y.bin`, `grid_z.bin`), in single or double precision (detected from the file size), or their HDF5 versions (`grid.h5`, ...; recognized by the extension); the number of cells is checked against the file. The domain lengths of the input and output grids must agree in every direction.

### `&bc`
Boundary conditions of `u`, `v`, `w`, `p` and the scalars, with the same names and conventions as in CaNS: `cbcvel(0:1,1:3,1:3)` and `bcvel(0:1,1:3,1:3)` for the velocity, `cbcpre(0:1,1:3)` and `bcpre(0:1,1:3)` for the pressure, and `cbcscal(0:1,1:3,iscal)` and `bcscal(0:1,1:3,iscal)` for scalar `iscal` (`s_<iscal>`); types `'P'` (periodic), `'D'` (Dirichlet) or `'N'` (Neumann), values default to zero. The conditions of every variable processed must be given; the boundary conditions refer to the input fields, and are used for the ghost points of the interpolation. The periodic directions must be the same for all variables.

### `&variable`
One group per variable that is not `u`, `v`, `w`, `p` or a scalar (e.g., a vorticity component, located on cell edges), with its name, location and boundary conditions:

| parameter | meaning |
|---|---|
| `name` | name of the variable, as in `variables` and in the file names |
| `is_staggered(1:3)` | location per direction: on the cell faces (`T`) or at the cell centres (`F`, default) |
| `cbc(0:1,1:3)`, `bc(0:1,1:3)` | boundary conditions, as in `&bc` |

The locations of the other variables follow from their names: `u` is staggered along x (`'scc'`), `v` along y (`'csc'`), `w` along z (`'ccs'`), and `p` and the scalars are cell-centred (`'ccc'`).

### `&operations`
Field operations (optional), per variable in the order of `variables` (with the default variables, `u`, `v`, `w`, `p` are 1 to 4):

| parameter | default | meaning |
|---|---|---|
| `is_remove_mean(:)` | `F` | subtract the volume-weighted mean: that of the input before the interpolation, in double precision (so that the interpolation works with O(1) values, e.g. for a pressure whose mean has drifted to large values in a long run, which single precision cannot resolve around its fluctuations), and then the small residual mean of the interpolated field (interpolation does not preserve the mean exactly on non-uniform grids), so that the output mean is zero to round-off |
| `offset(:)` | `0.` | constant subtracted from the input (exclusive with `is_remove_mean`) |
| `vel_frame(1:3)` | `0., 0., 0.` | uniform velocity added to `u`, `v` and `w` after the interpolation (change of reference frame) |

The mean is weighted with the control volumes of the variable's location: the cells along cell-centred directions, and the spacings between cell centres along staggered directions (at a bounded upper boundary, the half cell next to it; the lower boundary face is not stored in the files). The boundary conditions refer to the input fields: a subtracted mean or offset shifts the Dirichlet values alike (the ghost points are shifted with the field). After a change of reference frame, the CaNS run that uses the output must set its wall velocities and `velf` for the new frame: e.g., for a channel in the frame moving with the bulk velocity `ub` of a lab-frame run, `vel_frame(1) = -ub`, and the new run has `bcvel(0:1,3,1) = -ub, -ub` and `velf(1) = 0.`. The value subtracted or added is printed for each variable and operation.

### `&tiling`
Periodic tiling (optional): the input domain is made of `ntiles(1) x ntiles(2) x ntiles(3)` copies of the input grid (`&grid_in` describes one tile), along periodic directions only, each filled with a checkpoint (e.g., independent snapshots of a run, which avoids the artificial periodicity of identical copies) rolled by a whole number of cells; the output grid must span the tiled domain (`ntiles*l`). The tiled field is not written: the program reads the input blocks from the files of the tiles directly.

| parameter | default | meaning |
|---|---|---|
| `ntiles(1:3)` | `1, 1, 1` | number of tiles per direction |
| `tile_inputs(:)` | `input_prefix` (or `input_file`) | input of each tile, x fastest (tile `(tx,ty,tz)`, from 0, is number `1 + tx + ntiles(1)*(ty + ntiles(2)*tz)`): the prefix of its files, or its legacy file; all on the same grid and in the same precision |
| `tile_shifts(1:3,:)` | `0` | periodic shift of each tile, in input cells: point `i` of the tile holds point `i + shift` of its checkpoint (as `numpy.roll(f, -shift)`); whole cells keep the staggering, and require a uniform input grid along that direction |

With `is_remove_mean`, the mean of each tile is removed separately before the interpolation (different snapshots have different pressure means, which would leave jumps between the tiles). The time and time step of the output are those of the first tile. See the example [`channel-tiling`](examples/channel-tiling/interp.nml).

### Boundary conditions and ghost points
As in CaNS (`set_bc`): along a periodic direction, the field wraps around; along a bounded direction, a cell-centred variable gets ghost points mirrored about the boundary (Dirichlet: `f(0) = 2*bc - f(1)`; Neumann: `f(0) = f(1) - dzc(0)*bc`), and a staggered variable gets its values on the boundary faces (Dirichlet: `f(0) = bc` and `f(n) = bc`; Neumann: `f(0) = f(1) - dzf(0)*bc` and `f(n) = f(n-1) + dzf(n)*bc`). Fields that are linear in space and consistent with these conditions are therefore interpolated exactly, also near the boundaries. Immersed boundaries get no special treatment: solid regions are interpolated like the fluid. The boundary conditions only serve to set the ghost points for the interpolation, with one value per boundary: for spatially varying Dirichlet values (e.g., an inflow profile, whose boundary values are not stored in the checkpoints), a Neumann condition with zero gradient (`'N'`, `0.`) gives the closest ghost points, and the run restarted from the output imposes its own conditions.

### HDF5 files
As written and read by CaNS and CaNS-EIGEN (`io_field_hdf5`): the field in dataset `fields/<var>` (dimensions `ng`, x fastest, in single or double precision), the time and the time step in `meta/time` and `meta/istep`, and the coordinates of the variable's location, ghost points included, in `grid/x`, `grid/y` and `grid/z`; with `is_compress_output = T`, the field is stored in chunks of at most 2 MiB aligned with the blocks of the tasks, gzip-compressed (level 4, with the shuffle filter). The grid files `grid_<dir>.h5` hold the cell centres in `rc` and the cell faces in `rf`. All tasks read and write their blocks collectively, as hyperslabs.

## Examples
The folder [`examples`](examples) has input files for typical cases, each described at the top of the file: copy the `interp.nml` of a case into the run directory, and adapt the file names, grid sizes and process grid.

| example | case |
|---|---|
| [`channel-refine`](examples/channel-refine/interp.nml) | CaNS turbulent channel (uniform x and y, non-uniform z) onto a finer grid, from a single-precision legacy checkpoint, as an initial condition in a moving reference frame |
| [`duct`](examples/duct/interp.nml) | CaNS-EIGEN periodic duct (non-uniform along all directions) with a passive scalar |
| [`duct-hdf5`](examples/duct-hdf5/interp.nml) | the same, with HDF5 input and compressed HDF5 output |
| [`channel-tiling`](examples/channel-tiling/interp.nml) | CaNS turbulent channel: a box twice as long and wide, from four snapshots (periodic tiling) |

A variable on the cell edges, e.g. the z vorticity computed in post-processing, can be added to `variables` with a group:
```fortran
&variable
name = 'omz', is_staggered(1:3) = T, T, F
cbc(0:1,1:3) = 'P','P',  'D','D',  'N','N'
/
```

## Merging and splitting per-task checkpoints
`exe/task-checkpoints [namelist file]` (default: `tasks.nml`) converts between a checkpoint saved with one file per task and a single checkpoint: `mode = 'merge'` writes the per-task files into one checkpoint, e.g. to interpolate it, and `mode = 'split'` writes a checkpoint, e.g. the output of `interpit`, into per-task files. A per-task file holds the variables `task_variables`, one after the other, over the block of its task (without halos) in a `task_dims(1) x task_dims(2) x task_dims(3)` decomposition of the grid, followed by the time and the time step, in single or double precision (detected from the file size); a pencil decomposition has a 1 in `task_dims`, e.g. x-aligned pencils `1, py, pz`. Along each direction, the first `mod(ng,task_dims)` tasks get one more point, as in 2DECOMP&FFT. The tasks of the program process the files one round at a time (ranks `myid`, `myid+ntasks`, ...), so that their number is independent of `task_dims`. Namelist `&tasks`:

| parameter | default | meaning |
|---|---|---|
| `mode` | | `'merge'` or `'split'` |
| `ng(1:3)` | | number of grid cells |
| `task_dims(1:3)` | | decomposition of the per-task files |
| `task_order` | `'row-major'` | rank order of the per-task files: `'row-major'` (the last coordinate varies fastest, as for an MPI Cartesian topology and cuDecomp; for x-aligned pencils, rank `iy*pz + iz`) or `'column-major'` (the first one fastest; rank `iy + py*iz`) |
| `task_prefix`, `task_suffix`, `task_digits` | `''`, `'.bin'`, `5` | per-task file of rank `r`: `<task_prefix><r, with task_digits digits><task_suffix>` |
| `task_variables(:)` | `'u', 'v', 'w', 'p'` | variables of the per-task files, in order |
| `task_precision` | `'input'` | split: precision of the per-task files (`'single'`, `'double'`, or `'input'`: that of the checkpoint) |
| `checkpoint_format` | `'per-variable'` | the single checkpoint: `'per-variable'` (`<checkpoint_prefix>_<var>.bin`), `'legacy'` (`checkpoint_file`, with `task_variables` in order) or `'hdf5'` (`<checkpoint_prefix>_<var>.h5`, without coordinates) |
| `checkpoint_prefix`, `checkpoint_file` | | prefix of the per-variable files, or legacy file |
| `checkpoint_precision` | `'input'` | merge: precision of the checkpoint (`'single'`, `'double'`, or `'input'`: that of the per-task files) |
| `variables(:)` | all `task_variables` | merge: the variables written; with a subset, only these are written (into an existing, complete file for a legacy checkpoint), e.g. to merge one field again with another offset |
| `offset(:)` | `0.` | constant subtracted from each variable, in the order of `task_variables`, in double precision before the conversion (e.g. a large mean pressure, which single precision cannot resolve) |
| `is_compress` | `F` | merge: compression of an HDF5 checkpoint |
| `is_reset_time` | `F` | set the time and time step to zero |

When merging, the time step must be the same in all files; the time may not (it is read per task at each restart of a run saving one file per task), so an inconsistent time only gives a warning, and the time of rank 0 is kept. For example, merging a checkpoint of 1024 x-aligned pencils (cuDecomp rank order) into a single-precision per-variable checkpoint, with a large mean pressure subtracted, and then merging the pressure again with another offset:
```fortran
&tasks
mode = 'merge', ng(1:3) = 11200, 12288, 2374
task_dims(1:3) = 1, 16, 64, task_order = 'row-major'
task_prefix = 'data/fld_myid_', task_suffix = '_0005.bin'
checkpoint_prefix = 'data/fld_i', checkpoint_precision = 'single'
offset(4) = 949381.
/
```
```fortran
&tasks
mode = 'merge', ng(1:3) = 11200, 12288, 2374
task_dims(1:3) = 1, 16, 64, task_order = 'row-major'
task_prefix = 'data/fld_myid_', task_suffix = '_0005.bin'
checkpoint_prefix = 'data/fld_i', checkpoint_precision = 'single'
variables = 'p', offset(4) = 949381.25
/
```
and splitting the output of an interpolation into 4 x 4 x 2 blocks:
```fortran
&tasks
mode = 'split', ng(1:3) = 512, 384, 256
task_dims(1:3) = 4, 4, 2, task_order = 'column-major'
task_prefix = 'data/fld_', task_suffix = '_0000.bin'
checkpoint_prefix = 'data/fld_o'
/
```

## Tests
`tests/run_tests.py` (Python 3 with NumPy) runs the program on small synthetic cases, with several process grids: identity on the same grid, exactness for linear fields at all variable locations, bitwise process-grid independence and agreement with an independent NumPy implementation, precision conversion, legacy and per-variable input, field operations, periodic tiling, merging and splitting per-task checkpoints, HDF5 input and output (with [h5py](https://www.h5py.org), for a build with HDF5; otherwise, the error message), the defaults, and the input checks. Run it after building, e.g. `python3 tests/run_tests.py` (see `--help`; the MPI launcher can be set with the environment variable `MPIRUN`).
