# interpolate-fields
This program interpolates DNS data to a new grid using linear interpolation. It can run and store the interpolated data in a massively-parallel setting using MPI I/O. It handles the grids and file formats of [CaNS](https://github.com/CaNS-World/CaNS) and [CaNS-EIGEN](https://github.com/CaNS-World/CaNS-EIGEN): each of the input and output grids can be uniform, non-uniform along z (CaNS), or non-uniform along all directions (CaNS-EIGEN). Information about the boundary conditions of the system must be prescribed too. It is a cleaned up version of the routines used to create the initial conditions for particle-laden turbulent channel transport with 10^10 Eulerian grid points.

## Building
Build with `make` (compiler and flags in `src/configs/build.conf`); the executable is `exe/interpit`. The program can also be built with [fpm](https://fpm.fortran-lang.org) (`FPM_FC=mpifort fpm build`) or [meson](https://mesonbuild.com) (`meson setup build && ninja -C build`). There are no compile-time options for the precision or the grids: the precision of the input is detected from the files, and that of the output is set at runtime.

## Usage
```
mpirun -n <ntasks> ./exe/interpit [namelist file]
```
The case is set in a namelist file (default: `interp.nml` in the working directory), described below. The program:

* reads, for each variable, the input checkpoint, in single or double precision (detected from the file size);
* imposes the boundary conditions on the ghost points, as CaNS does;
* interpolates the field onto the output grid by trilinear interpolation, in double precision, with weights from the actual coordinates of the variable's location; output points that coincide with input points reproduce them exactly;
* writes the output in the current CaNS checkpoint format, one file per variable, `<output_prefix>_<var>.bin` (the field followed by the time and the time step), in the chosen precision.

Each MPI task owns a block of the output grid and reads exactly the block of the input field that it needs (plus, when the output grid is coarser, the input points between the blocks, for the statistics), wrapping around periodic directions, so that any process grid can be used, with z decomposed or not, and without divisibility constraints. The variables are processed one at a time: the memory footprint per task is that of one input block, plus one output block in double precision (and, for single-precision output, a conversion buffer of at most 64 MiB). The program prints, per variable, the minimum, maximum and volume-weighted mean of the input and output fields, and, at the end, the time spent reading, interpolating and writing, and the memory of the field arrays.

## Input namelist

### `&interp`
| parameter | default | meaning |
|---|---|---|
| `dims(1:3)` | `1, 0, 0` | MPI process grid; zeros are filled by `MPI_DIMS_CREATE` (the default gives x-aligned pencils, as the files are stored with x fastest) |
| `variables(:)` | `'u', 'v', 'w', 'p'` | variables to process, any subset or order: `u`, `v`, `w`, `p`, scalars `s_001`, `s_002`, ... as in CaNS, or other variables described in a `&variable` group |
| `input_format` | `'per-variable'` | `'per-variable'`: current CaNS checkpoints, `<input_prefix>_<var>.bin`; `'legacy'`: single-file checkpoint `input_file` with u, v, w, p and any scalars, one after the other, then the time and time step |
| `input_prefix` | | prefix of the per-variable input files, e.g. `'data/fld'` for `data/fld_u.bin`, ... |
| `input_file` | | legacy single-file checkpoint |
| `output_prefix` | | prefix of the output files, `<output_prefix>_<var>.bin` (must differ from the input) |
| `output_precision` | `'input'` | `'single'`, `'double'`, or `'input'` (that of each input file) |
| `is_reset_time` | `F` | set the time and time step of the output to zero (e.g., for an initial condition); otherwise, they are those of the input |

### `&grid_in` and `&grid_out`
The input and output grids, set independently, per direction:

| parameter | meaning |
|---|---|
| `ng(1:3)` | number of grid cells |
| `l(1:3)` | domain length; needed for uniform directions only (if given with a grid file, it is checked against it) |
| `grid_files(1:3)` | grid file per direction; `''` (default) for a uniform grid |

The grid files are those written by CaNS (`grid.bin`, along z) and CaNS-EIGEN (`grid_x.bin`, `grid_y.bin`, `grid_z.bin`), in single or double precision (detected from the file size); the number of cells is checked against the file. The domain lengths of the input and output grids must agree in every direction.

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

### Boundary conditions and ghost points
As in CaNS (`set_bc`): along a periodic direction, the field wraps around; along a bounded direction, a cell-centred variable gets ghost points mirrored about the boundary (Dirichlet: `f(0) = 2*bc - f(1)`; Neumann: `f(0) = f(1) - dzc(0)*bc`), and a staggered variable gets its values on the boundary faces (Dirichlet: `f(0) = bc` and `f(n) = bc`; Neumann: `f(0) = f(1) - dzf(0)*bc` and `f(n) = f(n-1) + dzf(n)*bc`). Fields that are linear in space and consistent with these conditions are therefore interpolated exactly, also near the boundaries. Immersed boundaries get no special treatment: solid regions are interpolated like the fluid.

## Examples

### CaNS turbulent channel (uniform x and y, non-uniform z)
Interpolation of a channel at Re_tau = 5200 onto a finer grid along y and z, from a single-precision legacy checkpoint into single-precision per-variable files, to be used as an initial condition:
```fortran
&interp
dims(1:3) = 1, 4, 30
input_format = 'legacy', input_file = 'data/fld_i.bin'
output_prefix = 'data/fld_o', output_precision = 'single'
is_reset_time = T
/
&grid_in
ng(1:3) = 3456, 2970, 976, l(1:3) = 6.283185307179586, 3.141592653589793, 2.
grid_files(1:3) = '', '', 'data/grid_i.bin'
/
&grid_out
ng(1:3) = 3456, 3840, 1648, l(1:3) = 6.283185307179586, 3.141592653589793, 2.
grid_files(1:3) = '', '', 'data/grid_o.bin'
/
&bc
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1)  = 0.,0.,  0.,0.,  -1.,-1.
/
```
The wall velocity `bcvel(:,3,1) = -1` is that of a run in a reference frame moving with the bulk velocity.

### CaNS-EIGEN periodic duct (non-uniform along all directions)
Interpolation of all fields, including a passive scalar, of a duct (periodic along x, walls along y and z), with the grid files written by CaNS-EIGEN for the input and output runs:
```fortran
&interp
variables = 'u', 'v', 'w', 'p', 's_001'
input_prefix = 'input/data/fld', output_prefix = 'output/data/fld'
/
&grid_in
ng(1:3) = 256, 192, 192
grid_files(1:3) = 'input/data/grid_x.bin', 'input/data/grid_y.bin', 'input/data/grid_z.bin'
/
&grid_out
ng(1:3) = 384, 288, 288
grid_files(1:3) = 'output/data/grid_x.bin', 'output/data/grid_y.bin', 'output/data/grid_z.bin'
/
&bc
cbcvel(0:1,1:3,1) = 'P','P',  'D','D',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'D','D',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'D','D',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'N','N',  'N','N'
cbcscal(0:1,1:3,1) = 'P','P',  'N','N',  'D','D'
bcscal(0:1,1:3,1)  = 0.,0.,  0.,0.,  0.,1.
/
```
The output grid files can be generated by running CaNS-EIGEN with the output grid for zero time steps (`nstep = 0`); with `output_prefix = 'output/data/fld'`, the output files are those that CaNS-EIGEN reads when restarting in `output/`. A variable on the cell edges, e.g. the z vorticity computed in post-processing, can be added to `variables` with:
```fortran
&variable
name = 'omz', is_staggered(1:3) = T, T, F
cbc(0:1,1:3) = 'P','P',  'D','D',  'N','N'
/
```

### Namelists of earlier versions
The namelist was redesigned, and earlier `interp.nml` files need updating: `l`, `ni` and `no` become `ng` and `l` in `&grid_in` and `&grid_out`; `input_grid_file` and `output_grid_file` become `grid_files(3)` there; the boundary conditions move to `&bc`; `input_file` (legacy single file) needs `input_format = 'legacy'`; `output_file_prefix` becomes `output_prefix`; `dims` now has three entries, any of which may be decomposed; and `is_reset_time` defaults to `F`.

## Tests
`tests/run_tests.py` (Python 3 with NumPy) runs the program on small synthetic cases, with several process grids: identity on the same grid, exactness for linear fields at all variable locations, bitwise process-grid independence and agreement with an independent NumPy implementation, precision conversion, legacy and per-variable input, the defaults, and the input checks. Run it after building, e.g. `python3 tests/run_tests.py` (see `--help`; the MPI launcher can be set with the environment variable `MPIRUN`).
