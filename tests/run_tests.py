#!/usr/bin/env python3
"""Tests of interpolate-fields on small synthetic cases.

usage: run_tests.py [--exe EXE] [--workdir DIR] [--keep] [test ...]

Each test writes CaNS-format grid and checkpoint files, runs the program with several MPI process grids,
and checks the output: identity on the same grid, exactness for linear fields, process-grid independence,
precision conversion, legacy and per-variable input, variables on edges and corners, a comparison with an
independent NumPy implementation, field operations (mean removal, offsets, frame change), HDF5 input and output (with h5py, if the program is built with HDF5;
otherwise, the error message), and the input checks. The MPI launcher is taken from the environment
variable MPIRUN (default: 'mpirun --oversubscribe').
"""
import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile

import numpy as np

try:
    import h5py
except ImportError:
    h5py = None

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MPIRUN = shlex.split(os.environ.get('MPIRUN', 'mpirun --oversubscribe'))
DEFAULT_LOC = {'u': 'scc', 'v': 'csc', 'w': 'ccs', 'p': 'ccc'}


# ----------------------------------------------------------------------------------------------------------
# grids and files in the CaNS formats
# ----------------------------------------------------------------------------------------------------------

class Grid1d:
    """faces xf[0..n] (xf[0] = 0, xf[n] = l) and centres xc[1..n] (xc[0] unused) of a 1D grid"""

    def __init__(self, xf, xc=None):
        self.xf = np.asarray(xf, dtype=np.float64)
        self.n = len(self.xf) - 1
        self.l = self.xf[-1]
        if xc is None:
            xc = np.concatenate([[np.nan], 0.5*(self.xf[1:] + self.xf[:-1])])
        self.xc = np.asarray(xc, dtype=np.float64)

    def with_ghosts(self, periodic):
        """coordinates of the faces and centres for the indices 0..n+1, as in the program"""
        n, l = self.n, self.l
        xf = np.empty(n + 2)
        xc = np.empty(n + 2)
        xf[:n+1] = self.xf
        xc[1:n+1] = self.xc[1:]
        if periodic:
            xc[0], xc[n+1], xf[n+1] = xc[n] - l, xc[1] + l, xf[1] + l
        else:
            xc[0], xc[n+1], xf[n+1] = -xc[1], 2*l - xc[n], 2*l - xf[max(n-1, 0)]
        return xf, xc


def uniform(n, l):
    return Grid1d(np.arange(n + 1)*l/n, np.concatenate([[np.nan], (np.arange(1, n + 1) - 0.5)*l/n]))


def stretched(n, l, gamma=1.5, periodic=False, seed=0):
    """clustered at both ends (walls), or smoothly and randomly varying spacing (periodic)"""
    s = np.arange(n + 1)/n
    if periodic:
        rng = np.random.default_rng(seed)
        a = rng.uniform(-0.25, 0.25, 3)
        xf = s + sum(a[k]*np.sin(2*np.pi*(k + 1)*s)/(2*np.pi*(k + 1))*0.9 for k in range(3))
    else:
        xf = 0.5*(1 + np.tanh(gamma*(2*s - 1))/np.tanh(gamma))
    xf = xf*l
    xf[0], xf[-1] = 0., l
    return Grid1d(xf)


def write_grid(fname, g, dtype=np.float64):
    """CaNS grid file: dzc(1:n), dzf(1:n), zc(1:n), zf(1:n)"""
    xf, xc = g.with_ghosts(False)
    dzf = np.diff(g.xf)
    dzc = xc[2:] - xc[1:-1]
    np.concatenate([dzc, dzf, g.xc[1:], g.xf[1:]]).astype(dtype).tofile(fname)


def write_field(fname, a, dtype=np.float64, time=1.5, step=100):
    with open(fname, 'wb') as f:
        f.write(np.asarray(a, dtype=dtype).tobytes(order='F'))
        f.write(np.array([time, step], dtype=dtype).tobytes())


def write_legacy(fname, fields, dtype=np.float64, time=1.5, step=100):
    with open(fname, 'wb') as f:
        for a in fields:
            f.write(np.asarray(a, dtype=dtype).tobytes(order='F'))
        f.write(np.array([time, step], dtype=dtype).tobytes())


def read_field(fname, shape, dtype=np.float64):
    raw = np.fromfile(fname, dtype=dtype)
    n = int(np.prod(shape))
    assert raw.size == n + 2, f'{fname}: unexpected size'
    return raw[:n].reshape(shape, order='F'), raw[n:]


def write_field_h5(fname, name, a, dtype=np.float64, time=1.5, step=100):
    """HDF5 checkpoint as written by CaNS: fields/<name> (C order, i.e. z, y, x), meta/time, meta/istep"""
    with h5py.File(fname, 'w') as f:
        f.create_dataset(f'fields/{name}', data=np.asarray(a, dtype=dtype).T)
        f.create_dataset('meta/time', data=np.array([time], dtype=dtype))
        f.create_dataset('meta/istep', data=np.array([step], dtype=np.int32))


def read_field_h5(fname, name):
    with h5py.File(fname, 'r') as f:
        d = f[f'fields/{name}']
        return d[()].T, np.array([f['meta/time'][0], f['meta/istep'][0]]), d.dtype, d.chunks, d.compression, \
            {k: f[f'grid/{k}'][()] for k in 'xyz'}


def write_grid_h5(fname, g, dtype=np.float64):
    """CaNS grid file in HDF5: rc, rf, drc, drf (n values each)"""
    xf, xc = g.with_ghosts(False)
    with h5py.File(fname, 'w') as f:
        f['rc'] = g.xc[1:].astype(dtype)
        f['rf'] = g.xf[1:].astype(dtype)
        f['drc'] = (xc[2:] - xc[1:-1]).astype(dtype)
        f['drf'] = np.diff(g.xf).astype(dtype)


# ----------------------------------------------------------------------------------------------------------
# namelist and runs
# ----------------------------------------------------------------------------------------------------------

def fmt(v):
    if isinstance(v, bool):
        return 'T' if v else 'F'
    if isinstance(v, str):
        return "'" + v + "'"
    if isinstance(v, (list, tuple, np.ndarray)):
        return ', '.join(fmt(x) for x in np.ravel(np.asarray(v, dtype=object), order='F'))
    return repr(float(v)) if isinstance(v, float) else str(v)


def namelist(groups):
    """groups: list of (name, dict); a key may hold an index range, e.g. 'cbcvel(0:1,1:3,1)'"""
    out = []
    for name, d in groups:
        out.append('&' + name)
        out += [f'{k} = {fmt(v)}' for k, v in d.items()]
        out.append('/')
    return '\n'.join(out) + '\n'


def run(workdir, nml, nproc=1, expect_fail=False):
    with open(os.path.join(workdir, 'interp.nml'), 'w') as f:
        f.write(nml)
    cmd = MPIRUN + ['-n', str(nproc), EXE]
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True)
    if expect_fail:
        assert r.returncode != 0, 'the run was expected to fail:\n' + r.stdout
    elif r.returncode != 0:
        raise RuntimeError(f'run failed ({" ".join(cmd)}):\n{r.stdout}\n{r.stderr}')
    return r.stdout


def bc_groups(cbc, bc):
    """&bc entries for u, v, w, p and scalars from per-variable (cbc, bc) of shape (2, 3)"""
    d = {}
    for name, (c, b) in cbc_bc_items(cbc, bc):
        if name in 'uvw':
            ivel = 'uvw'.index(name) + 1
            d[f'cbcvel(0:1,1:3,{ivel})'] = c
            d[f'bcvel(0:1,1:3,{ivel})'] = b
        elif name == 'p':
            d['cbcpre(0:1,1:3)'] = c
            d['bcpre(0:1,1:3)'] = b
        elif name.startswith('s_'):
            iscal = int(name[2:])
            d[f'cbcscal(0:1,1:3,{iscal})'] = c
            d[f'bcscal(0:1,1:3,{iscal})'] = b
    return d


def cbc_bc_items(cbc, bc):
    for name in cbc:
        yield name, (np.asarray(cbc[name]).T, np.asarray(bc[name], dtype=float).T)


def case_nml(gi, go, files_i, files_o, variables, cbc, bc, locs, dims=(1, 0, 0), extra=None, ops=None):
    """namelist of a case: gi, go: lists of Grid1d per direction; files_*: grid file names ('' = uniform);
    cbc, bc: per variable, arrays (3, 2) of [direction][side]; locs: locations of non-default variables;
    extra: other `&interp` entries; ops: `&operations` entries"""
    interp = {'dims(1:3)': list(dims), 'variables': variables, 'input_format': 'per-variable',
              'input_prefix': 'fld_i', 'output_prefix': 'fld_o'}
    interp.update(extra or {})
    groups = [('interp', interp),
              ('grid_in', {'ng(1:3)': [g.n for g in gi], 'l(1:3)': [g.l for g in gi], 'grid_files(1:3)': files_i}),
              ('grid_out', {'ng(1:3)': [g.n for g in go], 'l(1:3)': [g.l for g in go], 'grid_files(1:3)': files_o}),
              ('bc', bc_groups(cbc, bc))]
    for name in variables:
        if name in locs:
            groups.append(('variable', {'name': name, 'is_staggered(1:3)': [c == 's' for c in locs[name]],
                                        'cbc(0:1,1:3)': np.asarray(cbc[name]).T,
                                        'bc(0:1,1:3)': np.asarray(bc[name], dtype=float).T}))
    if ops:
        groups.append(('operations', ops))
    return namelist(groups)


# ----------------------------------------------------------------------------------------------------------
# independent reference: ghost points from the boundary conditions, then linear interpolation per axis
# ----------------------------------------------------------------------------------------------------------

def reference(a, gi, go, loc, cbc, bc, periodic):
    """interpolates a (ni1 x ni2 x ni3, x fastest) from grids gi to go; loc e.g. 'scc'; cbc, bc: (3, 2)"""
    ext = np.asarray(a, dtype=np.float64)
    coords = []
    for d in range(3):
        n = gi[d].n
        xf, xc = gi[d].with_ghosts(periodic[d])
        dxf = np.empty(n + 2)
        dxf[1:] = np.diff(xf)
        dxf[0] = dxf[n] if periodic[d] else dxf[1]
        dxc = np.diff(xc)
        stag = loc[d] == 's'
        e = np.moveaxis(ext, d, 0)
        new = np.zeros((n + 2,) + e.shape[1:])
        new[1:n+1] = e
        if periodic[d]:
            new[0], new[n+1] = e[n-1], e[0]
        else:
            (c0, c1), (b0, b1) = cbc[d], bc[d]
            if not stag:
                new[0] = 2*b0 - new[1] if c0 == 'D' else new[1] - dxc[0]*b0
                new[n+1] = 2*b1 - new[n] if c1 == 'D' else new[n] + dxc[n]*b1
            else:
                new[0] = b0 if c0 == 'D' else new[1] - dxf[0]*b0
                new[n] = b1 if c1 == 'D' else new[n-1] + dxf[n]*b1
                new[n+1] = new[n]
        ext = np.moveaxis(new, 0, d)
        coords.append(xf if stag else xc)
    out = ext
    for d in range(3):
        n = gi[d].n
        xi = coords[d]
        xfo, xco = go[d].with_ghosts(periodic[d])
        xo = (xfo[1:go[d].n+1] if loc[d] == 's' else xco[1:go[d].n+1])
        mmax = n - 1 if loc[d] == 's' else n
        m = np.clip(np.searchsorted(xi[:mmax+1], xo, side='right') - 1, 0, mmax)
        w = (xo - xi[m])/(xi[m+1] - xi[m])
        e = np.moveaxis(out, d, 0)
        shape = (-1,) + (1,)*(e.ndim - 1)
        out = np.moveaxis(e[m]*(1 - w.reshape(shape)) + e[m+1]*w.reshape(shape), 0, d)
    return out


def weighted_mean(a, g, loc, periodic):
    """volume-weighted mean, with the control volumes of the program (centres: the cells; faces: the spacings
    between the centres, and the half cell at a bounded upper boundary)"""
    w = []
    for d in range(3):
        n = g[d].n
        xf, xc = g[d].with_ghosts(periodic[d])
        if loc[d] == 's':
            wd = xc[2:n+2] - xc[1:n+1]
            if not periodic[d]:
                wd[-1] = xf[n] - xc[n]
        else:
            wd = np.diff(xf[:n+1])
        w.append(wd)
    vol = w[0][:, None, None]*w[1][None, :, None]*w[2][None, None, :]
    return (vol*a).sum()/vol.sum()


def coordinates(g, loc, periodic):
    """3D coordinate arrays of a variable's location (the interior points 1..n along each direction)"""
    x = []
    for d in range(3):
        xf, xc = g[d].with_ghosts(periodic[d])
        x.append(xf[1:g[d].n+1] if loc[d] == 's' else xc[1:g[d].n+1])
    return np.meshgrid(*x, indexing='ij')


def random_field(g, loc, periodic, cbc, bc, seed):
    """smooth plus random field, with the staggered boundary values consistent with Dirichlet conditions"""
    rng = np.random.default_rng(seed)
    x, y, z = coordinates(g, loc, periodic)
    a = np.sin(2*np.pi*x/g[0].l)*np.cos(2*np.pi*y/g[1].l) + z/g[2].l + 0.1*rng.standard_normal(x.shape)
    for d in range(3):
        if loc[d] == 's' and not periodic[d] and cbc[d][1] == 'D':
            idx = [slice(None)]*3
            idx[d] = -1
            a[tuple(idx)] = bc[d][1]
    return a


# ----------------------------------------------------------------------------------------------------------
# cases
# ----------------------------------------------------------------------------------------------------------

CHANNEL_CBC = {'u': [['P', 'P'], ['P', 'P'], ['D', 'D']], 'v': [['P', 'P'], ['P', 'P'], ['D', 'D']],
               'w': [['P', 'P'], ['P', 'P'], ['D', 'D']], 'p': [['P', 'P'], ['P', 'P'], ['N', 'N']],
               's_001': [['P', 'P'], ['P', 'P'], ['D', 'N']], 's_002': [['P', 'P'], ['P', 'P'], ['N', 'D']]}
CHANNEL_BC = {'u': [[0, 0], [0, 0], [-1, -1]], 'v': [[0, 0], [0, 0], [0, 0]], 'w': [[0, 0], [0, 0], [0, 0]],
              'p': [[0, 0], [0, 0], [0, 0]], 's_001': [[0, 0], [0, 0], [1, 0.5]],
              's_002': [[0, 0], [0, 0], [-0.3, 2]]}


def walls_bc(names, kinds):
    """all directions bounded; kinds: per direction, (lower, upper) condition type for all variables"""
    cbc = {n: [list(k) for k in kinds] for n in names}
    bc = {n: [[0.3, -0.2], [0.1, 0.4], [-0.5, 0.25]] for n in names}
    return cbc, bc


def grids(kind):
    """input and output grids (per direction) of the test cases, and their grid-file flags"""
    if kind == 'uniform':
        return [uniform(12, 1.5), uniform(10, 1.), uniform(8, 1.)], \
               [uniform(17, 1.5), uniform(23, 1.), uniform(13, 1.)], [False]*3
    if kind == 'channel':
        return [uniform(12, 2*np.pi), uniform(10, np.pi), stretched(14, 2., 1.8)], \
               [uniform(18, 2*np.pi), uniform(13, np.pi), stretched(21, 2., 2.2)], [False, False, True]
    if kind == 'eigen':
        return [stretched(12, 1.5, periodic=True, seed=1), stretched(10, 1., 1.4), stretched(9, 1., 1.7)], \
               [stretched(19, 1.5, periodic=True, seed=2), stretched(15, 1., 2.0), stretched(14, 1., 1.2)], \
               [True]*3
    raise ValueError(kind)


def setup(workdir, gi, go, files, dtype_grid=np.float64):
    """writes the grid files of the non-uniform directions; returns their names"""
    fi, fo = [], []
    for d, c in enumerate('xyz'):
        fi.append(f'grid_i_{c}.bin' if files[d] else '')
        fo.append(f'grid_o_{c}.bin' if files[d] else '')
        if files[d]:
            write_grid(os.path.join(workdir, fi[-1]), gi[d], dtype_grid)
            write_grid(os.path.join(workdir, fo[-1]), go[d], dtype_grid)
    return fi, fo


def shape(g):
    return tuple(gg.n for gg in g)


# ----------------------------------------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------------------------------------

def test_identity(workdir):
    """output grid = input grid: the input is reproduced exactly, in single and double precision"""
    for kind in ['uniform', 'channel', 'eigen']:
        gi, _, files = grids(kind)
        if kind == 'eigen':
            # no Neumann condition at an upper boundary face of a staggered variable, which sets its value
            cbc, bc = walls_bc('uvwp', [['P', 'P'], ['N', 'D'], ['N', 'D']])
        else:
            cbc, bc = CHANNEL_CBC, CHANNEL_BC
        periodic = [cbc['p'][d][0] == 'P' for d in range(3)]
        for dtype in [np.float64, np.float32]:
            fi, fo = setup(workdir, gi, gi, files, dtype)
            data = {}
            for iv, name in enumerate('uvwp'):
                data[name] = random_field(gi, DEFAULT_LOC[name], periodic, cbc[name], bc[name], iv).astype(dtype)
                write_field(os.path.join(workdir, f'fld_i_{name}.bin'), data[name], dtype)
            for nproc, dims in [(1, (1, 1, 1)), (4, (1, 2, 2)), (3, (3, 1, 1))]:
                nml = case_nml(gi, gi, fi, fi, list('uvwp'), cbc, bc, {}, dims=dims)
                run(workdir, nml, nproc)
                for name in 'uvwp':
                    out, info = read_field(os.path.join(workdir, f'fld_o_{name}.bin'), shape(gi), dtype)
                    assert np.array_equal(out, data[name]), f'{kind} {name} {dtype.__name__} {dims}: not identical'
                    assert np.array_equal(info, [1.5, 100]), f'{kind}: time and step not kept'
        print(f'  identity, {kind} grid: exact (double and single precision, 3 process grids)')


def linear_case(kind):
    """grids, conditions and linear fields consistent with them, for all variable locations"""
    gi, go, files = grids(kind)
    locs = {'u': 'scc', 'v': 'csc', 'w': 'ccs', 'p': 'ccc', 's_001': 'ccc', 'omz': 'ssc', 'q': 'sss'}
    if kind == 'eigen':
        # bounded in all directions: Neumann conditions with the gradient of the field
        grad = np.array([0.7, -1.3, 2.1])
        periodic = [False]*3
        cbc = {n: [['N', 'N']]*3 for n in locs}
        bc = {n: [[grad[d], grad[d]] for d in range(3)] for n in locs}
        cases = [(cbc, bc, lambda x, y, z: 0.4 + grad[0]*x + grad[1]*y + grad[2]*z)]
        # Dirichlet conditions along z (field linear in z only), on both sides or mixed
        f = lambda x, y, z: 0.4 + 2.1*z + 0*x + 0*y
        for kz, bz in [(['D', 'D'], [0.4, 2.5]), (['N', 'D'], [2.1, 2.5]), (['D', 'N'], [0.4, 2.1])]:
            cases.append(({n: [['N', 'N'], ['N', 'N'], kz] for n in locs},
                          {n: [[0., 0.], [0., 0.], bz] for n in locs}, f))
    else:
        # channel: linear in z, any function of the periodic directions
        periodic = [True, True, False]
        cbc = {n: [['P', 'P'], ['P', 'P'], ['D', 'N']] for n in locs}
        bc = {n: [[0, 0], [0, 0], [0.4, -0.6]] for n in locs}
        f = lambda x, y, z: 0.4 - 0.6*z + 0*x + 0*y
        cases = [(cbc, bc, f)]
        cbc2 = {n: [['P', 'P'], ['P', 'P'], ['N', 'N']] for n in locs}
        bc2 = {n: [[0, 0], [0, 0], [1.7, 1.7]] for n in locs}
        cases.append((cbc2, bc2, lambda x, y, z: -0.2 + 1.7*z + 0*x + 0*y))
    return gi, go, files, locs, periodic, cases


def test_linear(workdir):
    """linear fields consistent with the boundary conditions are interpolated exactly, at every location"""
    for kind in ['channel', 'eigen']:
        gi, go, files, locs, periodic, cases = linear_case(kind)
        fi, fo = setup(workdir, gi, go, files)
        names = list(locs)
        for cbc, bc, f in cases:
            for name in names:
                write_field(os.path.join(workdir, f'fld_i_{name}.bin'), f(*coordinates(gi, locs[name], periodic)))
            other = {n: locs[n] for n in ['omz', 'q']}
            nml = case_nml(gi, go, fi, fo, names, cbc, bc, other, dims=(2, 1, 2))
            run(workdir, nml, 4)
            err = 0.
            for name in names:
                out, _ = read_field(os.path.join(workdir, f'fld_o_{name}.bin'), shape(go))
                exact = f(*coordinates(go, locs[name], periodic))
                err = max(err, np.abs(out - exact).max())
            assert err < 1e-13, f'{kind}: linear fields not exact, max error {err}'
        print(f'  linear fields, {kind} grids: exact to round-off at locations '
              f'{", ".join(sorted(set(locs.values())))} (max error {err:.1e})')


def test_layouts(workdir):
    """identical output (bitwise) for several process grids, and agreement with the NumPy reference"""
    for kind in ['channel', 'eigen']:
        gi, go, files = grids(kind)
        if kind == 'eigen':
            cbc, bc = walls_bc(['u', 'v', 'w', 'p', 's_001'], [['P', 'P'], ['D', 'N'], ['N', 'D']])
        else:
            cbc, bc = CHANNEL_CBC, CHANNEL_BC
        periodic = [cbc['p'][d][0] == 'P' for d in range(3)]
        names = ['u', 'v', 'w', 'p', 's_001']
        fi, fo = setup(workdir, gi, go, files)
        ref = {}
        for iv, name in enumerate(names):
            loc = DEFAULT_LOC.get(name, 'ccc')
            a = random_field(gi, loc, periodic, cbc[name], bc[name], 10 + iv)
            write_field(os.path.join(workdir, f'fld_i_{name}.bin'), a)
            ref[name] = reference(a, gi, go, loc, cbc[name], bc[name], periodic)
        outs = {}
        layouts = [(1, (1, 1, 1)), (2, (2, 1, 1)), (4, (2, 2, 1)), (6, (1, 2, 3)), (5, (1, 1, 5)), (8, (2, 2, 2)),
                   (7, (7, 1, 1)), (6, (0, 0, 0))]
        for nproc, dims in layouts:
            nml = case_nml(gi, go, fi, fo, names, cbc, bc, {}, dims=dims)
            run(workdir, nml, nproc)
            for name in names:
                out, _ = read_field(os.path.join(workdir, f'fld_o_{name}.bin'), shape(go))
                if name not in outs:
                    outs[name] = out
                    err = np.abs(out - ref[name]).max()
                    assert err < 1e-12, f'{kind} {name}: differs from the reference by {err}'
                assert out.tobytes() == outs[name].tobytes(), f'{kind} {name}: output depends on dims = {dims}'
        print(f'  process grids, {kind} grids: bitwise identical for dims = '
              f'{", ".join(str(d) for _, d in layouts)}; equal to the NumPy reference to round-off')


def test_precision(workdir):
    """double in, single out = double result rounded; single in and out; grid files in either precision"""
    gi, go, files = grids('eigen')
    cbc, bc = walls_bc('uvwp', [['P', 'P'], ['D', 'N'], ['N', 'D']])
    periodic = [True, False, False]
    data = {n: random_field(gi, DEFAULT_LOC[n], periodic, cbc[n], bc[n], 20 + i) for i, n in enumerate('uvwp')}
    results = {}
    for dtype_grid in [np.float64, np.float32]:
        fi, fo = setup(workdir, gi, go, files, dtype_grid)
        for dtype_in, prec_out in [(np.float64, 'input'), (np.float64, 'single'), (np.float32, 'input'),
                                   (np.float32, 'double')]:
            for n in 'uvwp':
                write_field(os.path.join(workdir, f'fld_i_{n}.bin'), data[n], dtype_in)
            nml = case_nml(gi, go, fi, fo, list('uvwp'), cbc, bc, {}, dims=(1, 2, 2),
                           extra={'output_precision': prec_out})
            run(workdir, nml, 4)
            dtype_out = {'input': dtype_in, 'single': np.float32, 'double': np.float64}[prec_out]
            for n in 'uvwp':
                out, info = read_field(os.path.join(workdir, f'fld_o_{n}.bin'), shape(go), dtype_out)
                assert np.array_equal(info, [1.5, 100])
                results[(dtype_grid, dtype_in, prec_out, n)] = out
        for n in 'uvwp':
            dd = results[(dtype_grid, np.float64, 'input', n)]
            assert np.array_equal(results[(dtype_grid, np.float64, 'single', n)], dd.astype(np.float32)), \
                f'{n}: single-precision output is not the double-precision result rounded'
            sd = results[(dtype_grid, np.float32, 'double', n)]
            assert np.array_equal(results[(dtype_grid, np.float32, 'input', n)], sd.astype(np.float32)), \
                f'{n}: single in, single out differs from single in, double out rounded'
    for n in 'uvwp':
        err = np.abs(results[(np.float32, np.float64, 'input', n)] - results[(np.float64, np.float64, 'input', n)]).max()
        assert 0 < err < 1e-5, f'{n}: single- and double-precision grid files give {err}'
    print('  precision: double in, single out = double result rounded; single in and out = single in, double out '
          'rounded; single- and double-precision grid files agree to single-precision round-off')


def test_formats(workdir):
    """legacy single-file and per-variable input give identical output; scalars; a subset of the variables"""
    gi, go, files = grids('channel')
    cbc, bc = CHANNEL_CBC, CHANNEL_BC
    periodic = [True, True, False]
    names = ['u', 'v', 'w', 'p', 's_001', 's_002']
    fi, fo = setup(workdir, gi, go, files)
    for dtype in [np.float64, np.float32]:
        data = [random_field(gi, DEFAULT_LOC.get(n, 'ccc'), periodic, cbc[n], bc[n], 30 + i).astype(dtype)
                for i, n in enumerate(names)]
        for n, a in zip(names, data):
            write_field(os.path.join(workdir, f'fld_i_{n}.bin'), a, dtype)
        write_legacy(os.path.join(workdir, 'fld_legacy.bin'), data, dtype)
        nml = case_nml(gi, go, fi, fo, names, cbc, bc, {}, dims=(1, 2, 2))
        run(workdir, nml, 4)
        per_var = {n: read_field(os.path.join(workdir, f'fld_o_{n}.bin'), shape(go), dtype) for n in names}
        nml = case_nml(gi, go, fi, fo, names, cbc, bc, {}, dims=(2, 1, 2),
                       extra={'input_format': 'legacy', 'input_file': 'fld_legacy.bin'})
        run(workdir, nml, 4)
        for n in names:
            out, info = read_field(os.path.join(workdir, f'fld_o_{n}.bin'), shape(go), dtype)
            assert np.array_equal(out, per_var[n][0]) and np.array_equal(info, per_var[n][1]), \
                f'{n}: legacy and per-variable input differ'
        # subset: only p, from the legacy file (with two scalars); the other outputs are not written
        for n in names:
            os.remove(os.path.join(workdir, f'fld_o_{n}.bin'))
        nml = case_nml(gi, go, fi, fo, ['p'], cbc, bc, {}, dims=(1, 3, 1),
                       extra={'input_format': 'legacy', 'input_file': 'fld_legacy.bin', 'is_reset_time': True})
        run(workdir, nml, 3)
        out, info = read_field(os.path.join(workdir, 'fld_o_p.bin'), shape(go), dtype)
        assert np.array_equal(out, per_var['p'][0]) and np.array_equal(info, [0, 0])
        assert sorted(f for f in os.listdir(workdir) if f.startswith('fld_o')) == ['fld_o_p.bin']
    print('  formats: legacy (with 2 scalars) and per-variable input identical, in single and double precision; '
          'subset (p only) with the time reset')


def test_defaults(workdir):
    """default variables (u, v, w, p), default locations, and no `&variable` group needed"""
    gi, go, files = grids('channel')
    fi, fo = setup(workdir, gi, go, files)
    periodic = [True, True, False]
    for i, n in enumerate('uvwp'):
        a = random_field(gi, DEFAULT_LOC[n], periodic, CHANNEL_CBC[n], CHANNEL_BC[n], 40 + i)
        write_field(os.path.join(workdir, f'fld_i_{n}.bin'), a)
    groups = [('interp', {'input_prefix': 'fld_i', 'output_prefix': 'fld_o'}),
              ('grid_in', {'ng(1:3)': shape(gi), 'l(1:3)': [g.l for g in gi], 'grid_files(3)': fi[2]}),
              ('grid_out', {'ng(1:3)': shape(go), 'l(1:2)': [g.l for g in go[:2]], 'grid_files(3)': fo[2]}),
              ('bc', bc_groups({n: CHANNEL_CBC[n] for n in 'uvwp'}, {n: CHANNEL_BC[n] for n in 'uvwp'}))]
    out = run(workdir, namelist(groups), 2)
    for n in 'uvwp':
        assert f'Variable {n} ({DEFAULT_LOC[n]})' in out, out
        a, _ = read_field(os.path.join(workdir, f'fld_i_{n}.bin'), shape(gi))
        o, _ = read_field(os.path.join(workdir, f'fld_o_{n}.bin'), shape(go))
        ref = reference(a, gi, go, DEFAULT_LOC[n], np.asarray(CHANNEL_CBC[n]), np.asarray(CHANNEL_BC[n], float),
                        periodic)
        assert np.abs(o - ref).max() < 1e-12
    print('  defaults: u, v, w, p at their default locations, without `variables` or `&variable` groups')


def test_input_errors(workdir):
    """input errors stop the run with a message"""
    gi, go, files = grids('channel')
    fi, fo = setup(workdir, gi, go, files)
    periodic = [True, True, False]
    a = random_field(gi, 'ccc', periodic, CHANNEL_CBC['p'], CHANNEL_BC['p'], 0)
    write_field(os.path.join(workdir, 'fld_i_p.bin'), a)
    write_field(os.path.join(workdir, 'fld_i_q.bin'), a[:-1])
    base = dict(gi=gi, go=go, files_i=fi, files_o=fo, cbc=CHANNEL_CBC, bc=CHANNEL_BC)
    checks = [
        (case_nml(variables=['q'], locs={'q': 'ccc'}, **{**base, 'cbc': {'q': CHANNEL_CBC['p']},
                                                          'bc': {'q': CHANNEL_BC['p']}}), 'size of checkpoint'),
        (case_nml(variables=['omz'], locs={}, **base), 'must be prescribed in a `&variable` group'),
        (case_nml(variables=['p'], locs={}, **{**base, 'cbc': {}, 'bc': {}}), 'boundary conditions'),
        (case_nml(variables=['p'], locs={}, **base, extra={'output_prefix': 'fld_i'}), 'would overwrite its input'),
        (case_nml(variables=['p'], locs={}, **base, dims=(1, 3, 0)), 'inconsistent with the number of MPI tasks'),
        (case_nml(variables=['p'], locs={}, **{**base, 'go': [uniform(18, 6.), go[1], go[2]]}),
         'domain lengths along x differ'),
        (case_nml(variables=['p'], locs={}, **{**base, 'files_i': ['', '', 'grid_o_z.bin']}), 'grid file'),
        (case_nml(variables=['p', 'p'], locs={}, **base), 'listed twice'),
    ]
    for nml, msg in checks:
        out = run(workdir, nml, 2, expect_fail=True)
        assert msg in out, f'expected "{msg}" in:\n{out}'
    print(f'  input errors: {len(checks)} cases stop with the expected message')


def hdf5_enabled(workdir):
    """whether the program is built with HDF5 (otherwise, a request for HDF5 output stops with a message)"""
    gi = [uniform(4, 1.), uniform(4, 1.), uniform(4, 1.)]
    cbc, bc = {'p': [['P', 'P']]*3}, {'p': [[0, 0]]*3}
    write_field(os.path.join(workdir, 'fld_i_p.bin'), np.zeros(shape(gi)))
    nml = case_nml(gi, gi, ['']*3, ['']*3, ['p'], cbc, bc, {}, extra={'output_format': 'hdf5'})
    with open(os.path.join(workdir, 'interp.nml'), 'w') as f:
        f.write(nml)
    r = subprocess.run(MPIRUN + ['-n', '1', EXE], cwd=workdir, capture_output=True, text=True)
    if 'HDF5 support is not compiled in' in r.stdout:
        return False
    assert r.returncode == 0, r.stdout + r.stderr
    return True


def test_hdf5(workdir):
    """HDF5 checkpoints and grid files, in and out, against the MPI-IO files; compression; process grids"""
    gi, go, files = grids('eigen')
    names = ['u', 'v', 'w', 'p', 's_001']
    cbc, bc = walls_bc(names, [['P', 'P'], ['D', 'N'], ['N', 'D']])
    periodic = [True, False, False]
    locs = {n: DEFAULT_LOC.get(n, 'ccc') for n in names}
    if not hdf5_enabled(workdir):
        fi, fo = setup(workdir, gi, go, files)
        for n in names:
            write_field(os.path.join(workdir, f'fld_i_{n}.bin'), random_field(gi, locs[n], periodic, cbc[n], bc[n], 0))
        checks = [{'input_format': 'hdf5'}, {'output_format': 'hdf5'}, None]
        for extra in checks:
            fh = fi if extra else fi[:2] + ['grid_i_z.h5']
            out = run(workdir, case_nml(gi, go, fh, fo, names, cbc, bc, {}, extra=extra), 2, expect_fail=True)
            assert 'HDF5 support is not compiled in' in out, out
        print(f'  HDF5 not compiled in: {len(checks)} requests (input, output, grid file) stop with the message')
        return
    assert h5py is not None, 'the HDF5 tests need h5py'
    for dtype in [np.float64, np.float32]:
        fi, fo = setup(workdir, gi, go, files, dtype)
        fi_h5 = [f.replace('.bin', '.h5') for f in fi]
        fo_h5 = [f.replace('.bin', '.h5') for f in fo]
        for d in range(3):
            write_grid_h5(os.path.join(workdir, fi_h5[d]), gi[d], dtype)
            write_grid_h5(os.path.join(workdir, fo_h5[d]), go[d], dtype)
        for i, n in enumerate(names):
            a = random_field(gi, locs[n], periodic, cbc[n], bc[n], 50 + i).astype(dtype)
            write_field(os.path.join(workdir, f'fld_i_{n}.bin'), a, dtype)
            write_field_h5(os.path.join(workdir, f'fld_i_{n}.h5'), n, a, dtype)
        out = {}
        for fin, fout, grid_ext in [('per-variable', 'per-variable', 'bin'), ('per-variable', 'hdf5', 'bin'),
                                    ('hdf5', 'per-variable', 'bin'), ('hdf5', 'hdf5', 'bin'), ('hdf5', 'hdf5', 'h5')]:
            gfi, gfo = (fi, fo) if grid_ext == 'bin' else (fi_h5, fo_h5)
            nml = case_nml(gi, go, gfi, gfo, names, cbc, bc, {}, dims=(2, 1, 2),
                           extra={'input_format': fin, 'output_format': fout})
            run(workdir, nml, 4)
            for n in names:
                if fout == 'hdf5':
                    a, info, dt, chunks, comp, xyz = read_field_h5(os.path.join(workdir, f'fld_o_{n}.h5'), n)
                    assert dt == dtype and chunks is None and comp is None, (dt, chunks, comp)
                    for d, k in enumerate('xyz'):
                        xf, xc = go[d].with_ghosts(periodic[d])
                        x = xf if locs[n][d] == 's' else xc
                        assert xyz[k].dtype == dtype and np.abs(xyz[k] - x.astype(dtype)).max() < 1e-6, (n, k)
                else:
                    a, info = read_field(os.path.join(workdir, f'fld_o_{n}.bin'), shape(go), dtype)
                assert np.array_equal(info, [1.5, 100]), info
                out.setdefault(n, a)
                assert a.tobytes() == out[n].tobytes(), f'{n}: {fin} in, {fout} out, {grid_ext} grids differ'
    print('  HDF5: identical values for any combination of .bin and .h5 input and output, and of .bin and .h5 '
          'grid files, in single and double precision; time, step and coordinates written as CaNS')
    # compression on and off, several process grids (with wrap-around reads along x), double precision
    fi, fo = setup(workdir, gi, go, files)
    for i, n in enumerate(names):
        write_field_h5(os.path.join(workdir, f'fld_i_{n}.h5'), n,
                       random_field(gi, locs[n], periodic, cbc[n], bc[n], 60 + i))
    ref = {}
    for nproc, dims, comp in [(1, (1, 1, 1), False), (4, (2, 1, 2), True), (6, (3, 2, 1), True), (5, (1, 1, 5), False),
                              (8, (2, 2, 2), True)]:
        nml = case_nml(gi, go, fi, fo, names, cbc, bc, {}, dims=dims,
                       extra={'input_format': 'hdf5', 'output_format': 'hdf5', 'is_compress_output': comp})
        run(workdir, nml, nproc)
        for n in names:
            a, _, _, chunks, compression, _ = read_field_h5(os.path.join(workdir, f'fld_o_{n}.h5'), n)
            assert (compression == 'gzip' and chunks is not None) if comp else (compression is None), (comp, chunks)
            ref.setdefault(n, a)
            assert a.tobytes() == ref[n].tobytes(), f'{n}: HDF5 output depends on dims = {dims} or compression'
    print('  HDF5: compression on and off give the same values, bitwise identical for 5 process grids')
    # input checks: a missing dataset, and a field with the wrong dimensions
    write_field_h5(os.path.join(workdir, 'fld_i_q.h5'), 'p', np.zeros(shape(gi)))
    write_field_h5(os.path.join(workdir, 'fld_i_r.h5'), 'r', np.zeros(shape(go)))
    for name, msg in [('q', 'has no dataset `fields/q`'), ('r', 'is not of dimensions')]:
        nml = case_nml(gi, go, fi, fo, [name], {name: cbc['p']}, {name: bc['p']}, {name: 'ccc'},
                       extra={'input_format': 'hdf5'})
        out = run(workdir, nml, 2, expect_fail=True)
        assert msg in out, out
    print('  HDF5: a missing dataset and wrong dimensions stop with a message')


def test_operations(workdir):
    """mean removal, constant offsets and change of reference frame"""
    for kind in ['channel', 'eigen']:
        gi, go, files = grids(kind)
        names = ['u', 'v', 'w', 'p']
        if kind == 'eigen':
            cbc, bc = walls_bc(names, [['P', 'P'], ['D', 'N'], ['N', 'D']])
        else:
            cbc, bc = CHANNEL_CBC, CHANNEL_BC
        periodic = [cbc['p'][d][0] == 'P' for d in range(3)]
        fi, fo = setup(workdir, gi, go, files)
        data, ref = {}, {}
        for i, n in enumerate(names):
            data[n] = random_field(gi, DEFAULT_LOC[n], periodic, cbc[n], bc[n], 70 + i)
            ref[n] = reference(data[n], gi, go, DEFAULT_LOC[n], np.asarray(cbc[n]), np.asarray(bc[n], float), periodic)

        def interp(ops, big=0., dtype=np.float64, prec='double', dims=(2, 1, 2), nproc=4):
            # p is shifted by `big`, with its Dirichlet values (the boundary conditions refer to the input)
            bcs = {**bc, 'p': [[b + big*(c == 'D') for b, c in zip(bc['p'][d], cbc['p'][d])] for d in range(3)]}
            for n in names:
                write_field(os.path.join(workdir, f'fld_i_{n}.bin'), data[n] + (big if n == 'p' else 0.))
            nml = case_nml(gi, go, fi, fo, names, cbc, bcs, {}, dims=dims, extra={'output_precision': prec}, ops=ops)
            log = run(workdir, nml, nproc)
            return {n: read_field(os.path.join(workdir, f'fld_o_{n}.bin'), shape(go), dtype)[0] for n in names}, log
        # mean removal of all variables, with p shifted by 1e6 (as the drifting pressure of a long run)
        out, log = interp({'is_remove_mean(1:4)': [True]*4}, big=1.e6)
        for n in names:
            m = weighted_mean(out[n], go, DEFAULT_LOC[n], periodic)
            assert abs(m) < 1e-14, f'{kind} {n}: output mean {m}'
            exact = ref[n] - weighted_mean(ref[n], go, DEFAULT_LOC[n], periodic)
            err = np.abs(out[n] - exact).max()
            assert err < (1e-9 if n == 'p' else 1e-13), f'{kind} {n}: differs from interp - mean by {err}'
        assert 'p: mean of the input subtracted:' in log and 'p: residual mean of the output subtracted:' in log
        # single-precision output: the field shifted by 1e6 gives the same field as without the shift
        sp_big = interp({'is_remove_mean(4)': True}, big=1.e6, dtype=np.float32, prec='single')[0]['p']
        sp_ref = interp({'is_remove_mean(4)': True}, dtype=np.float32, prec='single')[0]['p']
        ulp = np.spacing(np.abs(sp_ref).max())
        assert np.abs(sp_big - sp_ref).max() <= ulp, f'{kind}: single-precision output depends on the mean'
        nsame = np.count_nonzero(sp_big == sp_ref)/sp_big.size
        sp_nomean = interp({}, big=1.e6, dtype=np.float32, prec='single')[0]['p'].astype(np.float64) - 1.e6
        err_nomean = np.abs(sp_nomean - ref['p']).max()
        assert err_nomean > 1e-2
        # offsets, also of a variable with Dirichlet conditions (u at the walls)
        base = interp({})[0]
        out = interp({'offset(1:4)': [3.7, 0., -1.25, 949381.]})[0]
        for n, c in zip(names, [3.7, 0., -1.25, 949381.]):
            err = np.abs(out[n] - (base[n] - c)).max()
            assert err < 1e-13*max(1., abs(c)), f'{kind} {n}: offset {c} gives an error {err}'
        # change of reference frame, combined with mean removal of p
        uf = [0.3, -0.2, 0.1]
        out, log = interp({'vel_frame(1:3)': uf, 'is_remove_mean(4)': True}, dims=(1, 3, 2), nproc=6)
        for n, U in zip('uvw', uf):
            assert np.array_equal(out[n], base[n] + U), f'{kind} {n}: frame velocity not added exactly'
            assert f'{n}: frame velocity added:' in log
        assert abs(weighted_mean(out['p'], go, 'ccc', periodic)) < 1e-14
        print(f'  operations, {kind} grids: zero output mean (all locations, 1e6 shift of p), same single-precision '
              f'p with and without the shift ({100*nsame:.1f}% bitwise, the rest within 1 ulp of the field; without mean removal: '
              f'error {err_nomean:.2f}); offsets; frame velocity added exactly')
    # input checks
    for ops, msg in [({'is_remove_mean(4)': True, 'offset(4)': 1.}, 'are exclusive'),
                     ({'is_remove_mean(5)': True}, 'beyond the number of `variables`')]:
        nml = case_nml(gi, go, fi, fo, names, cbc, bc, {}, ops=ops)
        out = run(workdir, nml, 2, expect_fail=True)
        assert msg in out, out
    print('  operations: conflicting or out-of-range options stop with a message')


TESTS = {'identity': test_identity, 'linear': test_linear, 'layouts': test_layouts, 'precision': test_precision,
         'formats': test_formats, 'defaults': test_defaults, 'hdf5': test_hdf5, 'operations': test_operations,
         'errors': test_input_errors}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--exe', default=os.path.join(ROOT, 'exe', 'interpit'))
    parser.add_argument('--workdir', default=None, help='directory for the test files (default: a temporary one)')
    parser.add_argument('--keep', action='store_true', help='keep the test files')
    parser.add_argument('tests', nargs='*', default=list(TESTS))
    args = parser.parse_args()
    EXE = os.path.abspath(args.exe)
    base = tempfile.mkdtemp(prefix='interpit-tests-', dir=args.workdir)
    nfail = 0
    try:
        for name in args.tests:
            wd = os.path.join(base, name)
            os.makedirs(wd)
            print(f'{name}:', flush=True)
            try:
                TESTS[name](wd)
            except (AssertionError, RuntimeError) as e:
                nfail += 1
                print(f'  FAILED: {e}')
    finally:
        if args.keep:
            print(f'test files kept in {base}')
        else:
            shutil.rmtree(base)
    print('all tests passed' if nfail == 0 else f'{nfail} test(s) failed')
    sys.exit(1 if nfail else 0)
