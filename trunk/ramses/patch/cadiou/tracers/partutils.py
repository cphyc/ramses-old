from scipy.io import FortranFile as FF
from glob import glob
import os
import numpy as np
from tqdm import tqdm

from numba import jit

def part2dens(particles):
    # Find levelmin and max
    lvlmin = None
    lvlmax = None
    lvl = 0

    P = [particles*2**lvl]
    Pf = P.astype(int)
    barr = [Pf == P[-1]]

    while not np.all(barr[-1]):
        _any = np.any(barr[-1])
        if lvlmin is None and _any:
            lvlmin = lvl

        lvl += 1
        P.append(particles*2**lvl)
        Pf = P[-1].astype(int)
        barr.append(Pf == P[-1])

    lvlmax = lvl
    yield (lvlmin, lvlmax)

    dens = np.zeros((2**lvlmax, 2**lvlmax))
    for lvl in range(lvlmin, lvlmax+1):
        pos = P[lvl] * 2 * (lvlmax - lvl + 1)
        mask = barr[lvl]



#################################
# Read helpers
#################################
def read_one_cpu(output):
    f = FF(output)
    ncpu = f.read_ints()
    ndim = f.read_ints()
    npart = f.read_ints()
    localseed = f.read_ints()
    nstart_tot = f.read_ints()
    mstar_tot = f.read_ints()
    mstart_lost = f.read_ints()
    nsink = f.read_ints()

    x = np.zeros((ndim, npart), dtype=float)
    v = np.zeros((ndim, npart), dtype=float)
    m = np.zeros((npart), dtype=float)
    ind = np.zeros((npart), dtype=int)
    for dim in range(ndim):
        x[dim] = f.read_reals()
    for dim in range(ndim):
        v[dim] = f.read_reals()

    m = f.read_reals()
    ind = f.read_ints()
    lvl = f.read_ints()
    try:
        tp = f.read_reals()
        Z = f.read_reals()
    except TypeError:
        tp = np.zeros((npart))
        Z = np.zeros((npart))

    return ind, ndim, npart, x, v, m, lvl, tp, Z


def read_output(path):
    paths = glob(os.path.join(path, 'part_*'))

    _pos = [None]*len(paths)
    _vel = [None]*len(paths)
    _mass = [None]*len(paths)
    _lvl = [None]*len(paths)
    _ind = [None]*len(paths)
    _npart = [None]*len(paths)
    _Z = [None]*len(paths)
    _tp = [None]*len(paths)
    _cpus = [None]*len(paths)
    npart = 0
    for i, output in enumerate(paths):
        (_ind[i], ndim, _npart[i], _pos[i], _vel[i], _mass[i],
         _lvl[i], _tp[i], _Z[i]) = read_one_cpu(output)

        cpu = int(output.split('.out')[-1])
        _cpus[i] = cpu
        npart += _npart[i]

    # # Prevent bug when there is no particle
    # if i == 0:
    #     ndim = 1
    pos = np.zeros((ndim, npart))
    vel = np.zeros((ndim, npart))
    mass = np.zeros((npart))
    ind = np.zeros((npart), dtype=int)
    lvl = np.zeros((npart), dtype=int)
    cpus = np.zeros((npart), dtype=int)

    n = 0
    for i in range(len(paths)):
        ids = _ind[i] - 1
        pos[:, ids] = _pos[i]
        vel[:, ids] = _vel[i]
        ind[n:n+_npart[i]] = ids

        mass[ids] = _mass[i]
        isTracer = ids[_mass[i] == 0]
        mass[isTracer] = _tp[i]

        lvl[ids] = _lvl[i]
        cpus[ids] = _cpus[i]

        n += _npart[i]

    return ind, pos, vel, mass, lvl, cpus
