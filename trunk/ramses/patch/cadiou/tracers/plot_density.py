# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from glob import glob
from tqdm import tqdm
import os
from os.path import join as pjoin
from datetime import datetime

import sys
sys.path.append('/home/cadiou/codes/pymses')
import pymses
from pymses.filters import CellsToPoints
from pymses.analysis.slicing import SliceMap
from pymses.analysis.operator import ScalarOperator
from pymses.analysis import Camera

cam = Camera(line_of_sight_axis='z', up_vector='y')

mpl.rcParams['figure.figsize'] = (16,9)
mpl.rcParams['figure.dpi'] = 120

import argparse
parser = argparse.ArgumentParser(description='Plot')
parser.add_argument('--glob', type=str, default='output_*', help='Input pattern')
parser.add_argument('--outdir', type=str, default='plots', help='Directory to store plots')
parser.add_argument('-f', '--format', type=str, default='png', help='Format of the outputs (png, pdf, …)')


args = parser.parse_args()
print(args)

outputs = sorted(glob(args.glob))
prefix = args.outdir
ramsesdir = os.path.split(args.glob)[0]
ramsesdir = '.' if ramsesdir == '' else ramsesdir
ext = args.format
now = datetime.now()

if not os.path.exists(prefix):
    os.mkdir(prefix)

prevpos = None

from scipy.io import FortranFile as FF
from glob import glob
import os

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

    m   = f.read_reals()
    ind = f.read_ints()
    lvl = f.read_ints()
    try:
        tp  = f.read_reals()
        Z   = f.read_reals()
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
    npart = 0
    for i, cpu in enumerate(paths):
        _ind[i], ndim, _npart[i], _pos[i], _vel[i], _mass[i], _lvl[i], _tp[i], _Z[i] = read_one_cpu(cpu)
        npart += _npart[i]

    pos = np.zeros((ndim, npart))
    vel = np.zeros((ndim, npart))
    mass = np.zeros((npart))
    ind = np.zeros((npart), dtype=int)
    lvl = np.zeros((npart), dtype=int)

    n = 0
    for i in range(len(paths)):
        ids = _ind[i] - 1
        pos[:, ids] = _pos[i]
        vel[:, ids] = _vel[i]
        ind[n:n+_npart[i]]  = ids

        mass[ids] = _mass[i]
        isTracer = ids[_mass[i] == 0]
        mass[isTracer] = _tp[i]

        lvl[ids]  = _lvl[i]

        n += _npart[i]

    return ind, pos, vel, mass, lvl

for output in tqdm(outputs[::10]):
    outputn = int(output.split('_')[-1])

    # Load ramses output
    r = pymses.RamsesOutput(ramsesdir, outputn, verbose=False)

    nbin = 2**7 #int(np.sqrt(map.map.shape[0]))

    def saveAndNotify(fname):
        plt.savefig(fname)#, dpi=120)
        print(fname)

    ##########################################
    # Particles
    ##########################################
    # Get them
    _, pos, vel, mass, lvl = read_output(output)
    if prevpos is None:
        prevpos = pos.copy()

    # Select 1024 random particles
    strain = pos.shape[1] // (1024)

    # Estimate displacement
    x, y = pos[0:2, ::strain]
    vx, vy = (pos[0:2, ::strain] - prevpos[0:2, ::strain])
    prevpos = pos.copy()

    # Projecting on grid
    strain2 = 1
    hist_pt, epx, epy = np.histogram2d(pos[0, ::strain2], pos[1, ::strain2],
                                       range=[[0, 1], [0, 1]],
                                       bins=nbin)

    # Plot everything
    plt.clf()
    vmin = 0
    vmax = 4
    percell = 10
    ax = plt.subplot(121)
    plt.title('Particles')
    plt.pcolormesh(epx, epy, hist_pt.T/percell, cmap='viridis', vmin=vmin, vmax=vmax)
    cb = plt.colorbar()
    cb.set_label('density')

    plt.quiver(x, y, vx, vy, scale_units='xy', scale=.5, angles='xy')

    plt.xlim(0, 1)
    plt.ylim(0, 1)

    ##########################################
    # Gas
    ##########################################
    # Get AMR field
    vx_op = ScalarOperator(lambda dset: dset["vel"][:, 0], r.info["unit_density"])
    vy_op = ScalarOperator(lambda dset: dset["vel"][:, 1], r.info["unit_density"])
    rho_op = ScalarOperator(lambda dset: dset["rho"], r.info["unit_density"])
    amr = r.amr_source(['rho', 'vel'])
    rhomap = SliceMap(amr, cam, rho_op, use_multiprocessing=False)
    vxmap = SliceMap(amr, cam, vx_op, use_multiprocessing=False)
    vymap = SliceMap(amr, cam, vy_op, use_multiprocessing=False)


    # Plot
    ax = plt.subplot(122)
    plt.title('Gas map')

    plt.imshow(rhomap.map.T,
               extent=(0, 1, 0, 1),
               aspect='auto',
               vmin=vmin, vmax=vmax, # density
               cmap='viridis'
    )
    cb = plt.colorbar()
    cb.set_label(u'Density [g.cm³]')

    # Velocity map
    xx = np.linspace(0, 1, vxmap.map.shape[0])
    yy = np.linspace(0, 1, vxmap.map.shape[1])
    xs = ys = 8
    plt.quiver(xx[::xs], yy[::ys], vxmap.map.T[::xs, ::ys], vymap.map.T[::xs, ::ys], angles='xy',
               scale=70, scale_units='xy')

    plt.suptitle('$t=%.3f$' % r.info['time'])

    ##########################################
    # Store images
    ##########################################
    outputn = int(output.split('_')[-1])
    fname = pjoin(prefix, 'density_{:0>5}_{}.{}'.format(
        outputn, now.strftime('%Hh%M_%d-%m-%y'), ext))

    saveAndNotify(fname)
