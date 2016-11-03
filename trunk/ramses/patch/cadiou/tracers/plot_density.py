# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from glob import glob
from tqdm import tqdm
import os
from os.path import join as pjoin
from datetime import datetime
import argparse

import pymses
from pymses.analysis.slicing import SliceMap
from pymses.analysis.operator import ScalarOperator, MaxLevelOperator
from pymses.analysis import Camera

from partutils import read_output
#################################
# Viewing config
#################################
cam = Camera(line_of_sight_axis='z', up_vector='y')

mpl.rcParams['figure.figsize'] = (16, 9)
mpl.rcParams['figure.dpi'] = 120

#################################
# Parse input
#################################
parser = argparse.ArgumentParser(description='Plot')
parser.add_argument('--glob', type=str, default='output_*',
                    help='Input pattern')
parser.add_argument('--outdir', type=str, default='plots',
                    help='Directory to store plots')
parser.add_argument('-f', '--format', type=str, default='png',
                    help='Format of the outputs (png, pdf, …)')


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


def oneOutput(output):
    global prevpos
    outputn = int(output.split('_')[-1])

    # Load ramses output
    r = pymses.RamsesOutput(ramsesdir, outputn, verbose=False)

    nbin = 2**7  # int(np.sqrt(map.map.shape[0]))
    percell = 10.
    vmin = 0
    vmax = 4

    def saveAndNotify(fname):
        plt.savefig(fname)  # , dpi=120)
        print(fname)

    plt.clf()
    ##########################################
    # Gas
    ##########################################
    # Get AMR field
    vx_op = ScalarOperator(lambda dset: dset["vel"][:, 0],
                           r.info["unit_density"])
    vy_op = ScalarOperator(lambda dset: dset["vel"][:, 1],
                           r.info["unit_density"])
    rho_op = ScalarOperator(lambda dset: dset["rho"], r.info["unit_density"])
    amr = r.amr_source(['rho', 'vel'])
    rhomap = SliceMap(amr, cam, rho_op, use_multiprocessing=False)
    vxmap = SliceMap(amr, cam, vx_op, use_multiprocessing=False)
    vymap = SliceMap(amr, cam, vy_op, use_multiprocessing=False)
    lvlmap = SliceMap(amr, cam, MaxLevelOperator(), use_multiprocessing=False)

    # Plot
    plt.subplot(122)
    plt.title('Gas map')

    plt.imshow(rhomap.map.T[::-1],
               extent=(0, 1, 0, 1),
               aspect='auto',
               vmin=vmin, vmax=vmax,
               cmap='viridis')
    cb = plt.colorbar()
    cb.set_label(u'Density [g.cm³]')

    # Velocity map
    xx = np.linspace(0, 1, vxmap.map.shape[0])
    yy = np.linspace(0, 1, vxmap.map.shape[1])
    xs = ys = 8
    plt.quiver(xx[::xs], yy[::ys],
               vxmap.map.T[::xs, ::ys], vymap.map.T[::xs, ::ys],
               angles='xy',
               scale=70, scale_units='xy')

    ##########################################
    # Particles
    ##########################################
    # Get them
    _, pos, vel, mass, lvl = read_output(output)
    if prevpos is None:
        prevpos = pos.copy()

    # Select 1024 random particles
    strain = max(pos.shape[1] // (1024), 1)

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

    plt.subplot(121)
    plt.title('Particles')
    plt.pcolormesh(epx, epy, hist_pt.T/percell,
                   cmap='viridis', vmin=vmin, vmax=vmax)
    cb = plt.colorbar()
    cb.set_label('density')

    plt.quiver(x, y, vx, vy,
               scale_units='xy', scale=1, angles='xy', color='white')

    # Plot lvl contours (if need be)
    ncontours = np.ptp(lvlmap.map)
    if ncontours > 0:
        plt.contour(lvlmap.map.T,
                    extent=(epx[0], epx[-1], epy[0], epy[-1]),
                    levels=np.linspace(lvlmap.map.min()+0.5,
                                       lvlmap.map.max()-0.5,
                                       ncontours),
                    alpha=1)
    plt.xlim(0, 1)
    plt.ylim(0, 1)

    plt.suptitle('$t=%.3f$' % r.info['time'])

    ##########################################
    # Store images
    ##########################################
    outputn = int(output.split('_')[-1])
    fname = pjoin(prefix, 'density_{:0>5}_{}.{}'.format(
        outputn, now.strftime('%Hh%M_%d-%m-%y'), ext))

    saveAndNotify(fname)

    return pos

for output in tqdm(outputs[:70:2]):
    oneOutput(output)
