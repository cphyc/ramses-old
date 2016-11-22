# -*- coding: utf-8 -*-
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from partutils import read_output
import os
# Setup plot
minorloc = MultipleLocator(base=1./2**5)
majorloc = MultipleLocator(base=1./2**6)
# ticks = np.linspace(0, 1, 2**4)

# ax.xaxis.set_major_locator(majorloc)
# ax.xaxis.set_minor_locator(minorloc)

# ax.yaxis.set_major_locator(majorloc)
# ax.yaxis.set_minor_locator(minorloc)


i = 0

prevPos = cb = None
limitSet = False
outputs = []
prefix = '.'

def setPrefix(p = '.'):
    global prefix
    if os.path.isdir(p):
        prefix = p
    else:
        raise Exception('%s is not a directory' % p)

def refreshOutputs():
    global prefix

    print(prefix)
    outputs = sorted(glob(os.path.join(prefix, 'output_*')))

    return outputs


def nextOutput(nexti=None):
    global i, prevPos, cb, limitSet

    if nexti is not None and nexti >= 0:
        i = nexti
    elif nexti is not None and nexti < 0:
        i = i + nexti
    else:
        i = i + 1

    outputs = refreshOutputs()
    output = outputs[i]
    print('reading %s' % output)
    ind, pos, vel, mass, lvl, cpu = read_output(output)

    if prevPos is None:
        prevPos = pos.copy()

    dx = pos - prevPos

    ticks = np.linspace(0, 1, 2**6 + 1)
    ax = plt.gca()

    if limitSet:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
    else:
        xmin, xmax = 0, 1
        ymin, ymax = 0, 1

    fig = plt.gcf()
    ax = fig.gca()
    ax.cla()

    # xmin, xmax = 0.45, 0.55
    # ymin, ymax = 0.4, 0.5
    mask = (((pos[0, :] >= xmin) * (pos[0, :] < xmax) *
             (pos[1, :] >= ymin) * (pos[1, :] < ymax)) +
            ((prevPos[0, :] >= xmin) * (prevPos[0, :] < xmax) *
             (prevPos[1, :] >= ymin) * (prevPos[1, :] < ymax)))

    nbin = 512 #256
    H, ex, ey = np.histogram2d(*(pos[:, :]), bins=nbin, range=[[0, 1], [0, 1]])
    H = np.ma.array(H, mask=(H == 0))
    plt.pcolormesh(ex, ey, H.T,
                   vmin=1, vmax=64, cmap='viridis')
    # ax.scatter(*pos[:, mask], c='blue')

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.grid('on')
    ax.quiver(pos[0, mask], pos[1, mask], dx[0, mask], dx[1, mask],
              scale_units='xy', angles='xy', scale=1,
              color='grey', alpha=.5,
              pivot='tip',
              linewidths=.5)

    plt.title('Output %s' % output)

    limitSet = True

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.grid('on')

    plt.draw_if_interactive()
    plt.setp(plt.xticks()[1], rotation=90)

    prevPos = pos.copy()

def nextOutputHist(nexti=None):
    global i, prevPos, cb, limitSet

    if nexti is not None and nexti >= 0:
        i = nexti
    elif nexti is not None and nexti < 0:
        i = i + nexti
    else:
        i = i + 1

    outputs = refreshOutputs()
    output = outputs[i]
    print('Reading %s' % output)
    ind, pos, vel, mass, lvl, cpus = read_output(output)

    if prevPos is None:
        prevPos = pos.copy()

    dx = pos - prevPos

    ticks = np.linspace(0, 1, 2**5 + 1)
    ax = plt.gca()

    if limitSet:
        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()
    else:
        xmin, xmax = 0, 1
        ymin, ymax = 0, 1

    fig = plt.gcf()
    ax = fig.gca()
    ax.cla()

    # xmin, xmax = 0.45, 0.55
    # ymin, ymax = 0.4, 0.5
    mask = ((pos[0, :] >= xmin) * (pos[0, :] < xmax) *
            (pos[1, :] >= ymin) * (pos[1, :] < ymax))

    nbin = 2**4
    H, ex, ey = np.histogram2d(*pos[:, :], bins=nbin, range=[[0, 1], [0, 1]])
    H = np.ma.array(H, mask=(H == 0))
    plt.pcolormesh(ex, ey, H.T,
                   #vmin=1, vmax=64,
                   cmap='viridis')
    # ax.scatter(*pos[:, mask], c='blue')

    f = 1
    cpuMap = np.zeros((f*nbin, f*nbin))
    # Getting cpu map
    for cpu in range(cpus.min(), cpus.max()+1):
        mask = (cpus == cpu)

        # get center of domain and annotate it
        cx, cy = np.mean(pos[:, mask], axis=1)
        plt.annotate(str(cpu), (cx, cy))
        # get map for cpu
        Hcpu, _, _ = np.histogram2d(*pos[:, mask], bins=f*cpuMap.shape[0], range=[[0, 1], [0, 1]])
        tmparr = cpu*(Hcpu.flatten() > 0)
        cpuMap += np.reshape(tmparr, Hcpu.shape)

    if cpus.min() < cpus.max():
        plt.contour(cpuMap.T, extent=(ex[0], ex[-1], ey[0], ey[-1]),
                    levels=list(range(cpus.min(), cpus.max())), alpha=0.5)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.grid('on')
    # ax.quiver(pos[0, mask], pos[1, mask], dx[0, mask], dx[1, mask],
    #           scale_units='xy', angles='xy', scale=1,
    #           color='grey', alpha=.5,
    #           pivot='tip',
    #           linewidths=.5)

    plt.title('Output %s' % output)

    limitSet = True

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.grid('on')

    plt.draw_if_interactive()
    plt.setp(plt.xticks()[1], rotation=90)
    prevPos = pos.copy()

def countFlux():
    outputs = refreshOutputs()
    prev = read_output(outputs[0])
    for i in range(len(outputs) - 1):
        curr = read_output(outputs[i+1])
        ppos = prev[1]
        cpos = curr[1]

        surf1 = 0.473 # coarse
        surf2 = 0.5   # between levels
        surf3 = 0.516 # fine
        phi1, phi2, phi3 = [((ppos < y) * (cpos > y)).sum()
                            for y in [surf1, surf2, surf3]]
        yield (phi1, phi2, phi3)
        prev = curr

def traceParticle(ids, begin=0, end=-1):
    outputs = refreshOutputs()

    ids = np.array(ids)
    if end == -1:
        end = len(outputs)
    pos = np.zeros((len(ids), 2, end-begin))
    for i, output in enumerate(outputs[begin:end]):
        print(output)
        _, X, _, _, _ = read_output(output)

        pos[:, 0:2, i] = X[0:2, ids-1].T

    c = ['red', 'green', 'blue', 'cyan', 'grey', 'yellow', 'black']
    for i, _ in enumerate(ids):
        plt.plot(*pos[i, :, :], color=c[i % len(c)])
        plt.plot(*pos[i, :, -1], marker='o', color=c[i % len(c)])
    return pos
