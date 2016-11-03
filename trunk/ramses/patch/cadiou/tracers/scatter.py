# -*- coding: utf-8 -*-
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import numpy as np
from glob import glob
from partutils import read_output
# Setup plot
minorloc = MultipleLocator(base=1./2**8)
majorloc = MultipleLocator(base=1./2**7)

# ax.xaxis.set_major_locator(majorloc)
# ax.xaxis.set_minor_locator(minorloc)

# ax.yaxis.set_major_locator(majorloc)
# ax.yaxis.set_minor_locator(minorloc)


outputs = sorted(glob('output_*'))

i = 0

prevPos = None


def nextOutput():
    global i, outputs, prevPos

    output = outputs[i]
    ind, pos, vel, mass, lvl = read_output(output)

    if prevPos is None:
        prevPos = pos.copy()

    dx = pos - prevPos

    ticks = np.linspace(0, 1, 2**8 + 1)
    fig = plt.gcf()
    ax = fig.gca()
    ax.cla()

    xmin, xmax = 0.45, 0.55
    ymin, ymax = 0.75, 0.85
    mask = ((pos[0, :] >= xmin) * (pos[0, :] < xmax) *
            (pos[1, :] >= ymin) * (pos[1, :] < ymax))
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.grid('on')
    ax.quiver(pos[0, mask], pos[1, mask], dx[0, mask], dx[1, mask],
              scale_units='xy', angles='xy', scale=1,
              color='grey', alpha=.5,
              pivot='tip',
              linewidths=.5)
    ax.scatter(*pos[:, mask], c='blue')

    plt.title('Output %s' % output)

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)

    plt.draw_if_interactive()
    i += 1
    prevPos = pos.copy()
