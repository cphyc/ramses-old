# -*- coding: utf-8 -*-
import argparse
import numpy as np
import itertools

parser = argparse.ArgumentParser(description='Generate initial conditions for the tracers on a uniform grid.')

parser.add_argument('-l', '--level', help='The level of refinement of the grid', type=int, required=True)
parser.add_argument('-np', '--nparticle', help='Number of particles per grid point', type=int, default=1)
parser.add_argument('-nd', '--ndim', help='Number of dimensions', default=3, type=int)
parser.add_argument('-o', '--output', help='The output file', required=True, type=str)
parser.add_argument('--disp', help='A displacement in unit of 1/2**level', default=0, type=float)

args = parser.parse_args()

ntotpart = args.nparticle*(2**args.level)**args.ndim
particles = np.zeros((3, ntotpart))

level = args.level
partpercell = args.nparticle

with open(args.output, 'w') as f:
    for i, j, k in itertools.product(range(2**(level-1)), range(2**(level-1)), range(2**(level-1))):
        dx = 0.5**(level)
        x = (i + args.disp) / 2.**(level-1) + 0.5**(level)
        y = (j + args.disp) / 2.**(level-1) + 0.5**(level)
        z = (k + args.disp) / 2.**(level-1) + 0.5**(level)

        if args.ndim == 2:
            z = 0
            if k > 0:
                continue

        for p in range(partpercell):
            repeat = 1 if x > .5 else 2
            for _ in range(repeat):
                f.write('%.12f %.12f %.12f %.12f %.12f %.12f %.12f\n' % (x, y, z, 0, 0, 0, 1))

