#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- Required Packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse
import h5py
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import ArtistAnimation
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from imgcat import imgcat

# -------------- Import TerminalPlot
from src import TerminalPlot

# -------------- Arguments
# Initialize argparse:
parser = argparse.ArgumentParser(description='VATPY Terminal Plot Script', usage='pl [options] snapshot', formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# Positional arguments:
parser.add_argument('snapshot', help='Snapshot to analyse')

# Optional arguments:
parser.add_argument('-vmin', '--vmin', action='store', help='Colorbar vmin value', default=None)
parser.add_argument('-vmax', '--vmax', action='store', help='Colorbar vmax value', default=None)
parser.add_argument('-xlim', '--xlim', action='store', help='Axis xlim', default=None, nargs=2, type=float)
parser.add_argument('-ylim', '--ylim', action='store', help='Axis ylim', default=None, nargs=2, type=float)
parser.add_argument('-box', '--box', action='store', help='Box zoom in', default=None, nargs=2, type=float)
parser.add_argument('-bins', '--numberofbins', action='store', help='Number of bins', default=100, type=int)
parser.add_argument('-savepath', '--savepath', action='store', help='Save path', default=os.getcwd())
parser.add_argument('-saveformat', '--saveformat', action='store', help='Save format', default=None)
parser.add_argument('-style', '--mplstyle', action='store', help='Matplotlib style sheet', 
                    default='/home/astro/jpeterss/VATPY/mpl/tplot.mplstyle')

parser.add_argument('-info', '--information', action='store_true', help='Information about snapshot')

parser.add_argument('-dens', '--density', action='store_true', help='2D gas density plot')
parser.add_argument('-axis', '--axis', action='store', help='Look down axis', default='z')
parser.add_argument('-unit', '--unit', action='store', help='Gas density unit', default='cgs')
parser.add_argument('-cut', '--cut', action='store', help='Gas density cut', default=None)
parser.add_argument('-col', '--column', action='store', help='Gas density column', default=None, nargs=2, type=float)

parser.add_argument('-temp', '--temperature', action='store_true', help='2D gas temperature plot')
parser.add_argument('-phase', '--phasediagram', action='store_true', help='Phase diagram')
parser.add_argument('-num', '--numberdensity', action='store_true', help='Number density instead of mass density')
parser.add_argument('-resolution', '--resolution', action='store_true', help='Resolution plot')
parser.add_argument('-movie', '--movie', action='store_true', help='Create a movie')


# Read arguments from the command line:
args = parser.parse_args()


# Run VATPY TerminalPlot:
print('\nWelcome to VATPY:')

if args.snapshot:
    print(f'  * Begin by reading data of {args.snapshot}')
    v = TerminalPlot(file=args.snapshot, savepath=args.savepath, vmin=args.vmin, vmax=args.vmax, xlim=args.xlim, ylim=args.ylim, style=args.mplstyle, saveformat=args.saveformat)

if args.information:
    v.info()

if args.density:
    print(f'  * Now generate figure of the gas density')
    v.density(bins=args.numberofbins, axis=args.axis, unit=args.unit, cut=args.cut, column=args.column, box=args.box)

if args.temperature:
    print(f'  * Now generate figure of the gas temperature')
    v.temperature(bins=args.numberofbins, axis=args.axis, cut=args.cut, column=args.column, box=args.box)

if args.phasediagram:
	v.phase_diagram(bins=args.numberofbins, num=args.numberdensity)

if args.resolution:
    v.resolution(bins=args.numberofbins)

if args.movie:
	v.movie()

print('  * The End\n')


