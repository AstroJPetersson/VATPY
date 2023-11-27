#!/home/users/j/jpeterss/anaconda3/bin/python


# -------------- Config
homedir  = '/home/users/j/jpeterss'
mplstyle = 'tplot_ver1'


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
from vatpy import TerminalPlot


# -------------- Arguments
# Initialize argparse:
parser = argparse.ArgumentParser(description='VATPY Terminal Plot Script', usage='pl [options] snapshot', 
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'


# Positional arguments:
parser.add_argument('snapshot', help='Snapshot to analyse')


# Optional arguments:
parser.add_argument('-vmin', '--vmin', action='store', default=None, type=float, help='Colorbar vmin value')
parser.add_argument('-vmax', '--vmax', action='store', default=None, type=float, help='Colorbar vmax value')
parser.add_argument('-xlim', '--xlim', action='store', default=None, nargs=2, type=float, help='Axis xlim')
parser.add_argument('-ylim', '--ylim', action='store', default=None, nargs=2, type=float, help='Axis ylim')
parser.add_argument('-bins', '--numberofbins', action='store', default=100, type=int, help='Number of bins')
parser.add_argument('-levels', '--numberoflevels', action='store', default=5, type=int, help='Number of levels')
parser.add_argument('-savepath', '--savepath', action='store', default=os.getcwd(), help='Path to save at')
parser.add_argument('-saveformat', '--saveformat', action='store', default=None, help='Format to save in')
parser.add_argument('-style', '--mplstyle', action='store', default=f'{homedir}/VATPY/mpl/{mplstyle}.mplstyle', 
                    help='Matplotlib style option')
parser.add_argument('-interactive', '--interactive', action='store_true', default=False, 
                    help="Interactive plot instead of imgcat")
parser.add_argument('-imgcatwidth', '--imgcatwidth', action='store', default=None, type=int, 
                    help="Width of imgcat figure [number of characters]")
parser.add_argument('-info', '--information', action='store_true', help='General information about the snapshot')
parser.add_argument('-dens', '--density', action='store_true', help='Gas density map (either as a slice or column)')
parser.add_argument('-qty', '--quantity', action='store', default='mass', help='Quantity to show in the gas density map'
                    + ' (mass/n, or if available: HII/H2/HI/CO/He/e)')
parser.add_argument('-axis', '--axis', action='store', default='z', help='Axis of rotation')
parser.add_argument('-rotate', '--rotate', action='store', default=0, type=float, 
                    help='Amount to rotate around the rotation axis [degrees]')
parser.add_argument('-xrange', '--xrange', action='store', default=None, nargs=2, type=float, help='Interpolation xrange')
parser.add_argument('-yrange', '--yrange', action='store', default=None, nargs=2, type=float, help='Interpolation yrange')
parser.add_argument('-zrange', '--zrange', action='store', default=None, nargs=2, type=float, help='Interpolation zrange')
parser.add_argument('-temp', '--temperature', action='store_true', help='Gas temperature map (either as a slice or column)')
parser.add_argument('-cut', '--cut', action='store', default=None, help='Cut in the gas density/temperature map')
parser.add_argument('-stellar', '--stellar', action='store_true', help='Stellar density map')
parser.add_argument('-dm', '--darkmatter', action='store_true', help='Dark matter density map')
parser.add_argument('-sf', '--starformation', action='store_true', help='Star formation surface density map')
parser.add_argument('-movie', '--movie', action='store_true', help='Create a movie')


# Read arguments from the command line:
args = parser.parse_args()


# Run VATPY TerminalPlot:
print('\nWelcome to VATPY:')

if args.snapshot:
    print(f'  * Reading data of {args.snapshot}')
    v = TerminalPlot(file=args.snapshot, savepath=args.savepath, vmin=args.vmin, vmax=args.vmax, xlim=args.xlim, ylim=args.ylim, 
                     style=args.mplstyle, saveformat=args.saveformat, interactive=args.interactive, imgcatwidth=args.imgcatwidth)

if args.information:
    v.info()

if args.density:
    print(f'  * Generating gas surface density map\n')
    v.density(axis=args.axis, rotate=args.rotate, quantity=args.quantity, bins=args.numberofbins, 
              xrange=args.xrange, yrange=args.yrange, zrange=args.zrange, cut=args.cut)

if args.temperature:
    print(f'  * Generating gas mean temperature map\n')
    v.temperature(axis=args.axis, rotate=args.rotate, bins=args.numberofbins, 
                  xrange=args.xrange, yrange=args.yrange, zrange=args.zrange, cut=args.cut)

if args.stellar:
    print(f'  * Generating stellar surface density map\n')
    v.stellar(axis=args.axis, rotate=args.rotate, bins=args.numberofbins, xrange=args.xrange, yrange=args.yrange, zrange=args.zrange)

if args.darkmatter:
    print(f'  * Generating dark matter surface density map\n')
    v.darkmatter(axis=args.axis, rotate=args.rotate, bins=args.numberofbins, xrange=args.xrange, yrange=args.yrange, zrange=args.zrange)

if args.starformation:
    print(f'  * Generating star formation surface density map\n')
    v.star_formation(axis=args.axis, rotate=args.rotate, bins=args.numberofbins, xrange=args.xrange, yrange=args.yrange, zrange=args.zrange)

if args.movie:
    print(f'  * Creating a movie')
    v.movie()

print('  * Run completed\n')


# -------------- End of file


