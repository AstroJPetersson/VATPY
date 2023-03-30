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
parser = argparse.ArgumentParser(description='VATPY Terminal Plot Script', usage='pl [options] filename', formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# Positional arguments:
parser.add_argument('snapshot', help='Name of the snapshot you would like to analyse')

# Optional arguments:
parser.add_argument('-vmin', '--vmin', action='store', help='Colorbar vmin value', default=None)
parser.add_argument('-vmax', '--vmax', action='store', help='Colorbar vmax value', default=None)
parser.add_argument('-xlim', '--xlim', action='store', help='Axis xlim', default=None, nargs=2, type=float)
parser.add_argument('-ylim', '--ylim', action='store', help='Axis ylim', default=None, nargs=2, type=float)
parser.add_argument('-bins', '--numberofbins', action='store', help='Number of bins', default=100, type=int)
parser.add_argument('-savepath', '--savepath', action='store', help='Path to save at', default=os.getcwd())
parser.add_argument('-saveformat', '--saveformat', action='store', help='Format to save as', default=None)
parser.add_argument('-style', '--mplstyle', action='store', help='Matplotlib style', 
                    default='/home/astro/jpeterss/VATPY/mpl/tplot.mplstyle')

parser.add_argument('-dens', '--density', action='store_true', help='2D gas density plot')
parser.add_argument('-axis', '--lookdownaxis', action='store', help='Look down axis', default='z')
parser.add_argument('-unit', '--unit', action='store', help='Unit of the gas density', default='cgs')
parser.add_argument('-zslice', '--zslice', action='store', help='Slice of the gas density', default=None)
parser.add_argument('-zcol', '--zcol', action='store', help='Column slice of the  gas density', default=None, 
                    nargs=2, type=float)

parser.add_argument('-temp', '--temperature', action='store_true', help='2D gas temperature plot')
parser.add_argument('-phase', '--phase', action='store_true', help='Phase diagram of the gas')
parser.add_argument('-resolution', '--resolution', action='store_true', help='Resolution plot of the gas')
parser.add_argument('-jeans', '--jeans', action='store_true', help='Jeans length extension in the resolution plot')
parser.add_argument('-info', '--information', action='store_true', help='Information about snapshot')
parser.add_argument('-movie', '--movie', action='store_true', help='Generate and save a movie')

# Read arguments from the command line:
args = parser.parse_args()

# Run VATPY TerminalPlot:
if args.snapshot:
	v = TerminalPlot(file=args.snapshot, savepath=args.savepath, vmin=args.vmin, vmax=args.vmax, xlim=args.xlim, ylim=args.ylim, style=args.mplstyle, saveformat=args.saveformat)

if args.density:
	v.density(lookDownAxis=args.lookdownaxis, unit=args.unit, zSlice=args.zslice, zCol=args.zcol, bins=args.numberofbins)

if args.temperature:
	v.temperature()

if args.phase:
	v.phase()

if args.resolution:
    v.resolution(bins=args.bins, jeans=args.jeans)

if args.information:
    v.info()

if args.movie:
	v.movie()


