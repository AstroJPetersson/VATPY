#!/home/astro/jpeterss/anaconda3/bin/python

# -------------- packages
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py
import os
from matplotlib.animation import ArtistAnimation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from imgcat import imgcat
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

# -------------- import VATPY
from src.vatpy import Vatpy

# -------------- arguments
# initialize argparse:
hpage = ('powered by...\n\n'+
		 '`..         `..      `.       `... `......`.......  `..      `.. \n'+
		 ' `..       `..      `. ..          `..    `..    `.. `..    `.. \n'+
		 '  `..     `..      `.  `..         `..    `..    `..  `.. `.. \n'+
		 '   `..   `..      `..   `..        `..    `.......      `.. \n'+
		 '    `.. `..      `...... `..       `..    `..           `.. \n'+
		 '     `....      `..       `..      `..    `..           `.. \n'+
		 '      `..      `..         `..     `..    `..           `.. \n')
parser = argparse.ArgumentParser(description=hpage, usage='plot [options] filename', formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# positional arguments:
parser.add_argument('filename', help='Name of the snapshot you would like to analyse')

# Optional arguments:
parser.add_argument('-dens', '--density', action='store_true', help='Gas density plot')
parser.add_argument('-unit', '--unit', action='store', help='Unit of the gas density plot (number density, chemical species, etc)', default='cgs')
parser.add_argument('-zslice', '--zslice', action='store', help='Slice (at given z) of the gas density plot', default=None)
parser.add_argument('-vmin', '--vmin', action='store', help='Colorbar vmin value', default=None)
parser.add_argument('-vmax', '--vmax', action='store', help='Colorbar vmax value', default=None)
parser.add_argument('-xlim', '--xlim', action='store', help='Axis xlim', default=None, nargs=2, type=float)
parser.add_argument('-ylim', '--ylim', action='store', help='Axis ylim', default=None, nargs=2, type=float)
parser.add_argument('-bins', '--numberofbins', action='store', help='Number of bins', default=100, type=int)
parser.add_argument('-save', '--savepath', action='store', help='Path where to save figures, animations, etc', default=os.getcwd())
parser.add_argument('-movie', '--movie', action='store_true', help='Movie of how the gas density evolves')
parser.add_argument('-dist', '--distribution', action='store_true', help='Various PDFs of the gas density, temperature, pressure, etc')
parser.add_argument('-temp', '--temperature', action='store_true', help='Gas temperature plot')
parser.add_argument('-phase', '--phasediagram', action='store_true', help='Gas-phase diagram')
parser.add_argument('-cell', '--cellmasssizerelation', action='store_true', help='Cell-mass-size relation')
parser.add_argument('-style', '--mplstyle', action='store', help='Mpl style sheet', default='/home/astro/jpeterss/VATPY/mpl/style.mplstyle')
parser.add_argument('-sink', '--sink', action='store_true', help='Gas density plot with sink particles included')

# read arguments from the command line:
args = parser.parse_args()

# run VATPY:
if args.filename:
	v = Vatpy(f=args.filename, save=args.savepath, vmin=args.vmin, vmax=args.vmax, xlim=args.xlim, ylim=args.ylim, style=args.mplstyle)

if args.density:
	v.density_plot(unit=args.unit, zslice=args.zslice, bins=args.numberofbins)

if args.sink:
	v.sink_plot(unit=args.unit, zslice=args.zslice, bins=args.numberofbins)

if args.movie:
	v.movie(unit=args.unit, zslice=args.zslice, bins=args.numberofbins)

if args.distribution:
	v.distribution_plot()

if args.temperature:
	v.temperature_plot(zslice=args.zslice, bins=args.numberofbins)

if args.phasediagram:
	v.phase_plot()

if args.cellmasssizerelation:
	v.cell_mass_size_relation()




