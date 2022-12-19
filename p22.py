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

# -------------- Import Project
from vatpy import Project2022

# -------------- Plot arguments
# Initialize argparse:
parser = argparse.ArgumentParser(description='Project 2022 Script', usage='proj22 [options]', formatter_class=argparse.RawDescriptionHelpFormatter)
parser._actions[0].help='Show this help message and exit'

# Arguments:
parser.add_argument('-test', '--testfunction', help='Test function', action='store_true')

# Read arguments from the command line:
args = parser.parse_args()

# Run Project:
if args.testfunction:
    print('Works!')


