# -------------- Required Packages
import argparse
import os

# -------------- Import TerminalPlot
from vatpy import TerminalPlot

# -------------- Config
import config

homedir = config.homedir
mplstyle = config.mplstyle

# -------------- Arguments
# Initialize argparse:
formatter_class = argparse.RawDescriptionHelpFormatter
parser = argparse.ArgumentParser(description='VATPY Plot Script',
                                 usage='pl [options] snapshot',
                                 formatter_class=formatter_class)
parser._actions[0].help = 'Show this help message'

# Positional arguments:
parser.add_argument('snapshot', help='Snapshot to analyse')

# Optional arguments:
parser.add_argument('-info', '--information', action='store_true',
                    help='''
                    Get some general information about the snapshot
                    ''')
parser.add_argument('-dens', '--density', action='store_true',
                    help='''
                    Generate a gas column density map
                    ''')
parser.add_argument('-temp', '--temperature', action='store_true',
                    help='''
                    Generate a gas column temperature map
                    ''')
parser.add_argument('-stellar', '--stellar', action='store_true',
                    help='''
                    Generate a stellar surface density map
                    ''')
parser.add_argument('-dm', '--darkmatter', action='store_true',
                    help='''
                    Generate a dark matter surface density map
                    ''')
parser.add_argument('-sf', '--starformation', action='store_true',
                    help='''
                    Generate a star formation surface density map
                    ''')
parser.add_argument('-sa', '--stellarage', action='store_true',
                    help='''
                    Generate a star particle stellar age map
                    ''')
parser.add_argument('-bhevol', '--blackholeevolution', action='store_true',
                    help='''
                    Time evolution of the black hole sink particle properties
                    ''')

parser.add_argument('-interpolation', '--interpolation', action='store',
                    default='kdtree', type=str,
                    help='Interpolation technique (default: kdtree)')
parser.add_argument('-xrange', '--xrange', action='store', default=None,
                    nargs=2, type=float,
                    help='Interpolation xrange (default: 0 to boxsize)')
parser.add_argument('-yrange', '--yrange', action='store', default=None,
                    nargs=2, type=float,
                    help='Interpolation yrange (default: 0 to boxsize)')
parser.add_argument('-zrange', '--zrange', action='store', default=None,
                    nargs=2, type=float,
                    help='Interpolation zrange (default: 0 to boxsize)')
parser.add_argument('-box', '--box', action='store', default=None, nargs=2,
                    type=float,
                    help='Interpolation box (default: 0 to boxsize)')
parser.add_argument('-bins', '--bins', action='store', default=100, type=int,
                    help='Number of bins in x/y/z (default: 100)')
parser.add_argument('-sfb', '--starformationbins', action='store', default=100,
                    type=int, help='''
                    Number of bins in x/y when calculating the star formation
                    rate surface density (default: 100)
                    ''')
parser.add_argument('-vcr', '--variablecircradius', action='store_true', 
                    default=False, help='''
                    Whether the flags for a variable circularisation radius is
                    on or not
                    ''')

parser.add_argument('-qty', '--quantity', action='store', default='mass',
                    type=str, help='''
                    Gas quantity [mass/n/HI/HII/H2/CO/He/e] (default: mass)
                    ''')
parser.add_argument('-ul', '--unitlength', action='store', default='kpc',
                    type=str, help='Unit length [kpc/pc] (default: kpc)')
parser.add_argument('-axis', '--axis', action='store', default='z', type=str,
                    help='Axis of rotation (default: z)')
parser.add_argument('-rotate', '--rotate', action='store', default=0,
                    type=float,
                    help='Amount of rotation [degrees] (default: 0)')
parser.add_argument('-cut', '--cut', action='store', default=None, type=float,
                    help='''
                    Cut in the gas column density/temperature map
                    [coordinate in z] (default: none)
                    ''')
parser.add_argument('-bf', '--blackholefocus', action='store_true',
                    default=False, help='Centre the data on the BH')
parser.add_argument('-age', '--maxstellarage', action='store', default=100, 
                    type=float, help='Maximum stellar age of star particles')

parser.add_argument('-vmin', '--vmin', action='store', default=None,
                    type=float, help='Colorbar vmin value (default: none)')
parser.add_argument('-vmax', '--vmax', action='store', default=None,
                    type=float, help='Colorbar vmax value (default: none)')
parser.add_argument('-xlim', '--xlim', action='store', default=None,
                    nargs=2, type=float, help='Axis xlim (default: xrange)')
parser.add_argument('-ylim', '--ylim', action='store', default=None,
                    nargs=2, type=float, help='Axis ylim (default: yrange)')

parser.add_argument('-path', '--savepath', action='store', default=os.getcwd(),
                    help='''
                    Path to save file at (default: current working directory)
                    ''')
parser.add_argument('-format', '--saveformat', action='store', default='png',
                    type=str, help='''
                    Format to save file as (default: png)
                    ''')
parser.add_argument('-style', '--mplstyle', action='store',
                    default=f'{homedir}/VATPY/mpl/{mplstyle}.mplstyle',
                    type=str, help='''
                    Matplotlib style sheet (see VATPY config file)
                    ''')
parser.add_argument('-interactive', '--interactive', action='store',
                    default=True, type=bool, help='''
                    Whether to use an interactive plot or imgcat
                    (default: True)
                    ''')
parser.add_argument('-width', '--width', action='store', default=None,
                    type=int, help='''
                    Width of imgcat [number of characters] (default: none)
                    ''')

# Read arguments from the command line:
args = parser.parse_args()

# Run VATPY TerminalPlot:
print('\nWelcome to VATPY:')

if args.snapshot:
    print(f'  * Reading data of {args.snapshot}')
    v = TerminalPlot(file=args.snapshot, style=args.mplstyle,
                     savepath=args.savepath, saveformat=args.saveformat,
                     vmin=args.vmin, vmax=args.vmax, xlim=args.xlim,
                     ylim=args.ylim, unitlength=args.unitlength,
                     interactive=args.interactive, width=args.width)

if args.information:
    v.info()

if args.density:
    print('  * Generating a gas column density map')
    v.density(axis=args.axis, rotate=args.rotate, quantity=args.quantity,
              bins=args.bins, blackholefocus=args.blackholefocus,
              xrange=args.xrange, yrange=args.yrange, zrange=args.zrange,
              box=args.box, cut=args.cut)

if args.temperature:
    print('  * Generating a gas collumn temperature map')
    v.temperature(axis=args.axis, rotate=args.rotate, bins=args.bins,
                  blackholefocus=args.blackholefocus, xrange=args.xrange, 
                  yrange=args.yrange, zrange=args.zrange, box=args.box, 
                  cut=args.cut)

if args.stellar:
    print('  * Generating a stellar surface density map')
    v.stellar(axis=args.axis, rotate=args.rotate, bins=args.bins,
              xrange=args.xrange, yrange=args.yrange, zrange=args.zrange)

if args.darkmatter:
    print('  * Generating a dark matter surface density map')
    v.darkmatter(axis=args.axis, rotate=args.rotate, bins=args.bins,
                 xrange=args.xrange, yrange=args.yrange, zrange=args.zrange)

if args.starformation:
    print('  * Generating a star formation surface density map')
    v.star_formation(axis=args.axis, rotate=args.rotate, bins=args.bins,
                     sfb=args.starformationbins,
                     blackholefocus=args.blackholefocus, xrange=args.xrange,
                     yrange=args.yrange, zrange=args.zrange, box=args.box,
                     cut=args.cut)

if args.stellarage:
    print('  * Generating a star particle stellar age map')
    v.stellar_age(axis=args.axis, rotate=args.rotate, bins=args.bins,
                  age=args.maxstellarage, blackholefocus=args.blackholefocus, 
                  xrange=args.xrange, yrange=args.yrange, zrange=args.zrange, 
                  box=args.box, cut=args.cut)

if args.blackholeevolution:
    print('  * Generating a time evolution plot of the black hole sink'
          + ' particle properties')
    v.black_hole_evolution(vcr=args.variablecircradius)

print('  * Run completed\n')

# -------------- End of file
