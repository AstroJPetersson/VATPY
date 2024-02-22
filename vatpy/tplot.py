'''
Description:

Last updated: 2023-09-27
'''


# -------------- Required packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
import h5py
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import ArtistAnimation
from scipy.interpolate import NearestNDInterpolator
from scipy.spatial.transform import Rotation
from scipy.ndimage import gaussian_filter
from imgcat import imgcat

from .read import read_hdf5
from .get_gas_property import number_density
from .get_gas_property import temperature
from .interpolation import interpolate_to_2d


# -------------- TerminalPlot
class TerminalPlot:
    '''
    TerminalPlot: 
    '''
    def __init__(self, file, savepath=os.getcwd(), vmin=None, vmax=None, 
                 xlim=None, ylim=None, style=None, saveformat=None, 
                 interactive=False, imgcatwidth=None):
        # Variables:
        self.file        = file
        self.savepath    = savepath
        self.vmin        = vmin
        self.vmax        = vmax
        self.xlim        = xlim
        self.ylim        = ylim
        self.saveformat  = saveformat
        self.interactive = interactive
        self.imgcatwidth = imgcatwidth

        # Constants:
        self.G    = 6.67259e-8  #[cm^3 g^-1 s^-2]
        self.mp   = 1.6726e-24  #[g] 
        self.kb   = 1.380658e-16  #[erg K^-1]   
        self.Msol = 1.9891e33  #[g]
        self.pc   = 3.08567758e18  #[cm]
        self.kpc  = 3.08567758e21  #[cm]
        self.Myr  = 1e6 * 365.25 * 24 * 60 * 60  #[s]

        # Mpl style:
        plt.style.use(style)

    ##########################################################################
    ##########################################################################
    def info(self):
        '''
        Description: Print some general information about the given snapshot.
        '''
        # Read the data:
        h, iu = read_hdf5(file=self.file)
        
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        numPart = h['Header'].attrs['NumPart_ThisFile']

        print(f'  * Snapshot information')
        print(f'  | Time: {np.round(time, 2)} Myr')
        print(f'  | BoxSize: {np.round(boxSize, 2)} kpc')
        print(f'  |')
        print(f'  | Number of particles')
        print(f'  | PartType0 (gas):   {numPart[0]}')
        print(f'  | PartType1 (halo):  {numPart[1]}')
        print(f'  | PartType2 (disk):  {numPart[2]}')
        print(f'  | PartType3 (bulge): {numPart[3]}')
        print(f'  | PartType4 (stars): {numPart[4]}')
        print(f'  | PartType5 (bndry): {numPart[5]}')
        
        return 0

    ##########################################################################
    ##########################################################################
    def density(self, axis='z', rotate=0, quantity='mass', bins=100, 
                xrange=None, yrange=None, zrange=None, cut=None):
        '''
        Description: Surface density map of the gas. The function uses a
        nearest neighbour interpolator to interpolate the gas density onto 
        a generated grid, and then sums everything up along a given axis via 
        a 2d histogram. 
        '''
        # Read the data:
        h, iu   = read_hdf5(file=self.file)
        pos     = h['PartType0']['Coordinates'][:] * iu['ulength'] / self.kpc
        dens    = h['PartType0']['Density'][:] * iu['udens']
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc

        # Coordinate ranges:
        if not (xrange):
            xrange=(0, boxSize)
        if not (yrange):
            yrange=(0, boxSize)
        if not (zrange):
            zrange=(0, boxSize)

        # Quantity selection:
        if (quantity != 'mass'):
            numDens = number_density(h, iu)
            dens = numDens[quantity]

        # Rotation:
        if rotate != 0:
            rotation = Rotation.from_euler(axis, rotate, degrees=True)
            pos = rotation.apply(pos - boxSize/2)
            pos += boxSize/2

        # Interpolate the data:
        interpDens = interpolate_to_2d(pos=pos, unit=self.kpc, values=dens, bins=bins, 
                                       xrange=xrange, yrange=yrange, zrange=zrange, cut=cut)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94, wspace=0, hspace=0)

        im = ax.imshow(np.log10(interpDens), extent=(xrange[0], xrange[1], yrange[0], yrange[1]), 
                       origin='lower', vmin=self.vmin, vmax=self.vmax, cmap='inferno')
        ax.text(0.95, 0.05, f'{time:.2f} Myr', color='k', ha='right', va='bottom', transform=ax.transAxes, 
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'})
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)

        # Colorbar:
        if not cut:
            dim = 2
        else:
            dim = 3
        cbar_label = {
            'HII'  : r'$\log_{10}(\Sigma_{\mathrm{HII}} \ [\mathrm{cm}^{-%d}])$' % dim,
            'H2'   : r'$\log_{10}(\Sigma_{\mathrm{H}_2} \ [\mathrm{cm}^{-%d}])$' % dim,
            'HI'   : r'$\log_{10}(\Sigma_{\mathrm{HI}} \ [\mathrm{cm}^{-%d}])$' % dim,
            'CO'   : r'$\log_{10}(\Sigma_{\mathrm{CO}} \ [\mathrm{cm}^{-%d}])$' % dim,
            'He'   : r'$\log_{10}(\Sigma_{\mathrm{He}} \ [\mathrm{cm}^{-%d}])$'% dim,
            'e'    : r'$\log_{10}(\Sigma_{e^{-}} \ [\mathrm{cm}^{-%d}])$' % dim,
            'n'    : r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{cm}^{-%d}])$' % dim,
            'mass' : r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-%d}])$' % dim
        }

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=cbar_label[quantity])

        # Save figure:
        figname = f'density_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        print('  * Figure generated and saved')

        # Display figure:
        if self.interactive != True:
            imgcat(fig, width=self.imgcatwidth)
        else:
            print('  * Interactive display of the figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def temperature(self, axis='z', rotate=0, bins=100, xrange=None, 
                    yrange=None, zrange=None, cut=None):
        '''
        Description: Surface temperature map of the gas. The function uses a
        nearest neighbour interpolator to interpolate the gas temperature onto 
        a generated grid, and then computes the density-weighted temperature 
        along a given axis via a 2d histogram. 
        '''
        # Read the data:
        h, iu   = read_hdf5(file=self.file)
        pos     = h['PartType0']['Coordinates'][:] * iu['ulength'] / self.kpc
        dens    = h['PartType0']['Density'][:] * iu['udens']
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        temp    = temperature(h, iu)

        # Coordinate ranges:
        if not (xrange):
            xrange=(0, boxSize)
        if not (yrange):
            yrange=(0, boxSize)
        if not (zrange):
            zrange=(0, boxSize)

        # Rotation:
        if rotate != 0:
            rotation = Rotation.from_euler(axis, rotate, degrees=True)
            pos = rotation.apply(pos - boxSize/2)
            pos += boxSize/2

        # Generate the grid:
        X, dX = np.linspace(xrange[0], xrange[1], bins, retstep=True)
        Y, dY = np.linspace(yrange[0], yrange[1], bins, retstep=True)
        if cut:
            Z   = cut
            dZ  = 1
        else:
            Z, dZ = np.linspace(zrange[0], zrange[1], bins, retstep=True)
        XX, YY, ZZ = np.meshgrid(X, Y, Z)
        coord = np.column_stack((XX.flatten(), YY.flatten(), ZZ.flatten()))

        # Interpolate the data:
        interp_dens = NearestNDInterpolator(pos, dens)
        interp_temp = NearestNDInterpolator(pos, temp)
        interpDens = interp_dens(coord)
        interpTemp = interp_temp(coord)
        interpTemp = (np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), weights=(interpTemp * interpDens))[0]
                      / np.histogram2d(coord[:,0], coord[:,1], bins=(X, Y), weights=(interpDens))[0])
        
        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94, wspace=0, hspace=0)

        im = ax.imshow(np.log10(interpTemp), extent=(xrange[0], xrange[1], yrange[0], yrange[1]), 
                       origin='lower', vmin=self.vmin, vmax=self.vmax, cmap='afmhot')
        ax.text(0.95, 0.05, f'{time:.2f} Myr', color='k', ha='right', va='bottom', transform=ax.transAxes, 
                bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'})
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)

        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}(T \ [\mathrm{K}])$')

        # Save figure:
        figname = f'temp_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        print('  * Figure generated and saved')
        
        # Display figure:
        if self.interactive != True:
            imgcat(fig, width=self.imgcatwidth)
        else:
            print('  * Interactive display of the figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def stellar(self, axis='z', rotate=0, bins=100, xrange=None, yrange=None, zrange=None):
        # Read the data:
        h, iu   = read_hdf5(file=self.file)
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        
        # Stellar material:
        pos_disk   = h['PartType2']['Coordinates'] * iu['ulength'] / self.kpc
        #pos_bulge  = h['PartType3']['Coordinates'] * iu['ulength'] / self.kpc
        mass_disk  = np.full(len(pos_disk), h['Header'].attrs['MassTable'][2] * iu['umass'] / self.Msol)
        #mass_bulge = np.full(len(pos_bulge), h['Header'].attrs['MassTable'][3] * iu['umass'] / self.Msol)
        #pos = np.append(pos_disk, pos_bulge, axis=0)
        #mass = np.append(mass_disk, mass_bulge)
        pos = pos_disk
        mass = mass_disk


        # Coordinate ranges:
        if not (xrange):
            xrange=(0, boxSize)
        if not (yrange):
            yrange=(0, boxSize)
        bins_x, dx = np.linspace(xrange[0], xrange[1], bins, retstep=True)
        bins_y, dy = np.linspace(yrange[0], yrange[1], bins, retstep=True)
        if zrange:
            pos = pos[(pos[:,2] > np.min(zrange)) * (pos[:,2] < np.max(zrange))]
        
        # Rotation:
        if rotate != 0:
            rotation = Rotation.from_euler(axis, rotate, degrees=True)
            pos = rotation.apply(pos - boxSize/2)
            pos += boxSize/2

        # Histogram 2D:
        H, xedges, yedges = np.histogram2d(pos[:,0], pos[:,1], bins=[bins_x, bins_y])
        with np.errstate(divide = 'ignore'):
            H = np.log10(H.T / (dx * dy))
    
        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6.4))
        fig.subplots_adjust(left=0.18, right=0.82, bottom=0.14, top=0.94, wspace=0, hspace=0)

        im = ax.imshow(H, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], 
                       vmin=self.vmin, vmax=self.vmax, cmap='bone')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 
                'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', 
                transform=ax.transAxes)
        bone = mpl.colormaps['bone']
        ax.set_facecolor(bone(0))
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}$($\Sigma_\star$ [M$_\odot$ kpc$^{-2}$])')

        # Save figure:
        figname = f'stellar_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        print('  * Figure generated and saved')

        # Display figure:
        if self.interactive != True:
            imgcat(fig, width=self.imgcatwidth)
        else:
            print('  * Interactive display of the figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def darkmatter(self, axis='z', rotate=0, bins=100, xrange=None, yrange=None, zrange=None):
        # Read the data:
        h, iu   = read_hdf5(file=self.file)
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        
        # Dark matter:
        pos  = h['PartType1']['Coordinates'] * iu['ulength'] / self.kpc
        mass = np.full(len(pos), h['Header'].attrs['MassTable'][1] * iu['umass'] / self.Msol)

        # Coordinate ranges:
        if not (xrange):
            xrange=(0, boxSize)
        if not (yrange):
            yrange=(0, boxSize)
        bins_x, dx = np.linspace(xrange[0], xrange[1], bins, retstep=True)
        bins_y, dy = np.linspace(yrange[0], yrange[1], bins, retstep=True)
        if zrange:
            pos = pos[(pos[:,2] > np.min(zrange)) * (pos[:,2] < np.max(zrange))]
        
        # Rotation:
        if rotate != 0:
            rotation = Rotation.from_euler(axis, rotate, degrees=True)
            pos = rotation.apply(pos - boxSize/2)
            pos += boxSize/2

        # Histogram 2D:
        H, xedges, yedges = np.histogram2d(pos[:,0], pos[:,1], bins=[bins_x, bins_y])
        with np.errstate(divide = 'ignore'):
            H = np.log10(H.T / (dx * dy))
    
        # Plot:
        fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
        im = ax.imshow(H, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], 
                       vmin=self.vmin, vmax=self.vmax, cmap='magma')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 
                'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', 
                transform=ax.transAxes)
        bone = mpl.colormaps['magma']
        ax.set_facecolor(bone(0))
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}$($\Sigma_\star$ [M$_\odot$ kpc$^{-2}$])')

        # Save figure:
        figname = f'dm_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        print('  * Figure generated and saved')

        # Display figure:
        if self.interactive != True:
            imgcat(fig, width=self.imgcatwidth)
        else:
            print('  * Interactive display of the figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def star_formation(self, axis='z', rotate=0, bins=100, xrange=None, yrange=None, zrange=None):
        # Read the data:
        h, iu     = read_hdf5(file=self.file)
        boxSize   = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        time      = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        pos_gas   = h['PartType0']['Coordinates'] * iu['ulength'] / self.kpc
        dens_gas  = h['PartType0']['Density'] * iu['udens']
        pos_star  = h['PartType4']['Coordinates'] * iu['ulength'] / self.kpc
        mass_star = h['PartType4']['InitialMass'] * iu['umass'] / self.Msol
        time_star = h['PartType4']['StellarFormationTime'] * iu['utime'] / self.Myr

        # Coordinate ranges:
        if not (xrange):
            xrange=(0, boxSize)
        if not (yrange):
            yrange=(0, boxSize)
        if not (zrange):
            zrange=(0, boxSize)

        # Interpolate the gas density:
        interpDens = interpolate_to_2d(pos=pos_gas, unit=self.kpc, values=dens_gas, bins=160, xrange=xrange, yrange=yrange, zrange=zrange, cut=None)

        # Star formation surface density:
        mask_sf = (time - time_star < 10)
        xbins, dx = np.linspace(np.min(xrange), np.max(xrange), bins, retstep=True)
        ybins, dy = np.linspace(np.min(yrange), np.max(yrange), bins, retstep=True)
        H, xedges, yedges = np.histogram2d(pos_star[:,0][mask_sf], pos_star[:,1][mask_sf], bins=[xbins, ybins], weights=mass_star[mask_sf])
        with np.errstate(divide = 'ignore'):
            H = np.log10(H.T/(dx * dy) / 10)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
        im_gas = ax.imshow(np.log10(interpDens), vmin=self.vmin, vmax=self.vmax, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), 
                       origin='lower', cmap='Greys')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'}, 
                color='k', ha='right', va='bottom', transform=ax.transAxes)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im_gas, cax=cax, label='$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-2}])$')
        
        im_sf = ax.imshow(H, extent=(xrange[0], xrange[1], yrange[0], yrange[1]), origin='lower', cmap='plasma')
        cax = ax.inset_axes([0.02, 0.94, 0.7, 0.04])
        cb = fig.colorbar(im_sf, cax=cax, orientation='horizontal', location='bottom')
        cb.set_label(label='$\log_{10}(\Sigma_\mathrm{SFR} \ [\mathrm{M}_\odot \ \mathrm{yr}^{-1} \ \mathrm{kpc}^{-2}])$', size=12)
        cb.ax.tick_params(labelsize=12)

        # Save figure:
        figname = f'sf_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        print('  * Figure generated and saved')
        
        # Display figure:
        if self.interactive != True:
            imgcat(fig, width=self.imgcatwidth)
        else:
            print('  * Interactive display of the figure is now running')
            plt.show()

        return None

    ##########################################################################
    ##########################################################################
    def movie():

        return None


# -------------- End of file


