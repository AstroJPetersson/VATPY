# -------------- Required packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import argparse
import h5py
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import ArtistAnimation
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from imgcat import imgcat
from labellines import labelLines

from .read import read_hdf5
from .frame import density_frame, temperature_frame, number_density, temperature

# -------------- TerminalPlot
class TerminalPlot:
    '''
    TerminalPlot: 
    '''
    def __init__(self, file, savepath=os.getcwd(), vmin=None, vmax=None, xlim=None, ylim=None, 
                 style=None, saveformat=None, imgcatoff=False, imgcatw=None, imgcath=None):
        # Variables:
        self.file       = file
        self.savepath   = savepath
        self.vmin       = vmin
        self.vmax       = vmax
        self.xlim       = xlim
        self.ylim       = ylim
        self.saveformat = saveformat
        self.imgcatoff  = imgcatoff
        self.imgcatw    = imgcatw
        self.imgcath    = imgcath

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

    # ------------------------------------------------- #
    def info(self):
        '''
        Description: Information about the snapshot.
        '''
        # Read the data:
        h, iu = read_hdf5(file=self.file)
        
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        numPart = h['Header'].attrs['NumPart_ThisFile']

        print(f'  * Snapshot information')
        print(f'  | Time: {np.round(time, 2)} Myr')
        print(f'  | BoxSize: {boxSize} kpc')
        print(f'  |')
        print(f'  | Number of particles')
        print(f'  | PartType0 (gas):   {numPart[0]}')
        print(f'  | PartType1 (halo):  {numPart[1]}')
        print(f'  | PartType2 (disk):  {numPart[2]}')
        print(f'  | PartType3 (bulge): {numPart[3]}')
        print(f'  | PartType4 (stars): {numPart[4]}')
        print(f'  | PartType5 (bndry): {numPart[5]}')
        
        return 0


    def density(self, bins=100, axis='z', unit='cgs', cut=None, column=None, box=None):
        '''
        Description: Gas density plot. 
        '''
        # Density frame:
        dens, cbarLabel, time, boxMin, boxMax = density_frame(file=self.file, bins=bins, axis=axis, unit=unit, 
                                                              cut=cut, column=column, box=box)

        # Plot:
        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(np.log10(dens), vmin=self.vmin, vmax=self.vmax, extent=(boxMin, boxMax, boxMin, boxMax), 
                       origin='lower', cmap='BuPu')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', transform=ax.transAxes)
        ax.set_aspect('equal')
        if axis == 'x':
            xlabel = '$y$ [kpc]'
            ylabel = '$z$ [kpc]'
        elif axis == 'y':
            xlabel = '$x$ [kpc]'
            ylabel = '$z$ [kpc]'
        else:
            xlabel = '$x$ [kpc]'
            ylabel = '$y$ [kpc]'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=cbarLabel)
        
        fig.subplots_adjust(left=0.07, bottom=0.1, right=0.9, top=0.93, wspace=0, hspace=0)
        figname = f'density_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        if self.imgcatoff == False:
            imgcat(fig, width=self.imgcatw, height=self.imgcath)
        
        return 0


    def temperature(self, bins=100, axis='z', unit='cgs', cut=None, column=None, box=None):
        '''
        Description: Gas temperature plot.
        '''
        # Temperature frame:
        temp, time, boxMin, boxMax = temperature_frame(file=self.file, bins=bins, axis=axis, cut=cut, column=column, box=box)

        # Plot:
        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(np.log10(temp), vmin=self.vmin, vmax=self.vmax, extent=(boxMin, boxMax, boxMin, boxMax), 
                       origin='lower', cmap='hot')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', transform=ax.transAxes)
        ax.set_aspect('equal')
        ax.set_xlabel('$x$ [kpc]')
        ax.set_ylabel('$y$ [kpc]')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}(T \ [\mathrm{K}])$')
        
        fig.subplots_adjust(left=0.07, bottom=0.1, right=0.9, top=0.93, wspace=0, hspace=0)
        figname = f'temperature_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)

        return 0


    def phase(self, bins=100, num=False):
        '''
        Description: Gas phase diagram. 
        '''
        # Read the data:
        h, iu   = read_hdf5(file=self.file)
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        mass    = h['PartType0']['Masses'] * iu['umass'] / self.Msol
        dens    = h['PartType0']['Density'] * iu['udens']
        numDens = number_density(h, iu)
        temp    = temperature(h, iu)
        
        if num == True:
            dens   = numDens['Total']
            xlabel = '$\log_{10}(n \ [\mathrm{cm}^{-3}])$'
        else:
            xlabel = r'$\log_{10}(\rho \ [\mathrm{g} \ \mathrm{cm}^{-3}])$'
        H, xedges, yedges = np.histogram2d(np.log10(dens), np.log10(temp), bins=bins, weights=mass)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(np.log10(H.T), vmin=self.vmin, vmax=self.vmax, extent=(xedges[0], xedges[-1], yedges[0], yedges[-1]), origin='lower', cmap='PuBuGn')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', transform=ax.transAxes)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('$\log_{10}(T \ [\mathrm{K}])$')
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label='$\log_{10}(M \ [\mathrm{M}_\odot])$')
        fig.tight_layout()

        figname = f'phase_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)

        return 0


    def resolution(self, bins=100, smooth=0, nlvls=5):
        '''
        Description: Resolution plot. 
        '''
        # Read the data:
        h, iu = read_hdf5(file=self.file)
        mass = h['PartType0']['Masses'] * iu['umass']
        dens = h['PartType0']['Density'] * iu['udens']
        radius = ((3*mass) / (4*np.pi*dens))**(1/3)

        # 2D Histograms:
        H0, xedges0, yedges0 = np.histogram2d(np.log10(dens), np.log10(radius / self.pc), bins=bins)
        X0, Y0 = np.meshgrid((xedges0[0:-1] + xedges0[1:])/2, (yedges0[0:-1] + yedges0[1:])/2)
        H1, xedges1, yedges1 = np.histogram2d(np.log10(dens), np.log10(mass / self.Msol), bins=bins)
        X1, Y1 = np.meshgrid((xedges1[0:-1] + xedges1[1:])/2, (yedges1[0:-1] + yedges1[1:])/2)
        
        # Gaussian filter:
        if smooth > 0:
            H0 = gaussian_filter(H0, sigma=smooth)
            H1 = gaussian_filter(H1, sigma=smooth)
        
        # Log scale:
        with np.errstate(divide = 'ignore'):
            H0 = np.log10(H0.T)
            H1 = np.log10(H1.T)

        # Levels:
        lvlmax = np.max([np.max(H0), np.max(H1)])
        levels = np.linspace(0, np.ceil(lvlmax), nlvls)

        # Plot:
        fig, ax = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

        # Density vs radius:
        cf = ax[0].contourf(H0, vmin=self.vmin, vmax=self.vmax, levels=levels, extent=(xedges0[0], xedges0[-1], yedges0[0], yedges0[-1]), 
                            origin='lower', cmap='Greens')
        div = make_axes_locatable(ax[0])
        cax = div.append_axes('top', size="5%", pad=0)
        cb = fig.colorbar(cf, cax=cax, orientation='horizontal', label='$\log_{10}(N)$')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.xaxis.set_label_position('top')
        cb.set_ticks(np.arange(0, np.ceil(lvlmax)+1, 1))
        ax[0].set_ylabel('$\log_{10}(r_\mathrm{cell} \ [\mathrm{pc}])$')

        # Density vs mass:
        ax[1].contourf(H1, vmin=self.vmin, vmax=self.vmax, levels=levels, extent=(xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]), 
                       origin='lower', cmap='Greens')
        ax[1].set_ylabel('$\log_{10}(M_\mathrm{cell} \ [\mathrm{M}_\odot])$')
        ax[1].set_xlabel(r'$\log_{10}(\rho_\mathrm{gas} \ [\mathrm{g} \ \mathrm{cm}^{-3}]$')    
        ax[1].set_xlim(self.xlim)

        # Jeans length:
        ax[0].autoscale(False)
        xlim = ax[0].get_xlim()
        xlim = np.array(xlim)

        mu = (1 + 4 * 0.1)
        jeans10  = np.sqrt((15 * self.kb * 10) / (4 * np.pi * self.G * mu * self.mp * 10**(xlim)))
        jeans100 = np.sqrt((15 * self.kb * 100) / (4 * np.pi * self.G * mu * self.mp * 10**(xlim)))

        ax[0].plot(xlim, np.log10(jeans10/self.pc), c='k', ls='--', lw=1, alpha=0.6, label='$\Lambda_J$($T=10$ K)')   
        ax[0].plot(xlim, np.log10(jeans100/self.pc), c='k', ls='--', lw=1, alpha=0.6, label='$\Lambda_J$($T=100$ K)')
        xval = -24
        labelLines(ax[0].get_lines(), xvals=[xval, xval], ha='right', fontsize=14)

        fig.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, wspace=0, hspace=0.05)
        figname = f'resolution_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)
        
        return 0

    
    def stellar(self, bins=100, axis='z', box=None):
        # Read the data:
        h, iu   = read_hdf5(file=self.file)
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        
        # Stellar material:
        pos_disk   = h['PartType2']['Coordinates'] * iu['ulength'] / self.kpc
        pos_bulge  = h['PartType3']['Coordinates'] * iu['ulength'] / self.kpc
        mass_disk  = np.full(len(pos_disk), h['Header'].attrs['MassTable'][2] * iu['umass'] / self.Msol)
        mass_bulge = np.full(len(pos_bulge), h['Header'].attrs['MassTable'][3] * iu['umass'] / self.Msol)
        pos = np.append(pos_disk, pos_bulge, axis=0)
        mass = np.append(mass_disk, mass_bulge)

        # Histogram 2D:
        if box:
            binArray, binSize = np.linspace(np.min(box), np.max(box), bins, retstep=True)
        else:
            binArray, binSize = np.linspace(0, boxSize, bins, retstep=True)

        if axis == 'x':
            H, xedges, yedges = np.histogram2d(pos[:,1], pos[:,2], bins=binArray)
            xlabel = '$y$ [kpc]'
            ylabel = '$z$ [kpc]'
        elif axis == 'y':
            H, xedges, yedges = np.histogram2d(pos[:,0], pos[:,2], bins=binArray)
            xlabel = '$x$ [kpc]'
            ylabel = '$z$ [kpc]'
        else:
            H, xedges, yedges = np.histogram2d(pos[:,0], pos[:,1], bins=binArray)
            xlabel = '$x$ [kpc]'
            ylabel = '$y$ [kpc]'
        with np.errstate(divide = 'ignore'):
            H = np.log10(H.T / (binSize * binSize))
    
        # Plot:
        bone = mpl.colormaps['bone']

        fig, ax = plt.subplots(figsize=(7, 6))
        im = ax.imshow(H, origin='lower', extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]], vmin=self.vmin, vmax=self.vmax, cmap='bone')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', transform=ax.transAxes)
        ax.set_facecolor(bone(0))
        ax.set_aspect('equal')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=r'$\log_{10}$($\Sigma_\star$ [M$_\odot$ kpc$^{-2}$])')

        fig.subplots_adjust(left=0.07, bottom=0.1, right=0.9, top=0.93, wspace=0, hspace=0)
        figname = f'stellar_{self.file[-8:-5]}.png'
        fig.savefig(f'{self.savepath}/{figname}')
        if self.imgcatoff == False:
            imgcat(fig, width=self.imgcatw, height=self.imgcath)

        return 0


    def movie():

        return 0


# -------------- End of file


