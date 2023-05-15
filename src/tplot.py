# -------------- Required packages
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

from .read import read_hdf5
from .frame import density_frame, temperature_frame, number_density, temperature

# -------------- TerminalPlot
class TerminalPlot:
    '''
    TerminalPlot: 
    '''
    def __init__(self, file, savepath=os.getcwd(), vmin=None, vmax=None, xlim=None, ylim=None, 
                 style=None, saveformat=None):
        # Variables:
        self.file       = file
        self.savepath   = savepath
        self.vmin       = vmin
        self.vmax       = vmax
        self.xlim       = xlim
        self.ylim       = ylim
        self.saveformat = saveformat
        
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
        Description: General information about the snapshot.
        '''
        # Read the data:
        h, iu = read_hdf5(file=self.file)
        
        time    = h['Header'].attrs['Time'] * iu['utime'] / self.Myr
        boxSize = h['Header'].attrs['BoxSize'] * iu['ulength'] / self.kpc
        N0 = len(h['PartType0']['ParticleIDs'])

        print('\n########################################################')
        print(f'##### Some general information about {self.file} #####')
        print('########################################################\n')

        print(f'Time: {time} Myr')
        print(f'BoxSize: {boxSize} kpc\n')
        print(f'Number of particles')
        print(f'PartType0 (gas): {N0}')

        print('\n########################################################\n')
        
        return 0


    def density(self, bins=100, axis='z', unit='cgs', cut=None, column=None, box=None):
        '''
        Description: 2D interpolation plot of the gas density. 
        '''
        # Density frame:
        dens, cbarLabel, time, boxMin, boxMax = density_frame(file=self.file, bins=bins, axis=axis, unit=unit, 
                                                              cut=cut, column=column, box=box)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(np.log10(dens), vmin=self.vmin, vmax=self.vmax, extent=(boxMin, boxMax, boxMin, boxMax), 
                       origin='lower', cmap='BuPu')
        ax.text(0.95, 0.05, f'{round(time, 2)} Myr', bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.5, 'boxstyle': 'round'}, color='k', ha='right', va='bottom', transform=ax.transAxes)
        ax.set_aspect('equal')
        if axis == 'x':
            xlabel = '$y$ [pc]'
            ylabel = '$z$ [pc]'
        elif axis == 'y':
            xlabel = '$x$ [pc]'
            ylabel = '$z$ [pc]'
        else:
            xlabel = '$x$ [pc]'
            ylabel = '$y$ [pc]'
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        div = make_axes_locatable(ax)
        cax = div.append_axes('right', size='5%', pad=0)
        fig.colorbar(im, cax=cax, label=cbarLabel)
        fig.tight_layout()
        
        figname = f'density_{self.file[5:8]}.png'
        print(f'  * Done! Save path: {self.savepath}/{figname}')
        print(f'  * View figure:')
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)
        
        return 0


    def temperature(self, bins=100, axis='z', unit='cgs', cut=None, column=None, box=None):
        '''
        Description: 2D interpolation plot of the gas temperature.
        '''
        # Temperature frame:
        temp, time, boxMin, boxMax = temperature_frame(file=self.file, bins=bins, axis=axis, cut=cut, column=column, box=box)

        # Plot:
        fig, ax = plt.subplots(figsize=(8, 6))
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
        fig.tight_layout()
        
        figname = f'temperature_{self.file[5:8]}.png'
        print(f'  * Done! Save path: {self.savepath}/{figname}')
        print(f'  * View figure:')
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)

        return 0


    def phase_diagram(self, bins, num=False):
        '''
        Description: Phase diagram of the gas. 
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

        figname = f'phase_{self.file[5:8]}.png'
        print(f'  * Done! Save path: {self.savepath}/{figname}')
        print(f'  * View figure:')
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)

        return 0


    def resolution(self, bins=100, jeans=False):
        '''
        Description: Resolution plot of the gas cells. 
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

        # Plot:
        fig, ax = plt.subplots(2, 1, figsize=(8, 10), sharex=True)

        # Density vs radius:
        cf = ax[0].contourf(np.log10(H0.T), vmin=self.vmin, vmax=self.vmax, extent=(xedges0[0], xedges0[-1], yedges0[0], yedges0[-1]), 
                            origin='lower', cmap='cividis')
        div = make_axes_locatable(ax[0])
        cax = div.append_axes('top', size="5%", pad=0.3)
        cb = fig.colorbar(cf, cax=cax, orientation='horizontal', label='$\log_{10}(N_\mathrm{cell})$')
        cb.ax.xaxis.set_ticks_position('top')
        cb.ax.xaxis.set_label_position('top')
        ax[0].set_ylabel('$\log_{10}(r_\mathrm{cell} \ [\mathrm{pc}])$')

        # Density vs mass:
        ax[1].contourf(np.log10(H1.T), vmin=self.vmin, vmax=self.vmax, extent=(xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]), 
                       origin='lower', cmap='cividis')
        ax[1].set_ylabel('$\log_{10}(M_\mathrm{cell} \ [\mathrm{M}_\odot])$')
        ax[1].set_xlabel(r'$\log_{10}(\rho_\mathrm{gas} \ [\mathrm{g} \ \mathrm{cm}^{-3}]$')    
        ax[1].set_xlim(self.xlim)

        # Jeans length:
        if jeans == True:
            rho = np.linspace(1.0e-30, 1.0e-15, 100)
            gamma = 5/3
            mu = (1 + 4 * 0.1)
            mJ_10K = (np.pi**(5/2) / 6) * ((gamma*self.kb*10) / (self.G*mu*self.mp))**(3/2) / np.sqrt(rho)
            mJ_100K = (np.pi**(5/2) / 6) * ((gamma*self.kb*100) / (self.G*mu*self.mp))**(3/2) / np.sqrt(rho)
            lJ_10K = np.sqrt((np.pi*gamma*self.kb*10) / (self.G*mu*self.mp*rho))
            lJ_100K = np.sqrt((np.pi*gamma*self.kb*100) / (self.G*mu*self.mp*rho))

            ax[0].plot(np.log10(rho), np.log10(lJ_10K/self.pc), c='k', ls=':', lw=1, alpha=0.8, label='10 K')   
            ax[0].plot(np.log10(rho), np.log10(lJ_100K/self.pc), c='k', ls='--', lw=1, alpha=0.8, label='100 K')
            ax[0].legend(title=r'$\lambda_\mathrm{J}(T)$', loc='lower left')

        fig.tight_layout()
        figname = f'resolution_{self.file[5:8]}.png'
        print(f'  * Done! Save path: {self.savepath}/{figname}')
        print(f'  * View figure:')
        fig.savefig(f'{self.savepath}/{figname}')
        imgcat(fig)
        
        return 0


    def movie():

        return 0


# -------------- End of file


