# -------------- Required packages
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

# -------------- VATPY
class Vatpy:
	'''
	VATPY : Visualisation of Arepo in the Terminal using PYthon
	'''
	def __init__(self, f, save=os.getcwd(), vmin=None, vmax=None, xlim=None, ylim=None, style=None):
		# Variables:
		self.f    = f
		self.save = save
		self.vmin = vmin
		self.vmax = vmax
		self.xlim = xlim
		self.ylim = ylim
		
		# Constants:
		self.G    = 6.67259e-8  #[cm^3 g^-1 s^-2]
		self.mp   = 1.6726e-24  #[g] 
		self.kb   = 1.380658e-16  #[erg K^-1]	
		self.Msol = 1.9891e33  #[g]
		self.pc   = 3.08567758e18  #[cm]

		# Mpl style:
		plt.style.use(style)

	def read_hdf5(self, file):
		'''
		Reads given HDF5-file using h5py, and returns corresponding h5py-object, 
		togheter with a pyhton dictonary of the internal units (in cgs). 
		'''
		# Read hdf5-file:
		h = h5py.File(file, 'r')

		# Internal units:	
		ulength = h['Header'].attrs['UnitLength_in_cm']
		umass   = h['Header'].attrs['UnitMass_in_g']
		uvel    = h['Header'].attrs['UnitVelocity_in_cm_per_s']
		utime   = ulength/uvel
		udens   = umass/(ulength**3)
		uinterg = uvel**2
		iu = {
			'ulength' : ulength,
			'umass'   : umass,
			'uvel'    : uvel,
			'utime'   : utime,
			'udens'   : udens,
			'uinterg' : uinterg
		}
		
		return h, iu

	def number_densities(self, h, iu):
		'''
		Converts the gas density into number densities of different species.
		'''
		# Gas density:	
		rho = h['PartType0']['Density'] * iu['udens']

		# Chemical abundances:
		X = h['PartType0']['ChemicalAbundances']
		xH2, xHII, xCO = X[:,0], X[:,1], X[:,2]
		
		# Number density of hydrogen:
		xHe = 0.1
		mu = (1 + 4 * xHe)
		nH = rho / (mu * self.mp)
		
		# More number densities:
		num = {
			'HII' : xHII * nH, 
			'H2' : xH2 * nH,
			'HI' : (1 - xHII - 2*xH2) * nH,
			'CO' : xCO * nH,
			'He' : xHe * nH,
			'e'  : xHII * nH
		}

		return num

	def density_frame(self, file, unit='cgs', zslice=None, bins=100):	
		'''
		Makes an interpolated 2D frame of the gas density in the x and y plane.
		'''
		# Read data:
		h, iu = self.read_hdf5(file=file)
		Coord = h['PartType0']['Coordinates'][:]
		Dens = h['PartType0']['Density'][:] * iu['udens']
		Time = h['Header'].attrs['Time'] * iu['utime'] / (1e6 * 365.25 * 24 * 60 * 60)

		# Colorbar labels:
		label = {
			'HII' : r'$\log_{10}(\Sigma_{\mathrm{HII}} \ [\mathrm{cm}^{-2}])$',
			'H2'  : r'$\log_{10}(\Sigma_{\mathrm{H}_2} \ [\mathrm{cm}^{-2}])$',
			'HI'  : r'$\log_{10}(\Sigma_{\mathrm{HI}} \ [\mathrm{cm}^{-2}])$',
			'CO'  : r'$\log_{10}(\Sigma_{\mathrm{CO}} \ [\mathrm{cm}^{-2}])$',	
			'He'  : r'$\log_{10}(\Sigma_{\mathrm{He}} \ [\mathrm{cm}^{-2}])$',	
			'e'   : r'$\log_{10}(\Sigma_{e^{-}} \ [\mathrm{cm}^{-2}])$',
			'num' : r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{cm}^{-2}])$',
			'cgs' : r'$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-2}])$'
		}
		
		# Slice at z:
		x, dx = np.linspace(0, 40, bins, retstep=True)
		if zslice is not None:
			xx, yy, zz = np.meshgrid(x, x, zslice)
			dz = 1
		else:
			xx, yy, zz = np.meshgrid(x, x, x)
			dz = dx
		
		# Unit:
		if unit != 'cgs':
			num = self.number_densities(h, iu)
			num['num'] = np.sum(np.array(list(num.values())), axis=0)
			I = num[unit]
			clabel = label[unit]
		else:
			I = Dens
			clabel = label['cgs']
		
		# Interpolation:
		interp = griddata(Coord, I, (xx, yy, zz), method='nearest')
		interp = np.sum(interp, axis=2) * dz

		return interp.T, clabel, Time

	def density_plot(self, unit='cgs', zslice=None, bins=100):
		'''
		Interpolation plot of the gas density in the x and y plane. 
		'''
		# Density frame:
		interp, clabel, Time = self.density_frame(file=self.f, unit=unit, zslice=zslice, bins=bins)

		# Plot:
		fig, ax = plt.subplots(figsize=(6, 5))
		im = ax.imshow(np.log10(interp), vmin=self.vmin, vmax=self.vmax, extent=(0, 40, 0, 40), origin='lower', cmap='magma')
		ax.text(0.95, 0.05, f'{round(Time, 2)} Myr', color='w', ha='right', va='bottom', transform=ax.transAxes)
		ax.set_aspect('equal')
		ax.set_xlabel('$x$ [pc]')
		ax.set_ylabel('$y$ [pc]')
		ax.set_xlim(self.xlim)
		ax.set_ylim(self.ylim)
		div = make_axes_locatable(ax)
		cax = div.append_axes('right', size='5%', pad=0.05)
		fig.colorbar(im, cax=cax, label=clabel)
		
		fig.tight_layout()
		imgcat(fig, height=35)
		figname = f'density_{self.f[5:8]}.png'
		fig.savefig(f'{self.save}/{figname}')
		print(f'(\'{figname}\' saved at \'{self.save}\')')
		
		return 0

	def sink_plot(self, unit='cgs', zslice=None, bins=100):
		'''
		Interpolation plot of the gas density in the x and y plane + sink partilces
		'''
		# Density frame:
		interp, clabel, Time = self.density_frame(file=self.f, unit=unit, zslice=zslice, bins=bins)

		# Plot:
		fig, ax = plt.subplots(figsize=(6, 5))
		im = ax.imshow(np.log10(interp), vmin=self.vmin, vmax=self.vmax, extent=(0, 40, 0, 40), origin='lower', cmap='magma')
		ax.text(0.95, 0.05, f'{round(Time, 2)} Myr', color='w', ha='right', va='bottom', transform=ax.transAxes)
		ax.set_aspect('equal')
		ax.set_xlabel('$x$ [pc]')
		ax.set_ylabel('$y$ [pc]')
		ax.set_xlim(self.xlim)
		ax.set_ylim(self.ylim)
		div = make_axes_locatable(ax)
		cax = div.append_axes('right', size='5%', pad=0.05)
		fig.colorbar(im, cax=cax, label=clabel)
		
		# Sink particles:
		h, iu = self.read_hdf5(file=self.f)
		sinkcoord = h['PartType5']['Coordinates'][:]
		sinkmass = h['PartType5']['Masses'] * iu['umass']
		gascoord = h['PartType0']['Coordinates'][:]
		gasdens = h['PartType0']['Density'][:]
		gasie = h['PartType0']['InternalEnergy'][:]
		gaspress = (2/3) * gasdens * gasie
		gassoundspeed = np.sqrt((5/3)*gaspress/gasdens) * iu['uvel']
		ax.scatter(sinkcoord[:,0], sinkcoord[:,1], marker='.', c='w')
		
		# Sound speed of gas:
		bondi = np.zeros(np.shape(sinkmass))
		r = 0.1
		for i in range(0, len(sinkmass)):
			mask = np.linalg.norm(gascoord - sinkcoord[i], axis=1) < r
			c = np.median(gassoundspeed[mask])
			print(c)
			bondi[i] = 2 * self.G * sinkmass[i] / (c**2)

		print(sinkmass/self.Msol)
		print(bondi/self.pc)

		fig.tight_layout()
		imgcat(fig, height=35)
		figname = f'density_{self.f[5:8]}.png'
		fig.savefig(f'{self.save}/{figname}')
		print(f'(\'{figname}\' saved at \'{self.save}\')')
		
		return 0

	def movie(self, unit='cgs', zslice=None, bins=100):
		'''
		Generate a movie of how the gas density / number density of chemical species evolve with time. 
		'''
		# Number of snapshots:
		N = int((self.f)[5:8].lstrip('0'))
		
		# Generate artist objects:
		artist_objects = []
		fig, ax = plt.subplots(figsize=(8, 6))
		for i in range(0, N+1):
			# Create density frame:
			snap = 'snap_' + '000'[:3-len(str(i))] + str(i) + '.hdf5'
			print(f'Generating artist object for {snap}')
			interp, clabel, Time = self.density_frame(file=snap, unit=unit, zslice=zslice, bins=bins)

			# Create artist object:
			im = ax.imshow(np.log10(interp), vmin=self.vmin, vmax=self.vmax, extent=(0, 40, 0, 40), origin='lower', cmap='magma')
			txt = ax.text(0.95, 0.05, f'{round(Time, 2)} Myr', color='w', ha='right', va='bottom', transform=ax.transAxes)
			artist_objects.append([im, txt])
		
		ax.set_aspect('equal')
		ax.set_xlabel('$x$ [pc]')
		ax.set_ylabel('$y$ [pc]')
		div = make_axes_locatable(ax)
		cax = div.append_axes('right', size='5%', pad=0.05)
		fig.colorbar(im, cax=cax, label=clabel)
		fig.tight_layout()	
		
		# Create and save animation:
		saveas = f'{self.save}/movie_{self.f[5:8]}.mp4'
		ani = ArtistAnimation(fig, artist_objects, blit=True)
		ani.save(f'{saveas}')	
		print(f'Animation saved at: {saveas}')

		return 0

	def distribution_plot(self):
		'''
		UNDER CONSTRUCTION
		'''
		# Read data:
		h, iu = self.read_hdf5(file=self.f)
		rho = h['PartType0']['Density'] * iu['udens']
		ie = h['PartType0']['InternalEnergy'] * iu['uinterg']
		m = h['PartType0']['Masses'] * iu['umass']
		n, ntot = self.number_density(h, iu)		
		mu = (n['HII'] * 1 + n['HI'] * 1 + n['H2'] * 2 + n['He'] * 4 + n['CO'] * 14) / ntot
		T = (2 * mu * self.mp * ie) / (3 * self.kb)

		# Density:
		hist_rho, bins_rho = np.histogram(np.log10(rho), bins=35, weights=m)
		hist_rho = np.log10(hist_rho/np.sum(m))
		baseline_rho = int(min(hist_rho[np.isfinite(hist_rho)])) - 1
		hist_rho[~np.isfinite(hist_rho)] = baseline_rho

		hist_n, bins_n = np.histogram(np.log10(ntot), bins=35, weights=m)
		hist_n = np.log10(hist_n/np.sum(m))
		baseline_n = int(min(hist_n[np.isfinite(hist_n)])) - 1
		hist_n[~np.isfinite(hist_n)] = baseline_n

		# Internal energy:
		hist_ie, bins_ie = np.histogram(np.log10(ie), bins=35, weights=m)
		hist_ie = np.log10(hist_ie/np.sum(m))
		baseline_ie = int(min(hist_ie[np.isfinite(hist_ie)])) - 1
		hist_ie[~np.isfinite(hist_ie)] = baseline_ie

		# Temperature:	
		hist_T, bins_T = np.histogram(np.log10(T), bins=35, weights=m)
		hist_T = np.log10(hist_T/np.sum(m))
		baseline_T = int(min(hist_T[np.isfinite(hist_T)])) - 1
		hist_T[~np.isfinite(hist_T)] = baseline_T

		# Figure:
		fig, ax = plt.subplots(2, 2, figsize=(8, 8))	
		ax[0,0].stairs(hist_rho, bins_rho, baseline=baseline_rho, lw=2, fill=True)
		ax[0,0].set_xlabel(r'$\log_{10}(\rho \ [\mathrm{g} \ \mathrm{cm}^{-3}])$')
		ax[0,0].set_ylabel(r'$\log_{10}$(Norm. Mass-Weighted PDF)')
		
		ax[0,1].stairs(hist_ie, bins_ie, baseline=baseline_ie, lw=2, fill=True)
		ax[0,1].set_xlabel(r'$e \ [\mathrm{cm}^{2} \ \mathrm{s}^{-2}]$')
		
		ax[1,0].stairs(hist_n, bins_n, baseline=baseline_n, lw=2, fill=True)
		ax[1,0].set_xlabel(r'$\log_{10}(n \ [\mathrm{cm}^{-3}])$')
		ax[1,0].set_ylabel(r'$\log_{10}$(Norm. Mass-Weighted PDF)')

		ax[1,1].stairs(hist_T, bins_T, baseline=baseline_T, lw=2, fill=True)
		ax[1,1].set_xlabel(r'$\log_{10}(T \ [\mathrm{K}])$')
		ax[1,1].set_ylabel(r'$\log_{10}$(Norm. Mass-Weighted PDF)')
		
		fig.tight_layout()
		imgcat(fig, height=45)

		return 0

	def temperature_plot(self, zslice=None, bins=100):
		'''
		Temperature plot.
		'''
		# Read data:
		h, iu = self.read_hdf5(file=self.f)	
		Coord = h['PartType0']['Coordinates'][:]
		Dens = h['PartType0']['Density'] * iu['udens']
		Mass = h['PartType0']['Masses'][:]
		IntErg = h['PartType0']['InternalEnergy'] * iu['uinterg']
		Time = h['Header'].attrs['Time'] * iu['utime'] / (1e6 * 365.25 * 24 * 60 * 60)
		num = self.number_densities(h, iu)
		ntot = np.sum(np.array(list(num.values())), axis=0)
		mu = (num['HII'] * 1 + num['HI'] * 1 + num['H2'] * 2 + num['He'] * 4 + num['CO'] * 14) / ntot
		Temp = (2 * mu * self.mp * IntErg) / (3 * self.kb)
		
		# Interpolation:	
		x, dx = np.linspace(0, 40, bins, retstep=True)
		if zslice is not None:
			xx, yy, zz = np.meshgrid(x, x, zslice)
		else:
			xx, yy, zz = np.meshgrid(x, x, x)
		interp_temp = griddata(Coord, Temp, (xx, yy, zz), method='nearest')
		interp_dens = griddata(Coord, Dens, (xx, yy, zz), method='nearest')
		interp = np.sum(interp_temp * interp_dens, axis=2)/np.sum(interp_dens, axis=2)

		# Plot:
		fig, ax = plt.subplots(figsize=(6, 5))
		im = ax.imshow(np.log10(interp.T), vmin=self.vmin, vmax=self.vmax, extent=(0, 40, 0, 40), origin='lower', cmap='afmhot')
		ax.text(0.95, 0.05, f'{round(Time, 2)} Myr', color='k', ha='right', va='bottom', transform=ax.transAxes)
		ax.set_aspect('equal')
		ax.set_xlabel('$x$ [pc]')
		ax.set_ylabel('$y$ [pc]')
		div = make_axes_locatable(ax)
		cax = div.append_axes('right', size='5%', pad=0.05)
		fig.colorbar(im, cax=cax, label='$\log_{10}(T \ [\mathrm{K}])$')
		fig.tight_layout()
		fig.savefig(f'temp_{self.f[5:8]}.png')
		imgcat(fig, height=35)

		return 0

	def phase_plot(self):
		'''
		UNDER CONSTRUCTION
		'''
		# Read data:
		h, iu = self.read_hdf5(file=self.f)
		Dens = h['PartType0']['Density'] * iu['udens']
		Mass = h['PartType0']['Masses']
		IntErg = h['PartType0']['InternalEnergy'] * iu['uinterg']
		Time = h['Header'].attrs['Time'] * iu['utime'] / (1e6 * 365.25 * 24 * 60 * 60)
		m = h['PartType0']['Masses'] * iu['umass']
		num = self.number_densities(h, iu)
		ntot = np.sum(np.array(list(num.values())), axis=0)
		mu = (num['HII'] * 1 + num['HI'] * 1 + num['H2'] * 2 + num['He'] * 4 + num['CO'] * 14) / ntot
		T = (2 * mu * self.mp * IntErg) / (3 * self.kb)

		H, xedges, yedges = np.histogram2d(np.log10(ntot), np.log10(T), bins=80, weights=Mass)

		# Plot:
		fig, ax = plt.subplots(figsize=(6, 6))
		im = ax.imshow(np.log10(H.T), vmin=self.vmin, vmax=self.vmax, extent=(xedges[0], xedges[-1], yedges[0], yedges[-1]), origin='lower', cmap='viridis')
		ax.text(0.05, 0.95, f'{round(Time, 2)} Myr', color='k', ha='left', va='top', transform=ax.transAxes)
		ax.set_xlabel('$\log_{10}(n \ [\mathrm{cm}^{-3}])$')
		ax.set_ylabel('$\log_{10}(T \ [\mathrm{K}])$')
		ax.set_xlim(self.xlim)
		ax.set_ylim(self.ylim)
		div = make_axes_locatable(ax)
		cax = div.append_axes('right', size='5%', pad=0.05)
		fig.colorbar(im, cax=cax, label='$\log_{10}(M \ [\mathrm{M}_\odot])$')
		
		# Jeans length:
		rho = np.linspace(1.0e-10, 1.0e10, 500)
		gamma = 5/3
		mu = (1 + 4 * 0.1)
		T_lJ10 = ((10*self.pc)**2 * self.G * mu * self.mp * (rho*mu*self.mp)) / (np.pi * gamma * self.kb)
		T_lJ1 = ((1*self.pc)**2 * self.G * mu * self.mp * (rho*mu*self.mp)) / (np.pi * gamma * self.kb)

		ax.plot(np.log10(rho), np.log10(T_lJ10), c='k', ls='--', lw=1, alpha=0.8, label='10 pc')	
		ax.plot(np.log10(rho), np.log10(T_lJ1), c='k', ls=':', lw=1, alpha=0.8, label='1 pc')	
		ax.legend(frameon=False, loc='lower left')

		fig.tight_layout()
		imgcat(fig, height=50)
		figname = f'phasediagram_{self.f[5:8]}.png'
		fig.savefig(f'{self.save}/{figname}')
		print(f'(figure \'{figname}\' saved at \'{self.save}\')')

		return 0

	def convert_numdens(self, n):
		h, iu = self.read_hdf5(file=self.f)	
		xHe = 0.1
		mu = (1 + 4 * xHe)
		n_in_iu = mu * self.mp * n / iu['udens']
		print(f'A number density of {n} is equivalent to [internal units]: {n_in_iu}')

		return 0
	
	def convert_cgsdens(self, rho):
		h, iu = self.read_hdf5(file=self.f)
		rho_in_iu = rho / iu['udens']
		print(f'A density of {rho} in cgs is equivalent to [internal units]: {rho_in_iu}')

		return 0

	def cell_mass_size_relation(self):
		'''
		UNDER CONSTRUCTION
		'''
		# Read data:
		h, iu = self.read_hdf5(file=self.f)
		Mass = h['PartType0']['Masses'] * iu['umass']
		Dens = h['PartType0']['Density'] * iu['udens']
		Radius = ((3*Mass) / (4*np.pi*Dens))**(1/3)

		# 2D Histograms:
		H0, xedges0, yedges0 = np.histogram2d(np.log10(Dens), np.log10(Radius/self.pc), bins=100)
		X0, Y0 = np.meshgrid((xedges0[0:-1] + xedges0[1:])/2, (yedges0[0:-1] + yedges0[1:])/2)
		H1, xedges1, yedges1 = np.histogram2d(np.log10(Dens), np.log10(Mass/self.Msol), bins=100)
		X1, Y1 = np.meshgrid((xedges1[0:-1] + xedges1[1:])/2, (yedges1[0:-1] + yedges1[1:])/2)

		# Plot:
		fig, ax = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
		cf = ax[0].contourf(np.log10(H0.T), vmin=self.vmin, vmax=self.vmax, extent=(xedges0[0], xedges0[-1], yedges0[0], yedges0[-1]), origin='lower', cmap='YlGnBu')
		div = make_axes_locatable(ax[0])
		cax = div.append_axes('top', size="4%", pad=0)
		cb = fig.colorbar(cf, cax=cax, orientation='horizontal', label='$\log_{10}(N_\mathrm{cell})$')
		cb.ax.xaxis.set_ticks_position('top')
		cb.ax.xaxis.set_label_position('top')
		ax[0].set_ylabel('$\log_{10}(r_\mathrm{cell} \ [\mathrm{pc}])$')
		ax[0].set_ylim(-1, 0.5)
		ax[1].contourf(np.log10(H1.T), vmin=self.vmin, vmax=self.vmax, extent=(xedges1[0], xedges1[-1], yedges1[0], yedges1[-1]), 
					   origin='lower', cmap='YlGnBu')
		ax[1].set_ylabel('$\log_{10}(M_\mathrm{cell} \ [\mathrm{M}_\odot])$')
		ax[1].set_xlabel(r'$\log_{10}(\rho_\mathrm{gas} \ [\mathrm{g} \ \mathrm{cm}^{-3}]$')	
		ax[1].set_xlim(self.xlim)
		ax[1].set_ylim(self.ylim)

		# Jeans length:
		rho = np.linspace(1.0e-30, 1.0e-15, 100)
		gamma = 5/3
		mu = (1 + 4 * 0.1)
		mJ_10K = (np.pi**(5/2) / 6) * ((gamma*self.kb*10) / (self.G*mu*self.mp))**(3/2) / np.sqrt(rho)
		mJ_100K = (np.pi**(5/2) / 6) * ((gamma*self.kb*100) / (self.G*mu*self.mp))**(3/2) / np.sqrt(rho)
		lJ_10K = np.sqrt((np.pi*gamma*self.kb*10) / (self.G*mu*self.mp*rho))
		lJ_100K = np.sqrt((np.pi*gamma*self.kb*100) / (self.G*mu*self.mp*rho))

		ax[0].plot(np.log10(rho), np.log10(lJ_10K/self.pc), c='k', ls=':', lw=1, alpha=0.8, label='10 K')	
		ax[0].plot(np.log10(rho), np.log10(lJ_100K/self.pc), c='k', ls='--', lw=1, alpha=0.8, label='100 K')
		ax[0].legend(title=r'$\lambda_\mathrm{J}(T)$', frameon=False, loc='lower left')

		fig.tight_layout()
		imgcat(fig, height=45)
		figname = f'cellMassSize_{self.f[5:8]}.png'
		fig.savefig(f'{self.save}/{figname}')
		print(f'(figure \'{figname}\' saved at \'{self.save}\')')
		
		return 0

	def experiment(self):
		h, iu = self.read_hdf5(file=self.f)
		print(len(h['PartType0']['Masses'][:]))

		return 0

# -------------- End of file


