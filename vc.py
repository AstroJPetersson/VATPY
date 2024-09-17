import cmd
import matplotlib.pyplot as plt
from imgcat import imgcat

from vatpy import read_hdf5
from vatpy import get_gas_density_image_cli, get_gas_temperature_image_cli

class VatpyCLI(cmd.Cmd):
    # ------- Cmd attributes
    intro  = '\nWelcome to VATPY (CLI version)\nAuthor: Jonathan Petersson\nLast updated: 2024-07-09\n'
    prompt = f'(vatpy) > '

    # ------- Variables
    snap     = '000'
    data     = None
    status   = 'not read'
    bins     = 100
    xrange   = None
    yrange   = None
    zrange   = None
    rot_axis = 'z'
    rot_ang  = 0
    ulength  = 'kpc'

    vmin     = None
    vmax     = None
    cmap_gas = 'viridis'

    # ------- Commands
    def get_info(self):
        if self.status == 'read':
            spaces = ' ' * 4
        else:
            spaces = ''
        info = ('Table of parameter values updated\n'
                '=======================================\n'
                '  Parameters       |  Values\n'
                '---------------------------------------\n'
               f'  snapshot         |  {self.snap} ({self.status}){spaces}\n'
               f'  bins             |  {self.bins}\n'
               f'  xrange           |  {self.xrange}\n'
               f'  yrange           |  {self.yrange}\n'
               f'  zrange           |  {self.zrange}\n'
               f'  rotation axis    |  {self.rot_axis}\n'
               f'  rotation angle   |  {self.rot_ang}\n'
               f'  vmin, vmax       |  {self.vmin}, {self.vmax}\n'
                '=======================================\n')
        print(info)

    def do_snap(self, line):
        snap_nr = '000'[len(str(line)):] + str(line)
        if (snap_nr == self.snap):
            if (self.status == 'read'):
                print('Snapshot already selected and read')
            else:
                print('Snapshot already selected but not read')
        else:
            print('Selecting a new snapshot to read')
            self.status = 'not read'
            self.data   = None
            self.snap = '000'[len(str(line)):] + str(line)
            self.get_info()

    def do_read(self, line):
        if (self.status == 'read'):
            print('Snapshot already read')
        else:
            print(f'Reading data of snap_{self.snap}.hdf5')
            self.data   = read_hdf5(f'snap_{self.snap}.hdf5')
            self.status = 'read'
            self.get_info()

    def do_bins(self, line):
        print('Updating the number of bins')
        self.bins = int(line)
        self.get_info()
    
    def do_xrange(self, line):
        print('Updating the xrange')
        self.xrange = tuple([float(i) for i in line.split(' ')])
        self.get_info()
    
    def do_yrange(self, line):
        print('Updating the yrange')
        self.yrange = tuple([float(i) for i in line.split(' ')])
        self.get_info()
    
    def do_zrange(self, line):
        print('Updating the zrange')
        self.zrange = tuple([float(i) for i in line.split(' ')])
        self.get_info()

    def do_vmin(self, line):
        print('Updating vmin')
        self.vmin = float(line)
        self.get_info()

    def do_vmax(self, line):
        print('Updating vmax')
        self.vmax = float(line)
        self.get_info()

    def do_axis(self, line):
        print('Selecting a new rotation axis')
        self.rot_axis = str(line)
        self.get_info()

    def do_rotate(self, line):
        print('Selecting a new rotation angle')
        self.rot_ang = float(line)
        self.get_info()

    def do_box(self, line):
        print('Updating the x/y/z ranges to make a box')
        self.xrange = tuple([float(i) for i in line.split(' ')])
        self.yrange = tuple([float(i) for i in line.split(' ')])
        self.zrange = tuple([float(i) for i in line.split(' ')])
        self.get_info()

    def do_plot(self, line):
        if (self.status == 'read'):
            if line == 'gas':
                H, X, Y, = get_gas_density_image_cli(h=self.data[0], iu=self.data[1], axis=self.rot_axis, 
                                                     rotate=self.rot_ang, bins=self.bins, xrange=self.xrange, 
                                                     yrange=self.yrange, zrange=self.zrange, ulength=self.ulength)
                fig, ax = plt.subplots(figsize=(6, 4.5), layout='constrained')
                im = ax.imshow(H, extent=(X[0], X[-1], Y[0], Y[-1]), origin='lower', 
                               vmin=self.vmin, vmax=self.vmax, cmap=self.cmap_gas)
                ax.set_xlabel('$x$' + f'[{self.ulength}]')
                ax.set_ylabel('$y$' + f'[{self.ulength}]')
                cax = ax.inset_axes([1, 0, 0.05, 1])
                cbar_label = '$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-2}])$'
                fig.colorbar(im, cax=cax, label=cbar_label)
                imgcat(fig)
        else:
            print('Error: please read the data before plotting!')

    def do_exit(self, line):
        return True

if __name__ == '__main__':
    VatpyCLI().cmdloop()

