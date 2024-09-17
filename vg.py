# -------------- Import necessary packages
import os
import numpy as np
import tkinter as tk
import threading
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,  
NavigationToolbar2Tk)
from PIL import Image, ImageTk

# -------------- Import GUI specific functions
from vatpy import get_gas_density_image, get_gas_temperature_image

# -------------- Configuration
import config
homedir = config.homedir

# -------------- Define functions
def thread_func():
    th = threading.Thread(target=generate_image) 
    th.start() 

def generate_image():
    # Get parameters:
    snap_var      = snap.get()
    ulength_var   = float(ulength.get())
    if xrange.get():
        xrange_var = xrange.get()
        xrange_var = xrange_var.replace('(', '')
        xrange_var = xrange_var.replace(')', '')
        xrange_var = tuple([float(i) for i in xrange_var.split(',')])
    else:
        xrange_var = None
    if yrange.get():
        yrange_var = yrange.get()
        yrange_var = yrange_var.replace('(', '')
        yrange_var = yrange_var.replace(')', '')
        yrange_var = tuple([float(i) for i in yrange_var.split(',')])
    else:
        yrange_var = None
    if zrange.get():
        zrange_var = zrange.get()
        zrange_var = zrange_var.replace('(', '')
        zrange_var = zrange_var.replace(')', '')
        zrange_var = tuple([float(i) for i in zrange_var.split(',')])
    else:
        zrange_var = None
    bins_var      = int(bins.get())
    phys_prop_var = phys_prop.get()
    if vmin.get():
        vmin_var = float(vmin.get())
    else:
        vmin_var = None
    if vmax.get():
        vmax_var = float(vmax.get())
    else:
        vmax_var = None
    
    # Enable textbox:
    textbox.config(state=tk.NORMAL)
    
    # Check if file exists:
    snap_nr = '000'[len(str(snap_var)):] + str(snap_var)
    file = f'{os.getcwd()}/snap_{snap_nr}.hdf5'
    if not os.path.isfile(file):
        textbox.insert(tk.END, f'> Error! No file with the name {file} detected!\n')
        textbox.config(state=tk.DISABLED)
        return 0
        
    # Gas density image:
    if phys_prop_var == 'gas_density':
        # Make sure axis is clean:
        textbox.insert(tk.END, '> Cleaning & preparing canvas for image output...')
        ax.cla()
        textbox.insert(tk.END, ' Done!\n')

        # Get the data for the image:
        textbox.insert(tk.END, '> Generating gas surface density map...')
        H, X, Y = get_gas_density_image(file=file, bins=bins_var, xrange=xrange_var, yrange=yrange_var, 
                                        zrange=zrange_var, ulength=ulength_var)
        textbox.insert(tk.END, ' Done!\n')
        
        # Plot:
        textbox.insert(tk.END, '> Now do some plotting...')
        im = ax.imshow(H, extent=(X[0], X[-1], Y[0], Y[-1]), origin='lower', vmin=vmin_var, vmax=vmax_var, 
                       cmap='magma')
        ax.set_xlabel('$x$ [i.u.]')
        ax.set_ylabel('$y$ [i.u.]')
        cax = ax.inset_axes([1, 0, 0.05, 1])
        fig.colorbar(im, cax=cax, label='$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-2}])$')
        textbox.insert(tk.END, ' Done!\n')
        
    # Gas temperature image:
    elif phys_prop_var == 'gas_temperature':
        # Make sure axis is clean:
        textbox.insert(tk.END, '> Cleaning & preparing canvas for image output...')
        ax.cla()
        textbox.insert(tk.END, ' Done!\n')

        # Get the data for the image:
        textbox.insert(tk.END, '> Generating gas temperature map...')
        H, X, Y = get_gas_temperature_image(file=file, bins=bins_var, xrange=xrange_var, yrange=yrange_var, 
                                            zrange=zrange_var, ulength=ulength_var)
        textbox.insert(tk.END, ' Done!\n')
        
        # Plot:
        textbox.insert(tk.END, '> Now do some plotting...')
        im = ax.imshow(H, extent=(X[0], X[-1], Y[0], Y[-1]), origin='lower', vmin=vmin_var, vmax=vmax_var, 
                       cmap='afmhot')
        ax.set_xlabel('$x$ [i.u.]')
        ax.set_ylabel('$y$ [i.u.]')
        cax = ax.inset_axes([1, 0, 0.05, 1])
        fig.colorbar(im, cax=cax, label='$\log_{10}(T \ [\mathrm{K}])$')
        textbox.insert(tk.END, ' Done!\n')
    
    # Error message:
    else:
        textbox.insert(tk.END, ' Error!\n> Something went wrong in the selection of physical property... please debug!')

    # Update the canvas & toolbar:
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=3, padx=10, pady=10, rowspan=14, sticky=tk.EW)
    
    toolbar.update()
    toolbar.grid(row=14, column=3, rowspan=3)
    
    # Final textbox message:
    textbox.insert(tk.END, '> Image successfully generated!\n')
    textbox.config(state=tk.DISABLED)
    textbox.see('end')

    return 1


def apply_image():
    # Get parameters:
    phys_prop_var = phys_prop.get()
    if vmin.get():
        vmin_var = float(vmin.get())
    else:
        vmin_var = None
    if vmax.get():
        vmax_var = float(vmax.get())
    else:
        vmax_var = None

    # Enable textbox:
    textbox.config(state=tk.NORMAL)
    
    # Retrive axis data:
    H = ax.get_images()[0]
    H = H.get_array()
    X = ax.get_xlim()
    Y = ax.get_ylim()

    # Decide which colormap & label to use:
    if phys_prop_var == 'gas_density':
        label = '$\log_{10}(\Sigma_\mathrm{Gas} \ [\mathrm{g} \ \mathrm{cm}^{-2}])$'
        cmap  = 'magma'
    elif phys_prop_var == 'gas_temperature':
        label = '$\log_{10}(T \ [\mathrm{K}])$'
        cmap = 'afmhot'
    else:
        textbox.insert(tk.END, ' Error!\n> Something went wrong in the selection of physical property... please debug!')

    # Make sure axis is clean:
    textbox.insert(tk.END, '> Cleaning & preparing canvas for image output...')
    ax.cla()
    textbox.insert(tk.END, ' Done!\n')

    # Update axis image with vmin/vmax:
    textbox.insert(tk.END, '> Updating with newest vmin & vmax values...')
    im = ax.imshow(H, extent=(X[0], X[-1], Y[0], Y[-1]), origin='lower', vmin=vmin_var, vmax=vmax_var, cmap=cmap)
    ax.set_xlabel('$x$ [i.u.]')
    ax.set_ylabel('$y$ [i.u.]')
    cax = ax.inset_axes([1, 0, 0.05, 1])
    fig.colorbar(im, cax=cax, label=label)
    textbox.insert(tk.END, ' Done!\n')

    # Update the canvas & toolbar:
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=3, padx=10, pady=10, rowspan=14, sticky=tk.EW)
    
    toolbar.update()
    toolbar.grid(row=14, column=3, rowspan=3)
    
    # Final textbox message:
    textbox.insert(tk.END, '> Image successfully updated!\n')
    textbox.config(state=tk.DISABLED)
    textbox.see('end')

    return 1

# -------------- Run a GUI via tkinter
# Window:
root = tk.Tk()
root.title('VATPY GUI')
root.configure(background='white')
root.minsize(800, 700)
root.maxsize(800, 700)
root.geometry('800x700+0+0')

# Logo display:
logo = Image.open(f'{homedir}/VATPY/logo/vatpy.png')
logo_resize = logo.resize((300, 200))
logo_img = ImageTk.PhotoImage(logo_resize)
tk.Label(root, image=logo_img, borderwidth=2, bg='black').grid(row=0, column=0, padx=10, pady=10, rowspan=3, columnspan=2)

# Snapshot selection:
tk.Label(root, text='Snapshot', font=('Courier', 12), bg='white', anchor='w').grid(row=5, column=0, padx=10, sticky=tk.W)
snap = tk.StringVar(root, '000')
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=snap).grid(row=5, column=1, padx=10, sticky=tk.W)

# Physical property selection:
tk.Label(root, text='Property', font=('Courier', 12), bg='white', anchor='w').grid(row=6, column=0, padx=10, sticky=tk.W)
phys_prop = tk.StringVar(root, 'gas_density')
tk.Radiobutton(root, text='Gas Density', font=('Courier', 10), bg='white', borderwidth=0, highlightthickness=0, 
               variable=phys_prop, value='gas_density').grid(row=7, column=0, padx=10, sticky=tk.W)
tk.Radiobutton(root, text="Stellar Density", font=('Courier', 10), bg='white', borderwidth=0, highlightthickness=0, 
               variable=phys_prop, value='stellar_density').grid(row=7, column=1, padx=0, sticky=tk.W)
tk.Radiobutton(root, text="Gas Temperature", font=('Courier', 10), bg='white', borderwidth=0, highlightthickness=0, 
               variable=phys_prop, value='gas_temperature').grid(row=8, column=1, padx=0, sticky=tk.W)

# x/y/z range selection:
tk.Label(root, text='Range in x', font=('Courier', 12), bg='white', anchor='w').grid(row=9, column=0, padx=10, sticky=tk.W)
xrange = tk.StringVar()
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=xrange).grid(row=9, column=1, padx=10, sticky=tk.W)

tk.Label(root, text='Range in y', font=('Courier', 12), bg='white', anchor='w').grid(row=10, column=0, padx=10, sticky=tk.W)
yrange = tk.StringVar()
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=yrange).grid(row=10, column=1, padx=10, sticky=tk.W)

tk.Label(root, text='Range in z', font=('Courier', 12), bg='white', anchor='w').grid(row=11, column=0, padx=10, sticky=tk.W)
zrange = tk.StringVar()
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=zrange).grid(row=11, column=1, padx=10, sticky=tk.W)

# Number of bins:
tk.Label(root, text='Number of Bins', font=('Courier', 12), bg='white', anchor='w').grid(row=12, column=0, padx=10, sticky=tk.W)
bins = tk.StringVar(root, '100')
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=bins).grid(row=12, column=1, padx=10, sticky=tk.W)

# Rotation:
tk.Label(root, text='Rotation Axis', font=('Courier', 12), bg='white', anchor='w').grid(row=13, column=0, padx=10, sticky=tk.W)
rot_axis = tk.StringVar(root, 'z')
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=rot_axis).grid(row=13, column=1, padx=10, sticky=tk.W)

tk.Label(root, text='Rotation Angle', font=('Courier', 12), bg='white', anchor='w').grid(row=14, column=0, padx=10, sticky=tk.W)
rot_angle = tk.StringVar(root, '0')
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=rot_angle).grid(row=14, column=1, padx=10, sticky=tk.W)

# Unit selection:
tk.Label(root, text='Unit Selection', font=('Courier', 12), bg='white', anchor='w').grid(row=15, column=0, padx=10, sticky=tk.W)
ulength = tk.StringVar(root, '3.08567758e21')
tk.Radiobutton(root, text='kpc', font=('Courier', 10), bg='white', borderwidth=0, highlightthickness=0, 
               variable=ulength, value='3.08567758e21').grid(row=16, column=0, padx=10, sticky=tk.W)
tk.Radiobutton(root, text="pc", font=('Courier', 10), bg='white', borderwidth=0, highlightthickness=0, 
               variable=ulength, value='3.08567758e18').grid(row=16, column=1, padx=0, sticky=tk.W)

# Figure object:
fig, ax = plt.subplots(figsize=(4, 3.5), layout='constrained')
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.set_aspect('equal')
ax.text(0.5, 0.5, 'Empty Canvas', ha='center', va='center')

# Canvas object:
canvas = FigureCanvasTkAgg(fig, master=root) 
canvas.get_tk_widget().grid(row=0, column=3, padx=10, pady=10, rowspan=14, sticky=tk.EW)
canvas.draw()

# Matplotlib toolbar:
toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar=False) 

# Textbox:
bot_right_frame = tk.Frame(root, width=400, height=120)
bot_right_frame.grid(row=17, column=3, padx=10, pady=10, rowspan=4, sticky=tk.N)

textbox = tk.Text(bot_right_frame, bd=2, width=60, height=7.5, font=('Courier', 8), bg='white', state=tk.DISABLED)
textbox.grid(row=0, column=0, sticky=tk.NS)

# Scrollbar widget:
scrollbar = tk.Scrollbar(bot_right_frame, orient='vertical', command=textbox.yview)
scrollbar.grid(row=0, column=1, sticky=tk.NS)
textbox['yscrollcommand'] = scrollbar.set

# Generate image:
tk.Button(root, text='Generate', font=('Courier', 12), command=thread_func).grid(row=17, column=0, padx=10, pady=5, columnspan=2, sticky=tk.EW)

# Colorbar vmin/vmax:
tk.Label(root, text='Colorbar vmin', font=('Courier', 12), bg='white', anchor='w').grid(row=18, column=0, padx=10, sticky=tk.W)
vmin = tk.StringVar()
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=vmin).grid(row=18, column=1, padx=10, sticky=tk.W)

tk.Label(root, text='Colorbar vmax', font=('Courier', 12), bg='white',  anchor='w').grid(row=19, column=0, padx=10, sticky=tk.W)
vmax = tk.StringVar()
tk.Entry(root, bd=2, width=10, font=('Courier', 10), textvariable=vmax).grid(row=19, column=1, padx=10, sticky=tk.W)

# Apply vmin/vmax:
tk.Button(root, text='Apply', font=('Courier', 12), command=apply_image).grid(row=20, column=0, padx=10, pady=5, columnspan=2, sticky=tk.EW)

# End of window:
root.mainloop()


