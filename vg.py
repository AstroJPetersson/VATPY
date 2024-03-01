import tkinter as tk
from PIL import Image, ImageTk

def generate_image():


    return None

root = tk.Tk()  # create a root widget
root.title('Tk Example')
root.configure(background='white')
root.minsize(800, 600)  # width, height
root.maxsize(800, 600)
root.geometry('800x600+0+0')  # width x height + x + y

# Create left and right frames
#left_frame = tk.Frame(root, width=180, height=, bg='grey')
#left_frame.grid(row=0, column=0, padx=10, pady=10)

#right_frame = tk.Frame(root, width=380, height=380, bg='grey')
#right_frame.grid(row=0, column=1, padx=10, pady=100)

# Logo
#logo_frame = tk.Frame(left_frame, width=180, height=120)
#logo_frame.grid(row=0, column=0, padx=0, pady=0)

# Logo
logo = Image.open('vatpy.png')
logo_resize = logo.resize((180, 120))
logo_img = ImageTk.PhotoImage(logo_resize)
tk.Label(root, image=logo_img, borderwidth=2, bg='black').grid(row=0, column=0, padx=10, pady=10)

# Settings
#settings = tk.Frame(left_frame, width=180, height=280)
#settings.grid(row=1, column=0, padx=0, pady=0, sticky=tk.W)

tk.Label(root, text='Settings', anchor='w').grid(row=1, column=0, padx=10, sticky=tk.W)

tk.Label(root, text='Snapshot', anchor='w').grid(row=2, column=0, padx=10, sticky=tk.W)
snap = tk.Entry(root, bd=3, width=15).grid(row=3, column=0, padx=10, sticky=tk.W)

tk.Label(root, text='Property', anchor='w').grid(row=4, column=0, padx=10, sticky=tk.W)
prop = tk.StringVar(root, 'Gas Density')
tk.Radiobutton(root, text='Gas Density', variable=prop, value='Gas Density').grid(row=5, column=0, padx=10, sticky=tk.W)
tk.Radiobutton(root, text="Gas Temperature", variable=prop, value='Gas Temperature').grid(row=5, column=1, padx=0, sticky=tk.W)

#tk.Label(root, text='xrange', anchor='w').grid(row=4, column=0, padx=10, sticky=tk.W)
#tk.Label(root, text='yrange', anchor='w').grid(row=5, column=0, padx=10, sticky=tk.W)
#tk.Label(root, text='zrange', anchor='w').grid(row=6, column=0, padx=10, sticky=tk.W)
#tk.Label(root, text='bins', anchor='w').grid(row=7, column=0, padx=10, sticky=tk.W)
#tk.Label(root, text='Rot. Axis', anchor='w').grid(row=8, column=0, padx=10, sticky=tk.W)
#tk.Label(root, text='Rot. Angle', anchor='w').grid(row=9, column=0, padx=10, sticky=tk.W)


#snap = tk.Entry(settings, bd=3, width=7)
#snap.grid(row=1, column=1, padx=0, pady=0)


# Property menu:
# Gender label and dropdown widgets
#prop = tk.Menubutton(settings, text='Property')
#prop.grid(row=2, column=1, padx=0, pady=0)
#prop.menu = tk.Menu(prop, tearoff=0)
#rop['menu'] = prop.menu

# choices in gender dropdown menu
#prop.menu.add_cascade(label='Gas Density')
#prop.menu.add_cascade(label='Gas Temperature')
#prop.menu.add_cascade(label='Stellar Density')
#prop.grid()

#prop_selection.add_cascade(label='Gas Density')
#prop_selection.add_cascade(label='Gas Temperature')


#
# Buttons
#buttons = tk.Frame(left_frame, width=180, height=40)
#buttons.grid(row=2, column=0, padx=0, pady=0)
#generate = tk.Button(buttons, text='Generate image', command=generate_image, bg="#6FAFE7", width=15)
#generate.grid(row=0, column=0)


root.mainloop()




