import sys
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

print('\n Welcome to VATPY Image Display (by Jonathan Petersson)')

# File:
file = str(sys.argv[1])
print(f'  * Displaying {file}')

# Figure:
img = np.asarray(Image.open(file))
ratio = np.shape(img)[0] / np.shape(img)[1]
plt.figure(figsize=(6/ratio, 6), layout='constrained')
plt.imshow(img)
plt.axis('off')
plt.show()

print('  * Done!\n')
