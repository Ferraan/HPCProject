#!/usr/bin/env python
#
# Visualize image with numpy and matplotlib.
#
# Last revision: 24/01/2022
from __future__ import print_function, division

import os, sys, numpy as np
import matplotlib.pyplot as plt
import time

imgfile = sys.argv[1]

# Read info file
info = np.genfromtxt(imgfile+'.info',comments='#')
ch   = int(info[0]) # Number of channels
w    = int(info[1]) # Width
h    = int(info[2]) # Heigth
while True:
    # Load image data - double array
    data = np.fromfile(imgfile+'.bin').reshape(h,w,ch)

    # Convert to image (uint8)
    img = np.clip(data*(1<<8-1),0,1<<8-1).astype(np.uint8)

    # Show image
    plt.imshow(img)
    plt.show()
    time.sleep(10)
