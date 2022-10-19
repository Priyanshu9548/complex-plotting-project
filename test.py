from __future__ import division
from numpy import arange, pi, cos, sin, tan
import numpy as np
import colorsys
import sys

np.set_printoptions(threshold=sys.maxsize)

resolution = 400 # Number of values of phi (and of theta) plotted. 400 is fast, 1000 is pretty

polarPlane = False # When True plane is drawn with more pixels closer to 0; when False, pixels are evenly distributed with cartesian coordinates.
trimMag = 5 # Disc extends to |z| = trimMag, or square to Re(z), Im(z) = +/-trimMag if polarPlane is False

phi, theta = np.mgrid[0:pi:(resolution+1)*1j, -pi:pi:(resolution+1)*1j] #complex step number just makes it inclusive

# Do stereographic projection. Theta is measured counter-clockwise from the negative real direction (-pi), phi down from the Northern half of the vertical axis (so pi/2 is the equator, pi the S pole)
z = 2*cos(0.5*(pi-phi))*np.exp(1j*theta)

# Define cartesian cordinates for spherical plots
cartX = sin(phi) * cos(theta)
cartY = sin(phi) * sin(theta)
cartZ = sin(phi)+1

# Define coordinates for planar plots
if polarPlane:
    # Define section of the polar grid to project onto planar plots
    trimPhi = np.arctan(trimMag/2)*(-2)+pi
    trimmedZ = z[int(len(z) * trimPhi/pi):]
    trimmedZeros = np.zeros_like(np.real(trimmedZ))
    trimmedTheta = theta[int(len(theta) * trimPhi/pi):]
else:
    reZ, imZ = np.mgrid[-trimMag:trimMag:(resolution+1)*1j, -trimMag:trimMag:(resolution+1)*1j] #complex step number just makes it inclusive
    trimmedZ = reZ + 1j*imZ
    trimmedZeros = np.zeros_like(np.real(trimmedZ))
    trimmedTheta = np.angle(trimmedZ)

    
# Define color look up table for contours
magLut = np.zeros((256,4), dtype='uint8')
magLut[:, 3] = np.linspace(255, 0, 256) # set alpha channel
#white contours:
#for i in range(3):
#    magLut[:,i] = np.linspace(255, 255, 256)

# Define color look up table for argument
argLut = np.zeros((256,4), dtype='uint8')

# set alpha channel to 255
argLut[:, 3] = np.linspace(255, 255, 256)
# color by hue: red - green - blue - not green (adj for angles from -pi to pi)
for i in range(128,192):
    rgb = colorsys.hls_to_rgb((i-128)/64/3, 0.5, 1.0)
    for j in range(3):
        argLut[i,j] = int(rgb[j]*255)
for i in range(192,256):
    rgb = colorsys.hls_to_rgb((1/3) + (1/3)*(i-192)/64, 0.5, 1.0)
    for j in range(3):
        argLut[i,j] = int(rgb[j]*255)
for i in range(0,64):
    rgb = colorsys.hls_to_rgb((2/3) + (1/6)*(i-0)/64, 0.5, 1.0)
    for j in range(3):
        argLut[i,j] = int(rgb[j]*255)
for i in range(64,128):
    rgb = colorsys.hls_to_rgb((5/6) + (1/6)*(i-64)/64, 0.5, 1.0)
    for j in range(3):
        argLut[i,j] = int(rgb[j]*255)

# print(f'Magnitude\n{magLut}\n')
# print(f'Argument\n{argLut}\n')

print("\n\nCase - 2\n\n")

argLut = np.zeros((256,4), dtype='uint8')
argLut[:, 3] = np.linspace(255, 0, 256) # set alpha channel
#white contours:
#for i in range(3):
#    magLut[:,i] = np.linspace(255, 255, 256)

# Define color look up table for argument
magLut = np.zeros((256,4), dtype='uint8')
# set alpha channel to 255
magLut[:, 3] = np.linspace(255, 255, 256)
# color by hue: red - green - blue - not green (adj for angles from -pi to pi)
for i in range(128,192):
    rgb = colorsys.hls_to_rgb((i-128)/64/3, 0.5, 1.0)
    for j in range(3):
        magLut[i,j] = int(rgb[j]*255)
for i in range(192,256):
    rgb = colorsys.hls_to_rgb((1/3) + (1/3)*(i-192)/64, 0.5, 1.0)
    for j in range(3):
        magLut[i,j] = int(rgb[j]*255)
for i in range(0,64):
    rgb = colorsys.hls_to_rgb((2/3) + (1/6)*(i-0)/64, 0.5, 1.0)
    for j in range(3):
        magLut[i,j] = int(rgb[j]*255)
for i in range(64,128):
    rgb = colorsys.hls_to_rgb((5/6) + (1/6)*(i-64)/64, 0.5, 1.0)
    for j in range(3):
        magLut[i,j] = int(rgb[j]*255)
print(f'Magnitude\n{magLut}\n')
# print(f'Argument\n{argLut}\n')
