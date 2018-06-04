import sys
import os
#import matplotlib.pyplot as plt
#import matplotlib
import pickle
from astropy.io import fits
#from matplotlib.patches import Circle
#from mpl_toolkits.axes_grid1 import make_axes_locatable
#from matplotlib import cm

"Main"

print("Arg 1: Map_x Arg2: Map_y : Compaeres Mapx with Mapy and plots the difference")

grid_x = []
grid_y = []

if len(sys.argv) == 3 :

	print ("Opening ", sys.argv[1])
	hdulist1 = fits.open(sys.argv[1])
	#hdulist.info()
	Map1 = hdulist1[0].data

	print ("Opening ", sys.argv[2])
	hdulist2 = fits.open(sys.argv[2])
	#hdulist.info()
	Map2 = hdulist2[0].data
	
	hdulist1.close()
	hdulist2.close()

	Mapdif = Map1 - Map2

	
	drive, path_and_file = os.path.splitdrive(sys.argv[1])
	path, file1 = os.path.split(path_and_file)
	drive, path_and_file = os.path.splitdrive(sys.argv[2])
	path, file2 = os.path.split(path_and_file)

	# Map 1
	fig,ax = plt.subplots(1)
	ax.set_aspect('equal')
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)

	im = ax.imshow(Map1,norm=matplotlib.colors.LogNorm())
	fig.colorbar(im,cax=cax,orientation='vertical')
	# Map 2
	fig,ax = plt.subplots(1)
	ax.set_aspect('equal')
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)

	im = ax.imshow(Map2,norm=matplotlib.colors.LogNorm())
	fig.colorbar(im,cax=cax,orientation='vertical')
	# Map 1 - Map 2
	fig,ax = plt.subplots(1)
	ax.set_aspect('equal')
	divider = make_axes_locatable(ax)
	cax = divider.append_axes('right', size='5%', pad=0.05)

	im = ax.imshow(Mapdif)
	fig.colorbar(im,cax=cax,orientation='vertical')

	plt.show()

	#pickle.dump(fig, open('Diff_'+file1+'_' + file2 + '.fig.pickle', 'wb')) 

	