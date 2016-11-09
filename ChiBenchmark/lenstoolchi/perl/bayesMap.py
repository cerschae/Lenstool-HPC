#!/usr/bin/env python
# Add option --nocube to compute mean and stddev on the fly (possible NaN is too noisy)
import glob
import re
import pyfits
import numpy
import os
import sys

cube=True
if len(sys.argv) > 1 and sys.argv[1] == '--nocube':
    cube=False

names = glob.glob('tmp/*fits')
root = re.findall('tmp/([A-Za-z]+)', names[0])[0]
n=len(names)
print 'Combine '+repr(n)+' FITS images in tmp/'

a=pyfits.open(names[0])
hdr=a[0].header
nx = a[0].header['NAXIS1']
ny = a[0].header['NAXIS2']
if cube:
    cube=numpy.zeros([n, nx, ny])
    for i in range(n-1):
        print "Read file %d / %d\r"%(i+1, n),
	try:
        	cube[i,:,:]=pyfits.getdata(names[i])
	except:
		print "Error reading file " + names[i]

    avg = numpy.sum(cube, 0) / n
    for i in range(n):
            cube[i,:,:] -= avg
    
    std = numpy.sqrt(numpy.sum(cube * cube,0) / n)
else:
    avg=numpy.zeros((nx,ny))
    std=numpy.zeros((nx,ny))
    for i in range(n-1):
        print "Read file %d / %d\r"%(i+1, n),
	try:
        	img =pyfits.getdata(names[i])
	except:
		print "Error reading file " + names[i]

        avg += img
        std += img * img

    avg /= n
    std = numpy.sqrt(std / n - avg * avg)

sn = avg / std
fmean='mean'+root+'.fits'
fstd='std'+root+'.fits'
fsn='sn'+root+'.fits'
if os.path.isfile(fmean): os.unlink(fmean)
if os.path.isfile(fstd): os.unlink(fstd)
if os.path.isfile(fsn): os.unlink(fsn)
pyfits.writeto(fmean, avg, hdr)
pyfits.writeto(fstd, std, hdr)
pyfits.writeto(fsn, sn, hdr)

