# -*- coding: utf-8 -*-
# author: Juho-Petteri Lesonen, 17.11.2017

import pandas as pd
from astropy.io import fits
from astropy import wcs
import os
import numpy as np

# Import the catalogue CSV
df = pd.read_csv('karsittu_tiedot_all_lista_csv.csv')

# Create the lists for each column in the catalogue
# Identifications and imaging info
ID = df['specobjid']
run = df['run']
rerun = df['rerun']
camcol = df['camcol']
field = df['field']
obj = df['obj']

# Object coordinates
RA = df['Column2']
DEC = df['Column3']

# Petrosian90% radius
R90_u = df['petroR90_u']
R90_g = df['petroR90_g']
R90_r = df['petroR90_r']
R90_i = df['petroR90_i']
R90_z = df['petroR90_z']

pixel_size = 0.396

# Create empty lists for the whole catalogue and for the URL-address catalogue
image_directories = []
object_directories = []

# Loop for the lists
for x in range(len(ID)):
	image_directories.append([ID[x], run[x], rerun[x], camcol[x], field[x], obj[x], RA[x], DEC[x], 
	R90_u[x], R90_g[x], R90_r[x], R90_i[x], R90_z[x]])
	object_directories.append([run[x], rerun[x], camcol[x], field[x]])

# Find the every unique combination of run, camcol and field parameters to create the download url
uniq_object_directories = list(map(list, set(map(tuple, object_directories))))

# Find the indices in the total catalogue for each unique run-camcol-field url-addresses
for ii in range(1):
	indices = [i for i, x in enumerate(object_directories) if x == uniq_object_directories[ii]]
	
	run0 = uniq_object_directories[ii][0]
	rerun0 = uniq_object_directories[ii][1]
	camcol0 = uniq_object_directories[ii][2]
	field0 = uniq_object_directories[ii][3]
	
	u_url = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%s/%s/frame-u-%06d-%s-%04d.fits.bz2' % (run0, camcol0, run0, camcol0, field0)
	g_url = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%s/%s/frame-g-%06d-%s-%04d.fits.bz2' % (run0, camcol0, run0, camcol0, field0)
	r_url = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%s/%s/frame-r-%06d-%s-%04d.fits.bz2' % (run0, camcol0, run0, camcol0, field0)
	i_url = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%s/%s/frame-i-%06d-%s-%04d.fits.bz2' % (run0, camcol0, run0, camcol0, field0)
	z_url = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%s/%s/frame-z-%06d-%s-%04d.fits.bz2' % (run0, camcol0, run0, camcol0, field0)
	
	
	# Download the images using wget 
	os.system('wget %s' % u_url)
	os.system('wget %s' % g_url)
	os.system('wget %s' % r_url)
	os.system('wget %s' % i_url)
	os.system('wget %s' % z_url)
	
	# Unzip downloaded files and remove compressed items
	os.system('bzip2 -dk *.bz2')
	os.system('rm -r *.bz2')
	
	# Define files, headers and word-coordinate-system calibration
	
	u_file = 'frame-u-%06d-%s-%04d.fits' % (run0, camcol0, field0)
	g_file = 'frame-g-%06d-%s-%04d.fits' % (run0, camcol0, field0)
	r_file = 'frame-r-%06d-%s-%04d.fits' % (run0, camcol0, field0)
	i_file = 'frame-i-%06d-%s-%04d.fits' % (run0, camcol0, field0)
	z_file = 'frame-z-%06d-%s-%04d.fits' % (run0, camcol0, field0)
	
	image_header_u = fits.getheader(u_file)
	image_header_g = fits.getheader(g_file)
	image_header_r = fits.getheader(r_file)
	image_header_i = fits.getheader(i_file)
	image_header_z = fits.getheader(z_file)
	
	image_u_wcs = wcs.WCS(image_header_u)
	image_g_wcs = wcs.WCS(image_header_g)
	image_r_wcs = wcs.WCS(image_header_r)
	image_i_wcs = wcs.WCS(image_header_i)
	image_z_wcs = wcs.WCS(image_header_z)

	orig_fits_u = fits.open(u_file)
	orig_fits_g = fits.open(g_file)
	orig_fits_r = fits.open(r_file)
	orig_fits_i = fits.open(i_file)
	orig_fits_z = fits.open(z_file)
	
	orig_image_u = orig_fits_u[0].data
	orig_image_g = orig_fits_g[0].data
	orig_image_r = orig_fits_r[0].data
	orig_image_i = orig_fits_i[0].data
	orig_image_z = orig_fits_z[0].data
	
	
	# Image scaled according to the Petrosian 90% radius.
	imagescale = 2.0
	
	# Loop for each object within the same field-image
	for iii in range(len(indices)):
		num = indices[iii]
		ra0,dec0 = [RA[num], DEC[num]]
		
		# Determine largest radius from each filter for the image
		size = int(max([R90_u[num], R90_g[num], R90_r[num], R90_i[num], R90_z[num]]) / pixel_size)
		
		# Create a "background" image which we use to make the images square shaped
		bg_dims = (int(2*imagescale*size),int(2*imagescale*size))
		bg_image = np.zeros(bg_dims)
		
		# Determine coordinates for the center of the galaxy in each filter
		(xpix_u,ypix_u), = image_u_wcs.wcs_world2pix([[ra0, dec0]],0)
		(xpix_g,ypix_g), = image_g_wcs.wcs_world2pix([[ra0, dec0]],0)
		(xpix_r,ypix_r), = image_r_wcs.wcs_world2pix([[ra0, dec0]],0)
		(xpix_i,ypix_i), = image_i_wcs.wcs_world2pix([[ra0, dec0]],0)
		(xpix_z,ypix_z), = image_z_wcs.wcs_world2pix([[ra0, dec0]],0)
		
		
		# Cropping limits for each filter
		sizex1r = int(xpix_r - imagescale * size)
		sizex2r = int(xpix_r + imagescale * size)
		sizey1r = int(ypix_r - imagescale * size)
		sizey2r = int(ypix_r + imagescale * size)
		
		sizex1u = int(xpix_u - imagescale * size)
		sizex2u = int(xpix_u + imagescale * size)
		sizey1u = int(ypix_u - imagescale * size)
		sizey2u = int(ypix_u + imagescale * size)
		
		sizex1g = int(xpix_g - imagescale * size)
		sizex2g = int(xpix_g + imagescale * size)
		sizey1g = int(ypix_g - imagescale * size)
		sizey2g = int(ypix_g + imagescale * size)
		
		sizex1i = int(xpix_i - imagescale * size)
		sizex2i = int(xpix_i + imagescale * size)
		sizey1i = int(ypix_i - imagescale * size)
		sizey2i = int(ypix_i + imagescale * size)
		
		sizex1z = int(xpix_z - imagescale * size)
		sizex2z = int(xpix_z + imagescale * size)
		sizey1z = int(ypix_z - imagescale * size)
		sizey2z = int(ypix_z + imagescale * size)
		
		# Check if the limits go beyond the border of the field and limit if neccesary
		if sizex1r <= 0:
			sizex1r = 0
		if sizex2r >= 2048:
			sizex2r = 2048
		if sizey1r <= 0:
			sizey1r = 0
		if sizey2r >= 1489:
			sizey2r = 1489
			
		if sizex1u <= 0:
			sizex1u = 0
		if sizex2u >= 2048:
			sizex2u = 2048
		if sizey1u <= 0:
			sizey1u = 0
		if sizey2u >= 1489:
			sizey2u = 1489
			
		if sizex1g <= 0:
			sizex1g = 0
		if sizex2g >= 2048:
			sizex2g = 2048
		if sizey1g <= 0:
			sizey1g = 0
		if sizey2g >= 1489:
			sizey2g = 1489
		
		if sizex1i <= 0:
			sizex1i = 0
		if sizex2i >= 2048:
			sizex2i = 2048
		if sizey1i <= 0:
			sizey1i = 0
		if sizey2i >= 1489:
			sizey2i = 1489
		
		if sizex1z <= 0:
			sizex1z = 0
		if sizex2z >= 2048:
			sizex2z = 2048
		if sizey1z <= 0:
			sizey1z = 0
		if sizey2z >= 1489:
			sizey2z = 1489
		
		# Crop the images
		crop_image_u = orig_image_u[sizey1u:sizey2u, sizex1u:sizex2u]
		crop_image_g = orig_image_g[sizey1g:sizey2g, sizex1g:sizex2g]
		crop_image_r = orig_image_r[sizey1r:sizey2r, sizex1r:sizex2r]
		crop_image_i = orig_image_i[sizey1i:sizey2i, sizex1i:sizex2i]
		crop_image_z = orig_image_z[sizey1z:sizey2z, sizex1z:sizex2z]
		
		# Check for the shapes
		print(bg_image.shape,crop_image_u.shape,crop_image_g.shape,crop_image_r.shape,crop_image_i.shape,crop_image_z.shape)
		
		# If cropping produced non square images, here it is corrected using the background image
		if len(crop_image_u) < len(bg_image):
			c_u = bg_image.copy()
			c_u[:len(crop_image_u)] += crop_image_u
		else:
			c_u = crop_image_u.copy()
			c_u[:len(bg_image)] += bg_image

		if len(crop_image_g) < len(bg_image):
			c_g = bg_image.copy()
			c_g[:len(crop_image_g)] += crop_image_g
		else:
			c_g = crop_image_g.copy()
			c_g[:len(bg_image)] += bg_image
			
		if len(crop_image_r) < len(bg_image):
			c_r = bg_image.copy()
			c_r[:len(crop_image_r)] += crop_image_r
		else:
			c_r = crop_image_r.copy()
			c_r[:len(bg_image)] += bg_image
			
		if len(crop_image_i) < len(bg_image):
			c_i = bg_image.copy()
			c_i[:len(crop_image_i)] += crop_image_i
		else:
			c_i = crop_image_i.copy()
			c_i[:len(bg_image)] += bg_image
			
		if len(crop_image_z) < len(bg_image):
			c_z = bg_image.copy()
			c_z[:len(crop_image_z)] += crop_image_z
		else:
			c_z = crop_image_z.copy()
			c_z[:len(bg_image)] += bg_image
		
		
		# Write fits files of each filter with the header of the original field image
		fits.writeto('u_%s_.fits' % ID[num], c_u, image_header_u, overwrite=True)
		fits.writeto('g_%s_.fits' % ID[num], c_g, image_header_g, overwrite=True)
		fits.writeto('r_%s_.fits' % ID[num], c_r, image_header_r, overwrite=True)
		fits.writeto('i_%s_.fits' % ID[num], c_i, image_header_i, overwrite=True)
		fits.writeto('z_%s_.fits' % ID[num], c_z, image_header_z, overwrite=True)
	
	# Remove original frame files to save space
	# os.system('rm -r frame-*.fits')
	
