##########################################
#           Crop a dicom serie
###########################################
import numpy as np
import pydicom

dir_ct_in = './CT/'
dir_ct_out = './CT_crop/'
xmin, xmax, ymin, ymax = 120, 385, 40, 380 # X,Y limits in voxels
zmin = -74.0 # in mm !

#####################

filelist = os.listdir(dir_ct_in)
i,N = 0,len(filelist)
print ''
for f in filelist:
	if f.endswith(".dcm"):
		### Read file infos
		ds = pydicom.read_file(dir_ct_in+f)
		im = ds.pixel_array
		origin = ds.ImagePositionPatient
		spacing = ds.PixelSpacing

		### filter z position
		if(origin[2]>zmin):
			### Crop image
			im = im[ymin:ymax,xmin:xmax]
			origin[0] = str(float(origin[0]) + float(spacing[0])*xmin)
			origin[1] = str(float(origin[1]) + float(spacing[1])*ymin)
			
			#print origin
			### Save file
			ds.PixelData = im.tostring()
			ds.Rows = np.shape(im)[0]
			ds.Columns = np.shape(im)[1]

			#origin[2] = str(float(origin[2]) + float(ds.SliceLocation))
			ds.ImagePositionPatient = origin
			ds.save_as(dir_ct_out + 'volume_{0}.dcm'.format(i)) # Save the file
			i = i+1
			print 'saving CT slices ... '+'█'*int(30*i/N) + '|'*(int(30*(N-i)/N)-1)+'\r', # progress bar
print '\ndone!'

