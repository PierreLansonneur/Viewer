##########################################
# Crop the array stored in a dicom file
##########################################
import numpy as np
import pydicom

dir_ct = '/media/sf_Linux_Shared/102771/SFUD/'
filename = 'merged'

### Read file
print('\nReading file '+dir_ct+filename+'.dcm')
ds = pydicom.read_file(dir_ct+filename+'.dcm')
vol = ds.pixel_array
origin = ds.ImagePositionPatient
spacing = ds.PixelSpacing

### Crop limits 102771
xmin, xmax = 129, 494
ymin, ymax = 44, 476
zmin, zmax = 257, 623
vol = vol[zmin:zmax, ymin:ymax, xmin:xmax]
"""
### Crop limits 12420
xmin, xmax = 140, 500
ymin, ymax = 100, 540
zmin, zmax = 275, 675
vol = vol[zmin:zmax, ymin:ymax, xmin:xmax]
"""
### Write cropped file
ds.Rows = np.shape(vol)[1]
ds.Columns = np.shape(vol)[2]
ds.NumberOfFrames = np.shape(vol)[0]
origin[0] = str( float(origin[0]) + float(spacing[0])*xmin )
origin[1] = str(float(origin[1]) + float(spacing[1])*ymin)
spacing_z = float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0])
origin[2] = str(float(origin[2]) + spacing_z*zmin)
ds.SliceLocation = ds.SliceLocation + spacing_z*zmin
ds.PixelData = vol.tostring()
ds.ImagePositionPatient = origin
ds.GridFrameOffsetVector = ds.GridFrameOffsetVector[0:np.shape(vol)[0]]

ds.save_as(dir_ct + filename + '_cropped.dcm')

print('File cropped successfully!')

