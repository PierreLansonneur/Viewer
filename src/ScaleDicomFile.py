##############################################
# change dose scaling factor of a DICOM file
##############################################
import pydicom
import numpy as np

directory = '/media/sf_Linux_Shared/102771/ctc4/'
filename =  'merged'
coeff = 100/0.05174145973019193

### Read file
ds = pydicom.read_file(directory+filename+'.dcm')
vol = np.array(ds.pixel_array)

'''
### Change data type to np.uint8
info = np.iinfo(vol.dtype)
vol = vol.astype(np.float64)/info.max
vol = 255*vol
vol = vol.astype(np.uint8)
coeff = (float(info.max)/255.)*coeff
ds.BitsAllocated = 8
ds.BitsStored = 8
ds.HighBit = 8
print 'new data type:', vol.dtype
'''

### change dose scaling factor
print 'new scaling:', coeff*ds.DoseGridScaling
ds.DoseGridScaling = coeff*ds.DoseGridScaling

### Write file
ds.PixelData = vol.tostring()
ds.save_as(directory + filename + '_Isogray.dcm')
print('\nFile saved successfully!')
