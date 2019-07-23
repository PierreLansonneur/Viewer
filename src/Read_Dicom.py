##########################################
# reads the info (tags) from a dicom file
#	    P. Lansonneur - 2019
###########################################

import pydicom
from pydicom.dataset import Dataset, FileDataset

##################### Read Dicom info #######################
filename = './MyDicomFile.dcm'

ds = pydicom.read_file(filename, force=True)
print ds

with open('tmp.txt', 'w') as f:     print >> f , ds  
