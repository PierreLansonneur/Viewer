##########################################
# reads the info (tags) from a dicom file
#	 P. Lansonneur - 2018
###########################################

import pydicom
from pydicom.dataset import Dataset, FileDataset

##################### Read Dicom info #######################
#filename = '/media/sf_Linux_Shared/MHD_files/dosimetry/TOPAS_MBRT/Dose_MBRT_1.dcm'
#filename = '/media/sf_Linux_Shared/MHD_files/CT_head_RAULT/CT.1.2.392.200036.9116.2.6.1.48.1211376057.1438919546.6618030.dcm'
#filename = '/media/sf_Linux_Shared/MHD_files/dosimetry/RAULT/RD.rault_AS3N.dcm'
#filename = '/media/sf_Linux_Shared/souris_CONV/volume0000.dcm'
#filename = '/media/sf_Linux_Shared/0000/CT/volume_263.dcm'
#ID = '102771'
#ID = '101534'
ID = '12420'
#ID = '104005'

#filename = '/media/sf_Linux_Shared/'+ID+'/pi/RD'+ID+'.dcm'
filename = '/media/sf_Linux_Shared/'+ID+'/SFUD/merged.dcm'

ds = pydicom.read_file(filename, force=True)
print ds
#print ds.SeriesDescription=='PatientLETScorer [MeV/mm/(g/cm3)]'

with open('tmp.txt', 'w') as f:	
    print >> f , ds  

'''
filename1 = '/media/sf_Linux_Shared/'+ID+'/pi/RD'+ID+'.dcm'
ds1 = pydicom.read_file(filename1, force=True)
print 'ok'
ds.DoseUnits='GY'
ds.DoseType='PHYSICAL'
ds.DoseSummationType='PLAN'
ds.ReferencedRTPlanSequence = ds1.ReferencedRTPlanSequence
ds.ReferencedStructureSetSequence = ds1.ReferencedStructureSetSequence
ds.save_as('/media/sf_Linux_Shared/'+ID+'/Dose_PS1N_ctc4_test.dcm') # Save the file
print 'ok'
'''
"""
ds.PatientPosition = 'HFS'
print ds.PatientPosition
ds.save_as('/media/sf_Linux_Shared/'+ID+'/PBS_22Fev2019/RD'+ID+'_test.dcm') # Save the file
"""
