# -*- coding: utf-8 -*-
###################################################
#   	     Open DICOM and .mhd files
# 		   	----
#           	P. Lansonneur 2019
###################################################

def OpenFile(filename=None): 
	"""
	Open a .dcm or a .mhd file
	"""
	global volume, spacing, dim_x, dim_y, dim_z, origin, CT_open, filename_CT, dir_ini
	ct_swapY, ct_swapZ = False, False

	if(filename==None):
		types = [('All files', '*.dcm *.mhd'), ('DCM files', '*.dcm'), ('MHD files', '*.mhd')]
		file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = types)
	else:	file_path = filename

	filename_CT = file_path
        dir_ini = str(file_path.rsplit('/', 1)[0])+'/'
	print 'Opening file...'

	### .dcm file ###
	if(file_path.endswith('.dcm')):
		ds = pydicom.read_file(file_path)
		ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian 
		volume = ds.pixel_array

                try:
		        spacing[0:1] = ds.PixelSpacing
	                origin[0:1] = ds.ImagePositionPatient
                except Exception:
		        spacing = ds.PixelSpacing
	                origin = ds.ImagePositionPatient

                if (ds.Modality == 'RTDOSE'):
		        if ("DoseGridScaling" in ds):	volume = float(ds.DoseGridScaling)*volume
                else:
	                if ("RescaleSlope" in ds):	volume = float(ds.RescaleSlope)*volume
	                if ("RescaleIntercept" in ds):	volume = volume + float(ds.RescaleIntercept)

	        if(len(np.shape(volume))==3):
		        spacing = [ float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]), float(spacing[1]),float(spacing[0])]
		        origin = [float(origin[2]),float(origin[1]),float(origin[0])]

                ct_swapZ =(ds.ImageOrientationPatient[0:3] == [1, 0, 0])
                ct_swapY =(ds.ImageOrientationPatient[3:6] == [0, 1, 0])

		#if ds.SeriesDescription=='PatientLETScorer [MeV/mm/(g/cm3)]':	SetIntensityRange(volume,0,15)

	### .mhd file ###
	if(file_path.endswith('.mhd')):	
    		itkimage = sitk.ReadImage(file_path)    			# Reads the image using SimpleITK
    		volume = sitk.GetArrayFromImage(itkimage)
		spacing = np.array(list(reversed(itkimage.GetSpacing())))    	# Read the spacing along each dimension
		origin = np.array(list(reversed((itkimage.GetOrigin()))))	# Read the origin of the ct_scan
		text_file = open(file_path, "r")
		tmp = text_file.readlines()
		ct_swap = (tmp[8][-4:-1] == 'RAI')

	if(len(np.shape(volume))==3):
        	dim_x, dim_y, dim_z = np.shape(volume)[0], np.shape(volume)[1], np.shape(volume)[2]

	# Dealing with image orientation
        print '  ct_swapY, ct_swapZ :', ct_swapY, ct_swapZ
	if(ct_swapY == True):
                volume = np.flip(volume,1) # flip volume, Y direction
                origin[1] = origin[1] + dim_y*spacing[1]                
        if(ct_swapZ == True):
                volume = np.flip(volume,2) # flip volume, Z direction
                origin[2] = origin[2] + dim_z*spacing[2]                
        if(ct_swapZ == True)and(ct_swapY == True):      spacing[1], spacing[2] = spacing[2], spacing[1]

	if(len(np.shape(volume))==2):	dim_x, dim_y, dim_z = np.shape(volume)[0], np.shape(volume)[1], 0

	Set_axes_lim_init()
	Set_scales()
	CT_open = True
	Update_all()

	print '  file successfully opened!'
	
def OpenDicomSerie(dirname=None):
	"""
	Open a dicom serie
	"""
	global volume, dim_x, dim_y, dim_z, spacing, origin, CT_open, filename_CT, dir_ini
        ct_swapY, ct_swapZ = False, False
        
	print 'Opening DICOM serie ... '

	# Opening file
	if(dirname==None):
		file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = [('DICOM files', '*.dcm')])
		filelist = os.listdir(os.path.dirname(file_path))
	else:
		filelist = os.listdir(dirname)
		file_path = dirname + filelist[0]

	filename_CT = file_path
        dir_ini = str(file_path.rsplit('/', 1)[0])+'/'

	# Getting dimensions
	ds = pydicom.read_file(file_path)
	sp = ds.PixelSpacing
	ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
	ct_swapZ =(ds.ImageOrientationPatient[0:3] == [1, 0, 0])
	ct_swapY =(ds.ImageOrientationPatient[3:6] == [0, 1, 0])

        dim_x = 0
        for f in filelist:
                if f.endswith(".dcm"): dim_x = dim_x + 1        

	dim_y, dim_z = np.shape(ds.pixel_array)[1], np.shape(ds.pixel_array)[0]
    
	volume = np.zeros((dim_x, dim_y,dim_z))
        slicelocation = np.zeros(dim_x)

	# creating volume
	for f,i in zip(filelist,range(dim_x)):
		if f.endswith(".dcm"):
			ds = pydicom.read_file(os.path.dirname(file_path)+'/'+f)
			ds.file_meta.transfersyntaxuid = pydicom.uid.ImplicitVRLittleEndian 
			volume[i,:,:] = ds.pixel_array
			if('slicelocation' in ds):	slicelocation[i] = ds.SliceLocation
			else:	slicelocation[i] = ds.ImagePositionPatient[2]
            
	order = np.argsort(slicelocation)
        slicelocation = slicelocation[order] # slicelocation is now sorted
    
	spacing = [float(slicelocation[1] - slicelocation[0]),float(sp[1]), float(sp[0])]
	origin = [float(slicelocation[0]),float(ds.ImagePositionPatient[1]),float(ds.ImagePositionPatient[0])]
	volume = volume[order,:,:] # volume is now sorted

	if ("RescaleSlope" in ds):	volume = float(ds.RescaleSlope)*volume
	if ("RescaleIntercept" in ds):	volume = volume + float(ds.RescaleIntercept)

	# Dealing with image orientation
        print '  ct_swapY, ct_swapZ :', ct_swapY, ct_swapZ
	if(ct_swapY == True):
                volume = np.flip(volume,1) # flip volume, Y direction
                origin[1] = origin[1] + dim_y*spacing[1]                
        if(ct_swapZ == True):
                volume = np.flip(volume,2) # flip volume, Z direction
                origin[2] = origin[2] + dim_z*spacing[2]                
        if(ct_swapZ == True)and(ct_swapY == True):      spacing[1], spacing[2] = spacing[2], spacing[1]

	Set_axes_lim_init()
	Set_scales()
	CT_open = True
	Update_all()

	print('  file successfully opened!')

def OpenDosi(filename=None): 
	"""
	Open a dosimetry file with a .dcm or .mhd extension
	"""
	global dosi, spacing_dosi, dim_x_dosi, dim_y_dosi, dim_z_dosi, dosi_open, isodose_show, origin_dosi, filename_dosi
	dosi_swapY,dosi_swapZ = False, False

	types = [('All files', '*.dcm *.mhd'), ('DCM files', '*.dcm'), ('MHD files', '*.mhd')]

	if(filename==None):	file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = types)
	else:	file_path = filename

	filename_dosi = file_path

	print('Opening RD file ...')

	### .dcm file ###
	if(file_path.endswith('.dcm')):
		ds = pydicom.read_file(file_path)
		ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian 
		scaling_dosi = float(ds.DoseGridScaling)
		dosi = scaling_dosi*ds.pixel_array
		sp = ds.PixelSpacing
		spacing_dosi = [ float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]), float(sp[1]),float(sp[0])]
		origin_dosi = ds.ImagePositionPatient
		origin_dosi = [float(origin_dosi[2]),float(origin_dosi[1]),float(origin_dosi[0])]
		dosi_swapZ =(ds.ImageOrientationPatient[0:3] == [1, 0, 0])
                dosi_swapY =(ds.ImageOrientationPatient[3:6] == [0, 1, 0])

		#if ds.SeriesDescription=='PatientLETScorer [MeV/mm/(g/cm3)]':	SetIntensityRange(dosi,0,15)

	### .mhd file ###
	if(file_path.endswith('.mhd')):	
    		itkimage = sitk.ReadImage(file_path)    				# Reads the image using SimpleITK
    		dosi = sitk.GetArrayFromImage(itkimage)
		spacing_dosi = np.array(list(reversed(itkimage.GetSpacing())))    	# Read the spacing along each dimension
		origin_dosi = np.array(list(reversed((itkimage.GetOrigin()))))		# Read the origin
		text_file = open(file_path, "r")
		tmp = text_file.readlines()
		dosi_swap = (tmp[8][-4:-1] == 'RAI')

	if(len(np.shape(volume))==3):	dim_x_dosi, dim_y_dosi, dim_z_dosi = np.shape(dosi)[0], np.shape(dosi)[1], np.shape(dosi)[2]

	if(len(np.shape(volume))==2):	dim_x_dosi, dim_y_dosi, dim_z_dosi = np.shape(dosi)[0], np.shape(dosi)[1], 1

	#print 'dosi type', dosi.dtype
	
	# Dealing with image orientation
	if(dosi_swapY == True):
		dosi = np.flip(dosi,1) # flip volume
		origin_dosi[1] = origin_dosi[1] + dim_y_dosi*spacing_dosi[1]		
	if(dosi_swapZ == True):
		dosi = np.flip(dosi,2) # flip volume
		origin_dosi[2] = origin_dosi[2] + dim_z_dosi*spacing_dosi[2]
	if(dosi_swapY == True)and(dosi_swapZ == True):
		spacing_dosi[1], spacing_dosi[2] = spacing_dosi[2], spacing_dosi[1]

        print '  dosi_swapY, dosi_swapZ :', dosi_swapY, dosi_swapZ

	dosi_open = True
	isodose_show = True
	check1.select()
	Update_all()

	print('  file successfully opened!')
