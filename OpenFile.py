# -*- coding: utf-8 -*-
def OpenFile(filename=None): 
	"""
	Open a .dcm or a .mhd file
	"""
	global volume, spacing, dim_x, dim_y, dim_z, info, origin, CT_open, filename_CT
	ct_swap = False

	if(filename==None):
		types = [('All files', '*.{dcm,mhd}'), ('DCM files', '*.dcm'), ('MHD files', '*.mhd')]
		file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = types)
	else:	file_path = filename

	filename_CT = file_path

	print 'Opening file...'

	### .dcm file ###
	if(file_path.endswith('.dcm')):
		ds = pydicom.read_file(file_path)#,defer_size = '2 MB', force=True)
		ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian 
		volume = ds.pixel_array
		sp = ds.PixelSpacing

		#if ("RescaleSlope" in ds):	volume = float(ds.RescaleSlope)*volume
		#if ("RescaleIntercept" in ds):	volume = volume - float(ds.RescaleIntercept)
		if ("DoseGridScaling" in ds):	volume = float(ds.DoseGridScaling)*volume

		spacing = [ float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]), float(sp[1]),float(sp[0])]
		origin = ds.ImagePositionPatient
		origin = [float(origin[2]),float(origin[1]),float(origin[0])]
		info = 'Name: {0}\nBits allocated: {1}'.format(ds.PatientName, ds.BitsAllocated)
		if ("PatientPosition" in ds):  ct_swap = (ds.PatientPosition == 'HFS')
		if ("ImageOrientationPatient" in ds):	ct_swap = (ds.ImageOrientationPatient == [1, 0, 0, 0, 1, 0])

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
		#for i in range(len(tmp)):	info = info + tmp[i]

	if(len(np.shape(volume))==3):	dim_x, dim_y, dim_z = np.shape(volume)[0], np.shape(volume)[1], np.shape(volume)[2]

	if(len(np.shape(volume))==2):	dim_x, dim_y, dim_z = np.shape(volume)[0], np.shape(volume)[1], 1

	#print 'ct type:', volume.dtype

	if(ct_swap == True):
		volume = np.flip(volume,1) # flip volume
		volume = np.flip(volume,2) # flip volume
		spacing[1], spacing[2] = spacing[2], spacing[1]
		origin[1] = origin[1] + dim_y*spacing[1]
		origin[2] = origin[2] + dim_z*spacing[2]

	print 'ct_swap:', ct_swap

	print 'File successfully opened!'

	Set_axes_lim_init()
	Set_scales()

	CT_open = True
	Update_all()

def OpenDicomSerie(dirname=None):
	"""
	Open a dicom serie
	"""
	global volume, dim_x, dim_y, dim_z, spacing, origin, CT_open, filename_CT
	ct_swap = False

	# Opening file
	if(dirname==None):
		file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = [('DICOM files', '*.dcm')])
		filelist = os.listdir(os.path.dirname(file_path))
	else:
		filelist = os.listdir(dirname)
		file_path = dirname + filelist[0]

	filename_CT = file_path

	# Getting dimensions
	ds = pydicom.read_file(file_path)
	sp = ds.PixelSpacing
	ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
	ct_swap = (ds.PatientPosition == 'HFS')
	dim_x, dim_y, dim_z = len(filelist), np.shape(ds.pixel_array)[1], np.shape(ds.pixel_array)[0]
	volume, volume_tmp = np.zeros((dim_x,dim_y,dim_z)), np.zeros((dim_x,dim_y,dim_z))
	print('')

	# Creating volume
	i=0
	order = [[]]*dim_x
	for f in filelist:
		if f.endswith(".dcm"):
			ds = pydicom.read_file(os.path.dirname(file_path)+'/'+f)
			ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian 
			volume_tmp[i,:,:] = ds.pixel_array
			if('SliceLocation' in ds):	order[i] = [i,ds.SliceLocation]
			else:	order[i] = [i,ds.ImagePositionPatient[2]]
			i=i+1
			print 'Creating CT volume ... '+'█'*int(30*i/dim_x) + '|'*int(30*(dim_x-i)/dim_x)+'\r', # progress bar
	order = sorted(order, key = itemgetter(1))
	print('')
	spacing = [float(order[1][1] - order[0][1]),float(sp[1]), float(sp[0])]
	origin = ds.ImagePositionPatient
	origin = [float(order[0][1]),float(origin[1]),float(origin[0])]

	# Sorting volume
	for j in range(dim_x-1):
		k = order[j][0]
		volume[j,:,:] = volume_tmp[k,:,:]
		print 'Sorting CT slices  ... '+'█'*int(30*j/(dim_x-2)) + '|'*int(30*(dim_x-j)/(dim_x-2))+'\r', # progress bar
	print('\nFile successfully opened!')

	if ("RescaleSlope" in ds):	volume = float(ds.RescaleSlope)*volume
	if ("RescaleIntercept" in ds):	volume = volume + float(ds.RescaleIntercept) 

	if(ct_swap == True):
		volume = np.flip(volume,1) # flip volume
		volume = np.flip(volume,2) # flip volume
		origin[1] = - origin[1]
		origin[2] = - origin[2]

	print 'ct_swap:', ct_swap

	Set_axes_lim_init()
	Set_scales()
	CT_open = True
	Update_all()

def OpenDosi(filename=None): 
	"""
	Open a dosimetry file with a .dcm or .mhd extension
	"""
	global dosi, spacing_dosi, dim_x_dosi, dim_y_dosi, dim_z_dosi, dosi_open, isodose_show, origin_dosi, filename_dosi
	dosi_swap = False

	types = [('All files', '*.{dcm,mhd}'), ('DCM files', '*.dcm'), ('MHD files', '*.mhd')]

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
		#dosi = (100./0.000280315)*scaling_dosi*ds.pixel_array
		#dosi = ds.pixel_array
		sp = ds.PixelSpacing
		spacing_dosi = [ float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]), float(sp[1]),float(sp[0])]
		origin_dosi = ds.ImagePositionPatient
		origin_dosi = [float(origin_dosi[2]),float(origin_dosi[1]),float(origin_dosi[0])]
		if ("PatientPosition" in ds):	dosi_swap = (ds.PatientPosition == 'HFS')
		if ("ImageOrientationPatient" in ds):	dosi_swap = (ds.ImageOrientationPatient == [1, 0, 0, 0, 1, 0])

		if ds.SeriesDescription=='PatientLETScorer [MeV/mm/(g/cm3)]':	SetIntensityRange(dosi,0,15)

	### .mhd file ###
	if(file_path.endswith('.mhd')):	
    		itkimage = sitk.ReadImage(file_path)    				# Reads the image using SimpleITK
    		dosi = sitk.GetArrayFromImage(itkimage)
		spacing_dosi = np.array(list(reversed(itkimage.GetSpacing())))    	# Read the spacing along each dimension
		origin_dosi = np.array(list(reversed((itkimage.GetOrigin()))))		# Read the origin
		text_file = open(file_path, "r")
		tmp = text_file.readlines()
		dosi_swap = (tmp[8][-4:-1] == 'RAI')
		#dosi_swap = True
		print tmp[8][-4:-1]

	print('File successfully opened!')

	if(len(np.shape(volume))==3):	dim_x_dosi, dim_y_dosi, dim_z_dosi = np.shape(dosi)[0], np.shape(dosi)[1], np.shape(dosi)[2]

	if(len(np.shape(volume))==2):	dim_x_dosi, dim_y_dosi, dim_z_dosi = np.shape(dosi)[0], np.shape(dosi)[1], 1

	#dosi = np.array(dosi, dtype = 'uint16')
	#dosi = np.array(dosi, dtype = 'float16')
	#print 'dosi type', dosi.dtype

	#if (ID=='104005')and(file_path!=dir_ini+ID+'/pi/RD'+ID+'.dcm'):
	if (file_path==dir_ini+'104005/SFUD/merged.dcm'):
		print '**** 104005 ******'
		dosi_swap=False
		#dosi = np.flip(dosi,1) # flip volume
		#dosi = np.flip(dosi,2) # flip volume
		origin_dosi[1] = origin_dosi[1] -327 #+ (dim_y_dosi*spacing_dosi[1]-dim_y*spacing[2]) +20 # -344
		origin_dosi[2] = origin_dosi[2] -187 #+ (dim_z_dosi*spacing_dosi[2]-dim_z*spacing[2]) +10 # -192 

	if(dosi_swap == True):
		dosi = np.flip(dosi,1) # flip volume
		dosi = np.flip(dosi,2) # flip volume
		spacing_dosi[1], spacing_dosi[2] = spacing_dosi[2], spacing_dosi[1]
		origin_dosi[1] = origin_dosi[1] + dim_y_dosi*spacing_dosi[1]
		origin_dosi[2] = origin_dosi[2] + dim_z_dosi*spacing_dosi[2]

	print 'dosi_swap:', dosi_swap

	dosi_open = True
	isodose_show = True
	check1.select()
	Update_all()

