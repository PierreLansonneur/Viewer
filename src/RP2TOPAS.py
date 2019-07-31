#########################################################
# 	Convert RT-Plan dicom to TOPAS file
# 		     (PBS & DS)
#	 	P. Lansonneur - 2019
#########################################################

def UpdateRPBoundingBox(ext7,ext8,ext9, rot1, rot2, rot3, ax_):
        """
        Update dosi bounding box for TOPAS calculation (red rectangle)
        """
        if(ax_==ax1):
	        width_x, width_y = (wy_max.get()-wy_min.get())*spacing[1], (wz_max.get()-wz_min.get()-1)*spacing[2]
                x_min, y_min = ext7[1] + wy_min.get()*spacing[1], ext7[2] + wz_min.get()*spacing[2]
	        if(rot1==True): 
                        y_min, x_min = -ext7[2] + wy_min.get()*spacing[1], ext7[1] + wz_min.get()*spacing[2] 
        	        width_x, width_y = width_y, width_x

        if(ax_==ax2):
	        width_x, width_y = (wx_max.get()-wx_min.get())*spacing[0], (wz_max.get()-wz_min.get()-1)*spacing[2]
                x_min, y_min = ext8[0] + wx_min.get()*spacing[0], ext8[2] + wz_min.get()*spacing[2]
	        if(rot2==True):
                        y_min, x_min = ext8[2] + wx_min.get()*spacing[0], -ext8[0] + wz_min.get()*spacing[2]
        	        width_x, width_y = width_y, width_x

        if(ax_==ax3):
	        width_x, width_y = (wx_max.get()-wx_min.get())*spacing[0], (wy_max.get()-wy_min.get()-1)*spacing[1]
	        x_min, y_min = ext9[0] + wx_min.get()*spacing[0], ext9[2] + wy_min.get()*spacing[1]
	        if(rot3==True): 
                        y_min, x_min = ext9[2] + wx_min.get()*spacing[0], -ext9[0] + wy_min.get()*spacing[1]
        	        width_x, width_y = width_y, width_x

        try:	del ax_.patches[-1]
        except Exception:	pass

        ax_.add_patch(patches.Rectangle((x_min,y_min),width_x,width_y,linewidth=1,edgecolor='r',facecolor='none'))

def OpenRP(filename=None):
	global filename_RP
	"""
	Open a RP file with a .dcm extension
	"""

	if (CT_open == False): 
		tkMessageBox.showerror('Error', 'You must open a CT first') #print 'You must open a CT first..'
		return

	if(filename==None):	file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = [('DCM files', '*.dcm')])
	else:	file_path = filename
	filename_RP = file_path

	RP_window()
        #ReadRP(createTopasfile = False, check = False)
        ReadRP(createTopasfile = False, check = True)

def RP_window():
	"""
	visualize RP file and create a TOPAS simulation file
	"""
	global RP_crop_show, wx_min, wx_max, wy_min, wy_max, wz_min, wz_max, window_RP, defaulttroughcolor

	window_RP = Toplevel(window)
	window_RP.geometry("{0}x{1}+{2}+{3}".format(int(0.47*w_px), 100, int(0.53*w_px), 0))
	window_RP.title('RP file: {0}'.format(filename_RP.rsplit('/', 1)[-1]))

        '''
	### Menu bar ##################
	menubar_RP = Menu(window_RP)
	window_RP.config(menu=menubar_RP)
	filemenu_RP = Menu(menubar_RP, tearoff=0)
	#menubar_RP.add_cascade(label="Open", menu=filemenu_RP)
	#filemenu_RP.add_command(label="select file..", command=OpenRP)
	menubar_RP.add_command(label="select file..", command=OpenRP)
        '''

        off_grid = 0

        ### TOPAS file settings
	wx_min_val = IntVar()
	wx_max_val = IntVar()
	wy_min_val = IntVar()
	wy_max_val = IntVar()
	wz_min_val = IntVar()
	wz_max_val = IntVar()

        scale_l = 250
	wx_min = Scale(window_RP, from_=1, to=dim_x, orient=HORIZONTAL, variable=wx_min_val, command=w1move, showvalue = 0, length=scale_l)
	wx_max = Scale(window_RP, from_=1, to=dim_x, orient=HORIZONTAL, variable=wx_max_val, command=w1move, showvalue = 0, length=scale_l)
	wy_min = Scale(window_RP, from_=1, to=dim_y, orient=HORIZONTAL, variable=wy_min_val, command=w1move, showvalue = 0, length=scale_l)
	wy_max = Scale(window_RP, from_=1, to=dim_y, orient=HORIZONTAL, variable=wy_max_val, command=w1move, showvalue = 0, length=scale_l)
	wz_min = Scale(window_RP, from_=1, to=dim_z, orient=HORIZONTAL, variable=wz_min_val, command=w1move, showvalue = 0, length=scale_l)
	wz_max = Scale(window_RP, from_=1, to=dim_z, orient=HORIZONTAL, variable=wz_max_val, command=w1move, showvalue = 0, length=scale_l)
        """
	wx_min = Scale(window_RP, from_=0, to=dim_x-1, orient=HORIZONTAL, variable=wx_min_val, command=w1move, showvalue = 0, length=scale_l)
	wx_max = Scale(window_RP, from_=0, to=dim_x-1, orient=HORIZONTAL, variable=wx_max_val, command=w1move, showvalue = 0, length=scale_l)
	wy_min = Scale(window_RP, from_=0, to=dim_y-1, orient=HORIZONTAL, variable=wy_min_val, command=w1move, showvalue = 0, length=scale_l)
	wy_max = Scale(window_RP, from_=0, to=dim_y-1, orient=HORIZONTAL, variable=wy_max_val, command=w1move, showvalue = 0, length=scale_l)
	wz_min = Scale(window_RP, from_=0, to=dim_z-1, orient=HORIZONTAL, variable=wz_min_val, command=w1move, showvalue = 0, length=scale_l)
	wz_max = Scale(window_RP, from_=0, to=dim_z-1, orient=HORIZONTAL, variable=wz_max_val, command=w1move, showvalue = 0, length=scale_l)

	wx_max.set(dim_x-1)
	wy_max.set(dim_y-1)
	wz_max.set(dim_z-1)
        """
	#s_grid = Spinbox(window_RP, from_=0.1, to=5, increment=0.1, width=3)

	wx_max.set(dim_x)
	wy_max.set(dim_y)
	wz_max.set(dim_z)

	txmin_val = Label(window_RP, textvariable = wx_min_val)
	txmax_val = Label(window_RP, textvariable = wx_max_val)
	tymin_val = Label(window_RP, textvariable = wy_min_val)
	tymax_val = Label(window_RP, textvariable = wy_max_val)
	tzmin_val = Label(window_RP, textvariable = wz_min_val)
	tzmax_val = Label(window_RP, textvariable = wz_max_val)

	txmin = Label(window_RP, text=" X from :")
	txmax = Label(window_RP, text=" to ")
	tymin = Label(window_RP, text=" Y from :")
	tymax = Label(window_RP, text=" to ")
	tzmin = Label(window_RP, text=" Z from :")
	tzmax = Label(window_RP, text=" to ")
	#tgrid = Label(window_RP, text=" grid : ")
	tspace = Label(window_RP, text="        ")
	tspace2 = Label(window_RP, text="        ")
	tcrop = Label(window_RP, text="keep CT dim.")

	wx_min.grid(row=off_grid+0,column=2)
	wx_max.grid(row=off_grid+0,column=4)
	wy_min.grid(row=off_grid+1,column=2)
	wy_max.grid(row=off_grid+1,column=4)
	wz_min.grid(row=off_grid+2,column=2)
	wz_max.grid(row=off_grid+2,column=4)
	#s_grid.grid(row=off_grid+3,column=1,sticky=W+N+S)

	txmin.grid(row=off_grid+0,column=0)
	txmax.grid(row=off_grid+0,column=3)
	tymin.grid(row=off_grid+1,column=0)
	tymax.grid(row=off_grid+1,column=3)
	tzmin.grid(row=off_grid+2,column=0)
	tzmax.grid(row=off_grid+2,column=3)
	tspace.grid(row=off_grid+3,column=1)
	tspace2.grid(row=off_grid+3,column=5)

	txmin_val.grid(row=off_grid+0,column=1)
	txmax_val.grid(row=off_grid+0,column=5)
	tymin_val.grid(row=off_grid+1,column=1)
	tymax_val.grid(row=off_grid+1,column=5)
	tzmin_val.grid(row=off_grid+2,column=1)
	tzmax_val.grid(row=off_grid+2,column=5)

	RP_bt1 = Button(window_RP, text="create TOPAS file(s)", command=CreateTopasFile)
	RP_bt1.grid(row=off_grid+1,column=6,rowspan=2,columnspan=2,sticky=W+N+S)
	checkRP1 = Checkbutton(window_RP, command = KeepCTDim)
	checkRP1.grid(row=off_grid+0,column=6,sticky=W+N+S)
        tcrop.grid(row=off_grid+0,column=7,sticky=W+N+S)

        defaulttroughcolor = wx_min.cget('troughcolor') # #b3b3b3

	RP_crop_show = True

	Update_all()

def KeepCTDim():
        global crop_RP, keepCTDim

	keepCTDim = not keepCTDim

        if (keepCTDim==True):
	        wx_min.set(1)
	        wy_min.set(1)
	        wz_min.set(1)
	        wx_max.set(dim_x)
	        wy_max.set(dim_y)
	        wz_max.set(dim_z)

                for w_ in [wx_min,wy_min,wz_min,wx_max,wy_max, wz_max]:
                        w_.config(state=DISABLED, troughcolor='lightgrey')

        if (keepCTDim==False):
                for w_ in [wx_min,wy_min,wz_min,wx_max,wy_max, wz_max]:
                        w_.config(state=NORMAL, troughcolor=defaulttroughcolor)
        Update_all()

def CreateTopasFile():
        """
        Create TOPAS MC file
        """
        global RP_crop_show
        ReadRP(createTopasfile = True, check = False)
        RP_crop_show = False
        Update_all()

def R0(E0):	
	"""
	return the range in mm for E0 expressed in MeV
	"""
	return	10*0.0022*(E0**1.77) 

def ReadRP(createTopasfile = False, check = True):
	global ID, window_RP, RP_crop_show
        global BeamName, N, GantryAngle, PatientSupportAngle, SnoutID, DCI, epRS, TransX, TransY, TransZ, BeamDose, BeamMeterset

	print('\nOpening file {0}'.format(filename_RP))
	ds = pydicom.read_file(filename_RP,force=True)

	if ("TargetPrescriptionDose" in ds.DoseReferenceSequence[0]):	
                print 'TargetPrescriptionDose :', ds.DoseReferenceSequence[0].TargetPrescriptionDose, 'Gy'
	        #disp0 = Label(window_RP, text=str(ds.DoseReferenceSequence[0].TargetPrescriptionDose))
                #disp0.grid(row= 0 , column=1, sticky=W)

	ibs = ds.IonBeamSequence

	if ID=='':	ID = ds.PatientID

	### Check RP modality (DS or PBS)
	RP_type = ''
	for k in range(len(ibs)):
		BeamName = str(ibs[k].BeamName )
		if(BeamName=='SETUP')or(BeamName=='setup'):	continue
		if(ibs[k].TreatmentDeliveryType=='SETUP'):	continue	
		if ibs[k].ScanMode=='NONE':			RP_type = 'DS'
		if ibs[k].ScanMode=='MODULATED':		RP_type = 'PBS'

	if RP_type=='DS':

		### Dose per fraction
		BeamDose = np.zeros(len(ibs))
		BeamMeterset = np.zeros(len(ibs))
		fgs = ds.FractionGroupSequence 
		for i in range(len(fgs)):
			rbs = fgs[i].ReferencedBeamSequence
			for j in range(len(rbs)):
				if(ibs[k].TreatmentDeliveryType=='SETUP'):      continue
                                try:	
				        BeamDose[rbs[j].ReferencedBeamNumber-1] = rbs[j].BeamDose		# Dose per fraction
				        BeamMeterset[rbs[j].ReferencedBeamNumber-1] = rbs[j].BeamMeterset	# UM per fraction
                                except Exception:       pass

		### Reads tags for each beams (orientation, isocenter, snout, etc...)
		for k in range(len(ibs)):

			if(ibs[k].TreatmentDeliveryType=='SETUP'):	continue

			icps = ibs[k].IonControlPointSequence
			BeamName = ibs[k].BeamName
			BeamName = BeamName.replace(' ','_')
			E = icps[0].NominalBeamEnergy

			### collimator
			#if('IonBlockSequence' in ibs[k]):
			iblocks = ibs[k].IonBlockSequence 
			BlockThickness = iblocks[0].BlockThickness
			BlockData = np.array(iblocks[0].BlockData)

		 	### compensator
			ircs = ibs[k].IonRangeCompensatorSequence 
			CRows,CColumns = ircs[0].CompensatorRows, ircs[0].CompensatorColumns
			CompensatorMillingToolDiameter = ircs[0].CompensatorMillingToolDiameter
			CompensatorThicknessData = np.array(ircs[0].CompensatorThicknessData).reshape((CRows,CColumns))
			CPixelSpacing = ircs[0].CompensatorPixelSpacing 
			CPosition = ircs[0].CompensatorPosition
			ext = [CPosition[0],CPosition[0]+CColumns*CPixelSpacing[0],CPosition[1]-CRows*CPixelSpacing[1],CPosition[1]]

			### range shifter, snout, angles, orientation
			epRS=0							
			if (ibs[k].RangeShifterSequence[0].RangeShifterNumber==1):	epRS=65
			if (ibs[k].RangeShifterSequence[0].RangeShifterNumber==2):	epRS=30

			GantryAngle = icps[0].GantryAngle  					# Gantry angle
			PatientSupportAngle = icps[0].PatientSupportAngle	# Table angle
			DCI =  	0.1*icps[0].SnoutPosition
			SnoutID = ibs[k].SnoutSequence[0].SnoutID
			dSnout =-6.5
			if (SnoutID==400):	dSnout=1.8
			retraction = DCI + dSnout

			# Isocenter and Patient translation
			#IsocenterPosition = [0,0,0]
			IsocenterPosition = icps[0].IsocenterPosition 	# in mm
                        if(keepCTDim==False):   crop_RP = [wz_min.get()-1, wz_max.get(), wy_min.get()-1, wy_max.get(), wx_min.get()-1, wx_max.get()]
                        else:                   crop_RP = [0, dim_z,0, dim_y,0, dim_x]

                        TransX, TransY, TransZ = PatientTranslation(IsocenterPosition, PatientSupportAngle, crop_RP)
			#if (pss[0].PatientPosition=='HFP'):	TransY, TransZ = -TransY, -TransZ

                        if(check == True):
			        print '\n------ Beam',ibs[k].BeamName , '------'
			        print 'Energy :		', E, 'MeV'
			        print "GantryAngle :		", GantryAngle, 'deg'
			        print "PatientSupportAngle :	", PatientSupportAngle, 'deg'
			        print "SnoutID :		", SnoutID
			        print "DCI :			{0:.1f} cm".format(DCI)
			        print "Range Shifter :		", epRS, 'mm'
			        #print "retraction :		", retraction, 'cm'
			        print "TransX :		{0:.1f} mm".format(TransX)
			        print "TransY :		{0:.1f} mm".format(TransY)
			        print "TransZ :		{0:.1f} mm".format(TransZ)
			        print 'Beam Dose :		', BeamDose[k], 'Gy'
			        print 'Beam Meterset :		', BeamMeterset[k], 'UM'

			        check_DS(BeamName, BlockData, CompensatorThicknessData, ext)

                        if(createTopasfile == True):
			        CreateApertureFileDS(BeamName,BlockData)
			        CreateCompensatorFileDS(BeamName,ircs)

	if RP_type=='PBS':

		### Number of fractions
		print 'Reference volume :\t', ds.DoseReferenceSequence[0].DoseReferenceDescription
		print 'Number of Fractions :\t', ds.FractionGroupSequence[0].NumberOfFractionsPlanned

		### Reads tags for each beams (orientation, isocenter, snout, etc...)
		pss = ds.PatientSetupSequence

		for k in range(len(ibs)):
			icps = ibs[k].IonControlPointSequence # one icps per beam
			BeamName = str(ibs[k].BeamName )
			BeamName = BeamName.replace(' ','_')

			if(BeamName=='SETUP')or(BeamName=='setup'):	continue 

			### Retrieve the spots energy, UM and position
			N,E,x,y,UM = 0,[],[],[],[]	# total number of spots, energy (MeV), x, y and UM for each spots

			for j in range(len(icps)):

				x = np.concatenate((x,icps[j].ScanSpotPositionMap[0::2]))
				y = np.concatenate((y,icps[j].ScanSpotPositionMap[1::2]))
				if(icps[j].NumberOfScanSpotPositions==1):	UM = np.concatenate((UM,np.array([icps[j].ScanSpotMetersetWeights])))
				else:	UM = np.concatenate((UM,icps[j].ScanSpotMetersetWeights))
				E = np.concatenate((E,icps[j].NominalBeamEnergy*np.ones((icps[j].NumberOfScanSpotPositions))))

			# remove spots with zero weights
			index = np.where(UM>0)
			E, x, y, UM = E[index], x[index], y[index], UM[index]
			N = len(UM)

			# range shifter thickness (in mm)
			epRS=0	
			#if (ibs[k].RangeShifterSequence[0].RangeShifterNumber==1):	epRS=65 # Isogray tags
			#if (ibs[k].RangeShifterSequence[0].RangeShifterNumber==2):	epRS=30 # Isogray tags
			if (ibs[k].RangeShifterSequence[0].RangeShifterID=='RS30'):	epRS=30 # ECLIPSE tags
			if (ibs[k].RangeShifterSequence[0].RangeShifterID=='RS65'):	epRS=65 # ECLIPSE tags

			GantryAngle = icps[0].GantryAngle  					# Gantry angle
			PatientSupportAngle = icps[0].PatientSupportAngle	# Table angle
			DCI = 0.1*icps[0].SnoutPosition						# DCI
			SnoutID = ibs[k].SnoutSequence[0].SnoutID			# snout
			dSnout = -6.5
			if (SnoutID==400):	dSnout=1.8
			retraction = DCI + dSnout

			# Isocenter and Patient translation
			IsocenterPosition = icps[0].IsocenterPosition 	# in mm
			#IsocenterPosition = [0,0,0]

                        if(keepCTDim==False):   crop_RP = [wz_min.get()-1, wz_max.get(), wy_min.get()-1, wy_max.get(), wx_min.get()-1, wx_max.get()]
                        else:                   crop_RP = [0, dim_z,0, dim_y,0, dim_x]

                        TransX, TransY, TransZ = PatientTranslation(IsocenterPosition, PatientSupportAngle, crop_RP)

			#PatientOrientation = 180
			if (pss[0].PatientPosition=='HFP'):	TransY, TransZ = -TransY, -TransZ

			# Dose per fraction
			BeamDose, BeamMeterset = 0,0
			if ("FractionGroupSequence" in ds):
				BeamDose = ds.FractionGroupSequence[0].ReferencedBeamSequence[k].BeamDose 				# Dose per fraction
				BeamMeterset = int(ds.FractionGroupSequence[0].ReferencedBeamSequence[k].BeamMeterset) 	# UM per fraction

			#E = 220*np.ones(N)

                        if(check == True):
			        print '\n------ Beam',BeamName , '------'
			        print 'Number of spots :	', N
			        print "GantryAngle :		", GantryAngle, 'deg'
			        print "PatientSupportAngle :	", PatientSupportAngle, 'deg'
			        print "SnoutID :		", SnoutID
			        print "DCI :			{0:.1f} cm".format(DCI)
			        print "Range Shifter :		", epRS, 'mm'
			        #print "retraction :		", retraction, 'cm'
			        print "TransX :		{0:.1f} mm".format(TransX)
			        print "TransY :		{0:.1f} mm".format(TransY)
			        print "TransZ :		{0:.1f} mm".format(TransZ)
			        print 'Beam Dose :		', BeamDose, 'Gy'
			        print 'Beam Meterset :		', BeamMeterset, 'UM'
			        #print 'IsocenterPosition:       ', IsocenterPosition

			        check_PBS(BeamName,E,x,y,UM)

                        if(createTopasfile == True):
        			grid_size = 1. # scorer binning in mm
        			CreateTopasFilePBS(ID, BeamName, SnoutID, epRS, DCI, retraction, GantryAngle, E ,x ,y , UM, N, PatientSupportAngle, TransX, TransY, TransZ, dim_x, dim_y, dim_z, spacing, grid_size)

	RP_crop_show=True
	#try:	window_RP.destroy()
	#except Exception:	pass
	Update_all()

def PatientTranslation(iso, PatientSupportAngle, crop_RP):#=None):
	"""
	return the patient translation coordinates
	iso: 			isocenter coordinates (in mm)
	PatientSupportAngle: 	table angle (in deg)
	crop_RP: 		limits of the dose matrix (in voxel unit).
				Default is CT limits.
	"""

        #print 'keepCTDim: ', keepCTDim, crop_RP

	transX = iso[2] - 0.5*(crop_RP[4]+crop_RP[5])*spacing[0] - origin[0] 
	transY = iso[0] - 0.5*(crop_RP[0]+crop_RP[1]-dim_z)*spacing[2]	 
	transZ = iso[1] - 0.5*(crop_RP[2]+crop_RP[3]-dim_y)*spacing[1]

	TransX = transX*np.cos(np.deg2rad(PatientSupportAngle)) + transY*np.sin(np.deg2rad(PatientSupportAngle))
	TransY = -transX*np.sin(np.deg2rad(PatientSupportAngle)) + transY*np.cos(np.deg2rad(PatientSupportAngle))
	TransZ = -transZ

	return TransX, TransY, TransZ

def CreateTopasFilePBS(ID, BeamName, SnoutID, epRS, DCI, retraction, GantryAngle, E ,x ,y , UM, N, PatientSupportAngle, TransX, TransY, TransZ, dim_x, dim_y, dim_z, spacing, grid_size):
	"""
	write the main.txt TOPAS file
	"""
	with open('./output/Main_{0}_{1}.txt'.format(ID,BeamName), 'w') as f:

		f.write('######################################\n')
		f.write('# TOPAS MC v3.2 simulation file\n')
		f.write('# Patient IPP: {0}\n'.format(ID))
		f.write('# Beam name: {0}\n'.format(BeamName))
		f.write('# https://github.com/PierreLansonneur\n')
		f.write('######################################\n\n')

		f.write('### GENERAL COMMANDS\n')
		f.write('b:Ge/CheckForOverlaps 		= "1"\n')
		f.write('b:Ge/CheckForUnusedComponents 	= "0"\n')
		f.write('i:Ts/NumberOfThreads 		= 0\n')
		f.write('b:Ts/ShowCPUTime 		= "1"\n')
		f.write('i:Ts/MaxInterruptedHistories	= 1000000000\n')
		f.write('i:Ts/ShowHistoryCountAtInterval = 0 # reduce size of slurm file\n\n')

		f.write('### WORLD VOLUME\n')
		f.write('d:Ge/World/HLX 			= 300 cm\n')
		f.write('d:Ge/World/HLY 			= 300 cm\n')
		f.write('d:Ge/World/HLZ 			= 300 cm\n\n')

		f.write('#### GEOMETRY\n')
		f.write('includeFile 			= Geometry/Nozzle.txt\n')
		f.write('includeFile 			= Geometry/Snout{0:.0f}.txt\n'.format(0.1*float(SnoutID)))
		f.write('#includeFile 			= Geometry/Aperture_mbrt.txt\n')

		# RS 65 mm, 30 mm or nothing
		if (epRS == 65 or epRS == 30):
			f.write('b:Ge/RS/Include 		= "1"\n')
			f.write('d:Ge/RS/HL 			= {0:.1f} mm\n'.format(epRS/2.) )
			f.write('d:Ge/RS/TransZ			= Ge/Snout/TransIso + Ge/RS/HL cm\n')
		else:	f.write('b:Ge/RS/Include 		= "0"\n')

		f.write('d:Ge/Snout/TransZ 		= -{0} cm + Ge/Snout/AntiTransIso\n'.format(DCI) ) 
		f.write('d:Ge/SnoutGroup/TransZ		= {0} cm\n'.format(10.-retraction))
		f.write('d:Ge/FM/HLZ	  		= {0:.3f} cm\n'.format(26.225+0.5*(10.-retraction)))
		f.write('d:Ge/FM/TransZ   		= {0:.3f} cm\n\n'.format(-87.026+0.5*(10.-retraction)))

		f.write('### SuperGroup for Gantry rotation\n')
		f.write('s:Ge/SuperGroup/Parent 		= "World"\n')
		f.write('s:Ge/SuperGroup/Type 		= "Group"\n')
		f.write('d:Ge/SuperGroup/TransX 		= 0 m\n')
		f.write('d:Ge/SuperGroup/TransY 		= 0 m\n')
		f.write('d:Ge/SuperGroup/TransZ 		= 0 m\n')
		f.write('d:Ge/SuperGroup/RotX 		= {0} deg # Gantry angle\n'.format(GantryAngle))	
		f.write('d:Ge/SuperGroup/RotY 		= 0 deg\n')
		f.write('d:Ge/SuperGroup/RotZ 		= 0 deg\n\n')

		f.write('### PATIENT\n')
		f.write('s:Ge/Patient/ImagingToMaterialConverter 	= "Schneider"\n')
		f.write('includeFile 					= HUtoMAT.txt # WARNING check the scanner model!\n\n') 

		f.write('### Patient group\n')
		f.write('s:Ge/Patient/Message 				= "Constructing Patient"\n')
		f.write('s:Ge/Patient/Type 				= "Group"\n')
		f.write('s:Ge/Patient/Parent 				= "World"\n')
		f.write('d:Ge/Patient/TransX 				= {0:.1f} mm\n'.format(TransX))
		f.write('d:Ge/Patient/TransY 				= {0:.1f} mm\n'.format(TransY))
		f.write('d:Ge/Patient/TransZ 				= {0:.1f} mm\n'.format(TransZ)) 
		f.write('d:Ge/Patient/RotX 				= 90 deg\n')
		f.write('d:Ge/Patient/RotY 				= {0} deg # 90 deg - table angle\n'.format(90-PatientSupportAngle))
		f.write('d:Ge/Patient/RotZ 				= 180 deg\n\n')
                """
		f.write('### Patient STL file (for visualization only)\n')
		f.write('s:Ge/CT_visu/Type 				= "TsCAD"\n')
		f.write('s:Ge/CT_visu/Parent 				= "Patient"\n')
		f.write('b:Ge/CT_visu/IsParallel 			= "True"\n')
		f.write('s:Ge/CT_visu/Material 				= "Vacuum"\n')
		f.write('d:Ge/CT_visu/TransX 				= {0:.1f} mm\n'.format(-0.5*dim_z*spacing[2]))
		f.write('d:Ge/CT_visu/TransY 				= {0:.1f} mm\n'.format(-0.5*dim_y*spacing[1]))
		f.write('d:Ge/CT_visu/TransZ 				= {0:.1f} mm\n'.format(-0.5*dim_x*spacing[0]))
		f.write('s:Ge/CT_visu/DrawingStyle 			= "Wireframe"\n')
		f.write('s:Ge/CT_visu/InputFile 				= "Geometry/{0}"\n'.format(ID))
		f.write('s:Ge/CT_visu/FileFormat 			= "stl"\n')
		f.write('d:Ge/CT_visu/Units 				= 1.0 mm\n\n')
                """
		f.write('### Patient modelisation\n')
		f.write('s:Ge/Patient/CT/Parent 				= "Patient"\n')
		f.write('s:Ge/Patient/CT/Type 				= "TsDicomPatient"\n')
		#f.write('b:Ge/Patient/CT/IsParallel 			= "{0}"\n'.format(keepCTDim))
		f.write('s:Ge/Patient/CT/Material 			= "G4_WATER"\n')
		f.write('s:Ge/Patient/CT/DicomDirectory 			= "./CT_{0}"\n'.format(ID))	###
		f.write('sv:Ge/Patient/CT/DicomModalityTags 		= 1 "CT"\n')
		f.write('d:Ge/Patient/CT/TransX 				= 0 mm\n')
		f.write('d:Ge/Patient/CT/TransY 				= 0 mm\n')
		f.write('d:Ge/Patient/CT/TransZ 				= 0 mm\n')
		f.write('d:Ge/Patient/CT/RotX 				= 0 deg\n')
		f.write('d:Ge/Patient/CT/RotY 				= 0 deg\n')
		f.write('d:Ge/Patient/CT/RotZ 				= 0 deg\n\n')
		#f.write('#sv:Ge/Patient/CT/ColorByRTStructNames 		= 2 "CONTOUR_EXTERNE" "PTV_1"\n')
		#f.write('#sv:Ge/Patient/CT/ColorByRTStructColors 		= 2 "green" "red"\n')
		#f.write('#b:Ge/Patient/CT/IgnoreInconsistentFrameOfReferenceUID = "True"\n\n')

	        f.write('### Crop CT dimensions\n')
	        f.write('i:Ge/Patient/CT/RestrictVoxelsXMin = {0}\n'.format(wz_min.get()))
	        f.write('i:Ge/Patient/CT/RestrictVoxelsXMax = {0}\n'.format(wz_max.get()))
	        f.write('i:Ge/Patient/CT/RestrictVoxelsYMin = {0}\n'.format(wy_min.get()))
	        f.write('i:Ge/Patient/CT/RestrictVoxelsYMax = {0}\n'.format(wy_max.get()))
	        f.write('i:Ge/Patient/CT/RestrictVoxelsZMin = {0}\n'.format(wx_min.get()))
	        f.write('i:Ge/Patient/CT/RestrictVoxelsZMax = {0}\n\n'.format(wx_max.get()))

		f.write('### Dose computation\n')
		f.write('s:Sc/PatientScorer/Component 			= "Patient/CT"\n')
		f.write('s:Sc/PatientScorer/Quantity 			= "DoseToWater"\n')
		f.write('s:Sc/PatientScorer/OutputFile 			= "./output/Dose_{0}_{1}"\n'.format(ID,BeamName))
		f.write('s:Sc/PatientScorer/OutputType 			= "DICOM"\n')
		f.write('b:Sc/PatientScorer/OutputToConsole 		= "0"\n')
		f.write('s:Sc/PatientScorer/IfOutputFileAlreadyExists 	= "Increment"\n')
		f.write('i:Sc/PatientScorer/XBins 			= {0:.0f}\n'.format(np.round((1.0/grid_size)*(1+wz_max.get()-wz_min.get())*spacing[2])))
		f.write('i:Sc/PatientScorer/YBins 			= {0:.0f}\n'.format(np.round((1.0/grid_size)*(1+wy_max.get()-wy_min.get())*spacing[1])))
		f.write('i:Sc/PatientScorer/ZBins 			= {0:.0f}\n\n'.format(np.round((1.0/grid_size)*(1+wx_max.get()-wx_min.get())*spacing[0])))
		#f.write('#sv:Sc/PatientScorer/OnlyIncludeIfInRTStructure 	= 1 "CONTOUR_EXTERNE"\n')
		#f.write('#b:Sc/PatientScorer/OutputAfterRun 		= True\n\n')

		f.write('### PHYSICS (TOPAS default settings)\n')
		f.write('s:Ph/ListName 			= "Default"\n')
		f.write('s:Ph/Default/Type 		= "Geant4_Modular"\n')
		f.write('sv:Ph/Default/Modules 		= 7 "g4em-standard_opt3" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping" "g4radioactivedecay"\n')
		f.write('d:Ph/Default/EMRangeMin 	= 100 eV\n')
		f.write('d:Ph/Default/EMRangeMax 	= 230 MeV\n')
		f.write('i:Ph/Default/EMBinsPerDecade 	= 100\n')
		f.write('d:Ph/Default/CutForAllParticles = 0.05 mm\n')
		f.write('d:Ge/SM1/Dipole/MaxStepSize 	= 1 mm\n')
		f.write('d:Ge/SM2/Dipole/MaxStepSize 	= 1 mm\n\n')

                if (0.1*float(SnoutID)==10):	LMGW = ["Snout/BrassCone2","Snout/BrassTube2","Snout/BrassCone3"]
                if (0.1*float(SnoutID)==18):	LMGW = ["Snout/BrassCone2","Snout/BrassCone3"]
                if (0.1*float(SnoutID)==25):	LMGW = ["Snout/AlBrassTube2", "Snout/BrassCone2", "Snout/BrassCone3", "Snout/BrassCone4"]
                if (0.1*float(SnoutID)==40):	LMGW = ["Snout/SquarePart2", "Snout/SquarePart3"]
                if (epRS==65 or epRS==30):      LMGW = np.append(LMGW,["RS"])
                f.write('sv:Ph/Default/LayeredMassGeometryWorlds = {0} '.format(len(LMGW)))
	        for item in LMGW:	f.write('"{0}" '.format(item))
	        f.write('\n\n')
                '''
		f.write('### SPEED-UP the simulation\n')
		f.write('#b:Vr/UseVarianceReduction 				= "True"\n')
		f.write('#b:Vr/KillOtherParticles/Active 			= "True"\n')
		f.write('#sv:Vr/KillOtherParticles/HaveNoEffectInComponentsNamed= 1 "Patient/CT"\n')
		f.write('#sv:Vr/KillOtherParticles/OnlyTrackParticlesNamed 	= 1 "proton"\n\n')
                '''

		f.write('### Time Features Parameters\n')
		f.write('i:Tf/NumberOfSequentialTimes 	= {0}\n'.format(N)) 
		f.write('d:Tf/TimelineStart 		= 0 s\n')
		f.write('d:Tf/TimelineEnd 		= {0} s\n'.format(N+1)) 
		f.write('#####\n\n')

		f.write('### Spot Energy\n')
		f.write('s:Tf/Energy/Function 		= "Step"\n')
		f.write('dv:Tf/Energy/Times 		= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('dv:Tf/Energy/Values 		= {0} '.format(N))
		for i in E:	f.write('{0:.2f} '.format(i))
		f.write('MeV\n')
		f.write('#####\n\n')

		f.write('##### Energy Spread\n')
		f.write('s:Tf/EnergySpread/Function 	= "Step"\n')
		f.write('dv:Tf/EnergySpread/Times 	= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('uv:Tf/EnergySpread/Values 	= {0} '.format(N))
		for i in poly_Espread(E):	f.write('{0:.2f} '.format(i))
		f.write('\n')
		f.write('#####\n\n')

		f.write('##### Spot Sigmas\n')
		f.write('s:Tf/SpotSize/Function 		= "Step"\n')
		f.write('dv:Tf/SpotSize/Times 		= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('dv:Tf/SpotSize/Values 		= {0} '.format(N))
		for i in poly_width(E):	f.write('{0:.2f} '.format(i))
		f.write('mm\n')
		f.write('#####\n\n')

		f.write('##### Spot Divergences\n')
		f.write('s:Tf/BeamDivergence/Function 	= "Step"\n')
		f.write('dv:Tf/BeamDivergence/Times 	= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('uv:Tf/BeamDivergence/Values 	= {0} '.format(N))
		for i in poly_div(E):	f.write('{0:.4f} '.format(i))
		f.write('\n')
		f.write('#####\n\n')

		f.write('##### Field for SM1 scanning magnet\n')
		f.write('s:Tf/SM1Field/Function 			= "Step"\n')
		f.write('dv:Tf/SM1Field/Times 			= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('dv:Tf/SM1Field/Values 			= {0} '.format(N))
		for i in 0.1*y/poly_SM1(E):	f.write('{0:.4f} '.format(i))
		f.write('tesla\n')
		f.write('d:Ge/SM1/Dipole/MagneticFieldStrength 	= Tf/SM1Field/Value tesla\n')
		f.write('u:Ge/SM1/Dipole/MagneticFieldDirectionX	= 0.0\n')
		f.write('u:Ge/SM1/Dipole/MagneticFieldDirectionY	= {0:.4f}\n'.format(np.cos(np.deg2rad(GantryAngle))))
		f.write('u:Ge/SM1/Dipole/MagneticFieldDirectionZ	= {0:.4f}\n'.format(-np.sin(np.deg2rad(GantryAngle))))
		f.write('#####\n\n')

		f.write('##### Field for SM2 scanning magnet\n')
		f.write('s:Tf/SM2Field/Function 			= "Step"\n')
		f.write('dv:Tf/SM2Field/Times 			= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('dv:Tf/SM2Field/Values 			= {0} '.format(N))
		for i in 0.1*x/poly_SM2(E):	f.write('{0:.4f} '.format(i))
		f.write('tesla\n')
		f.write('d:Ge/SM2/Dipole/MagneticFieldStrength 	= Tf/SM2Field/Value tesla\n')
		f.write('u:Ge/SM2/Dipole/MagneticFieldDirectionX	= 1.0\n')
		f.write('u:Ge/SM2/Dipole/MagneticFieldDirectionY	= 0.0\n')
		f.write('u:Ge/SM2/Dipole/MagneticFieldDirectionZ	= 0.0\n')
		f.write('#####\n\n')

		f.write('##### Number of Histories\n')
		f.write('s:Tf/Particles/Function 	= "Step"\n')
		f.write('dv:Tf/Particles/Times 		= {0} '.format(N))
		for i in range(1,N+1):	f.write('{0} '.format(i))
		f.write('s\n')
		f.write('iv:Tf/Particles/Values 		= {0} '.format(N))
		for i in (0.1*UM/poly_UM(E)).astype(int):	f.write('{0:.0f} '.format(i)) # histories in Gigaprotons

		f.write('\n')
		f.write('#####\n\n')

		f.write('### SOURCE\n')
		f.write('s:So/Default/Type 		= "emittance"\n')
		f.write('s:So/Default/Component 		= "VW"\n')
		f.write('s:So/Default/BeamParticle 	= "proton"\n')
		f.write('d:So/Default/BeamEnergy 	= Tf/Energy/Value MeV\n')
		f.write('u:So/Default/BeamEnergySpread 	= Tf/EnergySpread/Value\n')
		f.write('s:So/Default/Distribution 	= "BiGaussian"\n')
		f.write('d:So/Default/SigmaX 		= Tf/SpotSize/Value mm\n')
		f.write('u:So/Default/SigmaXprime 	= Tf/BeamDivergence/Value\n')
		f.write('u:So/Default/CorrelationX 	= -0.882\n')
		f.write('d:So/Default/SigmaY 		= Tf/SpotSize/Value mm\n')
		f.write('u:So/Default/SigmaYprime 	= Tf/BeamDivergence/Value\n')
		f.write('u:So/Default/CorrelationY 	= -0.882\n')
		f.write('i:So/Default/NumberOfHistoriesInRun = 1 * Tf/Particles/Value\n\n')

		f.write('### GRAPHICS\n')
		f.write('b:Gr/Enable 			= "0"\n')
		f.write('s:Gr/ViewA/Type 		= "VRML"\n')
		f.write('b:Ge/World/Invisible 		= "1"\n')
		#f.write('#b:Ge/Patient/CT/Invisible 	= "1"\n')
		#f.write('#b:Gr/ViewA/IncludeAxes 	= "0"\n')

	print 'file {0} saved succesfully!'.format('./output/Main_{0}_{1}.txt'.format(ID,BeamName))

def CreateApertureFileDS(BeamName,BlockData):
	"""
	write the .ap file for the collimator geometry
	"""
	BlockDataX = BlockData[0::2]
	BlockDataY = BlockData[1::2]

	with open('./output/'+BeamName+'.ap', 'w') as f:
		f.write('{0}\n'.format(len(BlockDataX)))
		for i in range(len(BlockDataX)):	f.write('{0}, {1}\n'.format(BlockDataX[i],BlockDataY[i]))

        print('collimator file {0} created'.format('./output/'+BeamName+'.ap'))

def CreateCompensatorFileDS(BeamName,ircs):
	"""
	write the .rc file for the compensator geometry
	"""
	CRows,CColumns = ircs[0].CompensatorRows, ircs[0].CompensatorColumns
	CompensatorMillingToolDiameter = ircs[0].CompensatorMillingToolDiameter
	CompensatorThicknessData = np.array(ircs[0].CompensatorThicknessData).reshape((CRows,CColumns))
	CPixelSpacing = ircs[0].CompensatorPixelSpacing 
	CPosition = ircs[0].CompensatorPosition

	CompensatorThicknessData_small = CompensatorThicknessData[0::10,0::10].T
	CRows, CColumns = np.shape(CompensatorThicknessData_small)[0], np.shape(CompensatorThicknessData_small)[1]

	with open('./output/'+BeamName+'.rc', 'w') as f:
		f.write('{0}\n'.format(CRows))
		f.write('{0}\n'.format(np.nanmax(CompensatorThicknessData_small)))
		f.write('{0}\n'.format(CompensatorMillingToolDiameter))

		for i in range(CRows):
			f.write('{0} {1} {2} {3}\n'.format(CColumns, -10*CPixelSpacing[1], CPosition[1],  -CPosition[0]-i*10*CPixelSpacing[0]))
			for j in range(CColumns):	f.write( '{0} '.format( np.nanmax(CompensatorThicknessData_small) - CompensatorThicknessData_small[CRows-i-1,j]) )
			f.write('\n')

        print('compensator file {0} created'.format('./output/'+BeamName+'.rc'))

def check_PBS(BeamName,E,x,y,UM):
	"""
	graphical check for PBS beams
	"""
	window_RPcheck = Toplevel(window_RP)
	window_RPcheck.title(BeamName)
	fig2 = P.figure(facecolor='lightgrey')
	graph2 = FigureCanvasTkAgg(fig2, master=window_RPcheck)
	canvas2 = graph2.get_tk_widget()
	canvas2.grid(row=0, columnspan=3, column=0)

	disp4 = Label(window_RPcheck,  text=' Number of spots:\t{0}'.format(N))
	disp5 = Label(window_RPcheck,  text=' Gantry Angle : \t{0} deg'.format(GantryAngle))
	disp6 = Label(window_RPcheck,  text=' Table Angle :  \t{0} deg'.format(PatientSupportAngle))
	disp7 = Label(window_RPcheck,  text=' Snout ID :     \t{0}'.format(SnoutID))
	disp8 = Label(window_RPcheck,  text=' DCI :          \t{0:.1f} cm'.format(DCI))
	disp9 = Label(window_RPcheck,  text=' Range Shifter :\t{0:.1f} mm'.format(epRS))
	disp10 = Label(window_RPcheck, text=' TransX :       \t{0:.1f} mm'.format(TransX))
	disp11 = Label(window_RPcheck, text=' TransY :       \t{0:.1f} mm'.format(TransY))
	disp12 = Label(window_RPcheck, text=' TransZ :       \t{0:.1f} mm'.format(TransZ))
	disp13 = Label(window_RPcheck, text=' Beam Dose :    \t{0} Gy'.format(BeamDose))

	disp4.grid(row= 1, column=1, sticky=W)
	disp5.grid(row= 2, column=1, sticky=W)
	disp6.grid(row= 3, column=1, sticky=W)
	disp7.grid(row= 4, column=1, sticky=W)
	disp8.grid(row= 5, column=1, sticky=W)
	disp9.grid(row= 1, column=0, sticky=W)
	disp10.grid(row=2, column=0, sticky=W)
	disp11.grid(row=3, column=0, sticky=W)
	disp12.grid(row=4, column=0, sticky=W)
	disp13.grid(row=5, column=0, sticky=W)

	ax7 = fig2.gca()
	ax7.axis('equal')
	sctr = ax7.scatter(x, y, c=E, cmap='RdYlGn')
	#sctr = ax7.scatter(x, y, c=UM/poly_UM(E), cmap='RdYlGn')
	ax7.set_xlabel('x (mm)')
	ax7.set_ylabel('y (mm)')
	ax7.set_title('Beam-Eye View')
	P.colorbar(sctr, ax=ax7, format='%d MeV')
	P.tight_layout()
	#P.savefig(BeamName+'.pdf',  bbox_inches='tight', dpi=300)
	graph2.draw()

def check_DS(BeamName, BlockData, CompensatorThicknessData, ext):
	"""
	graphical check for DS beams
	"""
	window_RP = Toplevel(window)
	window_RP.title(BeamName)
	fig2 = P.figure(facecolor='lightgrey')
	graph2 = FigureCanvasTkAgg(fig2, master=window_RP)
	canvas2 = graph2.get_tk_widget()
	canvas2.grid(row=0, column=0)
	ax7 = fig2.gca()
	ax7.axis('equal')

	### compensator
	#comp = ax7.imshow(CompensatorThicknessData[0::10,0::10], extent=ext)
	comp = ax7.imshow(CompensatorThicknessData, extent=ext, cmap='jet')
	ax7.set_xlabel('x (mm)')
	ax7.set_ylabel('y (mm)')
	ax7.set_title('Beam-Eye View')

	### collimator
	BlockDataX = BlockData[0::2]
	BlockDataY = BlockData[1::2]
	vertices = np.append(BlockDataX,BlockDataY).reshape(2,-1).T	# get the vertices coordinates of the contour,
	path = Path(vertices,closed=True)
	patch = patches.PathPatch(path, facecolor='none', edgecolor='r', lw = 1)
	ax7.add_patch(patch)

	P.colorbar(comp, ax=ax7, format='%d mm')
	P.tight_layout()
	#P.savefig(BeamName+'.pdf',  bbox_inches='tight', dpi=300)
	graph2.draw()
