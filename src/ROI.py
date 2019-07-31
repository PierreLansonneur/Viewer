# -*- coding: utf-8 -*-
###################################################
#   	      Import RT-STRUCT files
# 		   	----
#           	P. Lansonneur 2019
###################################################

def open_ROI(filename=None): 
	"""
	Open a RT-struct file with .dcm extension
	"""
	global ROI_open, ROI_show, ds_ROI, ROI_infos, N_ROI, filename_ROI

	if (CT_open == False): 
                tkMessageBox.showerror('Error', 'You must open a CT first')
                return

	else:
        	types = [('RS files', '*.dcm')]
		if(filename==None):	file_path = tkFileDialog.askopenfilename(initialdir = dir_ini, filetypes = types)
		else:	file_path = filename

		filename_ROI = file_path
			
		print 'Opening RS file ... '
		ds_ROI = pydicom.read_file(file_path)
		N_ROI = len(ds_ROI.StructureSetROISequence) 	# number of ROI
		ROI_infos = np.empty((N_ROI,7), dtype=np.object)

		for ROI_index in range(N_ROI):	
			Get_ROI_infos(ds_ROI,ROI_index) # retrieve the info for all ROI
			isPTV = ROI_infos[ROI_index,0].upper().startswith('PTV')
			isCTV = ROI_infos[ROI_index,0].upper().startswith('CTV')
			isGTV = ROI_infos[ROI_index,0].upper().startswith('GTV')
			if(isPTV or isCTV or isGTV):	
                                ROI_infos[ROI_index,6] = ROI_3D(ROI_infos, ROI_index)
                                tmp_vol = 0.001*np.sum(ROI_infos[ROI_index,6])*spacing[0]*spacing[1]*spacing[2] # volume in cm3
                               	print ' found target volume: {0}\t({1:.1f} cm3)'.format(ROI_infos[ROI_index,0],tmp_vol)

		#ROI_infos[0,6] = ROI_3D(ROI_infos, 0) # contour externe

		print 'File successfully opened!'
		ROI_open = True
		ROI_show = True
		check2.select()
		Update()

                roimenu.entryconfig("Compute DVH...", state="normal")
                roimenu.entryconfig("Crop dosimetry", state="normal")
                roimenu.entryconfig("Normalize dosi to PTV dose", state="normal")

def Get_ROI_infos(ds, i):
	"""
	Get informations for the ROI 'i'.
	ROI informations are saved into the array ROI_infos[i,:]
	"""
	global ROI_infos
	ctrs = ds.ROIContourSequence
	x, y, z, N_vert, N_vert_cumul = [],[],[],[],[]	 
	x, y, z = np.array([]).astype(float), np.array([]).astype(float), np.array([]).astype(float)

	if('ContourSequence' in ctrs[i]):
		for j in range(len(ctrs[i].ContourSequence)):
			x = np.append(x,np.asfarray(ctrs[i].ContourSequence[j].ContourData[0::3],float))
			y = np.append(y,np.asfarray(ctrs[i].ContourSequence[j].ContourData[1::3],float))
			z = np.append(z,np.asfarray(ctrs[i].ContourSequence[j].ContourData[2::3],float))
			N_vert = np.append(N_vert,ctrs[i].ContourSequence[j].NumberOfContourPoints) 	 # number of vertices for each contour
			N_vert_cumul = np.append(N_vert_cumul,np.sum(N_vert[0:j]))

	color = ds_ROI.ROIContourSequence[i].ROIDisplayColor
	color = np.array([color[0]/255.,color[1]/255.,color[2]/255.])				 

	N_vert_cumul = np.array(N_vert_cumul).astype(int)

	ROI_infos[i,0] = ds_ROI.StructureSetROISequence[i].ROIName	# ROI Name
	ROI_infos[i,1] = x						# vertices x coordinates in mm
	ROI_infos[i,2] = y						# vertices y coordinates in mm
	ROI_infos[i,3] = z - origin[0]			                # vertices z coordinates in mm
	ROI_infos[i,4] = N_vert_cumul			                # cumulated number of vertices for each contour
	ROI_infos[i,5] = color					        # color in RGB format
	ROI_infos[i,6] = 0						# 0-1 array with same dimension as the CT volume.

def UpdateROI(array, slice_ax1, slice_ax2, slice_ax3,ext8,ext9):
	"""
	display the ROI contours
	"""
	global ax1, ax2, ax3

	ROI_Name = array[0]
	x = array[1]
	y = array[2]
	z = array[3]
	N_vert_cumul = array[4]
	ROI_color = array[5]
	ROI_volume = array[6]

	isPTV = ROI_Name.upper().startswith('PTV')
	isCTV = ROI_Name.upper().startswith('CTV')
	isGTV = ROI_Name.upper().startswith('GTV')

	if (rot1==True): 	y,x = x,y

	z_s = (z/spacing[0]).astype(int) # scale vertices coordinates to voxel unit

	if slice_ax1 in z_s: 
		index = np.where(z_s == slice_ax1) # indices of the slice concerned

		if(len(N_vert_cumul)>1): ### surface contours
			for count in range(len(index[0])):
				if (index[0][count] in N_vert_cumul):	
					mini = N_vert_cumul[np.where(N_vert_cumul==index[0][count])[0][0]-1]
					maxi = N_vert_cumul[np.where(N_vert_cumul==index[0][count])[0][0]]	# get 1st and last indices of the contour,
					if (maxi>mini):
						vertices = np.append(y[mini:maxi],x[mini:maxi]).reshape(2,-1).T	# get the vertices coordinates of the contour,
						path = Path(vertices,closed=True)								# make them a path,
						if(isPTV or isCTV or isGTV):	patch = patches.PathPatch(path, facecolor='none', edgecolor=ROI_color, lw = 1,zorder=3)
						else:	patch = patches.PathPatch(path, facecolor='none', edgecolor=ROI_color, lw = 1)
						ax1.add_patch(patch)

		elif(len(N_vert_cumul)==1):	ax1.plot([y[0]],[x[0]],color=ROI_color,marker = '*',zorder=3) ### points

	if isPTV or isCTV or isGTV:
		roi2 = (ROI_volume)[:,dim_y-1-slice_ax2,:].T
		roi3 = (ROI_volume)[:,:,dim_z-1-slice_ax3].T
		if rot2:	roi2 = np.fliplr(np.flipud(np.rot90(roi2)))
		if rot3:	roi3 = np.fliplr(np.flipud(np.rot90(roi3)))
		ax2.contour(roi2, 0, colors=[ROI_color], linewidths=1, extent=ext8)
		ax3.contour(roi3, 0, colors=[ROI_color], linewidths=1, extent=ext9)

def ROI_DVH_analysis():
	"""
	display the cumulative DVH of every ROI loaded
	"""
	global ax8, graph3, fig3, DVHData

	window_ROI_DVH = Toplevel(window)
	window_ROI_DVH.geometry("{0}x{1}+{2}+{3}".format(900, 590, w_px, 0))
	window_ROI_DVH.title('DVH')

	fig3 = P.figure(facecolor='lightgrey')
	fig3.set_size_inches(9,5.5)
	graph3 = FigureCanvasTkAgg(fig3, master=window_ROI_DVH)
	canvas3 = graph3.get_tk_widget()
	toolbar_frame3 = Frame(window_ROI_DVH)
	toolbar3 = NavigationToolbar2TkAgg(graph3, toolbar_frame3)
	ax8 = fig3.gca()
        P.subplots_adjust(left=0.08,right=0.8, bottom=0.09, top=0.96)
	roi_bt = Button(window_ROI_DVH, text="Export..", command=ExportDVH)

	canvas3.grid(row=3, column=0, columnspan=5, sticky=W+E+N+S)
	toolbar_frame3.grid(row=4,column=0, columnspan=4, sticky=W+E+N+S)
	roi_bt.grid(row=4,column=4,sticky=E+N+S)

        DVHData= np.empty(N_ROI, dtype=np.object)

        print '\nCalculating DVH..'
        print 'Structure        \tmean dose (Gy)\tPVDR (error)\ts-index\tc-index\th-index'

	for ROI_index in range(1,N_ROI):

		ROIName = ROI_infos[ROI_index,0].decode('utf-8')
		isPTV = ROIName.upper().startswith('PTV')
		isCTV = ROIName.upper().startswith('CTV')
		isGTV = ROIName.upper().startswith('GTV')

		DVH = ROI_DVH(ROI_infos[ROI_index,:])

		if(len(DVH)>0)and(DVH.max()>0):

			ROI_color = ROI_infos[ROI_index,5]

			# cumulative DVH
			DVHData[ROI_index] = ax8.hist(DVH, bins='auto', range=(0,DVH.max()), histtype='step', cumulative=-1, density=True, label=ROIName, color=ROI_color) 

			# differential DVH
			#DVHData[ROI_index] = ax8.hist(DVH, bins='auto', range=(0,DVH.max()), histtype='step', density=True, label=ROIName, color=ROI_color)

			ROI_color = 255*ROI_infos[ROI_index,5].astype(float)
			ROI_color = ROI_color.astype(int)
                        s_index = c_index = h_index = np.nan

                        volume_ROI = ROI_3D(ROI_infos, ROI_index, tag='CROP')
                        PVDR_ = PVDR(dosi[volume_ROI], 0.5 )

	                if( isPTV or isCTV or isGTV):

		                DVH_norm = DVH/np.mean(DVH) # normalize mean dose to 1
                                s_index = np.std(DVH_norm)
                                c_index_tmp = len( DVH[np.where(DVH_norm>=0.95)])
                                c_index_tmp = 0.001*c_index_tmp*spacing_dosi[0]*spacing_dosi[1]*spacing_dosi[2]
                                c_index = c_index_tmp/(0.001*np.sum(ROI_infos[ROI_index,6])*spacing[0]*spacing[1]*spacing[2])
                                h_index = (np.percentile(DVH,98)-np.percentile(DVH,2))/np.percentile(DVH,50)

                        print '{0}\t{1:.1f}\t\t{2:.2f} ({3:.2f})\t{4:.2f}\t{5:.2f}\t{6:.2f}'.format(ROIName.ljust(17), np.mean(DVH), PVDR_[0], PVDR_[1], s_index, c_index, h_index)

	ax8.set_ylabel('Volume')
	ax8.set_xlabel('Dose (Gy)')
        ax8.legend(bbox_to_anchor=(1.02,1),loc='upper left', frameon=False)
	graph3.draw()
	toolbar3.draw()
	Update()

def ROI_DVH_PTV():

        print 'structure\t\tPVDR\t\ts_index (%)'

        for ROI_index in range(1,N_ROI):

	        ROIName = ROI_infos[ROI_index,0].decode('utf-8')
	        isPTV = ROIName.upper().startswith('PTV')
	        isCTV = ROIName.upper().startswith('CTV')
	        isGTV = ROIName.upper().startswith('GTV')

	        DVH = ROI_DVH(ROI_infos[ROI_index,:])
	        #if(len(DVH)>0)and(DVH.max()>0)and( isPTV ):
	        if(len(DVH)>0)and(DVH.max()>0)and( isGTV ):


                        volume_ROI = ROI_3D(ROI_infos, ROI_index, tag='CROP')
                        PVDR_ = PVDR(dosi[volume_ROI], 0.005*np.mean(DVH) )	                

                        DVH_norm = DVH/np.mean(DVH) # normalize mean dose to 1
                        s_index = np.std(DVH_norm)
                        c_index_tmp = len( DVH[np.where(DVH_norm>=0.95)])
                        c_index_tmp = 0.001*c_index_tmp*spacing_dosi[0]*spacing_dosi[1]*spacing_dosi[2]
                        c_index = c_index_tmp/(0.001*np.sum(ROI_infos[ROI_index,6])*spacing[0]*spacing[1]*spacing[2])
                        h_index = (np.percentile(DVH,98)-np.percentile(DVH,2))/np.percentile(DVH,50)

                        print '{0}\t{1:.2f} Â± {2:.2f}\t{3:.1f}'.format(ROIName.ljust(17),  PVDR_[0], PVDR_[1], 100*s_index)

def ExportDVH():
	"""
	Export DVH data
	"""
	types = [('text', '*.txt')]
	filename = tkFileDialog.asksaveasfilename(initialdir = dir_ini+'output', initialfile='DVH',filetypes=types)

        with open(filename, 'w') as f:
                f.write('Dose, Volume (%)\n')

                for ROI_index in range(1,N_ROI):
		        ROIName = ROI_infos[ROI_index,0].decode('utf-8')
                        try:
                                f.write('\nStructure: {0}\n'.format(ROIName))
                                for i,j in zip(DVHData[ROI_index][1], DVHData[ROI_index][0]):   f.write('{0:.3f}, {1:.2f}\n'.format(i,100*j))
                        except Exception:	pass

        print('DVH saved !')

def ROI_DVH(array,tag='dosi'): 
	"""
	return the dose of each voxels inside a given ROI
	the output is a np.array
	
	Usage:
	>>> ROI_DVH(ROI_infos[ROI_index,:])
	"""
	DVH = []
	ROIName = array[0].decode('utf-8')
	x = origin_dosi[2] - array[1]
	y = origin_dosi[1] - array[2]
	z = array[3]
	N_vert_cumul = array[4]
	ROI_color = array[5]

        coordinates = np.dstack(np.meshgrid(np.arange(dim_y_dosi), np.arange(dim_z_dosi))).reshape(-1,2)

        # scale vertices coordinates to voxel unit
        x_s = (x/spacing_dosi[2]).astype(int) 
        y_s = (y/spacing_dosi[1]).astype(int)
        z_s = (z/spacing[0]).astype(int)

        if(tag=='CT'):
	        coordinates = np.dstack(np.meshgrid(np.arange(dim_y), np.arange(dim_z))).reshape(-1,2)
	        # scale vertices coordinates to voxel unit
	        x_s = (x/spacing[2]).astype(int) 
	        y_s = (y/spacing[1]).astype(int)
	        z_s = (z/spacing[0]).astype(int)

	for slice_ in range(1,dim_x):
		slice_dosi = slice_ax1_dosi(slice_)
		if slice_ in z_s: 
			index = np.where(z_s == slice_) # indices of the slice concerned

			if(len(N_vert_cumul)>1): ### surface contours
				for count in range(len(index[0])):
					if (index[0][count] in N_vert_cumul):	
						mini = N_vert_cumul[np.where(N_vert_cumul==index[0][count])[0][0]-1]
						maxi = N_vert_cumul[np.where(N_vert_cumul==index[0][count])[0][0]]	# get 1st and last indices of the contour,
						if (maxi>mini):
							vertices = np.append(y_s[mini:maxi],x_s[mini:maxi]).reshape(2,-1).T	# get the vertices coordinates of the contour,
							path = Path(vertices,closed=True)				# make them a path,

							# dose-volume histogram
                                                        if(tag=='CT'):
							        within = coordinates[path.contains_points(coordinates)]
							        DVH = np.append(DVH,volume[slice_-1, within[:,0], within[:,1]])

							elif(tag=='dosi')and(slice_dosi != int((slice_-1)*spacing[0]/spacing_dosi[0]))and(slice_dosi<dim_x_dosi):
								within = coordinates[path.contains_points(coordinates)]
								DVH = np.append(DVH,dosi[slice_dosi-1, within[:,0], within[:,1]])
								#dosi[slice_dosi, within[:,0], within[:,1]] = 100 # check
	return np.array(DVH)

def ROI_3D(ROI_infos, index, tag=None): 
	"""
	return a 3D binary array with same dimensions than the CT volume
	The voxel is set to 1 if it belongs to the ROI, 0 otherwise.
	
	Usage:
	>>> ROI_3D(ROI_infos, 0)
	[[[0. 0. ... 1. 1. 0.]]]
	"""
	array = ROI_infos[index,:]

	ROIName = array[0].decode('utf-8')
	N_vert_cumul = array[4]
	ROI_color = array[6]

	# scale vertices coordinates to voxel unit
	x_s = ((array[1] + origin[2])/spacing[2]).astype(int) 
	y_s = ((array[2] + origin[1])/spacing[1]).astype(int)
	z_s = (array[3]/spacing[0]).astype(int)

	coordinates = np.dstack(np.meshgrid(np.arange(dim_y), np.arange(dim_z))).reshape(-1,2)
	volume_ROI = np.zeros_like(volume, dtype=bool)
	tmp_ = dim_x

	if tag=='CROP':
		x_s = ((-array[1] + origin_dosi[2])/spacing_dosi[2]).astype(int) 
		y_s = ((-array[2] + origin_dosi[1])/spacing_dosi[1]).astype(int)
		z_s = ((array[3] + origin[0] - origin_dosi[0])/spacing_dosi[0]).astype(int)
		coordinates = np.dstack(np.meshgrid(np.arange(dim_y_dosi), np.arange(dim_z_dosi))).reshape(-1,2)
		volume_ROI = np.zeros_like(dosi, dtype=bool)
		tmp_ = dim_x_dosi

	#print 'dim_x',dim_x
	for slice_ in range(1,tmp_):
		if slice_ in z_s: 
			slice_index = np.where(z_s == slice_) # indices of the slice concerned
			if(len(N_vert_cumul)>1): ### surface contours
				for count in range(len(slice_index[0])):
					if (slice_index[0][count] in N_vert_cumul):	
						mini = N_vert_cumul[np.where(N_vert_cumul==slice_index[0][count])[0][0]-1]
						maxi = N_vert_cumul[np.where(N_vert_cumul==slice_index[0][count])[0][0]]	# get 1st and last indices of the contour,
						if (maxi>mini):
							vertices = np.append(y_s[mini:maxi],x_s[mini:maxi]).reshape(2,-1).T	# get the vertices coordinates of the contour,
							path = Path(vertices,closed=True)					# make them a path,
							within = coordinates[path.contains_points(coordinates)]
							#dosi[slice_-1, within[:,0], within[:,1]] = 100
							volume_ROI[slice_-1, within[:,0], within[:,1]] = True

	# quick fix of empty slices (bug)
	for slice_ in range(2,tmp_-2):
		if (np.sum(volume_ROI[slice_+2,:,:])>0)and(np.sum(volume_ROI[slice_-2,:,:])>0):
			if ( 1.*np.sum(volume_ROI[slice_,:,:])/np.sum(volume_ROI[slice_-2,:,:])<=0.05 ):
				if( 1.*np.sum(volume_ROI[slice_,:,:])/np.sum(volume_ROI[slice_+2,:,:])<=0.05 ):	
					volume_ROI[slice_,:,:] = volume_ROI[slice_+2,:,:]
					#print slice_

	# quick fix of empty slices (bug)
	for slice_ in range(1,tmp_-1):
		if (np.sum(volume_ROI[slice_+1,:,:])>0)and(np.sum(volume_ROI[slice_-1,:,:])>0):
			if ( 1.*np.sum(volume_ROI[slice_,:,:])/np.sum(volume_ROI[slice_-1,:,:])<=0.05 ):
				if( 1.*np.sum(volume_ROI[slice_,:,:])/np.sum(volume_ROI[slice_+1,:,:])<=0.05 ):	
					volume_ROI[slice_,:,:] = volume_ROI[slice_+1,:,:]

	return volume_ROI

def CropDosi(ROI_index = 0): 
	"""
	Crop the dose matrix to the ROI volume.
	The function may take time for large matrices.
	"""
	global dosi
	dosi.setflags(write=1) # needed to override dosi array

	print 'Cropping dosimetry ...'

	#volume_ROI = ROI_3D(ROI_infos, ROI_index)
	volume_ROI = ROI_3D(ROI_infos, ROI_index, tag='CROP') # ROI_index = 0 for contour externe

	dosi[~volume_ROI] = 0 # ~ is the invert operator

	volume_ROI = None
	Update_all()

	print 'Dosimetry cropped'

