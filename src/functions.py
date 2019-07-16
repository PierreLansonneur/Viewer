# -*- coding: utf-8 -*-
###################################################
#   	      General Viewer functions
# 		   	----
#           	P. Lansonneur 2019
###################################################
def on_key(event):
	if event.key == 'control':
		if toolbar._active == 'PAN':	toolbar.release_pan(event)
		if toolbar._active == 'ZOOM':	toolbar.release_zoom(event)

def on_mousewheel(event):	
	global w1, w2, w3
	if event.inaxes==ax1:	w1.set(w1.get()-event.step)
	if event.inaxes==ax2:	w2.set(w2.get()-event.step)
	if event.inaxes==ax3:	w3.set(w3.get()-event.step)

def on_press(event):
	global xp1, yp1, profile_show, ax7, graph2

        on_key(event)

	if event.key == 'control':
		xp1, yp1, profile_show = event.xdata, event.ydata, True 
		event.inaxes.plot(xp1, yp1, 'bo')
		graph1.draw()

	elif event.dblclick:
		if event.inaxes==ax1:	rot_ax1()
		if event.inaxes==ax2:	rot_ax2()
		if event.inaxes==ax3:	rot_ax3()

	else:
		try:
			profile_show = False
			ax7.clear()
			ax7.set_xticks([])
			graph2.draw()

		except Exception:	pass

def on_motion(event):

        on_key(event)

	if event.key == 'control' and profile_show and event.button == 1:

		if len(event.inaxes.lines)>1:	del event.inaxes.lines[-1]
		del event.inaxes.lines[-1]
		event.inaxes.plot([xp1,event.xdata], [yp1,event.ydata], 'bo:')
		graph1.draw()

def on_release(event):	
	global xp1, yp1, xp2, yp2, profile_show, profile_ax

	if profile_show:
		profile_ax = event.inaxes
		xp2, yp2 = event.xdata, event.ydata
		Update_profile()

	Update_all()

def Window_Profile():
	global ax7, graph2, fig2

	window_profile = Toplevel(window)
	window_profile.geometry("{0}x{1}+{2}+{3}".format(500, 200, int(0.53*w_px), 0))
	#window_profile.geometry("{0}x{1}+{2}+{3}".format(450, 150, int(0.53*w_px), 0))
	window_profile.title('Profile (ctrl-click to select)')
	fig2 = P.figure(facecolor='lightgrey')
	fig2.set_size_inches(5,2)
	graph2 = FigureCanvasTkAgg(fig2, master=window_profile)
	canvas2 = graph2.get_tk_widget()
	canvas2.grid(row=0, column=0, sticky=W+E+N+S)
	ax7 = fig2.gca()
	#ax7.set_xticks([])
	P.tight_layout()

	try:    Update_profile()
	except Exception:	pass

        """
	toolbar_frame2 = Frame(window_profile)
	toolbar2 = NavigationToolbar2TkAgg( graph2, toolbar_frame2 )
	toolbar_frame2.grid(row=4,column=0, columnspan=4, sticky=W+E+N+S)
        """

def PVDR(array, prominence_):
        """
       	return PVDR and associated error
        prominence_ is the minimal prominence of peaks
        """

        array = array.astype(float)
	# find peaks and heights
	peaks, properties = find_peaks(array, prominence=prominence_)
	#peaks, properties = find_peaks(array, distance=100, prominence=prominence_)

	# PVDR and associated error
        PVDR = err_PVDR = PVDR_prezado = np.nan

        if (len(peaks)>1): 
                tmp = np.argmax(array[peaks])
                #PVDR_prezado = np.max(array[peaks])/(np.min(array[peaks[tmp]-50:peaks[tmp]+50]))

                np.delete(peaks, tmp)
                np.delete(properties['prominences'], tmp)

	        max_ = np.average(array[peaks])
	        min_ = np.average(array[peaks] - properties['prominences'])
	        PVDR = max_/min_
	        err_max = np.std(array[peaks])
	        err_min = np.std(array[peaks] - properties['prominences'])
	        err_PVDR = (max_/min_)*np.sqrt( (err_max/max_)**2 + (err_min/min_)**2 )
	        #err_PVDR = np.std(array[peaks]/(array[peaks] - properties['prominences']))

        return PVDR, err_PVDR, 1#PVDR_prezado

def Update_profile():

	if dosi_open and isodose_show:
		xp1_vox, yp1_vox = mm2voxel(xp1, yp1, profile_ax, tag='dosi')
		xp2_vox, yp2_vox = mm2voxel(xp2, yp2, profile_ax, tag='dosi')
		prof = AxisImage(profile_ax, tag='dosi')
	else:
		xp1_vox, yp1_vox = mm2voxel(xp1, yp1, profile_ax)
		xp2_vox, yp2_vox = mm2voxel(xp2, yp2, profile_ax)
		prof = AxisImage(profile_ax)

	ax7.clear()
	#ax7.set_xticks([])
	ax7.set_xlabel('x (mm)')
        ax7.xaxis.set_label_coords(1.1,-0.07)

	ax7.set_ylabel('Dose (a.u.)')
        P.tight_layout()
        #dist = np.sqrt( (yp2-yp1)**2 + (xp2-xp1)**2 )

	xp,yp = np.linspace(xp1_vox,xp2_vox,1000), np.linspace(yp1_vox,yp2_vox,1000)
        dist_ = np.linspace(0,np.sqrt( (yp2-yp1)**2 + (xp2-xp1)**2 ),1000)
        prof_ = map_coordinates(prof, np.vstack((xp,yp)))

        #p_seuil = np.nanmax(dosi)
        #print '*\t{0:.2f} +/- {1:.2f}\t{2:.2f}'.format(PVDR(prof_, 0.1*p_seuil)[0],PVDR(prof_, 0.1*p_seuil)[1],PVDR(prof_, 0.1*p_seuil)[2])
        #ax7.set_title('PVDR : {0:.1f} +/- {1:.1f} ({2:.2f})'.format(PVDR(prof_, 0.1*p_seuil)[0],PVDR(prof_, 0.1*p_seuil)[1],PVDR(prof_, 0.1*p_seuil)[2]))

	ax7.plot(dist_, prof_)
	graph2.draw()

	#with open('tmp_profile.txt', 'a') as f:
	with open('tmp_profile.txt', 'w') as f:
		for i,j in zip(dist_, prof_):
			f.write('{0:.3f}, {1}\n'.format(i,j))

def mm2voxel(x_mm,y_mm,ax,tag='CT'):
	"""
	convert coordinates in mm to coordinates in voxel units 
	"""
	if tag=='CT':	orig, sp = origin, spacing
	if tag=='dosi':	orig, sp = origin_dosi, spacing_dosi

	x_vox, y_vox = 0,0

	if ax==ax1:			x_vox, y_vox = int((orig[2]-y_mm)/sp[2]), int((orig[1]-x_mm)/sp[1])
	if (ax==ax1)and(rot1==True):	x_vox, y_vox = int((orig[2]-x_mm)/sp[2]), int((orig[1]-y_mm)/sp[1])
	if ax==ax2:			x_vox, y_vox = int((orig[2]-y_mm)/sp[2]), int((-orig[0]+x_mm)/sp[0])
	if (ax==ax2)and(rot2==True):	x_vox, y_vox = int((orig[2]-x_mm)/sp[2]), int((-orig[0]+y_mm)/sp[0])
	if ax==ax3:			x_vox, y_vox = int((orig[1]-y_mm)/sp[1]), int((-orig[0]+x_mm)/sp[0])
	if (ax==ax3)and(rot3==True):	x_vox, y_vox = int((orig[1]-x_mm)/sp[1]), int((-orig[0]+y_mm)/sp[0])

	return x_vox, y_vox

def AxisImage(ax,tag='CT'):
	"""
	return the array corresponding to the selected axis
	"""
	if tag=='CT':	a1, a2, a3 = im1, im2, im3
	if tag=='dosi':	a1, a2, a3 = dos1, dos2, dos3

	array = []

	if ax==ax1:			array = a1
	if (ax==ax1)and(rot1==True):	array = np.fliplr(a1.T)
	if ax==ax2:			array = a2
	if (ax==ax2)and(rot2==True):	array = np.fliplr(a2.T)
	if ax==ax3:			array = a3
	if (ax==ax3)and(rot3==True):	array = np.fliplr(a3.T)

	return array

def show_isodose():	
	"""
	display isodoses on each axis
	"""
	global isodose_show
	isodose_show = not isodose_show
        try:    cbaxes.set_visible(isodose_show)
        except Exception:       pass
	Update_all()

def show_ROI():	
	"""
	display RT structures on axis 1
	"""
	global ROI_show
	ROI_show = not ROI_show
	Update_all()

def file_save():
	"""
	Save CT volume as a .mhd file
	"""
	types = [('MHD files', '*.mhd')]
	f = tkFileDialog.asksaveasfilename(initialdir = dir_ini+'output', filetypes=types)
	#itkimage = sitk.GetImageFromArray(volume.astype(np.int16)) # for big CT images
	itkimage = sitk.GetImageFromArray(volume)
	itkimage.SetSpacing([spacing[2],spacing[1],spacing[0]])
	itkimage.SetOrigin([origin[2],origin[1],origin[0]])
	writer = sitk.ImageFileWriter()
   	writer.SetFileName(f)
  	writer.Execute(itkimage)
	print('File saved successfully!')

def about():
	"""
	Developer informations
	"""
	window_about = Toplevel(window)
	window_about.title('About') 
	display = Label(window_about, text=" P. Lansonneur, 2018 \n pierre.lansonneur@curie.fr ")
	display.grid(row=0, column=0, sticky=W, padx = 10, pady = 10)

def FileInfo():
	"""
	Display Files informations
	"""
	window_info = Toplevel(window)
	window_info.title('File properties')

	try:
		disp0 = Label(window_info, text=' CT : '+filename_CT.rsplit('/', 1)[-1])
		disp1 = Label(window_info, text='  size :\t{0}\t{1}\t{2}\tpix\t'.format(dim_x, dim_y, dim_z))
		disp2 = Label(window_info, text='  size :\t{0:.1f}\t{1:.1f}\t{2:.1f}\tmm\t'.format(dim_x*spacing[0], dim_y*spacing[1], dim_z*spacing[2]))
		disp3 = Label(window_info, text='  spacing :\t{0:.2f}\t{1:.2f}\t{2:.2f}\tmm\t'.format(spacing[0], spacing[1], spacing[2]))
		disp4 = Label(window_info, text='  origin :\t{0:.2f}\t{1:.2f}\t{2:.2f}\tmm\t'.format(origin[0], origin[1], origin[2]))

		disp0.grid(row=0, column=0, sticky=W)
		disp1.grid(row=1, column=0, sticky=W)
		disp2.grid(row=2, column=0, sticky=W)
		disp3.grid(row=3, column=0, sticky=W)
		disp4.grid(row=4, column=0, sticky=W)

		if dosi_open:
			disp0 = Label(window_info, text='\n Dose : '+filename_dosi.rsplit('/', 1)[-1])
			disp1 = Label(window_info, text='  size :\t{0}\t{1}\t{2}\tpix\t'.format(dim_x_dosi, dim_y_dosi, dim_z_dosi))
			disp2 = Label(window_info, text='  size :\t{0:.1f}\t{1:.1f}\t{2:.1f}\tmm\t'.format(dim_x_dosi*spacing_dosi[0], dim_y_dosi*spacing_dosi[1], dim_z_dosi*spacing_dosi[2]))
			disp3 = Label(window_info, text='  spacing :\t{0:.2f}\t{1:.2f}\t{2:.2f}\tmm\t'.format(spacing_dosi[0], spacing_dosi[1], spacing_dosi[2]))
			disp4 = Label(window_info, text='  origin :\t{0:.2f}\t{1:.2f}\t{2:.2f}\tmm\t'.format(origin_dosi[0], origin_dosi[1], origin_dosi[2]))

			disp0.grid(row=5, column=0, sticky=W)
			disp1.grid(row=6, column=0, sticky=W)
			disp2.grid(row=7, column=0, sticky=W)
			disp3.grid(row=8, column=0, sticky=W)
			disp4.grid(row=9, column=0, sticky=W)

		if ROI_open:
			disp0 = Label(window_info, text='\n RT-struct : '+filename_ROI.rsplit('/', 1)[-1])
			disp0.grid(row=10, column=0, sticky=W)

			for ROI_index in range(N_ROI):
				ROI_color = 255*ROI_infos[ROI_index,5].astype(float)
				ROI_color = ROI_color.astype(int)
				ROI_color = '#%02x%02x%02x' % (ROI_color[0], ROI_color[1], ROI_color[2])
				if(ROI_index<10):	disp1 = Label(window_info, text='  {0}   '.format(ROI_index)+ROI_infos[ROI_index,0].decode('utf-8'))
				else:	disp1 = Label(window_info, text='  {0} '.format(ROI_index)+ROI_infos[ROI_index,0].decode('utf-8'))
				disp2 = Label(window_info, text='      ', bg=ROI_color)
				if(ROI_index<=20):
					disp1.grid(row=11+ROI_index, column=0, sticky=W)
					disp2.grid(row=11+ROI_index, column=1, sticky=W)
				else:
					disp1.grid(row=-10+ROI_index, column=2, sticky=W)
					disp2.grid(row=-10+ROI_index, column=3, sticky=W)

	except Exception:	pass

def Yrevert_ax1():
	global revY1	
	revY1 = True
	Update_all()

def Xrevert_ax1():	
	global revX1	
	revX1 = True
	Update_all()

def Xrevert_ax2():	
	global revX2	
	revX2 = True
	Update_all()

def Yrevert_ax2():	
	global revY2	
	revY2 = True
	Update_all()

def Yrevert_ax3():	
	global revY3	
	revY3 = True
	Update_all()

def Xrevert_ax3():	
	global revX3	
	revX3 = True
	Update_all()

def rot_ax1():	
	global rot1, w1_moved	
	rot1 = not rot1
	Set_axes_lim_init()
	Update_all()

def rot_ax2():	
	global rot2, w2_moved, x2_lim, y2_lim
	rot2 = not rot2
	Set_axes_lim_init()
	Update_all()

def rot_ax3():	
	global rot3, w3_moved	
	rot3 = not rot3
	Set_axes_lim_init()
	Update_all()

def Inv_Scale():
	"""
	Reverse grayscale (CT volume)
	"""	
	global inv_scale
	inv_scale = not inv_scale
	if (inv_scale == True):
		P.set_cmap('Greys_r')
		for ax in [ax1, ax2, ax3]: ax.set_facecolor('0')
	else:
		P.set_cmap('Greys')
		for ax in [ax1, ax2, ax3]: ax.set_facecolor('1')
	Update_all()

def Low_Contrast():
	"""
	Lower the image contrast for the CT volume
	"""	
	global volume	
	min_volume, max_volume = np.nanmin(volume), np.nanmax(volume)
	volume = np.array(volume,dtype='float')
	volume = volume-np.amin(volume)				# normalize array to [0,1]
	volume = volume/np.nanmax(volume)			# normalize array to [0,1]	
	volume = np.power(volume,0.5)  				# gamma correction
	volume = (max_volume-min_volume)*volume+min_volume 	# rescale	
	Update_all()

def High_Contrast():
	"""
	Rise the image contrast for the CT volume
	"""
	global volume	
	min_volume, max_volume = np.nanmin(volume), np.nanmax(volume)
	volume = np.array(volume,dtype='float')
	volume = volume-np.amin(volume)				# normalize array to [0,1]
	volume = volume/np.nanmax(volume)			# normalize array to [0,1]	
	volume = np.power(volume,2)   				# gamma correction
	volume = (max_volume-min_volume)*volume+min_volume 	# rescale		
	Update_all()

def FilterDosi():
	"""
	"""
	global dosi	
	from scipy import ndimage
	dosi.setflags(write=1) # allows to override volume
	dosi = ndimage.median_filter(dosi,3)
	#dosi = ndimage.gaussian_filter(dosi,sigma = 1)
	Update_all()

def SetIntensityRange(array,minI, maxI):
	"""
	Voxels out of bounds are set to 0
	"""
	global volume, dosi	
	array.setflags(write=1) # allows to override volume
	index_minI = np.where(array < minI)
	index_maxI = np.where(array >= maxI)
	array[index_minI] = 0 # minI
	array[index_maxI] = 0 # maxI
	Update_all()

def SetDisplayRange(minI, maxI):
	""" 
	Set axes colorscale 

	>>>Usage:
	SetDisplayRange(-1000, 2000) #Hounsfield scale
	"""

        for ax in [ax1, ax2, ax3]:
        	if dosi_open and isodose_show:
                        if(len(ax.get_images())>0):     ax.get_images()[0].set_clim(minI, maxI)
                else:
	                for im in ax.get_images():	im.set_clim(minI, maxI)

        """
	for ax in [ax1, ax2, ax3]:
		for im in ax.get_images():	im.set_clim(minI, maxI)

        #from matplotlib import colors
        #ax.get_images()[0].set_norm(colors.PowerNorm(gamma=1./2.))
        #ax.get_images()[0].set_norm(colors.PowerNorm(gamma=2.))
        """
def VolumeHistogram():
	"""
	Display the volume histogram and set the volume intensity range
	"""
	window_histo = Toplevel(window)
	window_histo.geometry("{0}x{1}+{2}+{3}".format(600, 335, int(0.53*w_px), 0))
	window_histo.title('Histogram')

	fig2 = P.figure(facecolor='lightgrey')
	fig2.set_size_inches(6,3)
	graph2 = FigureCanvasTkAgg(fig2, master=window_histo)
	canvas2 = graph2.get_tk_widget()
	toolbar_frame2 = Frame(window_histo)
	toolbar2 = NavigationToolbar2TkAgg( graph2, toolbar_frame2 )
	ax7 = fig2.gca() 

	canvas2.grid(row=3, column=0, columnspan=5, sticky=W+E+N+S)
	toolbar_frame2.grid(row=4,column=0, columnspan=4, sticky=W+E+N+S)

	# volume histo
	ax7.hist(np.ravel(volume), bins='auto', range=(-1100,2000), histtype='step', density=True)
	#ax7.set_xlabel('value')
	ax7.set_yticks([])
	P.tight_layout()

	graph2.draw()
	Update()

def RaiseError():	tkMessageBox.showerror('Error', 'Not yet implemented!')

def Set_scales(): 
	"""
	Set the limits of sliders
	"""
	global w1, w2, w3
	w1 = Scale(window, from_=0, to=dim_x-1, orient=VERTICAL, command=w1move, showvalue = 0)
	w2 = Scale(window, from_=0, to=dim_y-1, orient=VERTICAL, command=w2move, showvalue = 0)
	w3 = Scale(window, from_=0, to=dim_z-1, orient=VERTICAL, command=w3move, showvalue = 0)
	w1.grid(row=0,column=2,sticky=W+N+S)
	w2.grid(row=1,column=2,sticky=W+N+S)
	w3.grid(row=2,column=2,sticky=W+N+S)
	w1.set(int(dim_x/2.))
	w2.set(int(dim_y/2.))
	w3.set(int(dim_z/2.))

def Get_axes_lim():
	global x1_lim, y1_lim, x2_lim, y2_lim, x3_lim, y3_lim
	x1_lim, y1_lim = ax1.get_xlim(), ax1.get_ylim()
	x2_lim, y2_lim = ax2.get_xlim(), ax2.get_ylim()
	x3_lim, y3_lim = ax3.get_xlim(), ax3.get_ylim()

def Set_axes_lim():
	global x1_lim, x2_lim, x3_lim, y1_lim, y2_lim, y3_lim
	global revX1, revX2, revX3, revY1, revY2, revY3

	if(revX1 == True):	x1_lim = np.flip(x1_lim,0)
	if(revX2 == True):	x2_lim = np.flip(x2_lim,0)
	if(revX3 == True):	x3_lim = np.flip(x3_lim,0)	
	if(revY1 == True):	y1_lim = np.flip(y1_lim,0)
	if(revY2 == True):	y2_lim = np.flip(y2_lim,0)
	if(revY3 == True):	y3_lim = np.flip(y3_lim,0)

	ax1.set_xlim(x1_lim)
	ax1.set_ylim(y1_lim)
	ax2.set_xlim(x2_lim)
	ax2.set_ylim(y2_lim)
	ax3.set_xlim(x3_lim)
	ax3.set_ylim(y3_lim)

	ax1.get_shared_y_axes().join(ax1, ax4)
	ax2.get_shared_y_axes().join(ax2, ax5)
	ax3.get_shared_y_axes().join(ax3, ax6)

	revX1 = revX2 = revX3 = revY1 = revY2 = revY3 = False

def Set_axes_lim_init():
	global x1_lim, x2_lim, x3_lim, y1_lim, y2_lim, y3_lim

	x2_lim = x3_lim = (origin[0],origin[0]+dim_x*spacing[0])
	x1_lim = y3_lim = (origin[1]-dim_y*spacing[1],origin[1])
	y1_lim = y2_lim = (origin[2]-dim_z*spacing[2],origin[2])
        """
        if ID=='12420':
	        x2_lim = x3_lim = (-90,55)
	        x1_lim = y3_lim = (-130,130)
	        y1_lim = y2_lim = (-87,87)

        if ID=='102771':
	        x2_lim = x3_lim = (-70,100)
	        x1_lim = y3_lim = (-150,100)
	        y1_lim = y2_lim = (-80,80)
        """
	if(rot1 == True):	x1_lim, y1_lim = y1_lim, x1_lim
	if(rot2 == True):	x2_lim, y2_lim = y2_lim, x2_lim
	if(rot3 == True):	x3_lim, y3_lim = y3_lim, x3_lim

	Set_axes_lim()

def slice_ax1_dosi(slice_ax1):
	"""
	return the slice along axis 1 in the dosimetry volume 
	corresponding to the slice "slice_ax1" in the CT volume
	"""
	if (origin[0] + slice_ax1*spacing[0] <= origin_dosi[0]):	return 0
	if (origin[0] + slice_ax1*spacing[0] >= origin_dosi[0] + dim_x_dosi*spacing_dosi[0]):	return -1
	else:	return int((slice_ax1-int((origin_dosi[0] - origin[0])/spacing[0]))*spacing[0]/spacing_dosi[0])

def slice_ax2_dosi(slice_ax2):
	"""
	return the slice along axis 2 in the dosimetry volume 
	corresponding to the slice "slice_ax2" in the CT volume
	"""
	if (origin[1] + slice_ax2*spacing[1] <= origin_dosi[1]):	return 0
	if (int((slice_ax2+int((origin_dosi[1] - origin[1])/spacing[1]))*spacing[1]/spacing_dosi[1]) >= dim_y_dosi):	return -1
	else:	return int((slice_ax2+int((origin_dosi[1] - origin[1])/spacing[1]))*spacing[1]/spacing_dosi[1])

def slice_ax3_dosi(slice_ax3):
	"""
	return the slice along axis 3 in the dosimetry volume 
	corresponding to the slice "slice_ax3" in the CT volume
	"""
	if (origin[2] + slice_ax3*spacing[2] <= origin_dosi[2]):	return 0
	if (int((slice_ax3+int((origin_dosi[2] - origin[2])/spacing[2]))*spacing[2]/spacing_dosi[2]) >= dim_z_dosi):	return -1
	else:	return int((slice_ax3+int((origin_dosi[2] - origin[2])/spacing[2]))*spacing[2]/spacing_dosi[2])

def HounsfieldScale():
	global c_scale
	c_scale='HOUNSFIELD'
	Update_all()

def AutoScale():
	global c_scale
	c_scale=None
	Update_all()

def UserScale():
	global c_scale
	c_scale='USER'
	Update_all()

def w1move(self):
	global w1_moved
	w1_moved = True
	Update() 
	w1_moved = False

def w2move(self):
	global w2_moved
	w2_moved = True
	Update()
	w2_moved = False

def w3move(self):
	global w3_moved
	w3_moved = True
	Update()
	w3_moved = False

def w4move(self):
	global w4_moved
	w4_moved = True
 	Update()
	w4_moved = False

def set_Levels():
	"""
	Set the Isodose levels
	"""
	global levels, window_levels, s1, s2, Dvar
	window_levels = Toplevel(window)
	window_levels.title('Isodose display')

	minvar, maxvar = IntVar(window_levels), IntVar(window_levels)
	minvar.set("20")        # isodose min
	maxvar.set("100")       # isodose max
        Dvar.set(1)

	s1 = Spinbox(window_levels, from_=0, to=100, width=4, textvariable=minvar, disabledbackground='lightgrey')
	s2 = Spinbox(window_levels, from_=0, to=100, width=4, textvariable=maxvar, disabledbackground='lightgrey')
	ts1 = Label(window_levels, text="From ")
	ts2 = Label(window_levels, text="% to ")
	ts3 = Label(window_levels, text="% of ")
        rb1 = Radiobutton(window_levels, text = "D max", variable=Dvar, value = 1, command=set_Display_norm)
        rb2 = Radiobutton(window_levels, text = "D PTV", variable=Dvar, value = 2, command=set_Display_norm)
        rb3 = Radiobutton(window_levels, text = "gradient", variable=Dvar, value = 3, command=set_Display_norm)
	bt1 = Button(window_levels, text="show!", command=set_Levels_validate)

	ts1.grid(row=0,column=0)
	s1.grid( row=0,column=1)
	ts2.grid(row=0,column=2)
	s2.grid( row=0,column=3)
	ts3.grid(row=0,column=4)
	rb1.grid(row=0,column=5, sticky=W)
	rb2.grid(row=0,column=6, sticky=W)
	rb3.grid(row=0,column=7, sticky=W)
	bt1.grid(row=0,column=8, sticky=W)

        if(ROI_open==False):      rb2.config(state=DISABLED)

	window_levels.bind('<Return>', set_Levels_validate)

def set_Display_norm():
        if(Dvar.get()==1 or 2):
              s1.config(state=NORMAL)
              s2.config(state=NORMAL)
        if(Dvar.get()==3):
              s1.config(state=DISABLED)
              s2.config(state=DISABLED)

def set_Levels_validate(event=None):
	global levels, coeff_D_PTV, cbaxes
	levels = np.linspace(0.01*float(s1.get()),0.01*float(s2.get()),9)

        if(Dvar.get()==2):
                for ROI_index in range(N_ROI):
		                isPTV = ROI_infos[ROI_index,0].upper().startswith('PTV')
		                if( isPTV ):	
			                DVH = ROI_DVH(ROI_infos[ROI_index,:],'dosi')
                                        if (len(DVH)>1):        coeff_D_PTV = np.mean(DVH)
                print 'PTV average dose: ',coeff_D_PTV

        if(Dvar.get()==1 or 2):
                try:    cbaxes.remove()
                except Exception:       pass

        if(Dvar.get()==3):      cbaxes = fig1.add_axes([0.78,0.67,0.02,0.31])

	window_levels.destroy()
	Update_all()

def Rescale(tag='dosi'): 
	"""
	Rescale the dose matrix such that the mean dose given to the PTV is 100 Gy
	"""
	global volume, dosi
	volume.setflags(write=1) # needed to override vol array
	dosi.setflags(write=1) # needed to override dosi array

	for ROI_index in range(N_ROI):	
		isPTV = ROI_infos[ROI_index,0].upper().startswith('PTV')
		isCTV = ROI_infos[ROI_index,0].upper().startswith('CTV')
		isGTV = ROI_infos[ROI_index,0].upper().startswith('GTV')

		#if( isPTV ):	
		if( isCTV ):	
			DVH = ROI_DVH(ROI_infos[ROI_index,:],tag)
                        if (len(DVH)>1):
                                coeff_ = np.mean(DVH)
			        print tag,'\t', ROI_infos[ROI_index,0].decode('utf-8'),'\t',np.mean(DVH), 'Gy'

        if (tag=='CT'):  volume = (100./coeff_)*volume
        else:   dosi = (100./coeff_)*dosi

	Update_all()

def Demo():
	"""
	Test function
	ID: patient ID
	"""
	global ID

	ID, w1_set, w2_set, w3_set = '12420', 243, 231, 297# 80, 91, 84
	#ID, w1_set, w2_set, w3_set = '102771', 70, 157, 127 #191, 323, 233

	dirname_ct = dir_ini+ID+'/CT/'
	file_path_dosi = dir_ini+ID+'/pi/RD'+ID+'.dcm'
	file_path_ROI = dir_ini+ID+'/pi/RS'+ID+'.dcm'
	file_path_RP = dir_ini+ID+'/pi/RP'+ID+'.dcm'

	OpenDicomSerie(dirname_ct)
	#OpenFile(file_path_dosi)
	#OpenFile(dir_ini+ID+'/CT'+ID+'.mhd')
	OpenDosi(file_path_dosi)
	#OpenDosi(dir_ini+ID+'/LET/1mm/merged.dcm')
	#OpenDosi(dir_ini+ID+'/LET/1mm/LET_PS1N.dcm')
	#OpenDosi(dir_ini+ID+'/test_dosi/Dose_Field_2_105.dcm')
	#OpenDosi(dir_ini+ID+'/crop/Dose_PS1N_crop.dcm')
	#OpenDosi(dir_ini+ID+'/ctc6/Dose_Field_123_ctc6.mhd')
	#OpenDosi(dir_ini+ID+'/ctc4/merged.dcm')
	#OpenDosi(dir_ini+ID+'/scan_mbrt/Dose_PS1N_ctc4_w700.dcm')
	#OpenDosi(dir_ini+ID+'/SFUD/merged_.dcm')
	#OpenDosi(dir_ini+ID+'/SFUD/Dose_Field_3.dcm')
	#OpenDosi(dir_ini+ID+'/Dose_104005_PD2N.dcm')
	open_ROI(file_path_ROI)
	#OpenRP(file_path_RP)

	w1.set(w1_set)
	w2.set(w2_set)
	w3.set(w3_set)

        #Set_axes_lim_init()
        #Rescale(tag='dosi')
        #Rescale(tag='CT')
        #CropDosi(0)
        #ROI_DVH_analysis()
	Update_all()      
        #Window_Gamma()
