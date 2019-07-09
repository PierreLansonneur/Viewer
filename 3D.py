###################################################
#   Display iso surface and volume from numpy array 
# 		   	----
#           	P. Lansonneur 2018
###################################################
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings('ignore', '.*output shape of zoom.*')

def VTK_iso():
	"""
	Select the value of isosurface to display in 3D 
	"""
	global window_VTK_iso, vtk_s1, isosurf_level
	window_VTK_iso = Toplevel(window)
	window_VTK_iso.title('Isosurface visualization')
	isosurf_level = StringVar(window_VTK_iso)
	isosurf_level.set("0")

	if ((dosi_open==True) and (isodose_show==True)):	vtk_s1 = Spinbox(window_VTK_iso, from_=0, to=100, textvariable=isosurf_level)
	else:	vtk_s1 = Spinbox(window_VTK_iso, from_=np.nanmin(volume), to=np.nanmax(volume), textvariable=isosurf_level)

	vtk_bt1 = Button(window_VTK_iso, text="show!", command=VTK_iso_show)
	vtk_t1 = Label(window_VTK_iso, text="value")
	vtk_t1.grid(row=0,column=0)
	vtk_s1.grid(row=0,column=1)
	vtk_bt1.grid(row=0,column=2)

	window_VTK_iso.bind('<Return>', VTK_iso_show)

def VTK_iso_show(event=None):
	"""
	Display isosurface (using python VTK)
	"""
	print '\nshowing iso-surface', float(vtk_s1.get())
	mini = float(vtk_s1.get())
	window_VTK_iso.destroy()

	if ((ROI_open==True)and(ROI_show == True)):	VTK_isosurface(ROI_3D(ROI_infos, 23), spacing, 0.5, 1., True, True)
	elif ((dosi_open==True)and(isodose_show == True)):	VTK_isosurface(dosi, spacing_dosi, 0.01*mini, 0.5, True, False)
	else:	VTK_isosurface(volume, spacing, mini, 0.3, False, True)

def VTK_vol():
	"""
	Display 3D volume (using python VTK)
	"""
	if ((dosi_open==True)and(isodose_show == True)):        VTK_volume(dosi, spacing_dosi, norm=True, scale='rainbow')
	else:	        VTK_volume(volume, spacing, norm=True, scale='grey')


def VTK_isosurface(array, spacing, value, zoom_factor=1, norm=True, save = False, color= np.array([0,0,0])):
	######## Plot the iso-value surface of a numpy array ########
	# array   : numpy array.
	# spacing : associated spacing. exple: [0.5,0.5,1].
	# value   : value of the iso surface to display. exple: 0.5.
	# norm    : Type True to normalize the data between 0 and 1.
	#############################################################

	#print array.shape[2], array.shape[1], array.shape[0]
	'''
	if(array.shape[2]*array.shape[1]*array.shape[0]>1E7): # Big images.. --> reduce dimensions
		from scipy.ndimage import zoom
		zoom_factor = 0.75
		array = zoom(array, (zoom_factor, zoom_factor, zoom_factor))
	'''
	#from scipy.ndimage import zoom
	#myzoom = [zoom_factor, zoom_factor, zoom_factor]
	#array = zoom(array, myzoom, order=0)

	extent_3D = [0, array.shape[2]-1, 0, array.shape[1]-1, 0, array.shape[0]-1]
	array = np.float64(array) 	# set data type

	#"""
	if (norm==True):	
		array = -np.amin(array)+array
		array = array/np.nanmax(array)	# normalize array
	#"""
	# import data
	dataImporter = vtk.vtkImageImport()
	dataImporter.SetImportVoidPointer(array)
	dataImporter.SetDataScalarTypeToDouble()
	dataImporter.SetNumberOfScalarComponents(1)
	dataImporter.SetDataExtent(extent_3D)
	dataImporter.SetWholeExtent(extent_3D)
	#dataImporter.SetDataSpacing([spacing[2],spacing[1],spacing[0]])
	dataImporter.SetDataSpacing((1./zoom_factor)*np.array([spacing[2],spacing[1],spacing[0]]))

	# Extract surface
	#Extractor = vtk.vtkMarchingCubes()
	Extractor = vtk.vtkContourFilter()
	Extractor.SetInputConnection(dataImporter.GetOutputPort())
	Extractor.SetValue(0,value)

	Mapper = vtk.vtkPolyDataMapper()
	Mapper.SetInputConnection(Extractor.GetOutputPort())
	Mapper.ScalarVisibilityOff()

	### Write the stl file to disk
	if(save==True):
		stlWriter = vtk.vtkSTLWriter()
		stlWriter.SetFileName('./output.stl')
		stlWriter.SetInputConnection(Extractor.GetOutputPort())
		stlWriter.Write()

		print 'file output.stl succesfully created!'
		dim = Mapper.GetBounds()

	iso = vtk.vtkActor()
	iso.SetMapper(Mapper)
	if (np.sum(color)>0):	iso.GetProperty().SetColor(color[0],color[2],color[2])

	#################################
	'''
	volumeMapper = vtk.vtkGPUVolumeRayCastMapper()
	volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

	volumeColor = vtk.vtkColorTransferFunction()
	#volumeColor.AddRGBPoint(0, 1.0, 1.0, 1.0)
	volumeColor.AddRGBPoint(0,   0.8, 0.8, 0.8)
	volumeColor.AddRGBPoint(0.5, 0.2, 0.2, 0.2)
	volumeColor.AddRGBPoint(1,   0.0, 0.0, 0.0)

	# The opacity transfer function is used to control the opacity
	# of different tissue types.
	volumeScalarOpacity = vtk.vtkPiecewiseFunction()
	volumeScalarOpacity.AddPoint(0.0,  0.0)
	volumeScalarOpacity.AddPoint(0.5,  0.7)
	volumeScalarOpacity.AddPoint(1.0,  1.0)

	volumeProperty = vtk.vtkVolumeProperty()
	volumeProperty.SetColor(volumeColor)
	volumeProperty.SetScalarOpacity(volumeScalarOpacity)
	volumeProperty.SetInterpolationTypeToLinear()

	# link VTK_vol to volumeMapper and volumeProperty
	VTK_vol = vtk.vtkVolume()
	VTK_vol.SetMapper(volumeMapper)
	VTK_vol.SetProperty(volumeProperty)
	'''
	# Create the renderer, the render window, and the interactor
	ren = vtk.vtkRenderer()			# renderer
	ren.AddActor(iso)
	#ren.AddActor(VTK_vol)
	#ren.SetActiveCamera(camera)
	ren.SetBackground(1, 1, 1) 		# white background
	ren.ResetCameraClippingRange()

	win = vtk.vtkRenderWindow()		# render window
	win.AddRenderer(ren)
	win.SetSize(500, 500)

	inter = vtk.vtkRenderWindowInteractor()	# interactor.
	#inter = vtk.vtkTkRenderWindowInteractor(toolbar_frame)#,  rw=self.renWin, width=400, height=400) 
	inter.SetRenderWindow(win)

	# Interact with the data.
	inter.Initialize()
	win.Render()
	inter.Start()

def VTK_volume(array, spacing, norm=True, scale='grey'):
	########### Display a 3D rendering of a numpy array #########
	# array   : numpy array.
	# spacing : associated spacing. exple: [0.5,0.5,1].
	# scale   : 'grey' or 'rainbow'.
	# norm    : Type True to normalize the data between 0 and 1.
	#############################################################
	'''
	from scipy.ndimage import zoom
	zoom_factor = 0.75
	array = zoom(array, (zoom_factor, zoom_factor, zoom_factor))
	'''
	extent_3D = [0, array.shape[2]-1, 0, array.shape[1]-1, 0, array.shape[0]-1]
	array = np.float64(array) 	# set data type

	if (norm==True):	
		array = -np.amin(array)+array
		array = array/np.nanmax(array)	# normalize array

	# import data
	dataImporter = vtk.vtkImageImport()
	dataImporter.SetImportVoidPointer(array)
	dataImporter.SetDataScalarTypeToDouble()
	dataImporter.SetNumberOfScalarComponents(1)
	dataImporter.SetDataExtent(extent_3D)
	dataImporter.SetWholeExtent(extent_3D)
	dataImporter.SetDataSpacing([spacing[2],spacing[1],spacing[0]])

	volumeMapper = vtk.vtkGPUVolumeRayCastMapper()
	volumeMapper.SetInputConnection(dataImporter.GetOutputPort())

	# Set volume properties ############################
	volumeColor = vtk.vtkColorTransferFunction()

	if (scale=='grey'):
		volumeColor.AddRGBPoint(0, 1.0, 1.0, 1.0)
		volumeColor.AddRGBPoint(0.35, 0.8, 0.8, 0.8)
		volumeColor.AddRGBPoint(1,   0.0, 0.0, 0.0)

	if (scale=='rainbow'):
		volumeColor.AddRGBPoint(0,     0.278431372549, 0.278431372549, 0.858823529412)
		volumeColor.AddRGBPoint(0.143, 0, 0, 0.360784313725)
		volumeColor.AddRGBPoint(0.285, 0, 1, 1)
		volumeColor.AddRGBPoint(0.429, 0, 0.501960784314, 0)
		volumeColor.AddRGBPoint(0.571, 1, 1, 0)
		volumeColor.AddRGBPoint(0.714, 1, 0.380392156863, 0)
		volumeColor.AddRGBPoint(0.857, 0.419607843137, 0, 0)
		volumeColor.AddRGBPoint(1,     0.8784313725489999, 0.301960784314, 0.301960784314)
        #"""
	# The opacity transfer function is used to control the opacity
	# of different tissue types.
	volumeScalarOpacity = vtk.vtkPiecewiseFunction()
	volumeScalarOpacity.AddPoint(0.0,  0.0)
	#volumeScalarOpacity.AddPoint(0.35,  0.0)
	#volumeScalarOpacity.AddPoint(0.5,  0.7)
	volumeScalarOpacity.AddPoint(1.0,  1.0)
        #"""
	volumeProperty = vtk.vtkVolumeProperty()
	volumeProperty.SetColor(volumeColor)
	volumeProperty.SetScalarOpacity(volumeScalarOpacity)
	volumeProperty.SetInterpolationTypeToLinear()

	# link VTK_vol to volumeMapper and volumeProperty
	VTK_vol = vtk.vtkVolume()
	VTK_vol.SetMapper(volumeMapper)
	VTK_vol.SetProperty(volumeProperty)

	# Create the renderer, the render window, and the interactor
	ren = vtk.vtkRenderer()			# renderer
	ren.AddActor(VTK_vol)
	#ren.SetActiveCamera(camera)
	ren.SetBackground(1, 1, 1) 		# white background
	ren.ResetCameraClippingRange()

	win = vtk.vtkRenderWindow()		# render window
	win.AddRenderer(ren)
	win.SetSize(500, 500)

	inter = vtk.vtkRenderWindowInteractor()	# interactor. 
	inter.SetRenderWindow(win)

	# Interact with the data.
	inter.Initialize()
	win.Render()
	inter.Start()

def STLReader(filename):
	reader = vtk.vtkSTLReader()
	reader.SetFileName(filename)

	mapper = vtk.vtkPolyDataMapper()
	mapper.SetInputConnection(reader.GetOutputPort())
	dim = mapper.GetBounds()

	print '\nextent X (mm): {0:.1f}'.format(dim[1]-dim[0])
	print 'extent Y (mm): {0:.1f}'.format(dim[3]-dim[2])
	print 'extent Z (mm): {0:.1f}'.format(dim[5]-dim[4])

	actor = vtk.vtkActor()
	actor.SetMapper(mapper)
	 
	# Create a rendering window and renderer
	ren = vtk.vtkRenderer()
	renWin = vtk.vtkRenderWindow()
	renWin.AddRenderer(ren)
	 
	# Create a renderwindowinteractor
	iren = vtk.vtkRenderWindowInteractor()
	iren.SetRenderWindow(renWin)
	 
	# Assign actor to the renderer
	ren.AddActor(actor)
	 
	# Enable user interface interactor
	iren.Initialize()
	renWin.Render()
	iren.Start()

#STLReader('/media/sf_Linux_Shared/12420/12420.stl')
