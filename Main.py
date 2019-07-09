# -*- coding: utf-8 -*-
############################################
# Visualize dosimetry files with this viewer
# 		   ---
#           P. Lansonneur 2018
############################################

import matplotlib.pyplot as P
import numpy as np
import SimpleITK as sitk
from Tkinter import *
import tkFileDialog
import tkMessageBox
#from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import pydicom, pydicom.uid
from pydicom.dataset import Dataset, FileDataset
import datetime, time
import os
import matplotlib.gridspec as gridspec
from operator import itemgetter, attrgetter
import vtk
#from vtk.tk.vtkTkRenderWindowInteractor import vtkTkRenderWindowInteractor # not working :(
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="No contour levels were found")
warnings.filterwarnings("ignore", message="converting a masked element to nan")
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.ndimage import rotate, zoom, map_coordinates
from fpdf import FPDF
from scipy.signal import find_peaks

window = Tk()
window.title("Viewer")
w_px, h_px = window.winfo_screenwidth(), window.winfo_screenheight()
w_mm, h_mm = window.winfo_screenmmwidth(), window.winfo_screenmmheight()
w_in, h_in = w_mm / 25.4, h_mm / 25.4
#w_dpi, h_dpi = w_px/w_in, h_px/h_in
w_dpi, h_dpi = 100, 100
window.geometry("{0}x{1}+0+0".format(int(0.485*w_px), h_px))

############ Parameters initialization ############
### General
info=''							# file infos
ID=''							# Patient ID
low_contrast = high_contrast = inv_scale = False	# lower/raise contrast, invert greyscale
revX1 = revX2 = revX3  = False				# reverse X axis
revY1 = revY2 = revY3 = False				# reverse Y axis
#rot1 = rot2 = rot3 = False				# rotate axes
rot1 = rot2 = rot3 = True
w1_moved = w2_moved = w3_moved = w4_moved = False	# did sliders move? 
profile_show = False					# show profile
profile_ax = None					# axis selected
#c_scale = 'HOUNSFIELD'					# Grey scale limits
c_scale = None					        # Grey scale limits

### CT
CT_open = False						# is CT opened?
dim_x = dim_y = dim_z = 100				# CT scan dimensions
x1_lim=x2_lim=x3_lim=y1_lim=y2_lim=y3_lim=(0,dim_x)	# CT projections range
volume = np.zeros((dim_x,dim_y,dim_z))			# CT volume
spacing = [1, 1, 1]					# CT spacing
origin = [0, 0, 0]					# CT origin
im1 = im2 = im3 = im = volume[0,:,:]			# CT scan projections
ext7 = ext8 = ext9 = extent = [0,100,0,100]		# CT scan projections dimensions
dir_ini = '/media/sf_Linux_Shared/'			# CT scans directory
filename_CT = None

### dosimetry
origin_dosi = [0, 0, 0]					# dosi origin
spacing_dosi = [1, 1, 1]				# dosi spacing
dosi_open = False					# is dosi opened?
isodose_show = False					# show isodoses
levels = np.array([0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]) 	# isodose levels
#levels = np.array([0.5,0.6,0.7,0.75,0.8,0.85,0.9,1]) 	# isodose levels
Dvar = IntVar(window)                                   # dose display method
Dvar.set(1)                                             # normalize to the max
coeff_D_PTV = np.nan                                    # PTV average dose
dosemap = P.get_cmap('jet')
filename_dosi = None

### ROI
N_ROI = 0						# number of ROI
ROI_infos = np.empty((1,7), dtype=np.object)		# ROI informations
ROI_open = False					# is ROI file opened?
ROI_show = False					# show ROIs
filename_ROI = None

### RP
#RP_open = False					# is RP file opened?
RP_crop_show = False					# show volume of calculation
filename_RP = None    
crop_RP = None
keepCTDim = False                                       # make matrix dimensions equal to CT dimensions
#keepCTDim = True                                       # make matrix dimensions equal to CT dimensions

poly_Espread 	= np.poly1d([-0.00000023, 0.00008559, -0.01205498, 1.23066909])		# energy spread
poly_width 	= np.poly1d([-0.00000143, 0.00086454, -0.19620932, 21.66378475])	# sigma (mm)
poly_div 	= np.poly1d([-0.000000002, 0.000001141, -0.00021799, 0.0156])		# angular spread
poly_SM1 	= np.poly1d([0.0005266, -0.2767, 60.29])				# SM1 field (T), Ludo
poly_SM2 	= np.poly1d([-0.0004444, 0.2328, -50.43])				# SM2 field (T), Ludo
poly_UM 	= np.poly1d([-2.92661e-11,  2.31840e-08, -6.63090e-06,  9.298374e-04])	# MU/Gp coefficient

### Load functions
execfile("./src/functions.py")
execfile("./src/3D.py")
execfile("./src/OpenFile.py")
execfile("./src/RP2TOPAS.py")
execfile("./src/GammaIndex.py")
execfile("./src/Update.py")	
execfile("./src/ROI.py")

############ Panel settings ############

P.set_cmap('Greys')
fig1 = P.figure(facecolor='lightgrey')
fig1.set_size_inches(0.45*w_in, 0.8*h_in)
#fig1.set_size_inches(0.3*w_in, 0.8*h_in)
gs1 = gridspec.GridSpec( nrows=3, ncols=2, width_ratios=[7.5,1], bottom=0.01, top=0.99, left = 0.01, right=0.99, hspace=0, wspace=0)
#gs1 = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[50,0.01], bottom=0.01, top=0.99, left = 0.01, right=0.99, hspace=0, wspace=0)
ax1 = fig1.add_subplot(gs1[0,0])
ax2 = fig1.add_subplot(gs1[1,0])
ax3 = fig1.add_subplot(gs1[2,0])
ax4 = fig1.add_subplot(gs1[0,1])
ax5 = fig1.add_subplot(gs1[1,1])
ax6 = fig1.add_subplot(gs1[2,1])

for ax in [ax1, ax2, ax3]:	ax.axis('equal')
for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:	ax.set_xticks([]) 
for ax in [ax1, ax2, ax3, ax4, ax5, ax6]:	ax.set_yticks([])

t1 = ax1.set_ylabel("",labelpad=-15, rotation = "horizontal", position=(0.,0.9), color='b')
t2 = ax2.set_ylabel("",labelpad=-15, rotation = "horizontal", position=(0.,0.9), color='b')
t3 = ax3.set_ylabel("",labelpad=-15, rotation = "horizontal", position=(0.,0.9), color='b')

t4 = Label(window, text="  isodoses :")
t6 = Label(window, text="structures :")

graph1 = FigureCanvasTkAgg(fig1, master=window)
canvas1 = graph1.get_tk_widget()

toolbar_frame = Frame(window)
toolbar = NavigationToolbar2TkAgg( graph1, toolbar_frame )
#toolbar = NavigationToolbar2Tk( graph1, toolbar_frame )

toolbar.pan() # default widget

check1 = Checkbutton(window, command = show_isodose)
check2 = Checkbutton(window, command = show_ROI)

w1 = Scale(window, from_=0, to=dim_x-1, orient=VERTICAL, command=w1move, showvalue = 0)
w2 = Scale(window, from_=0, to=dim_y-1, orient=VERTICAL, command=w2move, showvalue = 0)
w3 = Scale(window, from_=0, to=dim_z-1, orient=VERTICAL, command=w3move, showvalue = 0)
w4 = Scale(window, from_=0, to=1000, orient=HORIZONTAL, command=w4move, showvalue = 0, length=0.405*w_px)
#w4 = Scale(window, from_=0, to=1000, orient=HORIZONTAL, command=w4move, showvalue = 0, length=0.205*w_px)

### matplotlib events handler

fig1.canvas.mpl_connect('key_press_event', on_key)
fig1.canvas.mpl_connect('scroll_event', on_mousewheel)
fig1.canvas.mpl_connect('button_press_event', on_press)
fig1.canvas.mpl_connect('motion_notify_event', on_motion)
fig1.canvas.mpl_connect('button_release_event', on_release)
fig1.canvas.get_tk_widget().focus_force() # detect event.key when 'button_press_event' is active

############ Items locations ############

t4.grid(row=3,column=0, columnspan=2, sticky=E)
#t5.grid(row=3,column=0, columnspan=2, sticky=W, padx = 5)
t6.grid(row=4,column=1, sticky=E)
canvas1.grid(row=0, column=0, rowspan=3, columnspan=2, sticky=W+E+N+S)
toolbar_frame.grid(row=4,column=0, sticky=W+E+N+S)
check1.grid(column=2, row=3)
check2.grid(column=2, row=4,sticky=E+N+W)
w1.grid(row=0,column=2,sticky=W+N+S)
w2.grid(row=1,column=2,sticky=W+N+S)
w3.grid(row=2,column=2,sticky=W+N+S)
w4.grid(row=3,column=0, columnspan=2, sticky=W+N+S, padx = 7)

################ Menu bar ##################

menubar = Menu(window)
window.config(menu=menubar)

### File menu
filemenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="File", menu=filemenu)
filemenu.add_command(label="Open", command=OpenFile)
filemenu.add_command(label="Open a DICOM serie", command=OpenDicomSerie)
filemenu.add_command(label="Import dosimetry", command=OpenDosi)
filemenu.add_command(label="Open RP file...", command=OpenRP)
#filemenu.add_command(label="Demo", command=ID_demo)

filemenu.add_separator()
filemenu.add_command(label="Save as...", command=file_save)
filemenu.add_command(label="Exit", command=window.quit)

### Edit menu
editmenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="Edit", menu=editmenu)

sub_editmenu = Menu(editmenu, tearoff=0)
editmenu.add_cascade(label="Revert X axis", menu=sub_editmenu)
sub_editmenu.add_command(label="ax1", command=Xrevert_ax1)
sub_editmenu.add_command(label="ax2", command=Xrevert_ax2)
sub_editmenu.add_command(label="ax3", command=Xrevert_ax3)

sub_editmenu = Menu(editmenu, tearoff=0)
editmenu.add_cascade(label="Revert Y axis", menu=sub_editmenu)
sub_editmenu.add_command(label="ax1", command=Yrevert_ax1)
sub_editmenu.add_command(label="ax2", command=Yrevert_ax2)
sub_editmenu.add_command(label="ax3", command=Yrevert_ax3)

sub_editmenu = Menu(editmenu, tearoff=0)
editmenu.add_cascade(label="Rotate", menu=sub_editmenu)
sub_editmenu.add_command(label="ax1", command=rot_ax1)
sub_editmenu.add_command(label="ax2", command=rot_ax2)
sub_editmenu.add_command(label="ax3", command=rot_ax3)

sub_editmenu = Menu(editmenu, tearoff=0)
editmenu.add_cascade(label="Scale", menu=sub_editmenu)
sub_editmenu.add_command(label="Auto", command=AutoScale)
sub_editmenu.add_command(label="Hounsfield", command=HounsfieldScale)
sub_editmenu.add_command(label="User defined", command=UserScale)
sub_editmenu.add_separator()
sub_editmenu.add_command(label="Contrast +", command=High_Contrast)
sub_editmenu.add_command(label="Contrast -", command=Low_Contrast)
sub_editmenu.add_command(label="Invert scale", command=Inv_Scale)

### Tool menu
toolmenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="Tools", menu=toolmenu)
toolmenu.add_command(label="Isodoses...", command=set_Levels)
toolmenu.add_command(label="Profile...", command=Window_Profile)
toolmenu.add_command(label="Histogram...", command=VolumeHistogram)
toolmenu.add_command(label="Gamma...", command=Window_Gamma)
toolmenu.add_command(label="Iso-surface...", command=VTK_iso)
"""
sub3Dmenu = Menu(toolmenu, tearoff=0)
toolmenu.add_cascade(label="3D",menu=sub3Dmenu)
sub3Dmenu.add_command(label="Iso-surface...", command=VTK_iso)
sub3Dmenu.add_command(label="volume display", command=VTK_vol)
"""
### ROI menu
roimenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="Structures", menu=roimenu)
roimenu.add_command(label="Import structures...", command=open_ROI)
roimenu.add_command(label="Compute DVH...", command=ROI_DVH_analysis,state=DISABLED)
roimenu.add_command(label="Crop dosimetry", command=CropDosi,state=DISABLED)
roimenu.add_command(label="Normalize dosi to PTV dose", command=Rescale,state=DISABLED)

### More menu
helpmenu = Menu(menubar, tearoff=0)
menubar.add_cascade(label="More", menu=helpmenu)
helpmenu.add_command(label="File infos", command=FileInfo)
helpmenu.add_command(label="About...", command=about)

menubar.add_command(label="Demo", command=Demo)

print ''
Update()
window.mainloop()
