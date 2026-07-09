# -*- coding: utf-8 -*-
################################################
# Visualize dose and dose-rate with this viewer
# 		            ---
#           P. Lansonneur 2020
################################################

import os, sys
import numpy as np
import pydicom, pydicom.uid
from scipy.ndimage import map_coordinates
from functools import partial
import qdarkstyle
import argparse
from copy import deepcopy
from contextlib import suppress
from PySide6.QtWidgets import *
from PySide6 import QtCore, QtGui
from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as P
import matplotlib as mpl
from natsort import natsorted, ns
import multiprocessing as mp
# from RP2TOPAS_IBA import CreateTopasFilePBS
from RP2TOPAS_ProBeam import CreateTopasFilePBS
import cv2
import ctypes
import json

ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID('Varian.FLASH.Viewer.1.0')
mpl.rcParams['savefig.dpi'] = 200
mpl.rc('font', **{'family':'sans-serif','weight':'normal','size':10})
    
CURRENT_PATH = os.path.dirname(os.path.abspath(sys.argv[0])) if getattr(sys, 'frozen', False) else os.path.dirname(os.path.abspath(__file__))
ICONS_PATH = os.path.join(CURRENT_PATH, 'icons')

def AForBBetween(a, b, bmin=None, bmax=None):
    """ return array a with b between bmin and bmax """
    
    mask = np.ones_like(b, dtype=bool)
    if bmin: mask &= b > bmin
    if bmax: mask &= b < bmax
    return np.where(mask, a, 0.0)

def openPixmap(fname, color=QtGui.QColor(255, 255, 255)):

    base = QtGui.QPixmap(fname)
    result = QtGui.QPixmap(base.size())
    result.fill(QtCore.Qt.transparent)
    painter = QtGui.QPainter(result)
    painter.drawPixmap(0, 0, base)  # draw original image (for alpha mask)
    painter.setCompositionMode(QtGui.QPainter.CompositionMode_SourceIn)
    painter.fillRect(result.rect(), color)  # apply color while keeping alpha
    painter.end()
    return result

def openIcon(fname, color=QtGui.QColor(255, 255, 255)): return QtGui.QIcon(openPixmap(fname, color)) # default: white

class Cosmo():
    """ Cosmo image class """
    
    def __init__(self, parent, orientation, path):        
        
        self.parent = parent
        self.img = [P.imread(path + rf'\\CosmoImages\\{i}_{orientation}.PNG') for i in ['Sagittal','Frontal','Transversal']]
    
    def extent(self,ax):
        
        self.axis = self.parent.ax[ax]
        xlim, ylim = self.axis.properties().get('xlim'), self.axis.properties().get('ylim')
        img = self.img[ax]
        ratio = np.shape(img)[1] / np.shape(img)[0] # aspect ratio
        height = abs(ylim[1] - ylim[0]) # axis height
        h = (0.14 if ratio<=1 else 0.07) * height # cosmo height
        w = ratio * h # cosmo width
        off = 0.03 * height # offset
        x0, x1 = xlim[0] + off, xlim[0] + off + w
        y0, y1 = ylim[0] + off, ylim[0] + off + h

        # Handle reversed axes
        if xlim[1] < xlim[0]: x0, x1 = xlim[0] - off, xlim[0] - off - w
        if ylim[1] < ylim[0]: y0, y1 = ylim[0] - off, ylim[0] - off - h

        return [x0, x1, y0, y1]
        
class RTGeneralVolume():
    """ Generic class for volumes """
    
    def __init__(self):

        self.is_opened = False # is opened?
        self.origin = [0, 0, 0] # origin
        self.spacing = [1, 1, 1] # spacing
        self.dim = [1, 1, 1] # dimensions
        self.volume = np.zeros((self.dim)) # volume
        self.orientation = 'HFS' # orientation
        self.patientName = '' # patient's name
        self.label = '' # label (colormap title)
        self.cmin, self.cmax = None, None # colorbar range
    
    def i2mm(self, i, ax): return self.origin[ax] + i * self.spacing[ax]
    
    def mm2i(self, z, ax): return round((z - self.origin[ax]) / self.spacing[ax])
    
    def image(self, ax, idx): return np.take(self.volume, idx, axis=[2, 1, 0][ax])
    
    def extent(self, ax): return [ext for i in range(3) if i != ax for ext in [self.i2mm(j,i) - 0.5*self.spacing[i] for j in [0, self.dim[i]]]] # image extent [mm]
        
    def mm2voxel(self, x, y, ax): return [self.mm2i(z,i) for z,i in zip([y,x],[k for k in [2,1,0] if k != ax])] # mm to voxel
    
    def bbox2RSExtentmm(self, bbox, ax): return [ext for i in range(3) if i != ax for ext in [self.i2mm(j,i) - 0.5*self.spacing[i] for j in [min(bbox[i]), max(bbox[i])]]] # bounding box limits to mm
        
    def RSExtentmm2bbox(self, RS, ax=0):
        """ get the bounding box from Rectangle Selector extent (mm) """
        
        if ax==0: ymin, ymax, zmin, zmax = RS[0].extents; xmin, xmax, _, __ = RS[2].extents
        if ax==1: xmin, xmax, zmin, zmax = RS[1].extents; ymin, ymax, _, __ = RS[0].extents
        if ax==2: xmin, xmax, ymin, ymax = RS[2].extents; _, __, zmin, zmax = RS[1].extents

        return [[np.clip(self.mm2i(z,i), 0, self.dim[i]) for z in [[xmin, xmax], [ymin, ymax], [zmin, zmax]][i]] for i in range(3)]
           
class RTCTVolume(RTGeneralVolume):
    """ Generic class for CT  """
    
    def __init__(self):
        super().__init__()
        self.cmap = P.get_cmap('Greys_r') # colormap
        self.volume = -1000. * np.ones((self.dim)) # volume
        self.label = 'Hounsfield Unit (HU)'
        
class RTDoseVolume(RTGeneralVolume):
    """ Generic class for dose and dose-rate """
    
    def __init__(self, label=''):
        
        super().__init__()
        self.label = label
        self.cmap = P.get_cmap('jet') # colormap

    def updateColorbarRange(self, cmin=None, cmax=None):
        """ update colorbar range """
        
        self.cmax = np.percentile(self.volume, 99.999) if cmax is None else cmax
        self.cmin = 0.01 * self.cmax if cmin is None else cmin
        
    def slice_index(self, CT, ax, i):
        """ return the slice index along axis 'ax' corresponding to the i-th CT slice.
        One one hand, z = CT.origin + i * CT.spacing. On the other hand, z = D.origin  + j * D.spacing """
        
        return self.mm2i(CT.i2mm(i,ax), ax)

    def Open(self, filename): 
        """ Open a RT-Dose DICOM file """
        
        if not filename: return
        
        if filename.endswith('.dcm') or filename.endswith('.data'):

            ds = pydicom.dcmread(filename, force=True)
            ds.file_meta.TransferSyntaxUID = pydicom.uid.ImplicitVRLittleEndian
            self.spacing = [float(ds.PixelSpacing[1]), float(ds.PixelSpacing[0]), abs(float(ds.GridFrameOffsetVector[1] - ds.GridFrameOffsetVector[0]))]
            self.volume = float(ds.DoseGridScaling) * ds.pixel_array
            self.dim = list(self.volume.shape)[::-1]
            self.origin = ds.ImagePositionPatient
            self.patientName = ds.PatientID
            
            print(f'\nfilename:\t{filename}\ndimensions:\t{self.dim}\norigin (mm):\t{self.origin}\nspacing (mm):\t{self.spacing}')
            
            X = ds.ImageOrientationPatient[0:3]
            Y = ds.ImageOrientationPatient[3:6]
            R = np.array([X, Y, np.cross(X, Y)])
            if np.any(R[~np.eye(R.shape[0], dtype=bool)] > 1E-12): print(f'\nOrientation is not supported:\n{R}'); return
            
            self.is_opened = True
   
class RTPlan():
    """ A RTPlan contains several RTFields """

    def __init__(self, parent, fieldIDs=None):

        self.parent = parent
        self.is_opened = False # is opened?
        self.selectedFields = [] # list of selected fields
        self.fields = [RTField(self,i) for i in fieldIDs] if fieldIDs else [] # list of fields
        
    def OpenDicom(self, filepath=None):
        """ Open a RT-Plan DICOM file """
        
        filepath = self.parent.OpenDicom('Plan', filepath)
        if not filepath: return
        
        self.ds = pydicom.dcmread(filepath)
        
        self.fields = [] # clear list of fields
        for f in self.ds.IonBeamSequence:
        
            X = [ x for cp in f.IonControlPointSequence[::2] for x in cp.ScanSpotPositionMap[0::2] ]
            Y = [ y for cp in f.IonControlPointSequence[::2] for y in cp.ScanSpotPositionMap[1::2] ]
            E = [ e for cp in f.IonControlPointSequence[::2] for e in cp.NumberOfScanSpotPositions*[float(cp.NominalBeamEnergy)] ]
            W = [ w for cp in f.IonControlPointSequence[::2] for w in cp.ScanSpotMetersetWeights ]
            
            field = RTField(self, f.BeamName)
            cp = f.IonControlPointSequence[0]
            field.GantryAngle = float(cp.GantryAngle)
            field.PatientSupportAngle = float(cp.PatientSupportAngle)
            field.SnoutPosition = 0.1*float(cp.SnoutPosition)
            field.IsocenterPosition = [float(i) for i in cp.IsocenterPosition]
            if f.NumberOfRangeShifters > 0:
                field.RSWET = float(cp.RangeShifterSettingsSequence[0].RangeShifterWaterEquivalentThickness)
                field.RSID = f.RangeShifterSequence[0].RangeShifterID
            field.SnoutID = f.SnoutSequence[0].SnoutID
            field.spots = RTSpotlist(field, X, Y, E, W)
            self.fields.append(field)
            
        self.is_opened = True
        self.parent.toolmenu.actions()[5].setEnabled(True) # Export to TOPAS
        
    def exportTopasFilePBS(self):
        """ create TOPAS MC simulation file """
        
        filedir = str(QFileDialog.getExistingDirectory(self.parent, "Select Export Directory", self.parent.datadir))
        if not filedir: return
        self.parent.datadir = filedir
        for f in self.fields: CreateTopasFilePBS(f, self.parent.CT, self.parent.bbox, filedir)
        
    def ImportDoseAndDoseRate(self, filedir=None):
        """ open dose and dose-rate files from directory """
        
        if not filedir: filedir = str(QFileDialog.getExistingDirectory(self.parent, "Select Dose Rate Directory", self.parent.datadir))
        if not filedir: return
        self.parent.datadir = filedir
        
        for f in self.fields:
            for vol, name in zip([f.Dose, f.DoseRate], [f'\RD.Dose_{f.ID}', f'\RD.DoseRate_{f.ID}']):
                # filename = filedir + name + '.data'
                filename = filedir + name + '.dcm'
                if not os.path.isfile(filename): self.parent.ThrowError(f'{filename} was not found'); return
                vol.Open(filename)
                
        self.is_opened = True
        self.parent.tabWidget.ClearAllTabs()
        self.parent.Data = deepcopy(self.fields[0].Dose)
        self.parent.Data.volume = sum(f.Dose.volume for f in self.fields)
        for i in range(4): self.parent.displaymenu.actions()[i].setEnabled(True) # CT + dose/doserate/overlaps
        self.parent.dispopt_menu.actions()[1].setEnabled(True) # dose-rate range
        self.parent.fieldmenu.actions()[0].trigger() # All fields
        self.parent.displaymenu.actions()[1].trigger() # CT + dose
        self.parent.update_all_axes()
        
class RTField():
    """ A field contains one dose and one dose-rate RTDoseVolume """

    def __init__(self, parent = None, ID=' '):

        self.parent = parent # RTPlan
        self.is_opened = False # is opened?
        self.ID = ID # field ID
        self.Dose = RTDoseVolume('Dose (Gy)') # dose matrix
        self.DoseRate = RTDoseVolume('Dose Rate (Gy/s)') # dose-rate matrix
        self.spots = RTSpotlist(self) # spots list
        self.GantryAngle = 0 # Gantry angle (deg)
        self.PatientSupportAngle = 0 # Patient support angle (deg)
        self.SnoutPosition = 0 # Snout Position (mm)
        self.IsocenterPosition = [0,0,0] # Isocenter Position (mm)
        self.RSWET = 0 # Range Shifter Water Equivalent Thickness (mm)
        self.RSID = '' # Range Shifter ID
        self.SnoutID = '' # Snout ID
        
class RTSpotlist():
    """ spot list information """

    def __init__(self, parent, x=[], y=[], E=[], MU=[]):
        
        self.parent = parent
        
        self.x = x # x coordinates (mm)
        self.y = y # y coordinates (mm)
        self.E = E # Energy (MeV)
        self.MU = MU # Monitor Units (MU)
        self.N = len(x) # number of spots
        
def get_roi_binary_mask(roi, dim, origin, spacing, contours=None, dx=0., dy=0.):
    """ Efficiently calculate the binary mask for a ROI with possibly multiple contours. """
    
    idx = np.zeros(dim[::-1], dtype=bool)
    
    contours = contours or roi.contours
    
    if len(contours) > 1: # avoid markers

        for cs in contours:
            x = (np.asarray(cs.x) - origin[0]) / spacing[0] + 0.5 + dx
            y = (np.asarray(cs.y) - origin[1]) / spacing[1] + 0.5 + dy
            z = round( (cs.z - origin[2]) / spacing[2] )
            if 0 <= z < dim[2]: idx[z] ^= cv2.fillPoly(np.zeros((dim[1], dim[0]), dtype=np.uint8), [(np.c_[x,y]).astype(np.int32)], 1).astype(bool)
            
    return idx
    
def get_roi_fraction_mask(roi, dim, origin, spacing, contours=None):
    
    return sum([ get_roi_binary_mask(roi, dim, origin, spacing, contours, dx, dy) for dx in [-1,-0.5,0,0.5,1] for dy in [-1,-0.5,0,0.5,1]]).astype(float) / 25
    
class RTROISet():
    """ Structure Set """
    
    def __init__(self, parent=None):
        
        self.parent = parent # main window
        self.is_opened = False # is opened?
        self.ROIs = [] # ROIs
        self.selected = [] # list of ROI selected
        
    def OpenDicom(self, filepath=None):        
        """ Open a RT-struct file with .dcm extension """
        
        filepath = self.parent.OpenDicom('Plan', filepath)
        if not filepath: return
        print(f'\nOpening {filepath} ...\n')
        
        self.ROIs = [] # clear list of ROIs
        self.ds = pydicom.dcmread(filepath)
            
        for rs, rcs in zip(self.ds.StructureSetROISequence, self.ds.ROIContourSequence):
            print(f'- {rs.ROIName}')
            color = [c / 255. for c in rcs.ROIDisplayColor]
            contours = [ RTROIContour([float(xyz) for xyz in cs.ContourData]) for cs in rcs.ContourSequence ] if 'ContourSequence' in rcs else []
            self.ROIs.append( RTROI(rs.ROIName, color, contours) )
            
        if self.parent.CT.is_opened:
            fun = partial(get_roi_binary_mask, dim=self.parent.CT.dim, origin=self.parent.CT.origin, spacing=self.parent.CT.spacing, contours=None, dx=0., dy=0.)
            self.ROIs = sorted(self.ROIs, key=lambda roi: -len(roi.contours)) # treat large ROIs first
            with mp.Pool() as pool: masks = pool.map(fun, self.ROIs)
            for mask,roi in zip(masks, self.ROIs): roi.mask = mask
        
            # #####
            # ptv = next(roi.mask for roi in self.ROIs if roi.name == "PTV")
            # oar = next(roi.mask for roi in self.ROIs if roi.name == "Parotid Left")
            # inter = np.sum(np.logical_and(ptv, oar))
            # print(f'\nPTV volume:\t{np.sum(ptv)} voxels')
            # print(f'Parotid Left:\t{np.sum(oar)} voxels')
            # print(f'Intersection:\t{inter} voxels')
            # print(f'Intersection / PTV: {100.0 * inter / np.sum(ptv) :.2f} %')            
            
        self.is_opened = True
        self.ROIs = sorted(self.ROIs, key=lambda roi: roi.name) # alphabetical order            
        self.parent.fill_ROImenu()            
        self.selected = self.ROIs
        self.parent.update_all_axes()
    
class RTROI():
    """ Structure """
    
    def __init__(self, name, color, contours):
        
        self.name = name # name
        self.color = color # color in RGB
        self.contours = contours # list of contours
        self.dvh = dict() # DVH
        self.drvh = dict() # DRVH
        self.mask = None
        
class RTROIContour():

    def __init__(self, xyz=[0,0,0]): self.x, self.y, self.z = xyz[0::3] + [xyz[0]], xyz[1::3] + [xyz[1]], xyz[2] # vertices coordinates [mm]

class RTMainWindow(QMainWindow):
    """ Main window """
    
    def __init__(self, args=None):

        super().__init__()
        
        self.args = args
        # print(self.args)        
        self.plan = RTPlan(self, self.args.fields) # Plan
        self.CT = RTCTVolume() # CT
        self.Data = RTDoseVolume('Dose (Gy)') # Dose/Dose-rate/overlap
        self.myROISet = RTROISet(self) # Structure set
        self.profile_show = False # profile tool enabled
        self.ctrl = False # is control key pressed?
        self.drmin = None # min dose rate threshold [Gy/s]
        self.drmax = None # max dose rate threshold [Gy/s]
        self.Dmin = 0.1 # min dose threshold [Gy]
        self.fig = [] # list of figures
        self.ax  = [] # list of axes
        self.w = [] # list of sliders
        self.moved = [] # did the sliders move?
        self.x_lim = [(0,1),(0,1),(0,1)] # x axes limits
        self.y_lim = [(0,1),(0,1),(0,1)] # y axes limits
        self.RS = [] # list of rectangle selectors
        self.bbox = [] # calculation volume
        self.datadir = '' # last known data directory
        
        self.setWindowTitle('Dose Rate Viewer')
        self.createGridLayout()
        self.createMenuBar()
        self.statusBar()
        
    def createGridLayout(self):

        for i in range(3):
            
            ### axial views
            fig = FigureCanvas(Figure())
            ax = fig.figure.gca()
            ax.axis('equal')
            ax.xaxis.set_visible(False)
            ax.yaxis.set_visible(False)
            for feat in ['top','bottom','left','right']: ax.spines[feat].set_linewidth(0)
            fig.figure.subplots_adjust(left=-1E-3, bottom=0, right=1, top=1, wspace=0, hspace=0)

            fig.figure.canvas.mpl_connect('scroll_event', self.on_mousewheel)
            fig.figure.canvas.mpl_connect('button_press_event', self.on_click)
            fig.figure.canvas.mpl_connect('motion_notify_event', self.on_motion)
            fig.figure.canvas.mpl_connect('button_release_event', self.on_click_released)
            
            ### sliders
            w = QSlider(QtCore.Qt.Horizontal, minimum=0, maximum=self.CT.dim[i]-1)
            w.valueChanged.connect(lambda val, i=i: self.update_axis(i, val))
            rs = mpl.widgets.RectangleSelector(ax, self.onRSselect)
            
            self.fig.append(fig)
            self.ax.append(ax)
            self.w.append(w)
            self.moved.append(True)
            self.RS.append(rs)
        
        ### colorbar
        self.cb = FigureCanvas(Figure())
        self.cb.figure.subplots_adjust(left=0.25, bottom=0.1, right=0.35, top=0.9)
        
        ### tab widget
        self.tabWidget = MyTabWidget()
        self.tabWidget.dvh.figure.canvas.mpl_connect('button_press_event', self.on_click)
        self.tabWidget.drvh.figure.canvas.mpl_connect('button_press_event', self.on_click)
        self.tabWidget.ddrvh.figure.canvas.mpl_connect('button_press_event', self.on_click)
        self.tabWidget.profile.figure.canvas.mpl_connect('button_press_event', self.on_click)
        
        layout = QGridLayout()
        layout.setVerticalSpacing(5)
        layout.setHorizontalSpacing(0)
        layout.setColumnStretch(0,10) # (col index, stretch factor)
        layout.setColumnStretch(1,1)
        layout.setColumnStretch(2,10)
        layout.setRowStretch(1,1)
        layout.setRowStretch(2,1)
        
        layout.addWidget(self.fig[2],1,0,1,1) # (int fromRow, int fromColumn, int rowSpan, int columnSpan)
        layout.addWidget(self.fig[1],1,2,1,1)
        layout.addWidget(self.fig[0],2,0,1,2)
        layout.addWidget(self.cb,1,1,1,1)
        layout.addWidget(self.w[2],0,0,1,2)
        layout.addWidget(self.w[1],0,2,1,2)
        layout.addWidget(self.w[0],3,0,1,2)
        layout.addWidget(self.tabWidget,2,2)

        wid = QWidget(self)
        self.setCentralWidget(wid)
        wid.setLayout(layout)
        wid.setFocusPolicy( QtCore.Qt.StrongFocus )
        wid.setFocus()
        
    def createMenuBar(self):
        
        def icon(name): return openIcon(os.path.join(ICONS_PATH, f'{name}.png'))
        
        menubar = self.menuBar()
        
        ### File menu
        fileMenu = menubar.addMenu('File ')
        fileMenu.addAction(QtGui.QAction(icon('folder_open'), 'Open a DICOM Serie', self, shortcut="Ctrl+O", triggered=self.OpenDicomSerie))
        fileMenu.addAction(QtGui.QAction(icon('roi'), 'Import Structures', self, triggered=self.myROISet.OpenDicom))
        fileMenu.addAction(QtGui.QAction(icon('plan'), 'Import Plan', self, triggered=self.plan.OpenDicom))
        if self.args.fields: fileMenu.addAction(QtGui.QAction(icon('file_open'), 'Import Dose and Dose-Rate', self, triggered=self.plan.ImportDoseAndDoseRate))
        else: fileMenu.addAction(QtGui.QAction(icon('file_open'), 'Import Dose', self, triggered=self.OpenDose))
        fileMenu.addSeparator()
        fileMenu.addAction(QtGui.QAction(icon('logout'), 'Exit', self, shortcut="Ctrl+Q",triggered=QApplication.instance().quit))
                
        ### Tool menu
        self.toolmenu = menubar.addMenu(' Tools ')
        histograms = ['DVH','DRVH','DDRVH','DR-DVH'] if self.args.fields else ['DVH']
        for i in histograms: self.toolmenu.addAction(QtGui.QAction(icon('bar_chart'), i,  self, enabled=True, triggered=lambda c, s=i: self.showDVH(s)))
        if self.args.fields: self.toolmenu.addAction(QtGui.QAction(icon('overlap'), 'Overlap Volume', self, enabled=True, triggered=self.getOverlapVolume))
        self.toolmenu.addSeparator()
        self.toolmenu.addAction(QtGui.QAction(icon('export'), 'Export to TOPAS', self, enabled=False, triggered=self.export2TOPAS))

        ### structures menu
        self.roimenu = menubar.addMenu(' Structures ')
        self.roimenu.addAction(QtGui.QAction('Unselect All', self, checkable=True, triggered=self.ROI_unselectAll))
        self.roimenu.addSeparator()
        self.roimenu.triggered.connect(self.ROI_selection)

        ### Display option menu
        self.dispopt_menu = menubar.addMenu(' Options ')   
        self.ctscale_menu = QMenu('CT Windowing', self)
        ag1 = QtGui.QActionGroup(self)
        for i in ['Auto','Brain','Lungs','Abdomen','Bone']: self.ctscale_menu.addAction(ag1.addAction(QtGui.QAction(i, self, checkable=True, triggered=lambda c, s=i: self.SetCTWindowing(s))))
        self.dispopt_menu.addMenu(self.ctscale_menu)        
        if self.args.fields: self.dispopt_menu.addAction(QtGui.QAction('Dose-Rate Range...', self, triggered=self.update_doserate_bounds, enabled=False))
        self.dispopt_menu.addAction(QtGui.QAction('Colorbar Range ...', self, triggered=self.ChangeColorbarRange))
        
        ### View menu
        ag2 = QtGui.QActionGroup(self)
        self.displaymenu = menubar.addMenu(' View ')
        options = ['CT Only', 'CT + Dose', 'CT + Dose-Rate', 'CT + overlaps'] if self.args.fields else ['CT Only', 'CT + Dose']
        for i in options: isCT=(i=='CT Only'); self.displaymenu.addAction(ag2.addAction(QtGui.QAction(i, self, checkable=True, checked=isCT, enabled=isCT)))
        self.displaymenu.triggered.connect(self.ChangeDisplayOption)
        
        ### Fields menu
        if self.args.fields:
            self.fieldmenu = menubar.addMenu(' Fields ')
            ag3 = QtGui.QActionGroup(self)
            self.fieldmenu.addAction(ag3.addAction(QtGui.QAction('Select All', self, checkable=True)))
            for f in self.args.fields: self.fieldmenu.addAction(ag3.addAction(QtGui.QAction(f, self, checkable=True)))            
            self.fieldmenu.triggered.connect(self.Field_selection)

        ### Infos, Help and About actions
        menubar.addAction(QtGui.QAction(' Infos ', self, triggered=self.Infos))
        menubar.addAction(QtGui.QAction(' ? ', self, triggered=self.About))
        
        ### crop action
        self.cropAction = QtGui.QAction(icon('crop'), '', self, enabled=False, checkable=True, triggered=self.update_all_axes)
        menubar.addAction(self.cropAction)
        
        ### crosshair action
        self.crosshairAction = QtGui.QAction(icon('crosshair'), '', self, checkable=True, triggered=self.update_all_axes)
        menubar.addAction(self.crosshairAction)
        
    def display_option(self): return [a.isChecked() for a in self.displaymenu.actions()].index(True)
        
    def OpenDicomSerie(self, filedir=None):
        """ Open a dicom serie """
        
        if not filedir: filedir = str(QFileDialog.getExistingDirectory(self, "Select Directory", self.datadir))
        if not filedir: return
        self.datadir = filedir
        
        dicoms = [pydicom.dcmread(os.path.join(filedir,f)) for f in os.listdir(filedir) if f.endswith(".dcm")]
        dicoms = [ds for ds in dicoms if ds.Modality=='CT']
        dicoms = sorted(dicoms, key=lambda ds: float(ds.ImagePositionPatient[2]))
        
        self.CT.volume = np.asarray([ds.pixel_array for ds in dicoms], dtype=float)
        self.CT.dim = [*dicoms[0].pixel_array.shape, len(dicoms)]
        self.CT.spacing = [float(s) for s in [*dicoms[0].PixelSpacing, dicoms[0].SliceThickness]]
        self.CT.origin = [float(ipp) for ipp in dicoms[0].ImagePositionPatient]
        self.CT.patientName = dicoms[0].PatientID
        self.CT.orientation = dicoms[0].PatientPosition
        if "RescaleSlope" in dicoms[0]: self.CT.volume *= float(dicoms[0].RescaleSlope)
        if "RescaleIntercept" in dicoms[0]: self.CT.volume += float(dicoms[0].RescaleIntercept)
        
        if self.CT.volume.ndim != 3: self.ThrowError(f'Only 3D patient images are supported (image dimension: {self.CT.volume.ndim}).'); return
        
        if not self.CT.orientation in ["HFS", "HFP", "FFS", "FFP"]: self.ThrowError(f'Orientation is not supported\n{R}'); return
        
        print(f'\nCT serie:\t{filedir}\ndimensions:\t{self.CT.dim}\norigin (mm):\t{self.CT.origin}\nspacing (mm):\t{self.CT.spacing}\norientation:\t{self.CT.orientation}')
        
        self.CT.is_opened = True
        self.cosmo = Cosmo(self, self.CT.orientation, CURRENT_PATH)
        self.bbox = [[0, self.CT.dim[i]] for i in range(3)]
        self.x_lim = [(self.CT.extent(i)[:2]) for i in range(3)]
        self.y_lim = [(self.CT.extent(i)[2:]) for i in range(3)]
        self.ctscale_menu.actions()[0].trigger() # 'Auto' scale
        self.cropAction.setEnabled(True) # Edit calculation volume
        self.set_sliders_limits()
        
        match self.CT.orientation: # foot-to-head convention
            case 'HFS': self.invert_axis(2,'y');  self.invert_axis(0,'x')
            case 'HFP': self.invert_axis(2,'y');  self.invert_axis(0,'x')
            case 'FFS': self.invert_axis(1,'xy'); self.invert_axis(0,'x')
            case 'FFP': self.invert_axis(0,'xy'); self.invert_axis(1,'xy'); self.invert_axis(2,'xy')
            
        self.update_all_axes()

    def OpenDicom(self, what, filepath=None): 
        """ Open a dose file with a .dcm extension """
        
        with suppress(Exception):
            if not filepath: filepath, _ = QFileDialog.getOpenFileName(self, f"Import {what}", self.datadir, "DICOM (*.dcm)")
            if not filepath: return
            self.datadir = os.path.dirname(filepath)
        
        if not (filepath.endswith('.dcm') or filepath.endswith('.data')): self.ThrowError(f'Extension not supported'); return
        
        return filepath 
        
    def OpenDose(self, filepath=None): 
        """ Open a dose file with a .dcm extension """
        
        filepath = self.OpenDicom('Dose', filepath)
        if not filepath: return
        
        self.Data.Open(filepath)
        self.Data.updateColorbarRange()
        self.plan.is_opened = True
        self.displaymenu.actions()[1].setEnabled(True)
        self.displaymenu.actions()[1].trigger()
        self.tabWidget.ClearAllTabs()
        self.update_all_axes()
        
    def export2TOPAS(self):
        
        if self.CT.is_opened and self.plan.is_opened: self.plan.exportTopasFilePBS()

    def update_doserate_bounds(self): DRDialog(self).exec()
            
    def Field_selection(self):
        """ select which field to display """
        
        field_idx = [a.isChecked() for a in self.fieldmenu.actions()].index(True) - 1
        self.plan.selectedFields = self.plan.fields if field_idx < 0 else [self.plan.fields[field_idx]]

        match self.display_option():
            
            case 1: # CT + dose
                self.Data.volume = sum(AForBBetween(f.Dose.volume, f.DoseRate.volume, self.drmin, self.drmax) for f in self.plan.selectedFields)
                
                if not self.drmin and not self.drmax: self.Data.label = 'Dose (Gy)'
                elif not self.drmax: self.Data.label = f'Dose (Gy), DR > {self.drmin} Gy/s'
                elif not self.drmin: self.Data.label = f'Dose (Gy), DR < {self.drmax} Gy/s'
                else: self.Data.label = f'Dose (Gy), {self.drmin} < DR < {self.drmax} Gy/s'

            case 2: # CT + dose rate
                if len(self.plan.selectedFields) > 1:
                    self.ThrowError('Dose-rate cannot be displayed for multiple fields.')
                    self.fieldmenu.actions()[1].trigger() # select first field
                    return
                else:
                    f = self.plan.selectedFields[0]
                    self.Data = deepcopy(f.DoseRate)
                    self.Data.volume = AForBBetween(f.DoseRate.volume, f.Dose.volume, self.Dmin, None)
                    self.Data.label = f'Dose Rate (Gy/s), D > {self.Dmin:.2f} Gy'
                    
            case 3: # CT + overlaps
                self.Data.volume = sum(np.where(f.Dose.volume > self.Dmin, 1, 0) for f in self.plan.selectedFields)
                self.Data.label = f'overlap, D > {self.Dmin:.2f} Gy'
            
        self.tabWidget.ClearAllTabs()
        self.Data.updateColorbarRange()
        self.update_all_axes()
        
    def fill_ROImenu(self):
        """ fill the ROI menu """

        for act in self.roimenu.actions()[1:]: self.roimenu.removeAction(act) # clear the menu
        
        for roi in self.myROISet.ROIs:
            
            icon = QtGui.QIcon()
            icon.addPixmap(openPixmap(ICONS_PATH + "\\square_FILL1.png", QtGui.QColor.fromRgbF(*roi.color)), QtGui.QIcon.Normal, QtGui.QIcon.On)
            icon.addPixmap(openPixmap(ICONS_PATH + "\\square_FILL0.png", QtGui.QColor(125, 125, 125)), QtGui.QIcon.Normal, QtGui.QIcon.Off)
            action = self.roimenu.addAction(icon, roi.name)
            action.setCheckable(True)
            action.setChecked(True)            
            
    def ROI_selection(self, action):
        """ triggered when ROIs are selected from menu """
        
        self.myROISet.selected = [ roi for roi,action in zip(self.myROISet.ROIs, self.roimenu.actions()[1:]) if action.isChecked() ]
        self.update_all_axes()
        self.roimenu.exec()

    def ROI_unselectAll(self):
        """ (un)select all the ROI from the ROI menu """

        if self.roimenu.actions()[0].isChecked():
            for act in self.roimenu.actions()[1:]: act.setChecked(False)            
            self.myROISet.selected = []
            self.roimenu.actions()[0].setText('Select All')
        else:
            for act in self.roimenu.actions()[1:]: act.setChecked(True)
            self.myROISet.selected = self.myROISet.ROIs
            self.roimenu.actions()[0].setText('Unselect All')
            
        self.update_all_axes()

    def ChangeDisplayOption(self):
        """ """
        
        if (self.args.fields is not None) and self.plan.is_opened: self.Field_selection()
        self.update_all_axes()

    def on_click(self,event):
        """ run when mouse button is clicked """
        
        if event.button == 3: # right click
            popMenu = QMenu(self) # show the context menu at the cursor position
            popMenu.addAction(QtGui.QAction('Export As ...', self, triggered=partial(self.export, event.canvas) ))
            popMenu.popup(QtGui.QCursor().pos())
            
        elif self.ctrl and not self.cropAction.isChecked():
            self.xp1, self.yp1, self.profile_show, self.profile_ax = event.xdata, event.ydata, True, event.inaxes
            ax = self.ax.index(event.inaxes)
            event.inaxes.plot(self.xp1, self.yp1, 'C0o')
            self.set_axis_limit(ax)

            with suppress(Exception): self.fig[ax].draw()

        else:
            self.profile_show = False
            self.tabWidget.ClearTab(self.tabWidget.profile)
            if not self.cropAction.isChecked(): self.update_all_axes()
        
    def on_motion(self,event):
        """ run when mouse is moved """
        
        if event.inaxes is None: return

        ### update profile
        if self.ctrl and self.profile_show and event.button == 1:

            if len(event.inaxes.lines)>0: event.inaxes.lines[-1].remove()
            event.inaxes.plot([self.xp1,event.xdata], [self.yp1,event.ydata], 'C0o:')
            self.fig[self.ax.index(event.inaxes)].draw()
            self.xp2, self.yp2 = event.xdata, event.ydata
            self.updateProfile(event.inaxes)

        ### update statusbar
        self.statusBar().clearMessage()
        with suppress(Exception):
            
            vol = self.Data if self.display_option() > 0 else self.CT
            unit = {0:'HU', 1:'Gy', 2:'Gy/s', 3:''}[self.display_option()]
            
            ax = self.ax.index(event.inaxes)
            x, y = vol.mm2voxel(event.xdata, event.ydata, ax) # voxels coordinates
            pos = self.w[ax].value() if vol is self.CT else vol.slice_index(self.CT, ax, self.w[ax].value())
            if x >= 0 and y >= 0: self.statusBar().showMessage(f'x: {event.xdata:.1f}, y: {event.ydata:.1f} [{vol.image(ax, pos)[x, y]:.1f} {unit}]')

    def on_click_released(self,event):
        """ run when mouse click is released """
        
        if event.button == 3: return # right click
        
        if self.profile_show:
            self.xp2, self.yp2 = event.xdata, event.ydata
            self.updateProfile(event.inaxes)
            self.update_all_axes()

    def keyPressEvent(self,e): self.ctrl = (e.key() == QtCore.Qt.Key_Control)
            
    def keyReleaseEvent(self,e): self.ctrl = False
        
    def on_mousewheel(self,e): 
    
        if self.ctrl: self.ZoomOnAxis(e)
        else: self.w[self.ax.index(e.inaxes)].setValue(int(self.w[self.ax.index(e.inaxes)].value() + e.step))

    def ZoomOnAxis(self, event):
        
        ### adjust axes limit
        scale = 1 / (1 + 0.1*event.step) if event.step > 0 else 1 - 0.1*event.step
        x, y, ax = event.xdata, event.ydata, self.ax.index(event.inaxes)
        xL =  x - self.x_lim[ax][0] # distance to the left
        xR = -x + self.x_lim[ax][1] # distance to the right
        yB =  y - self.y_lim[ax][0] # distance to the bottom
        yT = -y + self.y_lim[ax][1] # distance to the top
        self.x_lim[ax] = (x - xL*scale, x + xR*scale)
        self.y_lim[ax] = (y - yB*scale, y + yT*scale)
        self.set_axis_limit(ax)
        
        ### redraw cosmo
        if event.inaxes.get_images()[-1].properties().get('zorder')==4:
            event.inaxes.get_images()[-1].set_visible(False)
            event.inaxes.imshow(self.cosmo.img[ax], extent=self.cosmo.extent(ax), zorder=4)

        self.fig[ax].draw()
        
    def export(self, canvas, checked=None):
                
        filename, _ = QFileDialog.getSaveFileName(self, 'Save File', self.datadir, "PNG (*.png)" if canvas in self.fig else "JSON (*.json);;PNG (*.png)")
        if not filename: return

        if filename.endswith('.png'): canvas.figure.savefig(filename)
        
        else:
            with open(filename, 'w') as file:
                match canvas:
                    case self.tabWidget.dvh: data = {roi.name: roi.dvh for roi in self.myROISet.selected}
                    case self.tabWidget.drvh: data = {roi.name: roi.drvh for roi in self.myROISet.selected}
                    case self.tabWidget.profile: data = {'x': self.profile[0].tolist(),'y': self.profile[1].tolist()}
                file.write(json.dumps(data, sort_keys=True))
                    
    def resizeEvent(self, event): self.update_all_axes()

    def set_sliders_limits(self):
        
        for i in range(3):
            self.w[i].setMaximum( self.CT.dim[i] - 1 )
            self.w[i].setValue( self.CT.dim[i] / 2 )

    def set_axis_limit(self, i): self.ax[i].set(xlim=self.x_lim[i], ylim=self.y_lim[i])
        
    def invert_axis(self, ax=0, what=''):
        """ invert axis x/y limits """
        
        if 'x' in what: self.x_lim[ax] = self.x_lim[ax][::-1]
        if 'y' in what: self.y_lim[ax] = self.y_lim[ax][::-1]
        self.set_axis_limit(ax)

    def update_all_axes(self):
        """ Update all axes and colorbar """
        
        for i in range(3): self.update_axis(i, self.w[i].value())
            
        self.cb.figure.gca().cla()
        vol = self.Data if self.display_option() > 0 else self.CT
        self.cb.figure.colorbar(mpl.cm.ScalarMappable(mpl.colors.Normalize(vol.cmin, vol.cmax), vol.cmap), cax=self.cb.figure.gca(), label=vol.label) # colorbar
        self.cb.figure.gca().yaxis.set_label_position('left')
        self.cb.figure.canvas.draw()
        
    def update_axis(self, i, value):
        """ update the i-th axis according to its slider position value """

        def draw_crosshair(ax, idx):
            """ Draw vertical and horizontal guide lines for a given axis index. """
            x_ax, y_ax = [ k for k in range(3) if k != idx ]
            ax.axvline(x=self.CT.i2mm(self.w[x_ax].value(),x_ax), ls=':', c='w', alpha=0.5)
            ax.axhline(y=self.CT.i2mm(self.w[y_ax].value(),y_ax), ls=':', c='w', alpha=0.5)


        self.ax[i].cla() # Clear axes
        
        self.set_axis_limit(i)
        
        ### 1. CT
        self.ax[i].imshow(self.CT.image(i, value), vmin=self.CT.cmin, vmax=self.CT.cmax, extent=self.CT.extent(i), cmap=self.CT.cmap, interpolation='none', origin='lower')

        ### 2. Dose / Dose-rate / Overlaps
        if self.plan.is_opened and (self.display_option() > 0):
            pos = self.Data.slice_index(self.CT, i, value)
            if 0 <= pos < self.Data.dim[i]:
                img = np.ma.masked_where(self.Data.image(i, pos) < self.Data.cmin, self.Data.image(i, pos)) # = self.Data.image(i, pos)
                self.ax[i].imshow(img, vmin=self.Data.cmin, vmax=self.Data.cmax, alpha=0.5, cmap=self.Data.cmap, extent=self.Data.extent(i), interpolation='none', origin='lower')

        ### 3. Profile tool
        if self.profile_show:
            if i == self.ax.index(self.profile_ax):
                self.profile_ax.plot([self.xp1, self.xp2], [self.yp1, self.yp2], 'C0o-')
                self.updateProfile(self.profile_ax)
        
        ### 4. Structures
        for roi in self.myROISet.selected:
            if i==2:
                for cs in roi.contours:
                    if self.CT.mm2i(cs.z, i) == value: self.ax[i].plot(cs.x, cs.y, c=roi.color, lw=1, marker=None if len(cs.x)>2 else '*')
            else:
                if roi.mask is not None:
                    img = np.take(roi.mask, value, axis=[2, 1, 0][i])
                    if True in img: self.ax[i].contour(img.astype(float), levels=[0.5], extent=self.CT.extent(i), colors=[roi.color], linewidths=1., origin='lower')
                
        if i==2: self.ax[i].text(0.05, 0.95, self.CT.patientName, va='top', transform=self.ax[0].transAxes, c='w') # Patient name
        
        self.x_lim[i], self.y_lim[i] = self.ax[i].get_xlim(), self.ax[i].get_ylim() # Get axis limit
        
        if self.CT.is_opened: self.ax[i].imshow(self.cosmo.img[i], extent=self.cosmo.extent(i), zorder=4) # Cosmo
            
        ### Calculation volume
        if self.CT.is_opened and self.cropAction.isChecked():
            color = 'dodgerblue'
            lp = dict(facecolor=None, edgecolor=color, fill=False, lw=1)
            hp = dict(mfc=color, mec=color, ls='none', alpha=1, ms=5.)
            self.RS[i] = mpl.widgets.RectangleSelector(self.ax[i], self.onRSselect, button=[1], props=lp, handle_props=hp, interactive=True, ignore_event_outside=True, grab_range=40)
            self.RS[i].extents = self.CT.bbox2RSExtentmm(self.bbox, i) # (xmin, xmax, ymin, ymax)
                    
        ### Crosshair
        if self.crosshairAction.isChecked():
            draw_crosshair(self.ax[i], i)
            for j in range(3):
                if j != i:
                    if len(self.ax[j].lines)>1: self.ax[j].lines[-1].remove(); self.ax[j].lines[-1].remove()
                    draw_crosshair(self.ax[j], j)
                    self.fig[j].figure.canvas.draw()
                
        self.fig[i].figure.canvas.draw() # Draw
            
    def onRSselect(self, eclick, erelease):
        if self.cropAction.isChecked():
            ax = self.ax.index(erelease.inaxes)
            self.bbox = self.CT.RSExtentmm2bbox(self.RS, ax)
            
            # from RP2TOPAS_ProBeam import PatientTranslation
            # print(f'self.bbox: {self.bbox}')
            # # print(f'PatientTranslation: {PatientTranslation([90.0,-360.0, 500.0], 0.0, self.bbox, self.CT)}') # AAMD2017
            # print(f'PatientTranslation: {PatientTranslation([-10.0,10.2, 15.2], 0.0, self.bbox, self.CT)}') # GBM80
            # vox_res = 1.0
            # XBins = int( (max(self.bbox[0]) - min(self.bbox[0]))*self.CT.spacing[0]/vox_res )
            # YBins = int( (max(self.bbox[1]) - min(self.bbox[1]))*self.CT.spacing[1]/vox_res )
            # ZBins = int( (max(self.bbox[2]) - min(self.bbox[2]))*self.CT.spacing[2]/vox_res )
            # print(f'XBins: {XBins}')
            # print(f'YBins: {YBins}')
            # print(f'ZBins: {ZBins}')
            
            self.update_all_axes()
        
    def SetCTWindowing(self, scale):
            
        self.CT.cmin = {'Brain':0, 'Lungs':-1000,'Abdomen':-150,'Bone':200, 'Auto':-1000}[scale]
        self.CT.cmax = {'Brain':80,'Lungs':1000, 'Abdomen':250, 'Bone':1000,'Auto':1000 }[scale]
        self.update_all_axes()

    def updateProfile(self, profile_ax):
        """ Display profile along a line """
        
        self.tabWidget.tabs.setCurrentIndex(0)

        with suppress(Exception):
            ax = self.ax.index(profile_ax)
            pos = self.w[ax].value()
            vol = self.Data if self.display_option() > 0 else self.CT
            ylabel = {0:'HU', 1:'Dose (Gy)', 2:'Dose Rate (Gy/s)', 3:'Overlap'}[self.display_option()]
            pos = vol.slice_index(self.CT, ax, pos) if self.display_option() > 0 else pos

            ### 1. Compute line coordinates
            xp1_vox, yp1_vox = vol.mm2voxel(self.xp1, self.yp1, ax)
            xp2_vox, yp2_vox = vol.mm2voxel(self.xp2, self.yp2, ax)            
            xp,yp = np.linspace(xp1_vox, xp2_vox, 1000), np.linspace(yp1_vox, yp2_vox, 1000)
            
            dist_ = np.linspace(0, np.sqrt( (self.yp2-self.yp1)**2 + (self.xp2-self.xp1)**2 ), 1000)
            prof_ = map_coordinates(vol.image(ax, pos), np.vstack((xp,yp)), order=1, cval=-1000. if self.display_option()==0 else 0.)
            self.profile = [dist_, prof_]
            
            ### 2. Display profile
            self.tabWidget.ClearTab(self.tabWidget.profile)
            self.tabWidget.profile.figure.gca().plot(dist_, prof_)
            self.tabWidget.profile.figure.gca().set_ylabel(ylabel)
            self.tabWidget.profile.draw()
    
    def roiIndices(self, roi, vol):
        """ return the indices of the voxels inside a given ROI """

        ### special case where dose spacing is large w.r.t CT spacing
        zs = np.unique([cs.z for cs in roi.contours])
        skip = max(1, int(vol.spacing[2] / self.CT.spacing[2]))
        contours = [cs for cs in roi.contours if cs.z in zs[::skip]] if skip > 1 else roi.contours
        
        return get_roi_binary_mask(roi, vol.dim, vol.origin, vol.spacing, contours)
        
    def showDVH(self, option='DVH'):
        """ Display the Dose Volume Histogram of each structure selected. """
        
        if not self.plan.is_opened: return
        
        def histogram(data): v = np.linspace(0,100,200); return [0, *np.percentile(data, v), data.max()], [100, *v[::-1], 0]
            
        if option == 'DRVH':
            if len(self.plan.selectedFields) > 1: self.ThrowError('Dose-rate cannot be displayed for multiple fields.'); return
            if not self.plan.selectedFields[0].DoseRate.is_opened: self.ThrowError('Open a dose-rate distribution first'); return
        
        match option:
            case 'DVH': canvas, tabidx = self.tabWidget.dvh, 1
            case 'DRVH': canvas, tabidx = self.tabWidget.drvh, 2
            case 'DDRVH': canvas, tabidx = self.tabWidget.ddrvh, 3
            case 'DR-DVH': canvas, tabidx = self.tabWidget.dvh, 1
            
        self.tabWidget.tabs.setCurrentIndex(tabidx)
        self.tabWidget.ClearTab(canvas)
        ax = canvas.figure.gca()
        
        vol = self.plan.selectedFields[0].Dose if self.args.fields else self.Data
        
        for c,roi in enumerate(self.myROISet.selected):
            
            if c > 15: continue
            
            mask = self.roiIndices(roi, vol)
            # self.Data.volume[mask] = 100. # debug
            
            if not True in mask: continue
            
            if len(roi.contours) > 1:
                
                match option:
                    case 'DVH' | 'DR-DVH':
                        dose = sum([f.Dose.volume for f in self.plan.selectedFields]) if self.args.fields else self.Data.volume
                        roi.dvh['bins'], roi.dvh['vals'] = histogram(dose[mask])
                        ax.plot(roi.dvh['bins'], roi.dvh['vals'], label=roi.name, color=roi.color)
                        
                        if option == 'DR-DVH':
                            drdose = sum(AForBBetween(f.Dose.volume, f.DoseRate.volume, self.drmin, self.drmax) for f in self.plan.selectedFields)
                            ax.plot(*histogram(drdose[mask]), color=roi.color, ls=':')

                    case 'DRVH':
                        doserate = self.plan.selectedFields[0].DoseRate.volume
                        roi.drvh['bins'], roi.drvh['vals'] = histogram(doserate[mask])
                        ax.plot(roi.drvh['bins'], roi.drvh['vals'], label=roi.name, color=roi.color)

                    case 'DDRVH':
                        if not all([f.DoseRate.is_opened for f in self.plan.selectedFields]): self.ThrowError('Open a dose-rate distribution first'); return
                        for i,f in enumerate(self.plan.selectedFields): ax.scatter(f.Dose.volume[mask], f.DoseRate.volume[mask], label=roi.name if i==0 else None, marker='.', color=roi.color)

        ax.legend(bbox_to_anchor=(1.02,1), loc='upper left', frameon=False)
        canvas.draw()
        
    def getOverlapVolume(self):

        self.displaymenu.actions()[3].trigger() # CT + overlaps
        vox_volume_cc = np.prod(self.Data.spacing)*1E-3
        for o in [1,2,3,4,5]: print( f'overlap >= {o}: {sum(np.where(self.Data.volume.ravel()>=o,1,0))*vox_volume_cc:.1f} cc' )

    def ChangeColorbarRange(self):

        if not self.plan.is_opened: return
        vol = self.Data if self.display_option() > 0 else self.CT
        unit = {0:'HU', 1:'Gy', 2:'Gy/s', 3:''}[self.display_option()]
        ColorbarRangeDialog(self, vol, unit).exec()
        
    def About(self):
        """ display the About window """
        
        InfoWindow = RTGeneralDialog(self,'About Dose-Rate Viewer')
        
        form = QFormLayout()
        form.setLabelAlignment(QtCore.Qt.AlignRight)
        for i, j in zip(['Version:','Author:'],['09.07.2026', 'pierre.lansonneur@varian.com']):
            tmp = QLabel(i, alignment = QtCore.Qt.AlignTop)
            tmp.setStyleSheet("color: #8F8F8F")
            form.addRow(tmp, QLabel(j))
        gbox = QGroupBox('Dose-Rate Viewer')
        gbox.setLayout(form)
        
        layout = QGridLayout()
        layout.addWidget(gbox,0,1)
        InfoWindow.setLayout(layout)     
        InfoWindow.exec()
        
    def Infos(self):
        """ display the Infos window """
            
        def array_format(a, r=0): return '\t'.join(f'{i}' for i in np.round(a,r))
        
        def make_label(text, color="#8F8F8F"): l = QLabel(text); l.setStyleSheet(f"color: {color}"); return l
        
        def make_form(): f = QFormLayout(); f.setVerticalSpacing(0); f.setLabelAlignment(QtCore.Qt.AlignRight); return f
        
        InfoWindow = RTGeneralDialog(self,'Infos')
        layout = QVBoxLayout()
        
        for txt, vol in zip(['CT SERIE','DOSE'],[self.CT, self.Data]):
            
            if not vol.is_opened: continue

            form = make_form()
            form.addRow(make_label('Patient ID:'), QLabel(vol.patientName))
            form.addRow(make_label('Origin (mm):'), QLabel(array_format(vol.origin,1)))
            form.addRow(make_label('Size (voxel):'), QLabel(array_format(vol.dim)))
            form.addRow(make_label('Resolution (mm):'), QLabel(array_format(vol.spacing,3)))
            gbox = QGroupBox(txt)
            gbox.setLayout(form)
            layout.addWidget(gbox)
            
        if self.myROISet.is_opened:
            
            form = make_form()
            for roi in self.myROISet.ROIs:
                rgb_color = ','.join(f'{int(255*i)}' for i in roi.color)
                form.addRow(make_label('\u2B24', f'rgb({rgb_color})'), QLabel(roi.name))
                
            sbox = QWidget()
            sbox.setLayout(form)
            sbox.setStyleSheet("border: 0px")
            scroll = QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
            scroll.setWidget(sbox)
            scroll.setStyleSheet("QScrollArea {border: 0px}")
            gbox = QGroupBox(f"STRUCTURES ({len(self.myROISet.ROIs)})")
            vlayout = QVBoxLayout(gbox)
            vlayout.addWidget(scroll)
            layout.addWidget(gbox)
        
        InfoWindow.setLayout(layout)      
        InfoWindow.exec()

    def ThrowError(self, message): RTErrorDialog(message).exec()
    
class RTGeneralDialog(QDialog):

    def __init__(self, parent=None, title=' '):
        
        super().__init__()
        self.parent = parent
        self.setWindowTitle(title)
        self.setWindowIcon(QtGui.QIcon(CURRENT_PATH + "\\Icon.ico"))
    
class RTDialog(RTGeneralDialog):

    def __init__(self, parent=None, title=' '):
        
        super().__init__(parent, title)
        
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        self.formlayout = QFormLayout()
        self.formGroupBox = QGroupBox(f'')
        self.formGroupBox.setLayout(self.formlayout)
        self.layout = QVBoxLayout()
        self.layout.addWidget(self.formGroupBox)
        self.setLayout(self.layout)
            
    def accept(self): self.done(QDialog.Accepted)
    
    def reject(self): self.done(QDialog.Rejected)
        
class RTErrorDialog(RTDialog):

    def __init__(self, message='', parent=None):
        
        super().__init__(parent, 'Error')
        
        icon = QLabel()
        icon.setPixmap(openPixmap(os.path.join(ICONS_PATH, 'error.png')))
        self.formlayout.addRow(icon, QLabel(message))
        
class ColorbarRangeDialog(RTDialog):

    def __init__(self, parent, vol, unit):
        
        super().__init__(parent)
        self.resize(300, 30)
        self.unit = unit
        self.vol = vol
        self.formGroupBox.setTitle(f'Select colorbar range ({unit})')
        self.cmin = QLineEdit(f'{self.vol.cmin:.1f}')
        self.cmax = QLineEdit(f'{self.vol.cmax:.1f}')
        self.formlayout.addRow(QLabel("from "), self.cmin)
        self.formlayout.addRow(QLabel("to "), self.cmax)
        self.layout.addWidget(self.buttonBox)

    def accept(self):
        
        cmin, cmax = float(self.cmin.text()), float(self.cmax.text())
            
        if cmin < cmax:
            self.vol.cmax = cmax
            self.vol.cmin = cmin if self.unit == 'HU' else max(cmin, 0.001)
            self.done(QDialog.Accepted)
            self.parent.update_all_axes()
    
class DRDialog(RTDialog):

    def __init__(self, parent):
        
        super().__init__(parent)
        self.formGroupBox.setTitle(f'Select dose-rate (Gy/s)')
        self.drmin = QLineEdit(f'{self.parent.drmin:.1f}' if self.parent.drmin else '')
        self.drmax = QLineEdit(f'{self.parent.drmax:.1f}' if self.parent.drmax else '')
        self.formlayout.addRow(QLabel("from "), self.drmin)
        self.formlayout.addRow(QLabel("to "), self.drmax)
        self.layout.addWidget(self.buttonBox)

    def accept(self):
        
        drmin, drmax = None, None
        with suppress(ValueError): drmin = max(0., float(self.drmin.text()))
        with suppress(ValueError): drmax = max(0., float(self.drmax.text()))
            
        if drmin and drmax:
            if drmin >= drmax: self.parent.ThrowError(f'min dose rate must be smaller that max dose rate'); return
            
        self.parent.tabWidget.ClearTab(self.parent.tabWidget.dvh)
        self.parent.drmin = drmin
        self.parent.drmax = drmax
        self.parent.Field_selection()
        self.done(QDialog.Accepted)
    
class MyTabWidget(QWidget):
    
    def __init__(self):

        super().__init__()
        self.layout = QVBoxLayout(self)
        self.tabs = QTabWidget()
        self.profile = FigureCanvas(Figure())
        self.dvh = FigureCanvas(Figure())
        self.drvh = FigureCanvas(Figure())
        self.ddrvh = FigureCanvas(Figure())
        self.tabs.addTab(self.profile," Profile ")
        self.tabs.addTab(self.dvh," DVH ")
        self.tabs.addTab(self.drvh," DRVH ")
        self.tabs.addTab(self.ddrvh," 2D - DDRVH ")
        self.layout.addWidget(self.tabs)

    def ClearAllTabs(self):

        for c in [self.profile,self.dvh,self.drvh,self.ddrvh]: self.ClearTab(c)

    def ClearTab(self, canvas):

        canvas.figure.gca().clear()
        match canvas:
            case self.profile: xlabel, ylabel = 'Depth (mm)', ''
            case self.dvh: xlabel, ylabel = 'Dose (Gy)', 'Volume (%)'
            case self.drvh: xlabel, ylabel = 'Dose Rate (Gy/s)', 'Volume (%)'
            case self.ddrvh: xlabel, ylabel = 'Dose (Gy)', 'Dose Rate (Gy/s)'
        canvas.figure.gca().set_xlabel(xlabel)
        canvas.figure.gca().set_ylabel(ylabel)
        canvas.figure.set_tight_layout(True)
        canvas.draw()

if __name__ == '__main__':
        
    app = QApplication(sys.argv)    
    app.setFont(QtGui.QFont('Segoe UI', 9))
    app.setStyleSheet(qdarkstyle.load_stylesheet(qt_api='pyside6')) # qdarkstyle dark theme
    P.style.use('dark_background') # matplotlib dark theme

    parser = argparse.ArgumentParser(description=' ', prefix_chars='-+/')
    parser.add_argument('-fields', nargs='+') # field IDs
    args = parser.parse_args()
    args.fields = natsorted(args.fields, alg=ns.IGNORECASE) if args.fields else None
    window = RTMainWindow(args)
    window.setWindowIcon(QtGui.QIcon(CURRENT_PATH + "\\Icon.ico"))
    window.setStyleSheet("""QLabel {background-color : #252525}""")
        
    # patientID = 'Lung004'
    # window.OpenDicomSerie(rf'D:\Patients\{patientID}\CT')
    # window.myROISet.OpenDicom(rf'D:\Patients\{patientID}\RS.{patientID}.dcm')
    # # window.OpenDose(rf'D:\Patients\{patientID}\MBRT\IM_1.0mm_200MeV_ctc4.0mm\RD.Optimized.dcm')
    # window.plan.ImportDoseAndDoseRate(rf'D:\Patients\{patientID}\coarseIM_2mm_Square\WO_3Fx')
    # # window.plan.OpenDicom(rf'D:\Patients\{patientID}\RN.Test.dcm')
    
    # window.OpenDicomSerie(rf'C:\dev\HN_TG244_DMAT_CONV_2A_ADS_ioe3-debug\2026-07-08_10-35-08\dicom')
    # window.myROISet.OpenDicom(rf'C:\dev\HN_TG244_DMAT_CONV_2A_ADS_ioe3-debug\2026-07-08_10-35-08\dicom\RS.1.2.246.352.800.5258531985563228124.102485075693992585.dcm')
    # window.OpenDose(r'C:\dev\HN_TG244_DMAT_CONV_2A_ADS.RD.dcm')
    
    window.showMaximized() # full screen
    exit_code = app.exec()
    sys.exit(exit_code)