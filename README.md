# Viewer

A visualization tool distributed for Linux that can be used to:

- Visualize CT DICOM files,
- Import RT-Dose files,
- Import RT-Struct files,
- Import RT-PLAN files,
- Create TOPAS MC files for PBS and DS.

<img src="https://github.com/PierreLansonneur/Viewer/blob/master/output/capture.jpg" width="400" />

#### launch the viewer: `python Main.py`

#### Dependencies
- Python 2.7
- Tkinter `sudo apt-get install python-tk`
- pyDICOM `pip install pydicom`
- SimpleITK `pip install SimpleITK`
- vtk `pip install vtk`
- pyFPDF `pip install fpdf`

### Image navigation
The file menu allows you to open Dicom files (.dcm), dicom series or a metaimages (.mhd). Once opened, three different views (axial, sagittal and coronal plane) are displayed on the axes. You can navigate through the different slices using either the mouse scroll or the three sliders on the viewer side.
To zoom/unzoom on the images, press right-click and move the mouse on the corresponding axis. To move the images, you may press left-click and move the mouse on the corresponding axis. The image orientation can be changed either by double clicking on the axis or via the Edit >> Rotate menu.
The colorscale limits can be set either to ‘Hounsfield’ (from -1000 HU to 2000 HU for CT images) or to ‘Auto’ (for any other type) via the Scale menu. 

### Import RT-dose files
Once a CT is opened, RT-Dose files can be imported via the File >> Import dosimetry menu. The isododose contours appear then on top of the image. To visualize the dose matrix as a color gradient, set instead the display option to ‘gradient’ in the Tool >> Isodoses menu.

### Import RT-struct files
Once a CT is opened, DICOM RT-Struct files can be imported via the Structures >> Import structures menu. The names and colors of the structures imported can be visualized in the More >> File infos window. Once a RT-Dose file and a DICOM RT-Struct file are opened, you may access to the following tools: 
-	compute DVH: Calculates the cumulative DVH for every structures
-	crop Dosimetry: Crop the dose matrix to the extent of the external contour.
-	normalize dosi to PTV dose: Normalize the dose matrix such that the average dose is equal to 100 in a structure whose the name contains ‘PTV’.
 
### Import RT-Plan files
Once a CT is opened, protontherapy DICOM RT-Plan files can be imported via the File >> Open RP file menu. The interface supports files intended for:
-	Passive Scattering: The collimator contour and compensator thickness are shown for each field (BEV).
-	Pencil Beam Scanning: The spots position and energy are shown for each field (BEV).
 
For each fields, several informations are retrieved such as the range shifter, gantry and table angles, number of spots, etc… For PBS files, you can moreover create a TOPAS MC simulation file for each field. The file is made to be interfaced with the ICPO gantry nozzle geometry (see https://github.com/PierreLansonneur/TOPAS_Gantry for public version of the code). The dose matrix has a default resolution of 0.5 x 0.5 x 0.5 mm3 and the size can be adjusted with the slider

### Miscellaneous tools
- Iso-surface: Users can visualize files with 3D rendering using the Tools >> iso surface window vtk tool. If only the CT is shown, the iso-suface corresponding to the value selected is displayed. If the isodoses are shown, the iso-surface value is interpreted as the percent of the maximum of the dose matrix. If the structures are shown, the external contour is displayed. 

-	Gamma index analysis: Dose matrices can be compared with the gamma index tool in Tools >> Gamma. To compare two matrices, load first one matrix with ‘Open File’ and the other via ‘Import dosimetry’. Default criterion is 3% / 3mm.
 
-	Profiles: 1D profiles along a line can be visualized using the ‘Tool >> Profile’ widget. Select the line for which you want to show the distribution by pressing the ‘ctrl’ key and click on the corresponding axis. 
