A visualization tool to:

- Display CT DICOM files,
- Display RT-Dose files,
- Display RT-Struct files,
- Display dose and dose-rate distributions,
- Convert protontherapy RT-PLAN to TOPAS MC files.

<img src="https://github.com/PierreLansonneur/Viewer/screenshot.png" width="400" />

### Getting Started
* Install python 3.10 (e.g. from: https://www.python.org/downloads/windows/).
* Clone the repo:
~~~ 
git clone ssh://tfs-app.vms.ad.varian.com:22/VMS/1-VMS/_git/TPMS.Proton.FlashTools
~~~
* Create the virtual environment and install necessary libraries:
~~~
python -m venv venv
.\venv\Scripts\activate
pip install -r requirements.txt
~~~
* Launch the viewer with:
~~~
python Viewer.py
~~~

### Import CT series
Open a dicom series (CTRL + O) by navigating through the folder containing the DICOM files. 
Once opened, three different views (axial, sagittal and coronal plane) are displayed on the axes. 
You can navigate through the different slices using either the mouse scroll or the three sliders on the viewer side.
To zoom/unzoom on the images, press right-click and move the mouse on the corresponding axis.
1D profiles along a line can be visualized by pressing the 'ctrl' key and click on the corresponding axis.

### Import RT-dose files
Once a CT is opened, RT-Dose files can be imported via File >> Import Dose.
The colorbar range can be set via the Options menu.

### Import RT-struct files
Once a CT is opened, DICOM RT-Struct files can be imported via Structures >> Import Structures.
Once a RT-Dose file and a DICOM RT-Struct file are opened, you can compute the DVHs for the ROIs selected in the Structures menu.

### Dose-Rate visualization tools
To unlock dose-rate related features, launch the viewer by specifying the field IDs:
~~~
python Viewer.py -fields "Field 1" "Field 2" "Field 3" "Field 4" "Field 5"
~~~
The dose and dose-rate distributions in DICOM format can then be imported via File >> Import Dose and Dose-Rate.
The dose (respectively dose-rate) distribution needs to be in DICOM format and named 'RD.Dose_ID.dcm' (resp. 'RD.DoseRate_ID.dcm'), where ID is the field ID.
For indivudal fields, one can then display the dose-rate distributions, or compute the Dose-Rate Volume Histograms (DRVH) for the ROIs selected in the Structures menu.
A minimum dose-rate can be specified in the Options menu to see only the dose delivered above this threshold.
If a minimum dose-rate is specified, one can plot the Dose-Rate Dose Volume Histograms (DR-DVH) together with the regular DVH via the Tools menu.
The 2D histogram of dose and dose-rate for all fields combined can moreover be shown via the Tools menu.
 
### Convert RT-Plan to TOPAS MC files
If the field IDs are specified, protontherapy DICOM RT-Plan files can be imported via File >> Open RP file. The interface only supports files intended for Pencil Beam Scanning. 
For each field, a TOPAS MC simulation file can be exported based on the informations in the file, such as the range shifter, gantry and table angles, number of spots, etc…