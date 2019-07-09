# Viewer
Some useful radiotherapy tools

![](https://github.com/PierreLansonneur/Viewer/blob/master/output/screenshot.png =200x)

#### launch the viewer 
`ipython Main.py &`

#### Features

- Visualize CT DICOM files
- Import RT-Dose files
- Import RT-Struct files
- Import RT-PLAN files
- Create TOPAS MC files for PBS and DS

#### Dependencies

- Python 2.7
- Tkinter `sudo apt-get install python-tk`
- pyDICOM `pip install pydicom`
- SimpleITK `pip install SimpleITK`
- vtk `pip install vtk`
- pyFPDF `pip install fpdf`
