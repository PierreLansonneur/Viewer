import numpy as np

p_Esigma = np.poly1d([-0.00000023, 0.00008559, -0.01205498, 1.23066909]) # energy spread
p_size = np.poly1d([-0.00000143, 0.00086454, -0.19620932, 21.66378475]) # beam sigma (mm)
p_div = np.poly1d([-0.000000002, 0.000001141, -0.00021799, 0.0156]) # Beam divergence
p_SM1 = np.poly1d([0.0005266, -0.2767, 60.29]) # SM1 field (T)
p_SM2 = np.poly1d([-0.0004444, 0.2328, -50.43]) # SM2 field (T)
p_MUGp = np.poly1d([-2.92661e-11,  2.31840e-08, -6.63090e-06,  9.298374e-04]) # MU/Gp coefficient

def PatientTranslation(iso, PatientSupportAngle, bbox, CT):
    """
    return the patient translation coordinates
    iso: isocenter coordinates (mm)
    PatientSupportAngle: table angle (deg)
    bbox: bounding box of the dose matrix (voxel)
    """

    # 1. isocenter translation
    tx = CT.origin[0] + 0.5*(bbox['x_min'] + bbox['x_max'])*CT.spacing[0] - iso[0]
    ty = CT.origin[1] + 0.5*(bbox['y_min'] + bbox['y_max'])*CT.spacing[1] - iso[1]
    tz = CT.origin[2] + 0.5*(bbox['z_min'] + bbox['z_max'])*CT.spacing[2] - iso[2]

    # 2. Table Rotation
    cos, sin = np.cos(np.deg2rad(PatientSupportAngle)), np.sin(np.deg2rad(PatientSupportAngle))
    
    return np.array([[cos,0,sin],[0,1,0],[-sin,0,cos]]) @ [tx,ty,tz] # rotation around Y axis

def CreateTopasFilePBS(field, CT, bbox, filedir, vox_res=1.):
    
    # print('\nfield: {vars(field)}\nCT: {vars(CT)}\nbbox: {bbox}')
    
    retraction = field.SnoutPosition + 1.8 if field.SnoutID=='400' else field.SnoutPosition - 6.5
    N, x, y, E, MU = field.spots.N, np.array(field.spots.x), np.array(field.spots.y), np.array(field.spots.E), np.array(field.spots.MU)
    TransX, TransY, TransZ = PatientTranslation(field.IsocenterPosition, field.PatientSupportAngle, bbox, CT)
    
    GRAPHIC_SUPPORT = True
    # GRAPHIC_SUPPORT = False
    
    filename = filedir + f'/Main_{field.parent.PatientID.replace(" ", "_")}_{field.ID.replace(" ", "_")}.txt'
    with open(filename, 'w') as f:
    
        f.write('################################\n')
        f.write('# TOPAS MC v3.2 simulation file\n')
        f.write(f'# Patient    : {field.parent.PatientID}\n')
        f.write(f'# Field ID   : {field.ID}\n')
        f.write(f'# N spots    : {field.spots.N}\n')
        f.write(f'# Gantry     : {field.GantryAngle:.1f} deg\n')
        f.write(f'# Table      : {field.PatientSupportAngle:.1f} deg\n')
        f.write(f'# Isocenter  : {field.IsocenterPosition[0]:.2f}, {field.IsocenterPosition[1]:.2f}, {field.IsocenterPosition[2]:.2f} mm\n')
        f.write(f'# RSWET      : {field.RSWET:.1f} mm\n')
        f.write(f'# Resolution : {vox_res:.2f} mm\n')
        f.write(f'# Machine    : Institut Curie IBA nozzle\n')
        f.write('################################\n\n')

        f.write('### GENERAL COMMANDS\n')
        f.write('b:Ge/CheckForOverlaps = "1"\n')
        f.write('b:Ge/CheckForUnusedComponents = "0"\n')
        f.write('i:Ts/NumberOfThreads = 0\n')
        f.write('b:Ts/ShowCPUTime = "1"\n')
        f.write('i:Ts/MaxInterruptedHistories = 1000000000\n')
        f.write('i:Ts/ShowHistoryCountAtInterval = 0 # reduce size of slurm file\n\n')

        f.write('### WORLD VOLUME\n')
        f.write('d:Ge/World/HLX = 300 cm\n')
        f.write('d:Ge/World/HLY = 300 cm\n')
        f.write('d:Ge/World/HLZ = 300 cm\n\n')

        f.write('#### GEOMETRY\n')
        f.write('includeFile = Geometry/Nozzle.txt\n')
        f.write(f'includeFile = Geometry/Snout{0.1*float(field.SnoutID):.0f}.txt\n')
        # f.write(f'includeFile = Geometry/Snout{field.SnoutID}.txt\n')

        # RS 65 mm, 30 mm or nothing
        # if field.RSWET > 0:
            # f.write('b:Ge/RS/Include = "1"\n')
            # f.write(f'd:Ge/RS/HL = {0.5*field.RSWET:.1f} mm\n')
            # f.write('d:Ge/RS/TransZ = Ge/Snout/TransIso + Ge/RS/HL cm\n')
        if (field.RSID == 'RS65' or field.RSID == 'RS30'):
            if(field.RSID == 'RS65'):       RST = 65
            elif(field.RSID == 'RS30'):     RST = 30
            f.write('b:Ge/RS/Include = "1"\n')
            f.write(f'd:Ge/RS/HL = {0.5*RST:.1f} mm\n')
            f.write('d:Ge/RS/TransZ = Ge/Snout/TransIso + Ge/RS/HL cm\n')
        else:
            f.write('b:Ge/RS/Include = "0"\n')

        f.write(f'd:Ge/Snout/TransZ = -{field.SnoutPosition:.1f} cm + Ge/Snout/AntiTransIso\n') 
        f.write(f'd:Ge/SnoutGroup/TransZ = {10.-retraction:.1f} cm\n')
        f.write(f'd:Ge/FM/HLZ = {26.225+0.5*(10.-retraction):.3f} cm\n')
        f.write(f'd:Ge/FM/TransZ = {-87.026+0.5*(10.-retraction):.3f} cm\n\n')

        f.write('### SuperGroup for Gantry rotation\n')
        f.write('s:Ge/SuperGroup/Parent = "World"\n')
        f.write('s:Ge/SuperGroup/Type = "Group"\n')
        f.write('d:Ge/SuperGroup/TransX = 0 m\n')
        f.write('d:Ge/SuperGroup/TransY = 0 m\n')
        f.write('d:Ge/SuperGroup/TransZ = 0 m\n')
        f.write('d:Ge/SuperGroup/RotX = 90 deg\n')	
        f.write(f'd:Ge/SuperGroup/RotY = {field.GantryAngle:.1f} deg # Gantry angle\n')
        f.write('d:Ge/SuperGroup/RotZ = 0 deg\n\n')

        f.write('### PATIENT\n')
        f.write('s:Ge/Patient/ImagingToMaterialConverter = "Schneider"\n')
        f.write('includeFile = HUtoMAT.txt # WARNING check the scanner model!\n\n') 

        f.write('### Patient group\n')
        f.write('s:Ge/Patient/Message = "Constructing Patient"\n')
        f.write('s:Ge/Patient/Type = "Group"\n')
        f.write('s:Ge/Patient/Parent = "World"\n')
        f.write(f'd:Ge/Patient/TransX = {TransX:.1f} mm\n')
        f.write(f'd:Ge/Patient/TransY = {TransY:.1f} mm\n')
        f.write(f'd:Ge/Patient/TransZ = {TransZ:.1f} mm\n') 
        f.write('d:Ge/Patient/RotX = 0 deg\n')
        f.write(f'd:Ge/Patient/RotY = {field.PatientSupportAngle:.1f} deg # table angle\n')
        f.write('d:Ge/Patient/RotZ = 0 deg\n\n')
        
        if GRAPHIC_SUPPORT:
            f.write('### Patient STL file (for visualization only)\n')
            f.write('s:Ge/CT_visu/Type = "TsCAD"\n')
            f.write('s:Ge/CT_visu/Parent = "Patient"\n')
            f.write('b:Ge/CT_visu/IsParallel = "True"\n')
            # f.write('s:Ge/CT_visu/Material = "Vacuum"\n')
            f.write(f'd:Ge/CT_visu/TransX = {-0.5*(bbox["x_min"] + bbox["x_max"])*CT.spacing[0] - CT.origin[0]:.1f} mm\n')
            f.write(f'd:Ge/CT_visu/TransY = {-0.5*(bbox["y_min"] + bbox["y_max"])*CT.spacing[1] - CT.origin[1]:.1f} mm\n')
            f.write(f'd:Ge/CT_visu/TransZ = {-0.5*(bbox["z_min"] + bbox["z_max"])*CT.spacing[2] - CT.origin[2]:.1f} mm\n')
            f.write('s:Ge/CT_visu/DrawingStyle = "Wireframe"\n')
            f.write(f's:Ge/CT_visu/InputFile = "Geometry/{field.parent.PatientID}"\n')
            f.write('s:Ge/CT_visu/FileFormat = "stl"\n')
            f.write('d:Ge/CT_visu/Units = 1.0 mm\n\n')

        f.write('### Patient modelisation\n')
        f.write('s:Ge/Patient/CT/Parent = "Patient"\n')
        f.write('s:Ge/Patient/CT/Type = "TsDicomPatient"\n')
        # f.write(f'b:Ge/Patient/CT/IsParallel = "{keepCTDim}"\n')
        f.write('s:Ge/Patient/CT/Material = "G4_WATER"\n')
        # f.write(f's:Ge/Patient/CT/DicomDirectory = "./CT_{field.parent.PatientID}"\n')
        f.write(f's:Ge/Patient/CT/DicomDirectory = "./CT"\n')
        f.write('sv:Ge/Patient/CT/DicomModalityTags = 1 "CT"\n')
        f.write('d:Ge/Patient/CT/TransX = 0 mm\n')
        f.write('d:Ge/Patient/CT/TransY = 0 mm\n')
        f.write('d:Ge/Patient/CT/TransZ = 0 mm\n')
        f.write('d:Ge/Patient/CT/RotX = 0 deg\n')
        f.write('d:Ge/Patient/CT/RotY = 0 deg\n')
        f.write('d:Ge/Patient/CT/RotZ = 0 deg\n\n')
        # f.write('#sv:Ge/Patient/CT/ColorByRTStructNames = 2 "CONTOUR_EXTERNE" "PTV_1"\n')
        # f.write('#sv:Ge/Patient/CT/ColorByRTStructColors = 2 "green" "red"\n')
        # f.write('#b:Ge/Patient/CT/IgnoreInconsistentFrameOfReferenceUID = "True"\n\n')

        f.write('### Crop CT dimensions\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsXMin = {bbox["x_min"]+1}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsXMax = {bbox["x_max"]}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsYMin = {bbox["y_min"]+1}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsYMax = {bbox["y_max"]}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsZMin = {bbox["z_min"]+1}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsZMax = {bbox["z_max"]}\n\n')

        f.write('### Dose computation\n')
        f.write('s:Sc/PatientScorer/Component = "Patient/CT"\n')
        f.write('s:Sc/PatientScorer/Quantity = "DoseToWater"\n')
        f.write(f's:Sc/PatientScorer/OutputFile = "./output/Dose_{field.parent.PatientID}_{field.ID}"\n')
        f.write('s:Sc/PatientScorer/OutputType = "DICOM"\n')
        f.write('b:Sc/PatientScorer/OutputToConsole = "0"\n')
        f.write('s:Sc/PatientScorer/IfOutputFileAlreadyExists = "Increment"\n')
        f.write(f'i:Sc/PatientScorer/XBins = {int((bbox["x_max"] - bbox["x_min"])*CT.spacing[0]/vox_res):.0f}\n')
        f.write(f'i:Sc/PatientScorer/YBins = {int((bbox["y_max"] - bbox["y_min"])*CT.spacing[1]/vox_res):.0f}\n')
        f.write(f'i:Sc/PatientScorer/ZBins = {int((bbox["z_max"] - bbox["z_min"])*CT.spacing[2]/vox_res):.0f}\n\n')
        # f.write('#sv:Sc/PatientScorer/OnlyIncludeIfInRTStructure = 1 "CONTOUR_EXTERNE"\n')
        # f.write('#b:Sc/PatientScorer/OutputAfterRun = True\n\n')

        f.write('### PHYSICS (TOPAS default settings)\n')
        f.write('s:Ph/ListName = "Default"\n')
        f.write('s:Ph/Default/Type = "Geant4_Modular"\n')
        f.write('sv:Ph/Default/Modules = 7 "g4em-standard_opt3" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping" "g4radioactivedecay"\n')
        f.write('d:Ph/Default/EMRangeMin = 100 eV\n')
        f.write('d:Ph/Default/EMRangeMax = 250 MeV\n')
        f.write('i:Ph/Default/EMBinsPerDecade = 100\n')
        f.write('d:Ph/Default/CutForAllParticles = 0.05 mm\n')
        f.write('d:Ge/SM1/Dipole/MaxStepSize = 1 mm\n')
        f.write('d:Ge/SM2/Dipole/MaxStepSize = 1 mm\n\n')

        # # if (0.1*float(field.SnoutID)==10):	LMGW = ["Snout/BrassCone2","Snout/BrassTube2","Snout/BrassCone3"]
        # # if (0.1*float(field.SnoutID)==18):	LMGW = ["Snout/BrassCone2","Snout/BrassCone3"]
        # # if (0.1*float(field.SnoutID)==25):	LMGW = ["Snout/AlBrassTube2", "Snout/BrassCone2", "Snout/BrassCone3", "Snout/BrassCone4"]
        # # if (0.1*float(field.SnoutID)==40):	LMGW = ["Snout/SquarePart2", "Snout/SquarePart3"]
        if   field.SnoutID=='100':	LMGW = ["Snout/BrassCone2", "Snout/BrassTube2", "Snout/BrassCone3"]
        elif field.SnoutID=='180':	LMGW = ["Snout/BrassCone2", "Snout/BrassCone3"]
        elif field.SnoutID=='250':	LMGW = ["Snout/AlBrassTube2", "Snout/BrassCone2", "Snout/BrassCone3", "Snout/BrassCone4"]
        elif field.SnoutID=='400':	LMGW = ["Snout/SquarePart2", "Snout/SquarePart3"]
        # if field.RSWET > 0: LMGW += ["RS"]
        if (field.RSID == 'RS65' or field.RSID == 'RS30'):  LMGW += ["RS"]
        if GRAPHIC_SUPPORT: LMGW += ["CT_visu"]
        f.write(f'sv:Ph/Default/LayeredMassGeometryWorlds = {len(LMGW)} ' + ' '.join(f'"{i}"' for i in LMGW) +'\n\n')

        # f.write('### SPEED-UP the simulation\n')
        # f.write('#b:Vr/UseVarianceReduction = "True"\n')
        # f.write('#b:Vr/KillOtherParticles/Active= "True"\n')
        # f.write('#sv:Vr/KillOtherParticles/HaveNoEffectInComponentsNamed = 1 "Patient/CT"\n')
        # f.write('#sv:Vr/KillOtherParticles/OnlyTrackParticlesNamed = 1 "proton"\n\n')

        f.write('### Time Features Parameters\n')
        f.write(f'i:Tf/NumberOfSequentialTimes = {N}\n') 
        f.write('d:Tf/TimelineStart = 0 s\n')
        f.write(f'd:Tf/TimelineEnd = {N+1} s\n\n')

        f.write('### Spot Energy\n')
        f.write('s:Tf/Energy/Function = "Step"\n')
        f.write(f'dv:Tf/Energy/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/Energy/Values = {N} ' + ' '.join(f'{i:.2f}' for i in E) + ' MeV\n\n')

        f.write('##### Energy Spread\n')
        f.write('s:Tf/EnergySpread/Function = "Step"\n')
        f.write(f'dv:Tf/EnergySpread/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'uv:Tf/EnergySpread/Values = {N} ' + ' '.join(f'{i:.2f}' for i in p_Esigma(E)) +'\n\n')

        f.write('##### Spot Sigmas\n')
        f.write('s:Tf/SpotSize/Function = "Step"\n')
        f.write(f'dv:Tf/SpotSize/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/SpotSize/Values = {N} ' + ' '.join(f'{i:.2f}' for i in p_size(E)) + ' mm\n\n')

        f.write('##### Spot Divergences\n')
        f.write('s:Tf/BeamDivergence/Function = "Step"\n')
        f.write(f'dv:Tf/BeamDivergence/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'uv:Tf/BeamDivergence/Values = {N} ' + ' '.join(f'{i:.4f}' for i in p_div(E)) +'\n\n')

        f.write('##### SM1 scanning magnetic field\n')
        f.write('s:Tf/SM1Field/Function = "Step"\n')
        f.write(f'dv:Tf/SM1Field/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/SM1Field/Values = {N} ' + ' '.join(f'{i:.4f}' for i in 0.1*y/p_SM1(E)) + ' tesla\n')
        f.write('d:Ge/SM1/Dipole/MagneticFieldStrength = Tf/SM1Field/Value tesla\n')
        f.write(f'u:Ge/SM1/Dipole/MagneticFieldDirectionX = {np.cos(np.deg2rad(field.GantryAngle)):.4f}\n')
        f.write(f'u:Ge/SM1/Dipole/MagneticFieldDirectionY = {-np.sin(np.deg2rad(field.GantryAngle)):.4f}\n')
        f.write(f'u:Ge/SM1/Dipole/MagneticFieldDirectionZ = 0.0\n\n')

        f.write('##### SM2 scanning magnetic field\n')
        f.write('s:Tf/SM2Field/Function = "Step"\n')
        f.write(f'dv:Tf/SM2Field/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/SM2Field/Values = {N} ' + ' '.join(f'{i:.4f}' for i in 0.1*x/p_SM2(E)) + ' tesla\n')
        f.write('d:Ge/SM2/Dipole/MagneticFieldStrength = Tf/SM2Field/Value tesla\n')
        f.write('u:Ge/SM2/Dipole/MagneticFieldDirectionX = 0.0\n')
        f.write('u:Ge/SM2/Dipole/MagneticFieldDirectionY = 0.0\n')
        f.write('u:Ge/SM2/Dipole/MagneticFieldDirectionZ = 1.0\n\n')

        f.write('##### Number of Histories\n')
        f.write('s:Tf/Particles/Function = "Step"\n')
        f.write(f'dv:Tf/Particles/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'iv:Tf/Particles/Values = {N} ' + ' '.join(f'{int(i):.0f}' for i in 0.1*MU/p_MUGp(E)) +'\n\n')

        f.write('### SOURCE\n')
        f.write('s:So/Default/Type = "emittance"\n')
        f.write('s:So/Default/Component = "VW"\n')
        f.write('s:So/Default/BeamParticle = "proton"\n')
        f.write('d:So/Default/BeamEnergy = Tf/Energy/Value MeV\n')
        f.write('u:So/Default/BeamEnergySpread = Tf/EnergySpread/Value\n')
        f.write('s:So/Default/Distribution = "BiGaussian"\n')
        f.write('d:So/Default/SigmaX = Tf/SpotSize/Value mm\n')
        f.write('u:So/Default/SigmaXprime = Tf/BeamDivergence/Value\n')
        f.write('u:So/Default/CorrelationX = -0.882\n')
        f.write('d:So/Default/SigmaY = Tf/SpotSize/Value mm\n')
        f.write('u:So/Default/SigmaYprime = Tf/BeamDivergence/Value\n')
        f.write('u:So/Default/CorrelationY = -0.882\n')
        f.write('i:So/Default/NumberOfHistoriesInRun = 1 * Tf/Particles/Value\n\n')

        f.write('### GRAPHICS\n')
        f.write(f'b:Gr/Enable = "{GRAPHIC_SUPPORT}"\n')
        f.write('s:Gr/ViewA/Type = "VRML"\n')
        f.write('b:Ge/World/Invisible = "1"\n')
        f.write('b:Gr/ViewA/IncludeAxes = "1"\n')
        f.write('d:Gr/ViewA/AxesSize = 0.5 m\n')
        f.write('#b:Ge/Patient/CT/Invisible = "1"\n')
        
    print(f'TOPAS file {filename} saved succesfully') 