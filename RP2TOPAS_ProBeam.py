import numpy as np

p_Emean = np.poly1d([0.9952, 0.6901]) # mean energy, Rahman Med Phys 47(12) 2020.
p_Esigma = np.poly1d([-1.02E-7, 4.05E-5, -9.99E-3, 1.79]) # energy spread [%], Rahman Med Phys 47(12) 2020.
p_size = np.poly1d([-4.71E-7, 2.6E-4, -4.64E-2, 5.84]) # beam sigma (mm), Rahman Med Phys 47(12) 2020.
p_div = np.poly1d([-8.84E-10, 5.39E-7, 1.23E-4, 1.29E-2]) # beam divergence
p_cor = np.poly1d([-8E-7, 0.0015, 0.5143]) # beam correlation
p_MpMU = np.poly1d([-2.99E-5, 2.66E-2, 4.88E-1]) # Mega (1E6) protons per MU coefficient

def PatientTranslation(iso, PatientSupportAngle, bbox, CT):
    """
    return the patient translation coordinates
    iso: isocenter coordinates (mm)
    PatientSupportAngle: table angle (deg)
    bbox: bounding box of the dose matrix (voxel)
    """

    # 1. isocenter translation
    tx = CT.origin[0] + 0.5*sum(bbox[0])*CT.spacing[0] - iso[0]
    ty = CT.origin[1] + 0.5*sum(bbox[1])*CT.spacing[1] - iso[1]
    tz = CT.origin[2] + 0.5*sum(bbox[2])*CT.spacing[2] - iso[2]

    # 2. Table Rotation
    cos, sin = np.cos(np.deg2rad(PatientSupportAngle)), np.sin(np.deg2rad(PatientSupportAngle))
    
    return np.array([[cos,0,sin],[0,1,0],[-sin,0,cos]]) @ [tx,ty,tz] # rotation around Y axis

def CreateTopasFilePBS(field, CT, bbox, filedir, vox_res=1.):
    """ ProBeam source model from https://github.com/mr3536/VarianProbeamTOPASModel """
    
    # print('\nfield: {vars(field)}\nCT: {vars(CT)}\nbbox: {bbox}')
    
    N, x, y, E, MU = field.spots.N, np.array(field.spots.x), np.array(field.spots.y), np.array(field.spots.E), np.array(field.spots.MU)
    TransX, TransY, TransZ = PatientTranslation(field.IsocenterPosition, field.PatientSupportAngle, bbox, CT)
    SADx, SADy = 2000, 2560 # ProBeam Source Axis Distance [mm]
    GRAPHIC_SUPPORT = True
    # GRAPHIC_SUPPORT = False
    
    filename = filedir + f'/Main_{field.parent.PatientID.replace(" ", "_")}_{field.ID.replace(" ", "_")}.txt'
    with open(filename, 'w') as f:
    
        f.write('################################\n')
        f.write('# TOPAS MC v3.9 simulation file\n')
        f.write(f'# Patient    : {field.parent.PatientID}\n')
        f.write(f'# Field ID   : {field.ID}\n')
        f.write(f'# N spots    : {field.spots.N}\n')
        f.write(f'# Gantry     : {field.GantryAngle:.1f} deg\n')
        f.write(f'# Table      : {field.PatientSupportAngle:.1f} deg\n')
        f.write(f'# Isocenter  : {field.IsocenterPosition[0]:.2f}, {field.IsocenterPosition[1]:.2f}, {field.IsocenterPosition[2]:.2f} mm\n')
        f.write(f'# RSWET      : {field.RSWET:.1f} mm\n')
        f.write(f'# Resolution : {vox_res:.2f} mm\n')
        f.write(f'# Machine    : Varian ProBeam\n')
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

        f.write('#### PROBEAM NOZZLE GEOMETRY\n')
        f.write('s:Ge/Nozzle/Message = "Constructing ProBeam Nozzle"\n')
        f.write('s:Ge/Nozzle/Parent = "World"\n')
        f.write('s:Ge/Nozzle/Type = "Group"\n')
        f.write('d:Ge/Nozzle/TransX = 0 m\n')
        f.write('d:Ge/Nozzle/TransY = 0 m\n')
        f.write('d:Ge/Nozzle/TransZ = 0 m\n')
        f.write('d:Ge/Nozzle/RotX = 90 deg\n')
        f.write(f'd:Ge/Nozzle/RotY = {field.GantryAngle:.1f} deg # Gantry angle\n')
        f.write('d:Ge/Nozzle/RotZ = 0 deg\n\n')

        f.write('### Scanning Magnets \n')
        f.write('includeFile = Nozzle_ScanningMagnets.txt\n')
        f.write('s:Ge/SM1/Parent = "Nozzle"\n')
        f.write(f'd:Ge/SM1/TransZ = -{SADy:.1f} mm\n')
        f.write('s:Ge/SM2/Parent = "Nozzle"\n')
        f.write(f'd:Ge/SM2/TransZ = -{SADx:.1f} mm\n\n')

        f.write('### Snout (source)\n')
        f.write('s:Ge/Snout/Parent = "Nozzle"\n')
        f.write('s:Ge/Snout/Type =  "TsBox"\n')
        f.write('s:Ge/Snout/Material = "Air"\n')
        f.write('d:Ge/Snout/HLX = 20 cm\n')
        f.write('d:Ge/Snout/HLY = 20 cm\n')
        f.write('d:Ge/Snout/HLZ = .05 cm\n')
        f.write('d:Ge/Snout/TransX = 0. m\n')
        f.write('d:Ge/Snout/TransY = 0. m\n')
        # f.write(f'd:Ge/Snout/TransZ = -{field.SnoutPosition:.1f} cm\n')
        f.write(f'd:Ge/Snout/TransZ = -420 mm # fixed position\n')
        f.write('d:Ge/Snout/RotX = 0. deg\n')
        f.write('d:Ge/Snout/RotY = 0. deg\n')
        f.write('d:Ge/Snout/RotZ = 0. deg\n\n')
        
        if field.RSWET>0:
            RSthickness = field.RSWET/1.157 # PMMA (WER=1.157)
            f.write('### Range Shifter\n')
            f.write('s:Ge/RS/Parent = "Nozzle"\n')
            f.write('s:Ge/RS/Type =  "TsBox"\n')
            f.write('s:Ge/RS/Material = "G4_WATER"\n')
            # f.write('s:Ge/RS/Material = "G4_PLEXIGLASS"\n') # PMMA
            f.write('d:Ge/RS/HLX = 10 cm\n')
            f.write('d:Ge/RS/HLY = 10 cm\n')
            f.write(f'd:Ge/RS/HLZ = {0.5*field.RSWET:.1f} mm\n')
            # f.write(f'd:Ge/RS/HLZ = {0.5*RSthickness:.1f} mm\n')
            f.write('d:Ge/RS/TransX = 0. m\n')
            f.write('d:Ge/RS/TransY = 0. m\n')
            f.write(f'd:Ge/RS/TransZ = -{field.SnoutPosition-0.1:.1f} cm + Ge/RS/HLZ # retraction\n')
            f.write('d:Ge/RS/RotX = 0. deg\n')
            f.write('d:Ge/RS/RotY = 0. deg\n')
            f.write('d:Ge/RS/RotZ = 0. deg\n\n')
            
        f.write('### Beam Position\n')
        f.write('s:Ge/BeamPosition/Parent = "Snout"\n')
        f.write('s:Ge/BeamPosition/Type = "Group"\n')
        f.write('d:Ge/BeamPosition/TransZ = 0. m\n')
        f.write('d:Ge/BeamPosition/RotZ = 0. deg\n\n')

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
            f.write(f'd:Ge/CT_visu/TransX = {-0.5*sum(bbox[0])*CT.spacing[0] - CT.origin[0]:.1f} mm\n')
            f.write(f'd:Ge/CT_visu/TransY = {-0.5*sum(bbox[1])*CT.spacing[1] - CT.origin[1]:.1f} mm\n')
            f.write(f'd:Ge/CT_visu/TransZ = {-0.5*sum(bbox[2])*CT.spacing[2] - CT.origin[2]:.1f} mm\n')
            f.write('s:Ge/CT_visu/DrawingStyle = "Wireframe"\n')
            f.write(f's:Ge/CT_visu/InputFile = "Geometry/{field.parent.PatientID}"\n')
            f.write('s:Ge/CT_visu/FileFormat = "stl"\n')
            f.write('d:Ge/CT_visu/Units = 1.0 mm\n\n')

        f.write('### Patient modelisation\n')
        f.write('s:Ge/Patient/CT/Parent = "Patient"\n')
        f.write('s:Ge/Patient/CT/Type = "TsDicomPatient"\n')
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

        f.write('### Crop CT dimensions\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsXMin = {min(bbox[0])+1}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsXMax = {max(bbox[0])}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsYMin = {min(bbox[1])+1}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsYMax = {max(bbox[1])}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsZMin = {min(bbox[2])+1}\n')
        f.write(f'i:Ge/Patient/CT/RestrictVoxelsZMax = {max(bbox[2])}\n\n')

        f.write('### Dose computation\n')
        f.write('s:Sc/PatientScorer/Component = "Patient/CT"\n')
        f.write('s:Sc/PatientScorer/Quantity = "DoseToWater"\n')
        f.write(f's:Sc/PatientScorer/OutputFile = "./output/Dose_{field.parent.PatientID}_{field.ID}"\n')
        f.write('s:Sc/PatientScorer/OutputType = "DICOM"\n')
        f.write('b:Sc/PatientScorer/OutputToConsole = "0"\n')
        f.write('s:Sc/PatientScorer/IfOutputFileAlreadyExists = "Increment" # "Overwrite"\n')
        f.write(f'i:Sc/PatientScorer/XBins = {int((max(bbox[0]) - min(bbox[0]))*CT.spacing[0]/vox_res):.0f}\n')
        f.write(f'i:Sc/PatientScorer/YBins = {int((max(bbox[1]) - min(bbox[1]))*CT.spacing[1]/vox_res):.0f}\n')
        f.write(f'i:Sc/PatientScorer/ZBins = {int((max(bbox[2]) - min(bbox[2]))*CT.spacing[2]/vox_res):.0f}\n\n')

        f.write('### PHYSICS (TOPAS default settings)\n')
        f.write('s:Ph/ListName = "Default"\n')
        f.write('s:Ph/Default/Type = "Geant4_Modular"\n')
        f.write('sv:Ph/Default/Modules = 7 "g4em-standard_opt3" "g4h-phy_QGSP_BIC_HP" "g4decay" "g4ion-binarycascade" "g4h-elastic_HP" "g4stopping" "g4radioactivedecay"\n')
        f.write('d:Ph/Default/EMRangeMin = 100 eV\n')
        f.write('d:Ph/Default/EMRangeMax = 250 MeV\n')
        f.write('i:Ph/Default/EMBinsPerDecade = 100\n')
        f.write('d:Ph/Default/CutForAllParticles = 0.05 mm\n\n')
        
        # LMGW = ["Snout"]
        LMGW = []
        if field.RSWET > 0: LMGW += ["RS"]
        if GRAPHIC_SUPPORT: LMGW += ["CT_visu"]
        if len(LMGW)>0: f.write(f'sv:Ph/Default/LayeredMassGeometryWorlds = {len(LMGW)} ' + ' '.join(f'"{i}"' for i in LMGW) +'\n\n')

        f.write('### Time Features Parameters\n')
        f.write(f'i:Tf/NumberOfSequentialTimes = {N}\n') 
        f.write('d:Tf/TimelineStart = 0 s\n')
        f.write(f'd:Tf/TimelineEnd = {N+1} s\n\n')

        f.write('### Spot Energy\n')
        f.write('s:Tf/Energy/Function = "Step"\n')
        f.write(f'dv:Tf/Energy/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/Energy/Values = {N} ' + ' '.join(f'{i:.2f}' for i in p_Emean(E)) + ' MeV\n\n')

        f.write('##### Energy Spread\n')
        f.write('s:Tf/EnergySpread/Function = "Step"\n')
        f.write(f'dv:Tf/EnergySpread/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'uv:Tf/EnergySpread/Values = {N} ' + ' '.join(f'{i:.2f}' for i in p_Esigma(E)) +'\n\n')

        f.write('##### Spot X Positions\n')
        f.write('s:Tf/XPos/Function = "Step"\n')
        f.write(f'dv:Tf/XPos/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/XPos/Values = {N} ' + ' '.join(f'{i:.2f}' for i in ((SADx-420)/SADx)*x) + ' mm\n\n')

        f.write('##### Spot Y Positions\n')
        f.write('s:Tf/YPos/Function = "Step"\n')
        f.write(f'dv:Tf/YPos/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/YPos/Values = {N} ' + ' '.join(f'{-i:.2f}' for i in ((SADy-420)/SADy)*y) + ' mm\n\n')

        f.write('##### Spot X Angle\n')
        f.write('s:Tf/XAngle/Function = "Step"\n')
        f.write(f'dv:Tf/XAngle/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/XAngle/Values = {N} ' + ' '.join(f'{-np.rad2deg(np.arctan(i/SADy)):.2f}' for i in y) + ' deg\n\n')

        f.write('##### Spot Y Angle\n')
        f.write('s:Tf/YAngle/Function = "Step"\n')
        f.write(f'dv:Tf/YAngle/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/YAngle/Values = {N} ' + ' '.join(f'{-np.rad2deg(np.arctan(i/SADx)):.2f}' for i in x) + ' deg\n\n')

        f.write('##### Spot Sigmas\n')
        f.write('s:Tf/SpotSize/Function = "Step"\n')
        f.write(f'dv:Tf/SpotSize/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'dv:Tf/SpotSize/Values = {N} ' + ' '.join(f'{i:.2f}' for i in p_size(E)) + ' mm\n\n')

        f.write('##### Spot Divergences\n')
        f.write('s:Tf/BeamDivergence/Function = "Step"\n')
        f.write(f'dv:Tf/BeamDivergence/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'uv:Tf/BeamDivergence/Values = {N} ' + ' '.join(f'{i:.4f}' for i in p_div(E)) +'\n\n')

        f.write('##### Spot correlation\n')
        f.write('s:Tf/BeamCorrelation/Function = "Step"\n')
        f.write(f'dv:Tf/BeamCorrelation/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'uv:Tf/BeamCorrelation/Values = {N} ' + ' '.join(f'{i:.4f}' for i in p_cor(E)) +'\n\n')

        f.write('##### Number of Histories\n')
        f.write('s:Tf/Particles/Function = "Step"\n')
        f.write(f'dv:Tf/Particles/Times = {N} ' + ' '.join(f'{i}' for i in range(1,N+1)) + ' s\n')
        f.write(f'iv:Tf/Particles/Values = {N} ' + ' '.join(f'{int(i):.0f}' for i in 1000*MU*p_MpMU(E)) +'\n\n') 

        f.write('### SOURCE\n')
        f.write(f'd:Ge/BeamPosition/TransX = Tf/XPos/Value mm\n')
        f.write(f'd:Ge/BeamPosition/TransY = Tf/YPos/Value mm\n')
        f.write(f'd:Ge/BeamPosition/RotX = Tf/XAngle/Value deg\n')
        f.write(f'd:Ge/BeamPosition/RotY = Tf/YAngle/Value deg\n')        
        f.write('s:So/Default/Type = "emittance"\n')
        f.write('s:So/Default/Component = "BeamPosition"\n')
        f.write('s:So/Default/BeamParticle = "proton"\n')
        f.write('d:So/Default/BeamEnergy = Tf/Energy/Value MeV\n')
        f.write('u:So/Default/BeamEnergySpread = Tf/EnergySpread/Value\n')
        f.write('s:So/Default/Distribution = "BiGaussian"\n')
        f.write('d:So/Default/SigmaX = Tf/SpotSize/Value mm\n')
        f.write('u:So/Default/SigmaXprime = Tf/BeamDivergence/Value\n')
        f.write('u:So/Default/CorrelationX = Tf/BeamCorrelation/Value\n')
        f.write('d:So/Default/SigmaY = Tf/SpotSize/Value mm\n')
        f.write('u:So/Default/SigmaYprime = Tf/BeamDivergence/Value\n')
        f.write('u:So/Default/CorrelationY = Tf/BeamCorrelation/Value\n')
        f.write('i:So/Default/NumberOfHistoriesInRun = 1 * Tf/Particles/Value\n\n')

        f.write('### GRAPHICS\n')
        f.write(f'b:Gr/Enable = "{GRAPHIC_SUPPORT}"\n')
        f.write('s:Gr/ViewA/Type = "VRML"\n')
        f.write('b:Ge/World/Invisible = "1"\n')
        f.write('b:Gr/ViewA/IncludeAxes = "1"\n')
        f.write('d:Gr/ViewA/AxesSize = 0.5 m\n')
        f.write('#b:Ge/Patient/CT/Invisible = "1"\n')
        
    print(f'TOPAS file {filename} saved succesfully') 