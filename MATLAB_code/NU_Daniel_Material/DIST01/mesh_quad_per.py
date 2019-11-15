# -*- coding: mbcs -*-
#
# Abaqus/CAE script
# Created by A. R. Melro and M.A. Bessa on 28-Mar-2019 17:34:00
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=196.9, 
    height=171.7)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
Mdb()
#
# Fibres Sketches
#
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize= 1.200000000000000)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.EllipseByCenterPerimeter(center=( 0.200000000000000, 0.100000000000000), axisPoint1=(
     0.420000000000000,  0.100000000000000)
   ,axisPoint2=(
     0.200000000000000,  0.210000000000000))
mdb.models['Model-1'].ConstrainedSketch(name='Fibre_Sketch_Prov', objectToCopy=s)
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='Fibre_Sketch_Prov')
s.unsetPrimaryObject()
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=1.20000)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0.0, 0.0), point2=(0.40000, 0.20000))
s.rectangle(point1=(-0.44000, -0.22000), point2=(0.84000, 0.42000))
mdb.models['Model-1'].sketches.changeKey(fromName='__profile__', 
    toName='Fibre_Sketch_Trim')
s.unsetPrimaryObject()
#
# Fibres Part (Provisional)
#
s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=1.20000)
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=STANDALONE)
s1.sketchOptions.setValues(gridOrigin=(0.0, 0.0))
s1.retrieveSketch(sketch=mdb.models['Model-1'].sketches['Fibre_Sketch_Prov'])
p = mdb.models['Model-1'].Part(name='Fibre_Part_Prov', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Fibre_Part_Prov']
p.BaseSolidExtrude(sketch=s1, depth=0.40000)
s1.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Fibre_Part_Prov']
del mdb.models['Model-1'].sketches['__profile__']
del mdb.models['Model-1'].sketches['Fibre_Sketch_Prov']
#
# Fibres Trimmer Part
#
s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
    sheetSize=1.20000)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.sketchOptions.setValues(gridOrigin=(0.0, 0.0))
s.retrieveSketch(sketch=mdb.models['Model-1'].sketches['Fibre_Sketch_Trim'])
p = mdb.models['Model-1'].Part(name='Fibre_Part_Trim', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
p = mdb.models['Model-1'].parts['Fibre_Part_Trim']
p.BaseSolidExtrude(sketch=s, depth=0.40000)
s.unsetPrimaryObject()
p = mdb.models['Model-1'].parts['Fibre_Part_Trim']
del mdb.models['Model-1'].sketches['__profile__']
del mdb.models['Model-1'].sketches['Fibre_Sketch_Trim']
#
# Fibre Part
#
a = mdb.models['Model-1'].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models['Model-1'].parts['Fibre_Part_Prov']
a.Instance(name='Fibre_Instance_Prov', part=p, dependent=OFF)
p = mdb.models['Model-1'].parts['Fibre_Part_Trim']
a.Instance(name='Fibre_Instance_Trim', part=p, dependent=OFF)
a.InstanceFromBooleanCut(name='Fibre_Final_Part', 
    instanceToBeCut=mdb.models['Model-1'].rootAssembly.instances['Fibre_Instance_Prov'], 
    cuttingInstances=(a.instances['Fibre_Instance_Trim'], ), originalInstances=DELETE)
a.makeIndependent(instances=(a.instances['Fibre_Final_Part-1'], ))
mdb.models['Model-1'].rootAssembly.features.changeKey(
    fromName='Fibre_Final_Part-1', toName='Fibre_Final_Instance')
del mdb.models['Model-1'].parts['Fibre_Part_Prov']
del mdb.models['Model-1'].parts['Fibre_Part_Trim']
#
# Assembly and Meshing
a.setMeshControls(regions=allCells, technique=SWEEP, allowMapped=False, 
    elemShape=HEX_DOMINATED)
mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1,minSizeFactor=0.1, regions=(
mdb.models['Model-1'].rootAssembly.instances['Fibre_Final_Instance'], ), size=0.01)
mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
mdb.models['Model-1'].rootAssembly.instances['Fibre_Final_Instance'], ))
#
# Viewports
#
session.viewports['Viewport: 1'].assemblyDisplay.setValues(mesh=ON)
session.viewports['Viewport: 1'].view.fitView()
session.viewports['Viewport: 1'].view.setProjection(projection=PARALLEL)
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.meshOptions.setValues(
    meshTechnique=ON)
#
# Job creation
#
mdb.Job(name='mesh_prov', model='Model-1', type=ANALYSIS, explicitPrecision=SINGLE, 
    nodalOutputPrecision=SINGLE, description='', 
    parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT, 
    numDomains=1, userSubroutine='', numCpus=1, memory=90, 
    memoryUnits=PERCENTAGE, scratch='', echoPrint=OFF, modelPrint=OFF, 
    contactPrint=OFF, historyPrint=OFF)
import os
os.chdir(r'/home/frederic/Desktop/For_Frederic/For_Frederic/MATLAB_code/NU_Daniel_Material/DIST01')
mdb.jobs['mesh_prov'].writeInput(consistencyChecking=OFF)
