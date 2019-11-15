# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2018 replay file
# Internal Version: 2017_11_07-18.21.41 127140
# Run by frederic on Thu Mar 28 17:34:02 2019
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(1.36719, 1.36719), width=201.25, 
    height=135.625)
session.viewports['Viewport: 1'].makeCurrent()
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
execfile(
    '/home/frederic/Desktop/For_Frederic/For_Frederic/MATLAB_code/NU_Daniel_Material/DIST01/mesh_quad_per.py', 
    __main__.__dict__)
#: A new model database has been created.
#: The model "Model-1" has been created.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
#* NameError: name 'allCells' is not defined
#* File 
#* "/home/frederic/Desktop/For_Frederic/For_Frederic/MATLAB_code/NU_Daniel_Material/DIST01/mesh_quad_per.py", 
#* line 97, in <module>
#*     a.setMeshControls(regions=allCells, technique=SWEEP, allowMapped=False,
