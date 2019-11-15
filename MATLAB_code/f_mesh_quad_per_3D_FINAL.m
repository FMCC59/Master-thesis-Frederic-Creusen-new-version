%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Generate a Mesh in 3D                                 %
%                                                                         %
%  This file is part of the FEA_Automatic.m sequence of files             %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v3.0                                   %
%                                                                         %
%                Antonio Rui Melro - antonio.melro@fe.up.pt               %
%                             February 2009                               %
%                                                                         %
%                Miguel Anibal Bessa - mbessa@u.northwestern.edu          %
%                             September 2015                              %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [status_mesh] = f_mesh_quad_per_3D_FINAL
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global a b c R Fibre_pos N_fibre mat_name dir_name dummy1 dummy2 dummy3 ...
    S_base Em vm Xtm Xcm GcTm Ef1 vf12 GcTf vf23 Xtf vf21 abaqus_path ...
    cohesive_choice explicit_option hardening_curves N R_A R_B;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Material Properties From Text File                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(cd,'/Input/',mat_name,'.mat');
fid = fopen(file_name,'r');
%
for i=1:1:17, fgetl(fid); end
Ef1 = str2num(fgetl(fid)); Ef2 = str2num(fgetl(fid)); fgetl(fid);
vf12 = str2num(fgetl(fid)); fgetl(fid);
vf21 = vf12*Ef2/Ef1;
Gf12 = str2num(fgetl(fid)); Gf23 = str2num(fgetl(fid)); fgetl(fid);
vf23 = Ef2/(2*Gf23)-1;
alphaf1 = str2num(fgetl(fid)); alphaf2 = str2num(fgetl(fid)); fgetl(fid);
GcTf = str2num(fgetl(fid)); fgetl(fid);
Xtf = str2num(fgetl(fid)); Xcf = str2num(fgetl(fid)); fgetl(fid);
Densf = str2num(fgetl(fid));
%
for i=1:1:9, fgetl(fid); end
Em = str2num(fgetl(fid)); fgetl(fid); vm = str2num(fgetl(fid));
fgetl(fid); alpham = str2num(fgetl(fid));
fgetl(fid); T0 = str2num(fgetl(fid));
fgetl(fid); niup = str2num(fgetl(fid)); fgetl(fid);
GcTm = str2num(fgetl(fid)); fgetl(fid);
Xtm = str2num(fgetl(fid)); Xcm = str2num(fgetl(fid));
fgetl(fid); Densm = str2num(fgetl(fid));
fgetl(fid); Tref = str2num(fgetl(fid));
fgetl(fid); EPRref = str2num(fgetl(fid));
fgetl(fid); alpha_c = str2num(fgetl(fid));
fgetl(fid); alpha_t = str2num(fgetl(fid));
fgetl(fid); alpha_s = str2num(fgetl(fid));
fgetl(fid); beta_c = str2num(fgetl(fid));
fgetl(fid); beta_t = str2num(fgetl(fid));
fgetl(fid); beta_s = str2num(fgetl(fid));
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize variables                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
node=zeros(6e6,3);    % Maximum size of list of nodes
elem6=zeros(5e6,7);   % Maximum size of list of wedge elements
elem8=zeros(5e6,9);   % Maximum size of list of hex elements
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Python File for Abaqus CAE with Fibre Distribution               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Started creation of Python file');
%
file_name = strcat(dir_name,'/mesh_quad_per.py');
fid = fopen(file_name,'wt');
%
% Heading
%
fprintf(fid,'# -*- coding: mbcs -*-\n');
fprintf(fid,'#\n');
fprintf(fid,'# Abaqus/CAE script\n');
% fprintf(fid,'# Internal Version: 2008_11_03-14.47.54 87371\n');
fprintf(fid,'# Created by A. R. Melro and M.A. Bessa on %s\n',datestr(now));
fprintf(fid,'#\n');
fprintf(fid,'\n');
fprintf(fid,'# from driverUtils import executeOnCaeGraphicsStartup\n');
fprintf(fid,'# executeOnCaeGraphicsStartup()\n');
fprintf(fid,'#: Executing "onCaeGraphicsStartup()" in the site directory ...\n');
fprintf(fid,'from abaqus import *\n');
fprintf(fid,'from abaqusConstants import *\n');
fprintf(fid,'session.Viewport(name=''Viewport: 1'', origin=(0.0, 0.0), width=196.9, \n');
fprintf(fid,'    height=171.7)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].makeCurrent()\n');
fprintf(fid,'session.viewports[''Viewport: 1''].maximize()\n');
fprintf(fid,'from caeModules import *\n');
fprintf(fid,'from driverUtils import executeOnCaeStartup\n');
fprintf(fid,'executeOnCaeStartup()\n');
fprintf(fid,'Mdb()\n');
%
% Fibre Sketches
%
fprintf(fid,'#\n');
fprintf(fid,'# Fibres Sketches\n');
fprintf(fid,'#\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%18.15f)\n',3*a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
for i=1:1:N %N_fibre
    XC = Fibre_pos(i,2); YC = Fibre_pos(i,3);
    fprintf(fid,'s.EllipseByCenterPerimeter(center=(%18.15f,%18.15f), axisPoint1=(\n',XC,YC);
    fprintf(fid,'    %18.15f, %18.15f)\n',XC+R_A,YC);
    fprintf(fid,'   ,axisPoint2=(\n',XC,YC);
    fprintf(fid,'    %18.15f, %18.15f))\n',XC,YC+R_B);
end
fprintf(fid,'mdb.models[''Model-1''].ConstrainedSketch(name=''Fibre_Sketch_Prov'', objectToCopy=s)\n');
fprintf(fid,'mdb.models[''Model-1''].sketches.changeKey(fromName=''__profile__'', \n');
fprintf(fid,'    toName=''Fibre_Sketch_Prov'')\n');
fprintf(fid,'s.unsetPrimaryObject()\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',3*a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.rectangle(point1=(0.0, 0.0), point2=(%7.5f, %7.5f))\n',a,b);
fprintf(fid,'s.rectangle(point1=(%7.5f, %7.5f), point2=(%7.5f, %7.5f))\n',-2*R_A,-2*R_B,a+2*R_A,b+2*R_B);
fprintf(fid,'mdb.models[''Model-1''].sketches.changeKey(fromName=''__profile__'', \n');
fprintf(fid,'    toName=''Fibre_Sketch_Trim'')\n');
fprintf(fid,'s.unsetPrimaryObject()\n');
%
% Matrix Sketch
%
% fprintf(fid,'#\n');
% fprintf(fid,'# Matrix Sketch\n');
% fprintf(fid,'#\n');
% fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
% fprintf(fid,'    sheetSize=%7.5f)\n',3*a);
% fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
% fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
% fprintf(fid,'s.rectangle(point1=(0.0, 0.0), point2=(%7.5f, %7.5f))\n',a,b);
% fprintf(fid,'mdb.models[''Model-1''].sketches.changeKey(fromName=''__profile__'', \n');
% fprintf(fid,'    toName=''Matrix_Sketch'')\n');
% fprintf(fid,'s.unsetPrimaryObject()\n');
%
% Provisional Fibre Part
%
fprintf(fid,'#\n');
fprintf(fid,'# Fibres Part (Provisional)\n');
fprintf(fid,'#\n');
fprintf(fid,'s1 = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',3*a);
fprintf(fid,'g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints\n');
fprintf(fid,'s1.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s1.sketchOptions.setValues(gridOrigin=(0.0, 0.0))\n');
fprintf(fid,'s1.retrieveSketch(sketch=mdb.models[''Model-1''].sketches[''Fibre_Sketch_Prov''])\n');
fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''Fibre_Part_Prov'', dimensionality=THREE_D, \n');
fprintf(fid,'    type=DEFORMABLE_BODY)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'p.BaseSolidExtrude(sketch=s1, depth=%7.5f)\n',c);
fprintf(fid,'s1.unsetPrimaryObject()\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''Fibre_Sketch_Prov'']\n');
%
% Fibre Part Trimmer Only
%
fprintf(fid,'#\n');
fprintf(fid,'# Fibres Trimmer Part\n');
fprintf(fid,'#\n');
fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
fprintf(fid,'    sheetSize=%7.5f)\n',3*a);
fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
fprintf(fid,'s.sketchOptions.setValues(gridOrigin=(0.0, 0.0))\n');
fprintf(fid,'s.retrieveSketch(sketch=mdb.models[''Model-1''].sketches[''Fibre_Sketch_Trim''])\n');
fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''Fibre_Part_Trim'', dimensionality=THREE_D, \n');
fprintf(fid,'    type=DEFORMABLE_BODY)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
fprintf(fid,'p.BaseSolidExtrude(sketch=s, depth=%7.5f)\n',c);
fprintf(fid,'s.unsetPrimaryObject()\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
fprintf(fid,'del mdb.models[''Model-1''].sketches[''Fibre_Sketch_Trim'']\n');
%
% Matrix Part in Cube Shape Only
% %
% fprintf(fid,'#\n');
% fprintf(fid,'# Matrix Part (Cube only)\n');
% fprintf(fid,'#\n');
% fprintf(fid,'s = mdb.models[''Model-1''].ConstrainedSketch(name=''__profile__'', \n');
% fprintf(fid,'    sheetSize=%7.5f)\n',3*a);
% fprintf(fid,'g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints\n');
% fprintf(fid,'s.setPrimaryObject(option=STANDALONE)\n');
% fprintf(fid,'s.sketchOptions.setValues(gridOrigin=(0.0, 0.0))\n');
% fprintf(fid,'s.retrieveSketch(sketch=mdb.models[''Model-1''].sketches[''Matrix_Sketch''])\n');
% fprintf(fid,'session.viewports[''Viewport: 1''].view.fitView()\n');
% fprintf(fid,'p = mdb.models[''Model-1''].Part(name=''Matrix_Part'', dimensionality=THREE_D, \n');
% fprintf(fid,'    type=DEFORMABLE_BODY)\n');
% fprintf(fid,'p = mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
% fprintf(fid,'p.BaseSolidExtrude(sketch=s, depth=%7.5f)\n',c);
% fprintf(fid,'s.unsetPrimaryObject()\n');
% fprintf(fid,'p = mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
% fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=p)\n');
% fprintf(fid,'del mdb.models[''Model-1''].sketches[''__profile__'']\n');
% fprintf(fid,'del mdb.models[''Model-1''].sketches[''Matrix_Sketch'']\n');
% %
% Creation of Real Fibre Part
%
fprintf(fid,'#\n');
fprintf(fid,'# Fibre Part\n');
fprintf(fid,'#\n');
fprintf(fid,'a = mdb.models[''Model-1''].rootAssembly\n');
fprintf(fid,'a.DatumCsysByDefault(CARTESIAN)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'a.Instance(name=''Fibre_Instance_Prov'', part=p, dependent=OFF)\n');
fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
fprintf(fid,'a.Instance(name=''Fibre_Instance_Trim'', part=p, dependent=OFF)\n');
fprintf(fid,'a.InstanceFromBooleanCut(name=''Fibre_Final_Part'', \n');
fprintf(fid,'    instanceToBeCut=mdb.models[''Model-1''].rootAssembly.instances[''Fibre_Instance_Prov''], \n');
fprintf(fid,'    cuttingInstances=(a.instances[''Fibre_Instance_Trim''], ), originalInstances=DELETE)\n');
fprintf(fid,'a.makeIndependent(instances=(a.instances[''Fibre_Final_Part-1''], ))\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.features.changeKey(\n');
fprintf(fid,'    fromName=''Fibre_Final_Part-1'', toName=''Fibre_Final_Instance'')\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Fibre_Part_Prov'']\n');
fprintf(fid,'del mdb.models[''Model-1''].parts[''Fibre_Part_Trim'']\n');
%
% Creation of Real Matrix Part
%
% fprintf(fid,'#\n');
% fprintf(fid,'# Matrix Part\n');
% fprintf(fid,'#\n');
% fprintf(fid,'a = mdb.models[''Model-1''].rootAssembly\n');
% fprintf(fid,'p = mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
% fprintf(fid,'a.Instance(name=''Matrix_Instance'', part=p, dependent=OFF)\n');
% fprintf(fid,'a.InstanceFromBooleanCut(name=''Matrix_Final_Part'', \n');
% fprintf(fid,'    instanceToBeCut=mdb.models[''Model-1''].rootAssembly.instances[''Matrix_Instance''], \n');
% fprintf(fid,'    cuttingInstances=(a.instances[''Fibre_Final_Instance''], ), originalInstances=DELETE)\n');
% fprintf(fid,'a.makeIndependent(instances=(a.instances[''Matrix_Final_Part-1''], ))\n');
% fprintf(fid,'mdb.models[''Model-1''].rootAssembly.features.changeKey(\n');
% fprintf(fid,'    fromName=''Matrix_Final_Part-1'', toName=''Matrix_Final_Instance'')\n');
% fprintf(fid,'del mdb.models[''Model-1''].parts[''Matrix_Part'']\n');
% %
% Creation of Final Geometry
%
fprintf(fid,'#\n');
fprintf(fid,'# Assembly and Meshing\n');
fprintf(fid,'a.setMeshControls(regions=allCells, technique=SWEEP, allowMapped=False, \n');
fprintf(fid,'    elemShape=HEX_DOMINATED)\n')
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.seedPartInstance(deviationFactor=0.1,minSizeFactor=0.1, regions=(\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.instances[''Fibre_Final_Instance''], ), size=0.01)\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.generateMesh(regions=(\n');
fprintf(fid,'mdb.models[''Model-1''].rootAssembly.instances[''Fibre_Final_Instance''], ))\n');

% fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Final_Part'']\n');
% fprintf(fid,'a.Instance(name=''Fibre_Final_Instance'', part=p, dependent=OFF)\n');
% fprintf(fid,'p = mdb.models[''Model-1''].parts[''Fibre_Final_Part'']\n');
% fprintf(fid,'a.Instance(name=''Final_Stuff-1'', part=p, dependent=OFF)\n');
% %fprintf(fid,'a.InstanceFromBooleanMerge(name=''Final_Stuff'', instances=(\n');
% % a.instances[''Matrix_Final_Instance''], 
% %fprintf(fid,'a.instances[''Fibre_Final_Instance''], ), \n');
% %fprintf(fid,'keepIntersections=ON, originalInstances=DELETE, domain=GEOMETRY)\n');
% %fprintf(fid,'del mdb.models[''Model-1''].parts[''Matrix_Final_Part'']\n');
% fprintf(fid,'del mdb.models[''Model-1''].parts[''Fibre_Final_Part'']\n');
% fprintf(fid,'a.makeIndependent(instances=(a.instances[''Final_Stuff-1''], ))\n');
% fprintf(fid,'partInstances =(a.instances[''Final_Stuff-1''], )\n');
% fprintf(fid,'a.seedPartInstance(regions=partInstances, size=%7.5f, deviationFactor=0.1)\n',S_base);
% fprintf(fid,'elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)\n');
% fprintf(fid,'elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)\n');
% fprintf(fid,'elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)\n');
% fprintf(fid,'c1 = a.instances[''Final_Stuff-1''].cells\n');
% fprintf(fid,'pickedRegions =(c1, )\n');
% fprintf(fid,'a.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, \n');
% fprintf(fid,'    elemType3))\n');
% fprintf(fid,'allCells = a.instances[''Final_Stuff-1''].cells\n');
% fprintf(fid,'a.setMeshControls(regions=allCells, technique=SWEEP, allowMapped=False, \n');
% fprintf(fid,'    elemShape=HEX_DOMINATED)\n');
% fprintf(fid,'a.generateMesh(regions=partInstances, seedConstraintOverride=OFF,\n');
% fprintf(fid,'    meshTechniqueOverride=OFF)\n');
%
% View settings... Not really necessary...
%
fprintf(fid,'#\n');
fprintf(fid,'# Viewports\n');
fprintf(fid,'#\n');
fprintf(fid,'session.viewports[''Viewport: 1''].assemblyDisplay.setValues(mesh=ON)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].view.fitView()\n');
fprintf(fid,'session.viewports[''Viewport: 1''].view.setProjection(projection=PARALLEL)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=a)\n');
fprintf(fid,'session.viewports[''Viewport: 1''].assemblyDisplay.meshOptions.setValues(\n');
fprintf(fid,'    meshTechnique=ON)\n');
%
% Export mesh to .inp file
%
fprintf(fid,'#\n');
fprintf(fid,'# Job creation\n');
fprintf(fid,'#\n');
fprintf(fid,'mdb.Job(name=''mesh_prov'', model=''Model-1'', type=ANALYSIS, explicitPrecision=SINGLE, \n');
fprintf(fid,'    nodalOutputPrecision=SINGLE, description='''', \n');
fprintf(fid,'    parallelizationMethodExplicit=DOMAIN, multiprocessingMode=DEFAULT, \n');
fprintf(fid,'    numDomains=1, userSubroutine='''', numCpus=1, memory=90, \n');
fprintf(fid,'    memoryUnits=PERCENTAGE, scratch='''', echoPrint=OFF, modelPrint=OFF, \n');
fprintf(fid,'    contactPrint=OFF, historyPrint=OFF)\n');
fprintf(fid,'import os\n');
fprintf(fid,'os.chdir(r''%s'')\n',dir_name);
fprintf(fid,'mdb.jobs[''mesh_prov''].writeInput(consistencyChecking=OFF)\n');
%
fclose(fid);
%
disp(' ');
disp('Creation of Python file COMPLETED');
disp('Elapsed Time [min]: ');
disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Execute Python Script in ABAQUS CAE                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Started generation of temporary mesh (without cohesive elements)'); disp(' ');
unix(strcat(abaqus_path,' cae noGUI=',dir_name,'/mesh_quad_per.py'));
disp(' ');
disp('Generation of mesh COMPLETED');
disp('Elapsed Time [min]: ');
disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read mesh.inp file created by Abaqus CAE                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Started creation of .inp files');
%
file_name = strcat(dir_name,'/mesh_prov.inp');
fid = fopen(file_name,'r');
%
for i=1:1:17
    fgetl(fid);
end
%
% Reads nodal coordinates
%
status = 1; ncount = 0;
while status == 1
    tline=fgetl(fid);
    [data status] = str2num(tline);
    if status == 0, break, end
    ncount = ncount + 1;
    node(ncount,2) = data(1,2);
    node(ncount,3) = data(1,3);
    node(ncount,1) = data(1,4);
end
%
% Reads element data for elements
%
status = 1;
if str2num(tline(19)) == 6   % C3D6 elements
    el6_counter = 0;
    while status == 1
        tline=fgetl(fid);
        [data status] = str2num(tline);
        if status == 0, break, end
        el6_counter = el6_counter + 1;
        for i=1:1:7
            elem6(el6_counter,i) = data(1,i);
        end
    end
else                         % C3D8R elements
    el8_counter = 0;
    while status == 1
        tline=fgetl(fid);
        [data status] = str2num(tline);
        if status == 0, break, end
        el8_counter = el8_counter + 1;
        for i=1:1:9
            elem8(el8_counter,i) = data(1,i);
        end
    end
end
%
status = 1; 
if str2num(tline(19)) == 6   % C3D6 elements
    el6_counter = 0;
    while status == 1
        tline=fgetl(fid);
        [data status] = str2num(tline);
        if status == 0, break, end
        el6_counter = el6_counter + 1;
        for i=1:1:7
            elem6(el6_counter,i) = data(1,i);
        end
    end
else                         % C3D8R elements
    el8_counter = 0;
    while status == 1
        tline=fgetl(fid);
        [data status] = str2num(tline);
        if status == 0, break, end
        el8_counter = el8_counter + 1;
        for i=1:1:9
            elem8(el8_counter,i) = data(1,i);
        end
    end
end
%
fclose(fid);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distinguish fibre nodes from matrix nodes                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
node(:,4) = 0;
for i=1:1:ncount
    x1 = node(i,2); y1 = node(i,3);
    for j=1:1:N
        xx = Fibre_pos(j,2); yy = Fibre_pos(j,3);
        DIST_TMP = sqrt((xx-x1)^2 + (yy-y1)^2);
        if DIST_TMP <= 1.001*R
            node(i,4) = j; % it belongs to fibre j
            break
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Search, for each fibre, the nodes in the interface                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
node(:,5) = 0;
for i=1:1:ncount
    x1 = node(i,2); y1 = node(i,3);
    for j=1:1:N
        xx = Fibre_pos(j,2); yy = Fibre_pos(j,3);
        DIST_TMP = sqrt((xx-x1)^2 + (yy-y1)^2);
        if DIST_TMP <= 1.001*R && DIST_TMP >= 0.999*R
            node(i,5) = j; % it belongs to interface of fibre j
            break
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distinguish fibre elements from matrix elements                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
elem6(:,8) = 0;
for i=1:1:el6_counter
    k = 0;
    for j=2:1:7
        if j == 2
            if node(elem6(i,j),4) > 0
                k = k + 1;
            end
        else
            if node(elem6(i,j),4) > 0 && node(elem6(i,j),4)==node(elem6(i,j-1),4)
                k = k + 1;
            end
        end
    end
    if k == 6
        elem6(i,8) = 1;
    end
end
%
elem8(:,10) = 0;
for i=1:1:el8_counter
    k = 0;
    for j=2:1:9
        if j == 2
            if node(elem8(i,j),4) > 0
                k = k + 1;
            end
        else
            if node(elem8(i,j),4) > 0 && node(elem8(i,j),4)==node(elem8(i,j-1),4)
                k = k + 1;
            end
        end
    end
    if k == 8
        elem8(i,10) = 1;
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create cohesive element, if requested                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cohesive_choice == 1
    %
    % N_rows = ceil(c/S_base) + 1;
    N_rows = 1;
    for i=2:1:ncount
        if round(1e8*node(i,2))/1e8 == round(1e8*node(1,2))/1e8 && ...
                round(1e8*node(i,3))/1e8 == round(1e8*node(1,3))/1e8
            N_rows = N_rows + 1;
        end
    end
    %
    coh_thick = 1e-6;
    %
    ncoh = 0; nelcoh = 0;
    for i=1:1:N
        temp_count = 0; nlink_count = 0;
        clear temp_fibre cohe_nlink
        %
        % Collect the nodes in one interface
        %
        for j=1:1:ncount
            if node(j,5) == i
                temp_count = temp_count + 1;
                temp_fibre(temp_count,1) = j;
                temp_fibre(temp_count,2) = node(j,1);
                temp_fibre(temp_count,3) = node(j,2);
                temp_fibre(temp_count,4) = node(j,3);
                temp_fibre(temp_count,5) = node(j,5);
            end
        end
        %
        % Determine their polar coordinates relative to fibre centre
        %
        xc = Fibre_pos(i,2); yc = Fibre_pos(i,3);
        for j=1:1:temp_count
            xx = temp_fibre(j,3); yy = temp_fibre(j,4);
            temp_fibre(j,6) = sqrt((xc-xx)^2 + (yc-yy)^2);
            if Fibre_pos(i,4) == 0
                if yy-yc>0 && xx-xc>0
                    temp_fibre(j,7) = atan((yy-yc)/(xx-xc));
                elseif yy-yc>0 && xx-xc<0
                    temp_fibre(j,7) = pi - atan((yy-yc)/abs(xx-xc));
                elseif yy-yc<0 && xx-xc<0
                    temp_fibre(j,7) = pi + atan(abs(yy-yc)/abs(xx-xc));
                else
                    temp_fibre(j,7) = 2*pi - atan(abs(yy-yc)/(xx-xc));
                end
            else
                if xc > a - R
                    if yy-yc>0 && xx-xc>0
                        temp_fibre(j,7) = atan((yy-yc)/(xx-xc));
                    elseif yy-yc>0 && xx-xc<0
                        temp_fibre(j,7) = pi - atan((yy-yc)/abs(xx-xc));
                    elseif yy-yc<0 && xx-xc<0
                        temp_fibre(j,7) = pi + atan(abs(yy-yc)/abs(xx-xc));
                    else
                        temp_fibre(j,7) = 2*pi - atan(abs(yy-yc)/(xx-xc));
                    end
                elseif xc < R
                    if yy-yc>0 && xx-xc>0
                        temp_fibre(j,7) = pi + atan((yy-yc)/(xx-xc));
                    elseif yy-yc>0 && xx-xc<0
                        temp_fibre(j,7) = 2*pi - atan((yy-yc)/abs(xx-xc));
                    elseif yy-yc<0 && xx-xc<0
                        temp_fibre(j,7) = atan(abs(yy-yc)/abs(xx-xc));
                    else
                        temp_fibre(j,7) = pi - atan(abs(yy-yc)/(xx-xc));
                    end
                elseif yc > b-R
                    if yy-yc>0 && xx-xc>0
                        temp_fibre(j,7) = 2*pi - atan((xx-xc)/(yy-yc));
                    elseif yy-yc>0 && xx-xc<0
                        temp_fibre(j,7) = atan(abs(xx-xc)/(yy-yc));
                    elseif yy-yc<0 && xx-xc<0
                        temp_fibre(j,7) = pi - atan(abs(xx-xc)/abs(yy-yc));
                    else
                        temp_fibre(j,7) = pi + atan((xx-xc)/abs(yy-yc));
                    end
                elseif yc < R
                    if yy-yc>0 && xx-xc>0
                        temp_fibre(j,7) = pi - atan((xx-xc)/(yy-yc));
                    elseif yy-yc>0 && xx-xc<0
                        temp_fibre(j,7) = pi + atan(abs(xx-xc)/(yy-yc));
                    elseif yy-yc<0 && xx-xc<0
                        temp_fibre(j,7) = 2*pi - atan(abs(xx-xc)/abs(yy-yc));
                    else
                        temp_fibre(j,7) = atan((xx-xc)/abs(yy-yc));
                    end
                end
            end
        end
        %
        % Order the nodes based in angular position
        %
        if exist('temp_fibre','var') > 0
            temp_fibre = sortrows(temp_fibre,[2 7]);
        else
            continue
        end
        %
        % Move the nodes to the interior of the fibre
        % And create new nodes in the exterior of the fibre
        %
        for j=1:1:temp_count
            if xc < R
                node(temp_fibre(j,1),2) = node(temp_fibre(j,1),2) + coh_thick*cos(temp_fibre(j,7));
                node(temp_fibre(j,1),3) = node(temp_fibre(j,1),3) + coh_thick*sin(temp_fibre(j,7));
            elseif yc < R && xc < a-R
                node(temp_fibre(j,1),2) = node(temp_fibre(j,1),2) - coh_thick*sin(temp_fibre(j,7));
                node(temp_fibre(j,1),3) = node(temp_fibre(j,1),3) + coh_thick*cos(temp_fibre(j,7));
            elseif yc > b-R && xc < a-R
                node(temp_fibre(j,1),2) = node(temp_fibre(j,1),2) + coh_thick*sin(temp_fibre(j,7));
                node(temp_fibre(j,1),3) = node(temp_fibre(j,1),3) - coh_thick*cos(temp_fibre(j,7));
            else
                node(temp_fibre(j,1),2) = node(temp_fibre(j,1),2) - coh_thick*cos(temp_fibre(j,7));
                node(temp_fibre(j,1),3) = node(temp_fibre(j,1),3) - coh_thick*sin(temp_fibre(j,7));
            end
            if node(temp_fibre(j,1),2) < coh_thick
                node(temp_fibre(j,1),2) = 0;
            end
            if node(temp_fibre(j,1),2) > a-coh_thick
                node(temp_fibre(j,1),2) = a;
            end
            if node(temp_fibre(j,1),3) < coh_thick
                node(temp_fibre(j,1),3) = 0;
            end
            if node(temp_fibre(j,1),3) > b-coh_thick
                node(temp_fibre(j,1),3) = b;
            end
            %
            ncoh = ncoh + 1;
            nlink_count = nlink_count + 1;
            ncohxy(ncoh,1) = temp_fibre(j,2);
            if xc < R
                ncohxy(ncoh,2) = temp_fibre(j,3) - coh_thick*cos(temp_fibre(j,7));
                ncohxy(ncoh,3) = temp_fibre(j,4) - coh_thick*sin(temp_fibre(j,7));
            elseif yc < R && xc < a-R
                ncohxy(ncoh,2) = temp_fibre(j,3) + coh_thick*sin(temp_fibre(j,7));
                ncohxy(ncoh,3) = temp_fibre(j,4) - coh_thick*cos(temp_fibre(j,7));
            elseif yc > b-R && xc < a-R
                ncohxy(ncoh,2) = temp_fibre(j,3) - coh_thick*sin(temp_fibre(j,7));
                ncohxy(ncoh,3) = temp_fibre(j,4) + coh_thick*cos(temp_fibre(j,7));
            else
                ncohxy(ncoh,2) = temp_fibre(j,3) + coh_thick*cos(temp_fibre(j,7));
                ncohxy(ncoh,3) = temp_fibre(j,4) + coh_thick*sin(temp_fibre(j,7));
            end
            if ncohxy(ncoh,2) < coh_thick, ncohxy(ncoh,2) = 0; end
            if ncohxy(ncoh,2) > a-coh_thick, ncohxy(ncoh,2) = a; end
            if ncohxy(ncoh,3) < coh_thick, ncohxy(ncoh,3) = 0; end
            if ncohxy(ncoh,3) > b-coh_thick, ncohxy(ncoh,3) = b; end
            cohe_nlink(nlink_count,1) = temp_fibre(j,1);
            cohe_nlink(nlink_count,2) = ncoh + ncount;
        end
        %
        % Change the associativity of the matrix elements
        %
        for l=1:1:max(length(cohe_nlink))
            for j=1:1:el6_counter
                if elem6(j,8) == 0
                    for k=2:1:7
                        if elem6(j,k) == cohe_nlink(l,1)
                            elem6(j,k) = cohe_nlink(l,2);
                        end
                    end
                end
            end
            for j=1:1:el8_counter
                if elem8(j,10) == 0
                    for k=2:1:9
                        if elem8(j,k) == cohe_nlink(l,1)
                            elem8(j,k) = cohe_nlink(l,2);
                        end
                    end
                end
            end
        end
        %
        % Create Cohesive elements around each fibre
        %
        for l=1:1:N_rows-1
            for j=1:1:temp_count/N_rows-1
                nelcoh = nelcoh + 1;
                elcoh(nelcoh,1) = cohe_nlink(j+temp_count/N_rows*(l-1),1);
                elcoh(nelcoh,2) = cohe_nlink(j+temp_count/N_rows*(l-1)+1,1);
                elcoh(nelcoh,3) = cohe_nlink(j+temp_count/N_rows*l+1,1);
                elcoh(nelcoh,4) = cohe_nlink(j+temp_count/N_rows*l,1);
                elcoh(nelcoh,5) = cohe_nlink(j+temp_count/N_rows*(l-1),2);
                elcoh(nelcoh,6) = cohe_nlink(j+temp_count/N_rows*(l-1)+1,2);
                elcoh(nelcoh,7) = cohe_nlink(j+temp_count/N_rows*l+1,2);
                elcoh(nelcoh,8) = cohe_nlink(j+temp_count/N_rows*l,2);
            end
            if Fibre_pos(i,4) == 0
                j = j+1;
                nelcoh = nelcoh + 1;
                elcoh(nelcoh,1) = cohe_nlink(j+temp_count/N_rows*(l-1),1);
                elcoh(nelcoh,2) = cohe_nlink(temp_count/N_rows*(l-1)+1,1);
                elcoh(nelcoh,3) = cohe_nlink(temp_count/N_rows*l+1,1);
                elcoh(nelcoh,4) = cohe_nlink(j+temp_count/N_rows*l,1);
                elcoh(nelcoh,5) = cohe_nlink(j+temp_count/N_rows*(l-1),2);
                elcoh(nelcoh,6) = cohe_nlink(temp_count/N_rows*(l-1)+1,2);
                elcoh(nelcoh,7) = cohe_nlink(temp_count/N_rows*l+1,2);
                elcoh(nelcoh,8) = cohe_nlink(j+temp_count/N_rows*l,2);
            end
        end
    end
% 
else % if no cohesive elements are requested
    ncoh = 0; nelcoh = 0;
end % end if cohesive_choice == 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate mesh.inp file with the final mesh in new coordinates           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
file_name = strcat(dir_name,'/mesh.inp');
fid = fopen(file_name,'w');
%
% Heading
%
fprintf(fid,'** ***************************************************************************\n');
fprintf(fid,'**   Created by A.R. Melro and M.A. Bessa on %s\n',datestr(now));
fprintf(fid,'**\n');
fprintf(fid,'** ***************************************************************************\n');
fprintf(fid,'**\n');
%
% Nodes
%
fprintf(fid,'*NODE, NSET=GLOBAL\n');
for i=1:1:ncount
    fprintf(fid,' %9d, %9.8f, %9.8f, %9.8f\n',i,node(i,1),node(i,2),node(i,3));
end
%
if cohesive_choice == 1
    for i=1:1:ncoh
        fprintf(fid,' %9d, %9.8f, %9.8f, %9.8f\n',i+ncount,ncohxy(i,1),ncohxy(i,2),ncohxy(i,3));
    end
end
%
% Elements
%
fprintf(fid,'*ELEMENT, TYPE=C3D6, ELSET=MATRIX\n');
for i=1:1:el6_counter
    if elem6(i,8) == 0
        fprintf(fid,' %9d, %9d, %9d, %9d, %9d, %9d, %9d\n',elem6(i,1),elem6(i,2),...
            elem6(i,3),elem6(i,4),elem6(i,5),elem6(i,6),elem6(i,7));
    end
end
%
fprintf(fid,'*ELEMENT, TYPE=C3D8R, ELSET=MATRIX\n');
for i=1:1:el8_counter
    if elem8(i,10) == 0
        fprintf(fid,' %9d, %9d, %9d, %9d, %9d, %9d, %9d, %9d, %9d\n',elem8(i,1),elem8(i,2),...
            elem8(i,3),elem8(i,4),elem8(i,5),elem8(i,6),elem8(i,7),elem8(i,8),elem8(i,9));
    end
end
%
fprintf(fid,'*ELEMENT, TYPE=C3D6, ELSET=FIBRE\n');
for i=1:1:el6_counter
    if elem6(i,8) ~= 0
        fprintf(fid,' %9d, %9d, %9d, %9d, %9d, %9d, %9d\n',elem6(i,1),elem6(i,2),...
            elem6(i,3),elem6(i,4),elem6(i,5),elem6(i,6),elem6(i,7));
    end
end
%
fprintf(fid,'*ELEMENT, TYPE=C3D8R, ELSET=FIBRE\n');
for i=1:1:el8_counter
    if elem8(i,10) ~= 0
        fprintf(fid,' %9d, %9d, %9d, %9d, %9d, %9d, %9d, %9d, %9d\n',elem8(i,1),elem8(i,2),...
            elem8(i,3),elem8(i,4),elem8(i,5),elem8(i,6),elem8(i,7),elem8(i,8),elem8(i,9));
    end
end
%
if cohesive_choice == 1
    fprintf(fid,'*ELEMENT, TYPE=COH3D8, ELSET=COHESI\n');
    for i=1:1:nelcoh
        fprintf(fid,' %9d, %9d, %9d, %9d, %9d, %9d, %9d, %9d, %9d\n',...
            el6_counter+el8_counter+i,elcoh(i,1),elcoh(i,2),elcoh(i,3),...
            elcoh(i,4),elcoh(i,5),elcoh(i,6),elcoh(i,7),elcoh(i,8));
    end
end
%
% Dummy Nodes
%
fprintf(fid,'*NODE, NSET=NDUMMY\n');
dummy1 = ncount+ncoh+1;
dummy2 = ncount+ncoh+2;
dummy3 = ncount+ncoh+3;
fprintf(fid,' %9d, %9.8f, %9.8f, %9.8f\n',dummy1,2*c,0.0,0.0); % Epsilon_1i
fprintf(fid,' %9d, %9.8f, %9.8f, %9.8f\n',dummy2,0.0,2*a,0.0); % Epsilon_2j
fprintf(fid,' %9d, %9.8f, %9.8f, %9.8f\n',dummy3,0.0,0.0,2*b); % Epsilon_3k
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Which Nodes are in Vertices, Edges and Faces of the RVE          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(dir_name,'/pbcs.inp');
fid = fopen(file_name,'w');
%
vert=zeros(1,8);
edge1count=0; edge2count=0; edge3count=0; edge4count=0; edge5count=0;
edge6count=0; edge7count=0; edge8count=0; edge9count=0; edge10count=0;
edge11count=0; edge12count=0;
face1count=0; face2count=0; face3count=0; face4count=0; face5count=0; face6count=0;
%
for i=1:1:ncount
    xx = round(1e8*node(i,1))/1e8;
    yy = round(1e8*node(i,2))/1e8;
    zz = round(1e8*node(i,3))/1e8;
    aa = round(1e8*a)/1e8;
    bb = round(1e8*b)/1e8;
    cc = round(1e8*c)/1e8;
    if xx>=cc && yy==0 && zz==0; vert(1)=i; continue; end
    if xx>=cc && yy>=aa && zz==0; vert(2)=i; continue; end
    if xx>=cc && yy>=aa && zz>=bb; vert(3)=i; continue; end
    if xx>=cc && yy==0 && zz>=bb; vert(4)=i; continue; end
    if xx==0 && yy==0 && zz==0; vert(5)=i; continue; end
    if xx==0 && yy>=aa && zz==0; vert(6)=i; continue; end
    if xx==0 && yy>=aa && zz>=bb; vert(7)=i; continue; end
    if xx==0 && yy==0 && zz>=bb; vert(8)=i; continue; end
    if xx==cc && yy==0 , edge1count = edge1count+1; edge1(edge1count)=i; continue; end
    if xx==cc && yy==aa, edge2count = edge2count+1; edge2(edge2count)=i; continue; end
    if xx==0  && yy==aa, edge3count = edge3count+1; edge3(edge3count)=i; continue; end
    if xx==0  && yy==0 , edge4count = edge4count+1; edge4(edge4count)=i; continue; end
    if xx==cc && zz==0 , edge5count = edge5count+1; edge5(edge5count)=i; continue; end
    if xx==cc && zz==bb, edge6count = edge6count+1; edge6(edge6count)=i; continue; end
    if xx==0  && zz==bb, edge7count = edge7count+1; edge7(edge7count)=i; continue; end
    if xx==0  && zz==0 , edge8count = edge8count+1; edge8(edge8count)=i; continue; end
    if yy==0  && zz==0 , edge9count = edge9count+1; edge9(edge9count)=i; continue; end
    if yy==aa && zz==0 , edge10count = edge10count+1; edge10(edge10count)=i; continue; end
    if yy==aa && zz==bb, edge11count = edge11count+1; edge11(edge11count)=i; continue; end
    if yy==0  && zz==bb, edge12count = edge12count+1; edge12(edge12count)=i; continue; end
    if xx == cc, face1count = face1count + 1; face1(face1count) = i; continue; end
    if yy == aa, face2count = face2count + 1; face2(face2count) = i; continue; end    
    if xx == 0 , face3count = face3count + 1; face3(face3count) = i; continue; end
    if yy == 0 , face4count = face4count + 1; face4(face4count) = i; continue; end
    if zz == 0 , face5count = face5count + 1; face5(face5count) = i; continue; end
    if zz == bb, face6count = face6count + 1; face6(face6count) = i; continue; end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if PBCs were properly applied                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for i=1:1:8, if vert(i) == 0, disp('PBC Application failed. Vertices');...
            status_mesh = 0; return, end, end
if edge1count ~= edge3count, disp('PBC Application failed. Edge 1-3');...
        status_mesh = 0; return, end
if edge2count ~= edge4count, disp('PBC Application failed. Edge 2-4');...
        status_mesh = 0; return, end
if edge6count ~= edge8count, disp('PBC Application failed. Edge 6-8');...
        status_mesh = 0; return, end
if edge5count ~= edge7count, disp('PBC Application failed. Edge 5-7');...
        status_mesh = 0; return, end
if edge11count ~= edge9count, disp('PBC Application failed. Edge 11-9');...
    status_mesh = 0; return, end
if edge10count ~= edge12count, disp('PBC Application failed. Edge 10-12');...
        status_mesh = 0; return, end
if face2count~=face4count, disp('PBC Application failed. Face 2-4');...
        status_mesh = 0; return, end
if face1count~=face3count, disp('PBC Application failed. Face 1-3');...
        status_mesh = 0; return, end
if face5count~=face6count, disp('PBC Application failed. Face 5-6');...
        status_mesh = 0; return, end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define 3D Periodic Boundary Conditions for Vertices                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% Vertices 3-5
%
for ndof=1:3 % Loop over the degrees of freedom
    fprintf(fid,'*EQUATION\n');
    fprintf(fid,' 5\n');
    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
        vert(3),ndof,1,vert(5),ndof,-1,dummy1,ndof,-c,dummy2,ndof,-a);
    fprintf(fid,' %9d, %1d, %9.8f\n',dummy3,ndof,-b);
end
%
% Vertices 2-8
%
for ndof=1:3 % Loop over the degrees of freedom
    fprintf(fid,'*EQUATION\n');
    fprintf(fid,' 5\n');
    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
        vert(2),ndof,1,vert(8),ndof,-1,dummy1,ndof,-c,dummy2,ndof,-a);
    fprintf(fid,' %9d, %1d, %9.8f\n',dummy3,ndof,b);
end
%
% Vertices 7-1
%
for ndof=1:3 % Loop over the degrees of freedom
    fprintf(fid,'*EQUATION\n');
    fprintf(fid,' 5\n');
    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
        vert(7),ndof,1,vert(1),ndof,-1,dummy1,ndof,c,dummy2,ndof,-a);
    fprintf(fid,' %9d, %1d, %9.8f\n',dummy3,ndof,-b);
end
%
% Vertices 4-6
%
for ndof=1:3 % Loop over the degrees of freedom
    fprintf(fid,'*EQUATION\n');
    fprintf(fid,' 5\n');
    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
        vert(4),ndof,1,vert(6),ndof,-1,dummy1,ndof,-c,dummy2,ndof,a);
    fprintf(fid,' %9d, %1d, %9.8f\n',dummy3,ndof,-b);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define 3D Periodic Boundary Conditions for Edges                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Edges 2-4
%
for i=1:1:edge2count
    for j=1:1:edge4count
        if node(edge2(i),3) == node(edge4(j),3)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 4\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    edge2(i),ndof,1,edge4(j),ndof,-1,dummy1,ndof,-c,dummy2,ndof,-a);
            end
            break
        end
    end
end
%
% Edges 1-3
%
for i=1:1:edge1count
    for j=1:1:edge3count
        if node(edge1(i),3) == node(edge3(j),3)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 4\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    edge1(i),ndof,1,edge3(j),ndof,-1,dummy1,ndof,-c,dummy2,ndof,a);
            end
            break
        end
    end
end
%
% Edges 6-8
%
for i=1:1:edge6count
    for j=1:1:edge8count
        if node(edge6(i),2) == node(edge8(j),2)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 4\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    edge6(i),ndof,1,edge8(j),ndof,-1,dummy1,ndof,-c,dummy3,ndof,-b);
            end
            break
        end
    end
end
%
% Edges 5-7
%
for i=1:1:edge5count
    for j=1:1:edge7count
        if node(edge5(i),2) == node(edge7(j),2)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 4\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    edge5(i),ndof,1,edge7(j),ndof,-1,dummy1,ndof,-c,dummy3,ndof,b);
            end
            break
        end
    end
end
%
% Edges 11-9
%
for i=1:1:edge11count
    for j=1:1:edge9count
        if node(edge11(i),1) == node(edge9(j),1)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 4\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    edge11(i),ndof,1,edge9(j),ndof,-1,dummy2,ndof,-a,dummy3,ndof,-b);
            end
            break
        end
    end
end
%
% Edges 10-12
%
for i=1:1:edge10count
    for j=1:1:edge12count
        if node(edge10(i),1) == node(edge12(j),1)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 4\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    edge10(i),ndof,1,edge12(j),ndof,-1,dummy2,ndof,-a,dummy3,ndof,b);
            end
            break
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define 3D Periodic Boundary Conditions for Faces                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Faces 1-3
%
for i=1:1:face1count
    for j=1:1:face3count
        if node(face1(i),2)==node(face3(j),2) && node(face1(i),3)==node(face3(j),3)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 3\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    face1(i),ndof,1,face3(j),ndof,-1,dummy1,ndof,-c);
            end
            break
        end
    end
end
%
% Faces 2-4
%
for i=1:1:face2count
    for j=1:1:face4count
        if node(face2(i),1)==node(face4(j),1) && node(face2(i),3)==node(face4(j),3)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 3\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    face2(i),ndof,1,face4(j),ndof,-1,dummy2,ndof,-a);
            end
            break
        end
    end
end
%
% Faces 6-5
%
for i=1:1:face6count
    for j=1:1:face5count
        if node(face6(i),1)==node(face5(j),1) && node(face6(i),2)==node(face5(j),2)
            for ndof=1:3 % Loop over the degrees of freedom
                fprintf(fid,'*EQUATION\n');
                fprintf(fid,' 3\n');
                fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                    face6(i),ndof,1,face5(j),ndof,-1,dummy3,ndof,-b);
            end
            break
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Repeat PBCs search and implementation for nodes on ncohxy               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cohesive_choice == 1
    %
    cedge1count=0; cedge2count=0; cedge3count=0; cedge4count=0; cedge5count=0;
    cedge6count=0; cedge7count=0; cedge8count=0; cedge9count=0; cedge10count=0;
    cedge11count=0; cedge12count=0;
    cface1count=0; cface2count=0; cface3count=0; cface4count=0; cface5count=0; cface6count=0;
    %
    for i=ncount+1:1:ncoh+ncount
        xx = round(1e8*ncohxy(i-ncount,1))/1e8;
        yy = round(1e8*ncohxy(i-ncount,2))/1e8;
        zz = round(1e8*ncohxy(i-ncount,3))/1e8;
        aa = round(1e8*a)/1e8;
        bb = round(1e8*b)/1e8;
        cc = round(1e8*c)/1e8;
        if xx==cc && yy==0 , cedge1count = cedge1count+1; cedge1(cedge1count)=i; continue; end
        if xx==cc && yy==aa, cedge2count = cedge2count+1; cedge2(cedge2count)=i; continue; end
        if xx==0  && yy==aa, cedge3count = cedge3count+1; cedge3(cedge3count)=i; continue; end
        if xx==0  && yy==0 , cedge4count = cedge4count+1; cedge4(cedge4count)=i; continue; end
        if xx==cc && zz==0 , cedge5count = cedge5count+1; cedge5(cedge5count)=i; continue; end
        if xx==cc && zz==bb, cedge6count = cedge6count+1; cedge6(cedge6count)=i; continue; end
        if xx==0  && zz==bb, cedge7count = cedge7count+1; cedge7(cedge7count)=i; continue; end
        if xx==0  && zz==0 , cedge8count = cedge8count+1; cedge8(cedge8count)=i; continue; end
        % added my MAB
        if yy==aa && zz==0 , cedge10count = cedge10count+1; cedge10(cedge10count)=i; continue; end
        if yy==aa && zz==bb, cedge11count = cedge11count+1; cedge11(cedge11count)=i; continue; end
        if yy==0  && zz==bb, cedge12count = cedge12count+1; cedge12(cedge12count)=i; continue; end
        % end of ammendment
        if xx == cc, cface1count = cface1count + 1; cface1(cface1count) = i; continue; end
        if yy == aa, cface2count = cface2count + 1; cface2(cface2count) = i; continue; end
        if xx == 0 , cface3count = cface3count + 1; cface3(cface3count) = i; continue; end
        if yy == 0 , cface4count = cface4count + 1; cface4(cface4count) = i; continue; end
        if zz == 0 , cface5count = cface5count + 1; cface5(cface5count) = i; continue; end
        if zz == bb, cface6count = cface6count + 1; cface6(cface6count) = i; continue; end
    end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if  were properly applied                                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    if cedge1count ~= cedge3count, disp('PBC Application failed. cedge 1-3');...
            status_mesh = 0; return, end
    if cedge2count ~= cedge4count, disp('PBC Application failed. cedge 2-4');...
            status_mesh = 0; return, end
    if cedge6count ~= cedge8count, disp('PBC Application failed. cedge 6-8');...
            status_mesh = 0; return, end
    if cedge5count ~= cedge7count, disp('PBC Application failed. cedge 5-7');...
            status_mesh = 0; return, end
    % added by MAB
    if cedge11count ~= cedge9count, disp('PBC Application failed. cedge 11-9');...
        status_mesh = 0; return, end
    if cedge10count ~= cedge12count, disp('PBC Application failed. cedge 10-12');...
            status_mesh = 0; return, end
    % end of ammendment
    if cface2count~=cface4count, disp('PBC Application failed. cface 2-4');...
            status_mesh = 0; return, end
    if cface1count~=cface3count, disp('PBC Application failed. cface 1-3');...
            status_mesh = 0; return, end
    if cface5count~=cface6count, disp('PBC Application failed. cface 5-6');...
            status_mesh = 0; return, end
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define 3D Periodic Boundary Conditions for cedges                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % cedges 2-4
    %
    for i=1:1:cedge2count
        for j=1:1:cedge4count
            if ncohxy(cedge2(i)-ncount,3) == ncohxy(cedge4(j)-ncount,3)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 4\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cedge2(i),ndof,1,cedge4(j),ndof,-1,dummy1,ndof,-c,dummy2,ndof,-a);
                end
                break
            end
        end
    end
    %
    % cedges 1-3
    %
    
    for i=1:1:cedge1count
        for j=1:1:cedge3count
            if ncohxy(cedge1(i)-ncount,3) == ncohxy(cedge3(j)-ncount,3)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 4\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cedge1(i),ndof,1,cedge3(j),ndof,-1,dummy1,ndof,-c,dummy2,ndof,a);
                end
                break
            end
        end
    end
    %
    % cedges 6-8
    %
    
    for i=1:1:cedge6count
        for j=1:1:cedge8count
            if ncohxy(cedge6(i)-ncount,2) == ncohxy(cedge8(j)-ncount,2)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 4\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cedge6(i),ndof,1,cedge8(j),ndof,-1,dummy1,ndof,-c,dummy3,ndof,-b);
                end
                break
            end
        end
    end
    %
    % cedges 5-7
    %
    for i=1:1:cedge5count
        for j=1:1:cedge7count
            if ncohxy(cedge5(i)-ncount,2) == ncohxy(cedge7(j)-ncount,2)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 4\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cedge5(i),ndof,1,cedge7(j),ndof,-1,dummy1,ndof,-c,dummy3,ndof,b);
                end
                break
            end
        end
    end
    %
    % added by MAB:
    %
    % cedges 11-9
    %
    for i=1:1:cedge11count
        for j=1:1:cedge9count
            if ncohxy(cedge11(i)-ncount,1) == ncohxy(cedge9(j)-ncount,1)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 4\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cedge11(i),ndof,1,cedge9(j),ndof,-1,dummy2,ndof,-a,dummy3,ndof,-b);
                end
                break
            end
        end
    end
    %
    % cedges 10-12
    %
    for i=1:1:cedge10count
        for j=1:1:cedge12count
            if ncohxy(cedge10(i)-ncount,1) == ncohxy(cedge12(j)-ncount,1)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 4\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cedge10(i),ndof,1,cedge12(j),ndof,-1,dummy2,ndof,-a,dummy3,ndof,b);
                end
                break
            end
        end
    end
    % end of ammendment
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define 3D Periodic Boundary Conditions for Faces                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % cfaces 1-3
    %
    for i=1:1:cface1count
        for j=1:1:cface3count
            if ncohxy(cface1(i)-ncount,2)==ncohxy(cface3(j)-ncount,2) &&...
                        ncohxy(cface1(i)-ncount,3)==ncohxy(cface3(j)-ncount,3)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 3\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cface1(i),ndof,1,cface3(j),ndof,-1,dummy1,ndof,-c);
                end
                break
            end
        end
    end
    %
    % cfaces 2-4
    %
    for i=1:1:cface2count
        for j=1:1:cface4count
            if ncohxy(cface2(i)-ncount,1)==ncohxy(cface4(j)-ncount,1) &&...
                        ncohxy(cface2(i)-ncount,3)==ncohxy(cface4(j)-ncount,3)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 3\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cface2(i),ndof,1,cface4(j),ndof,-1,dummy2,ndof,-a);
                end
                break
            end
        end
    end
    %
    % cfaces 6-5
    %
    for i=1:1:cface6count
        for j=1:1:cface5count
            if ncohxy(cface6(i)-ncount,1)==ncohxy(cface5(j)-ncount,1) &&...
                        ncohxy(cface6(i)-ncount,2)==ncohxy(cface5(j)-ncount,2)
                for ndof=1:3 % Loop over the degrees of freedom
                    fprintf(fid,'*EQUATION\n');
                    fprintf(fid,' 3\n');
                    fprintf(fid,' %9d, %1d, %9.8f, %9d, %1d, %9.8f, %9d, %1d, %9.8f\n',...
                        cface6(i),ndof,1,cface5(j),ndof,-1,dummy3,ndof,-b);
                end
                break
            end
        end
    end
    %
    %
end % end if cohesive_choice == 1
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate file with material properties for elastic analysis             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(dir_name,'/matprops_e.inp');
fid = fopen(file_name,'w');
%
% Properties
%
fprintf(fid,'** Property 1 : Matrix Property\n');
fprintf(fid,'*SOLID SECTION, ELSET=MATRIX, MATERIAL=MATR, CONTROLS=SECT1\n');
fprintf(fid,' 1\n');
if explicit_option == 0 % If running implicit
    fprintf(fid,'*HOURGLASS STIFFNESS\n');
    fprintf(fid,' %6.4E\n',0.005*Em/2/(1+vm));
end
fprintf(fid,'** Property 2 : Fibre Property\n');
fprintf(fid,'*ORIENTATION, NAME=ORIENT_FIBRE\n');
fprintf(fid,' 1, 0, 0, 0, 1, 0\n');
fprintf(fid,' 3, 0.\n');
fprintf(fid,'*SOLID SECTION, ELSET=FIBRE, MATERIAL=FIBR, ORIENTATION=ORIENT_FIBRE, CONTROLS=SECT2\n');
if explicit_option == 0 % If running implicit
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT1, Hourglass=STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT2, Hourglass=STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
else % if running explicit
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT1, Hourglass=RELAX STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT2, Hourglass=RELAX STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
end

%
% Materials
%
fprintf(fid,'**\n');
fprintf(fid,'** Material 1 : Matrix\n');
fprintf(fid,'*MATERIAL, NAME=MATR\n');
fprintf(fid,'*ELASTIC\n');
fprintf(fid,' %6g, %6.4f, 0.\n',Em,vm);
fprintf(fid,'*EXPANSION, ZERO=%4g\n',T0);
fprintf(fid,' %6.5E, 0.\n',alpham);
fprintf(fid,' %6.5E, 200.\n',alpham);
fprintf(fid,'*DENSITY\n');
fprintf(fid,' %6.5E\n',Densm);
fprintf(fid,'**\n');
fprintf(fid,'** Material 2 : Fibre\n');
fprintf(fid,'*MATERIAL, NAME=FIBR\n');
fprintf(fid,'*ELASTIC, TYPE=ENGINEERING CONSTANTS\n');
fprintf(fid,' %7g, %7g, %7g, %6.4f, %6.4f, %6.4f, %7g, %7g,\n',Ef1,Ef2,Ef2,vf12,vf12,vf23,Gf12,Gf12);
fprintf(fid,' %7g , 0.\n',Gf23);
fprintf(fid,'*EXPANSION, TYPE=ORTHO, ZERO=%4g\n',T0);
fprintf(fid,' %6.5E, %6.5E, %6.5E, 0.\n',alphaf1,alphaf2,alphaf2);
fprintf(fid,' %6.5E, %6.5E, %6.5E, 200.\n',alphaf1,alphaf2,alphaf2);
fprintf(fid,'*DENSITY\n');
fprintf(fid,' %6.5E\n',Densf);
fprintf(fid,'**\n');
fprintf(fid,'**Initial Conditions\n');
fprintf(fid,'*INITIAL CONDITIONS, TYPE=TEMPERATURE\n');
fprintf(fid,'GLOBAL, %3g\n',T0);
if cohesive_choice == 1
    fprintf(fid,'**\n');
    fprintf(fid,'** Material 3 : Interface\n');
    fprintf(fid,'*COHESIVE SECTION, RESPONSE=TRACTION SEPARATION, MATERIAL=COHMAT, ELSET=COHESI\n');
    fprintf(fid,' 1.000\n');
    fprintf(fid,'*MATERIAL, NAME=COHMAT\n');
    fprintf(fid,'*DENSITY\n');
    fprintf(fid,' %6.5E\n',Densm);
    fprintf(fid,'*ELASTIC, TYPE=TRACTION\n');
    fprintf(fid,' 1E8, 1E8, 1E8\n');
end
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate file with material properties for elasto-plastic with damage   %
% analyses                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(dir_name,'/matprops_p.inp');
fid = fopen(file_name,'w');
%
% Properties
%

fprintf(fid,'** Property 1 : Matrix Property\n');
fprintf(fid,'*SOLID SECTION, ELSET=MATRIX, MATERIAL=MATR, CONTROLS=SECT1\n');
fprintf(fid,' 1\n');
if explicit_option == 0 % If running implicit
    fprintf(fid,'*HOURGLASS STIFFNESS\n');
    fprintf(fid,' %6.4E\n',0.005*Em/2/(1+vm));
end
fprintf(fid,'** Property 2 : Fibre Property\n');
fprintf(fid,'*ORIENTATION, NAME=ORIENT_FIBRE\n');
fprintf(fid,' 1, 0, 0, 0, 1, 0\n');
fprintf(fid,' 3, 0.\n');
fprintf(fid,'*SOLID SECTION, ELSET=FIBRE, MATERIAL=FIBR, ORIENTATION=ORIENT_FIBRE, CONTROLS=SECT2\n');
if explicit_option == 0 % If running implicit
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT1, Hourglass=STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT2, Hourglass=STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
    fprintf(fid,'*HOURGLASS STIFFNESS\n');
    fprintf(fid,' %6.4E\n',0.005*Gf12);
else % if running explicit
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT1, Hourglass=RELAX STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
    fprintf(fid,'*SECTION CONTROLS, NAME=SECT2, Hourglass=RELAX STIFFNESS\n');
    fprintf(fid,'1., 1., 1.\n');
end
%
% Materials
%
fprintf(fid,'**\n');
% Matrix properties
fprintf(fid,'** Material 1 : Matrix\n');
fprintf(fid,'*MATERIAL, NAME=MATR\n');
fprintf(fid,'*DEPVAR\n');
fprintf(fid,' 22\n');
fprintf(fid,'*DENSITY\n');
fprintf(fid,' %6.5E\n',Densm);
fprintf(fid,'*USER MATERIAL, CONSTANTS=9\n');
fprintf(fid,'** E, NIU, ALPHA, NIUp, XT, XC, Gc, EPRref\n');
fprintf(fid,' %6g, %6.4f, %6.5e, %6.4f, %6g, %6g, %6.4f, %6.4f\n',...
    Em,vm,alpham,niup,Xtm,Xcm,GcTm,EPRref);
fprintf(fid,'** Tref, ALPHA_C, ALPHA_T, ALPHA_S, BETA_C, BETA_T, BETA_S\n');
fprintf(fid,' %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f, %6.4f\n',...
    Tref,alpha_c,alpha_t,alpha_s,beta_c,beta_t,beta_s);
fprintf(fid,'**\n');
%
% Fibre properties
fprintf(fid,'** Material 2 : Fibre\n');
fprintf(fid,'*MATERIAL, NAME=FIBR\n');
fprintf(fid,'*ELASTIC, TYPE=ENGINEERING CONSTANTS\n');
fprintf(fid,' %7g, %7g, %7g, %6.4f, %6.4f, %6.4f, %7g, %7g,\n',Ef1,Ef2,Ef2,vf12,vf12,vf23,Gf12,Gf12);
fprintf(fid,' %7g , 0.\n',Gf23);
fprintf(fid,'*EXPANSION, TYPE=ORTHO, ZERO=%4g\n',T0);
fprintf(fid,' %6.5E, %6.5E, %6.5E, 0.\n',alphaf1,alphaf2,alphaf2);
fprintf(fid,' %6.5E, %6.5E, %6.5E, 200.\n',alphaf1,alphaf2,alphaf2);
fprintf(fid,'*DENSITY\n');
fprintf(fid,' %6.5E\n',Densm);
fprintf(fid,'**\n');
%
% Initial conditions (Temperature)
fprintf(fid,'**Initial Conditions\n');
fprintf(fid,'*INITIAL CONDITIONS, TYPE=TEMPERATURE\n');
fprintf(fid,'GLOBAL, %3g\n',T0);
fprintf(fid,'**\n');
%
if cohesive_choice == 1
    % Interface properties (cohesive elements between fibres and matrix)
    fprintf(fid,'** Material 3 : Interface\n');
    fprintf(fid,'*COHESIVE SECTION, RESPONSE=TRACTION SEPARATION, MATERIAL=COHMAT, ELSET=COHESI\n');
    fprintf(fid,' 1.000\n');
    fprintf(fid,'*MATERIAL, NAME=COHMAT\n');
    fprintf(fid,'*DENSITY\n');
    fprintf(fid,' %6.5E\n',Densf);
    fprintf(fid,'*ELASTIC, TYPE=TRACTION\n');
    fprintf(fid,' 1E8, 1E8, 1E8\n');
    fprintf(fid,'*DAMAGE INITIATION, CRITERION=MAXS\n');
    fprintf(fid,' 50e3, 75e3, 75e3\n');
    fprintf(fid,'*DAMAGE EVOLUTION, TYPE=ENERGY, MIXED MODE BEHAVIOR=BK, MODE MIX RATIO=ENERGY, POWER=1.45, SOFTENING=LINEAR\n');
    fprintf(fid,' 0.002, 0.006, 0.006\n');
    fprintf(fid,'*DAMAGE STABILIZATION\n');
    fprintf(fid,' 0.005\n');
end
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate file with node sets                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(dir_name,'/spnsets.inp');
fid = fopen(file_name,'w');
%
% Vertices
%
fprintf(fid,'*NSET, NSET=VERTICES\n');
for i=1:1:7
    fprintf(fid,' %9d,',vert(i));
end
fprintf(fid,' %9d\n',vert(8));
%
% Individual Vertices
%
for i=1:1:8
    vertex_string = strcat('VERTEX0',int2str(i));
    fprintf(fid,'*NSET, NSET=%s\n',vertex_string);
    fprintf(fid,' %9d\n',vert(i));
end
%
% Edges
%
fprintf(fid,'*NSET, NSET=EDGE01\n');
count = 0;
for i=1:1:edge1count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge1(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE02\n');
count = 0;
for i=1:1:edge2count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge2(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE03\n');
count = 0;
for i=1:1:edge3count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge3(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE04\n');
count = 0;
for i=1:1:edge4count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge4(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE05\n');
count = 0;
for i=1:1:edge5count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge5(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE06\n');
count = 0;
for i=1:1:edge6count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge6(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE07\n');
count = 0;
for i=1:1:edge7count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge7(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE08\n');
count = 0;
for i=1:1:edge8count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge8(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE09\n');
count = 0;
for i=1:1:edge9count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge9(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE10\n');
count = 0;
for i=1:1:edge10count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge10(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE11\n');
count = 0;
for i=1:1:edge11count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge11(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=EDGE12\n');
count = 0;
for i=1:1:edge12count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',edge12(i));
end
fprintf(fid,'\n');
%
% Faces
%
fprintf(fid,'*NSET, NSET=FACE01\n');
count = 0;
for i=1:1:face1count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',face1(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=FACE02\n');
count = 0;
for i=1:1:face2count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',face2(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=FACE03\n');
count = 0;
for i=1:1:face3count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',face3(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=FACE04\n');
count = 0;
for i=1:1:face4count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',face4(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=FACE05\n');
count = 0;
for i=1:1:face5count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',face5(i));
end
fprintf(fid,'\n');
%
fprintf(fid,'*NSET, NSET=FACE06\n');
count = 0;
for i=1:1:face6count
    count = count + 1;
    if count == 17
        fprintf(fid,'\n');
        count = 1;
    end
    if count > 1
        fprintf(fid,',');
    end
    fprintf(fid,' %6d',face6(i));
end
fprintf(fid,'\n');
%
% EVERY node in each FACE (including vertices and edges)
fprintf(fid,'*NSET, NSET=ALLnodes_Face01\n');
fprintf(fid,'Face01, EDGE01, EDGE02, EDGE05, EDGE06, VERTEX01, VERTEX02, VERTEX03, VERTEX04\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_Face03\n');
fprintf(fid,'Face03, EDGE04, EDGE03, EDGE08, EDGE07, VERTEX05, VERTEX06, VERTEX07, VERTEX08\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_Face02\n');
fprintf(fid,'Face02, EDGE02, EDGE03, EDGE10, EDGE11, VERTEX02, VERTEX03, VERTEX06, VERTEX07\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_Face04\n');
fprintf(fid,'Face04, EDGE01, EDGE04, EDGE09, EDGE12, VERTEX01, VERTEX04, VERTEX05, VERTEX08\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_Face05\n');
fprintf(fid,'Face05, EDGE05, EDGE08, EDGE09, EDGE10, VERTEX01, VERTEX02, VERTEX05, VERTEX06\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_Face06\n');
fprintf(fid,'Face06, EDGE06, EDGE07, EDGE12, EDGE11, VERTEX04, VERTEX03, VERTEX08, VERTEX07\n');
%
% EVERY node in each EDGE (including vertices)
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE01\n');
fprintf(fid,'EDGE01, VERTEX01, VERTEX04\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE02\n');
fprintf(fid,'EDGE02, VERTEX02, VERTEX03\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE03\n');
fprintf(fid,'EDGE03, VERTEX06, VERTEX07\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE04\n');
fprintf(fid,'EDGE04, VERTEX05, VERTEX08\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE05\n');
fprintf(fid,'EDGE05, VERTEX01, VERTEX02\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE06\n');
fprintf(fid,'EDGE06, VERTEX04, VERTEX03\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE07\n');
fprintf(fid,'EDGE07, VERTEX08, VERTEX07\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE08\n');
fprintf(fid,'EDGE08, VERTEX05, VERTEX06\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE09\n');
fprintf(fid,'EDGE09, VERTEX05, VERTEX01\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE10\n');
fprintf(fid,'EDGE10, VERTEX06, VERTEX02\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE11\n');
fprintf(fid,'EDGE11, VERTEX07, VERTEX03\n');
%
fprintf(fid,'*NSET, NSET=ALLnodes_EDGE12\n');
fprintf(fid,'EDGE12, VERTEX08, VERTEX04\n');
%
fprintf(fid,'**\n');
%
% Finally, define surfaces that are needed for PSEUDO-Periodic Boundary Condtions
fprintf(fid,'*** Define surfaces from node sets:\n');
fprintf(fid,'*Surface, type=NODE, name=Entire_FACE06\n');
fprintf(fid,'ALLnodes_FACE06, 1.\n');
%
fprintf(fid,'*Surface, type=NODE, name=Entire_FACE05\n');
fprintf(fid,'ALLnodes_FACE05, 1.\n');
%
fprintf(fid,'*Surface, type=NODE, name=Entire_FACE04\n');
fprintf(fid,'ALLnodes_FACE04, 1.\n');
%
fprintf(fid,'*Surface, type=NODE, name=Entire_FACE03\n');
fprintf(fid,'ALLnodes_FACE03, 1.\n');
%
fprintf(fid,'*Surface, type=NODE, name=Entire_FACE02\n');
fprintf(fid,'ALLnodes_FACE02, 1.\n');
%
fprintf(fid,'*Surface, type=NODE, name=Entire_FACE01\n');
fprintf(fid,'ALLnodes_FACE01, 1.\n');
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Hardening Curves From Text File                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(cd,'/Input/',hardening_curves,'.txt');
fid = fopen(file_name,'r');
%
for i=1:1:19, fgetl(fid); end
% Compression Hardening curve
NCC = str2num(fgetl(fid)); fgetl(fid);
YCC = zeros(NCC,2);
% Get equivalent plastic strain:
for i=1:1:NCC
    YCC(i,1)=str2num(fgetl(fid));
end
fgetl(fid); fgetl(fid);
% Get stresses:
for i=1:1:NCC
    YCC(i,2)=str2num(fgetl(fid));
end
%
%
% Tension Hardening Curve
for i=1:1:4, fgetl(fid); end
NCT = str2num(fgetl(fid)); fgetl(fid);
YCT = zeros(NCT,2);
% Get equivalent plastic strain:
for i=1:1:NCT
    YCT(i,1)=str2num(fgetl(fid));
end
fgetl(fid); fgetl(fid);
% Get stresses:
for i=1:1:NCT
    YCT(i,2)=str2num(fgetl(fid));
end
%
%
% Tension Hardening Curve
for i=1:1:4, fgetl(fid); end
NCS = str2num(fgetl(fid)); fgetl(fid);
YCS = zeros(NCS,2);
% Get equivalent plastic strain:
for i=1:1:NCS
    YCS(i,1)=str2num(fgetl(fid));
end
fgetl(fid); fgetl(fid);
% Get stresses:
for i=1:1:NCS
    YCS(i,2)=str2num(fgetl(fid));
end
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write Hardening Curves to a file to be included in the USER MATERIAL    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
file_name = strcat(dir_name,'/hardening_laws_for_vumat.f');
fid = fopen(file_name,'w');
%
fprintf(fid,'C Hardening Curves for VUMAT\n');
fprintf(fid,'C ==================================================================\n');
fprintf(fid,'C\n');
fprintf(fid,'      PARAMETER (NCC=%d,NCT=%d,NCS=%d)\n',round(NCC),round(NCT),round(NCS));
fprintf(fid,'      DIMENSION YCC(NCC,2),YCT(NCT,2),YCS(NCS,2)\n');
%
% Write Compression hardening curve
fprintf(fid,'C\n');
fprintf(fid,'C     Compression hardening curve\n');
fprintf(fid,'      data (YCC(i,1), i=1,NCC) /\n');
for i=1:1:(NCC-1)
    fprintf(fid,'     *  %6.9f,\n',YCC(i,1));
end
fprintf(fid,'     *  %6.9f /\n',YCC(NCC,1));
%
fprintf(fid,'      data (YCC(i,2), i=1,NCC) /\n');
for i=1:1:(NCC-1)
    fprintf(fid,'     *  %6.9f,\n',YCC(i,2));
end
fprintf(fid,'     *  %6.9f /\n',YCC(NCC,2));
%
% Write Tension hardening curve
fprintf(fid,'C\n');
fprintf(fid,'C     Tension hardening curve\n');
fprintf(fid,'      data (YCT(i,1), i=1,NCT) /\n');
for i=1:1:(NCT-1)
    fprintf(fid,'     *  %6.9f,\n',YCT(i,1));
end
fprintf(fid,'     *  %6.9f /\n',YCT(NCT,1));
%
fprintf(fid,'      data (YCT(i,2), i=1,NCT) /\n');
for i=1:1:(NCT-1)
    fprintf(fid,'     *  %6.9f,\n',YCT(i,2));
end
fprintf(fid,'     *  %6.9f /\n',YCT(NCT,2));
%
%
% Write Shear hardening curve
fprintf(fid,'C\n');
fprintf(fid,'C     Shear hardening curve\n');
fprintf(fid,'      data (YCS(i,1), i=1,NCS) /\n');
for i=1:1:(NCS-1)
    fprintf(fid,'     *  %6.9f,\n',YCS(i,1));
end
fprintf(fid,'     *  %6.9f /\n',YCS(NCS,1));
%
fprintf(fid,'      data (YCS(i,2), i=1,NCS) /\n');
for i=1:1:(NCS-1)
    fprintf(fid,'     *  %6.9f,\n',YCS(i,2));
end
fprintf(fid,'     *  %6.9f /\n',YCS(NCS,2));
%
fprintf(fid,'C\n');
%
fclose(fid);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display Miscellaneous Information                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
disp(' ');
disp('Generation of .inp files COMPLETED');
disp('Elapsed Time [min]: ');
disp(toc/60);
status_mesh = 1;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
