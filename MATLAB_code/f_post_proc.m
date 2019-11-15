%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Post-Process FE Analyses                              %
%                                                                         %
%  This file is part of the FEA_Automatic.m sequence of files             %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v4.0                                   %
%                                                                         %
%                Antonio Rui Melro - antonio.melro@fe.up.pt               %
%                             February 2010                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function f_post_proc(step_ID)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global dir_fea_name s11 s22 s33 s12 s13 s23 a b c;
%
disp(' ');
disp('Started post-processing results');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-Process Results for 1 <= StepID <= 6                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if step_ID >= 1 && step_ID <=6
    if step_ID == 1
        disp('Post-Processing Results to determine 1st Column of Constitutive Tensor');
    elseif step_ID == 2
        disp('Post-Processing Results to determine 2nd Column of Constitutive Tensor');
    elseif step_ID == 3
        disp('Post-Processing Results to determine 3rd Column of Constitutive Tensor');
    elseif step_ID == 4
        disp('Post-Processing Results to determine 4th Column of Constitutive Tensor');
    elseif step_ID == 5
        disp('Post-Processing Results to determine 5th Column of Constitutive Tensor');
    elseif step_ID == 6
        disp('Post-Processing Results to determine 6th Column of Constitutive Tensor');
    end
    disp(' ');
    %
    file_name = strcat(dir_fea_name,'/post_proc.py');
    fid = fopen(file_name,'wt');
    %
    % Heading
    %
    fprintf(fid,'#\n');
    fprintf(fid,'# Abaqus/Viewer Version 6.7-5 script\n');
    fprintf(fid,'# Internal Version: 2008_01_21-15.23.13 80075\n');
    fprintf(fid,'# Created by Antonio Rui Melro on %s\n',datestr(now));
    fprintf(fid,'#\n');
    fprintf(fid,'\n');
    fprintf(fid,'# from driverUtils import executeOnCaeGraphicsStartup\n');
    fprintf(fid,'# executeOnCaeGraphicsStartup()\n');
    fprintf(fid,'#: Executing "onCaeGraphicsStartup()" in the site directory ...\n');
    fprintf(fid,'from abaqus import *\n');
    fprintf(fid,'from abaqusConstants import *\n');
    fprintf(fid,'session.Viewport(name=''Viewport: 1'', origin=(0.0, 0.0), width=176.1, \n');
    fprintf(fid,'    height=191.2)\n');
    fprintf(fid,'session.viewports[''Viewport: 1''].makeCurrent()\n');
    fprintf(fid,'session.viewports[''Viewport: 1''].maximize()\n');
    fprintf(fid,'from viewerModules import *\n');
    fprintf(fid,'from driverUtils import executeOnCaeStartup\n');
    fprintf(fid,'executeOnCaeStartup()\n');
    %
    % Open results file and print stress totals
    %
    fprintf(fid,'o1 = session.openOdb(name=''%s/run.odb'')\n',dir_fea_name);
    fprintf(fid,'session.viewports[''Viewport: 1''].setValues(displayedObject=o1)\n');
    fprintf(fid,'odb = session.odbs[''%s/run.odb'']\n',dir_fea_name);
    fprintf(fid,'session.fieldReportOptions.setValues(printXYData=ON, printMinMax=OFF, printTotal=OFF)\n');
    fprintf(fid,'session.writeFieldReport(\n');
    fprintf(fid,'    fileName=''%s/stresses.rpt'', append=OFF, \n',dir_fea_name);
    fprintf(fid,'    sortItem=''Element Label'', odb=odb, step=0, frame=1, \n');
    fprintf(fid,'    outputPosition=INTEGRATION_POINT, variable=((''IVOL'', INTEGRATION_POINT), (\n');
    fprintf(fid,'    ''S'', INTEGRATION_POINT, ((COMPONENT, ''S11''), (COMPONENT, ''S22''), (\n');
    fprintf(fid,'    COMPONENT, ''S33''), (COMPONENT, ''S12''), (COMPONENT, ''S13''), (COMPONENT, \n');
    fprintf(fid,'    ''S23''), )), ))\n');
    fprintf(fid,'session.odbs[''%s/run.odb''].close()\n',dir_fea_name);
    %
    fclose(fid);
    %
    % Runs the Post-processing utility
    %
    tmp=pwd;
    cd(dir_fea_name);
    unix(strcat('/home/mab031/abaqus/Commands/abq6112 vi noGUI=post_proc.py'));
    cd(tmp);
    %
    % Reads from the generated report file
    %
    file_name = strcat(dir_fea_name,'/stresses.rpt');
    fid = fopen(file_name,'r');
    %
    for i=1:1:19
        fgetl(fid);
    end
    %
    s11 = 0; s22 = 0; s33 = 0; s12 = 0; s13 = 0; s23 = 0;
    for i=1:1:2
        status = 1;
        while status == 1
            tline=fgetl(fid);
            [data status] = str2num(tline);
            if status == 0, break, end
            s11 = s11 + data(1,3)*data(1,4);
            s22 = s22 + data(1,3)*data(1,5);
            s33 = s33 + data(1,3)*data(1,6);
            s12 = s12 + data(1,3)*data(1,7);
            s13 = s13 + data(1,3)*data(1,8);
            s23 = s23 + data(1,3)*data(1,9);
        end
        if i == 2, break, end
        for j=1:1:6
            fgetl(fid);
        end
    end
    %
    vol = a*b*c;
    s11 = s11/vol; s22 = s22/vol; s33 = s33/vol;
    s12 = s12/vol; s13 = s13/vol; s23 = s23/vol;
    %
    fclose(fid);
end
%
disp(' ');
disp('Post-Processing COMPLETED');
disp('Elapsed Time [min]: ');
disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
