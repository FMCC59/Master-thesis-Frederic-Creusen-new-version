%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Generate and Run FE Analyses                          %
%                                                                         %
%  This file is part of the FEA_Automatic.m sequence of files             %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v4.0                                   %
%                                                                         %
%                Antonio Rui Melro - antonio.melro@fe.up.pt               %
%                              October 2010                               %
%                                                                         %
%                Miguel Anibal Bessa - mbessa@u.northwestern.edu          %
%                              August 2015                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   LIST OF STEP IDENTIFIERS                                              %
%                                                                         %
%        11 --> Positive Longitudinal Normal Load                         %
%        12 --> Positive Transverse Normal Load (22)                      %
%        13 --> Positive Transverse Normal Load (33)                      %
%        14 --> Positive Longitudinal Shear Load (12)                     %
%        15 --> Positive Longitudinal Shear Load (13)                     %
%        16 --> Positive Transverse Shear Load                            %
%        21 --> Negative Longitudinal Normal Load                         %
%        22 --> Negative Transverse Normal Load (22)                      %
%        23 --> Negative Transverse Normal Load (33)                      %
%        24 --> Negative Longitudinal Shear Load (12)                     %
%        25 --> Negative Longitudinal Shear Load (13)                     %
%        26 --> Negative Transverse Shear Load                            %
%       236 --> Transverse Compression + Transverse Shear                 %
%       123 --> Biaxial Stress State S22 + S33                            %
%       124 --> Transv. Tensile + Long. Shear;                            %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function f_fea_3D_explicit(stepID)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global dir_name dummy1 dummy2 dummy3 dir_fea_name abaqus_path ...
    cohesive_choice explicit_option plastic_option pbcs_option a b c;
%
disp(' '); disp(' ');
disp('Started writing the main input file: run.inp');
%
if strcmp(stepID, 'p1')
    disp('Analysis with a Positive Longitudinal Normal Load');
    disp('LoadID_p1');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p1');
    %
elseif strcmp(stepID, 'p2')
    disp('Analysis with a Positive Transverse Normal Load (22)');
    disp('LoadID_p2');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p2');
    %
elseif strcmp(stepID, 'p3')
    disp('Analysis with a Positive Transverse Normal Load (33)');
    disp('LoadID_p3');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p3');
    %
elseif strcmp(stepID, 'p4')
    disp('Analysis with a Positive Longitudinal Shear Load (12)');
    disp('LoadID_p4');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p4');
    %
elseif strcmp(stepID, 'p5')
    disp('Analysis with a Positive Longitudinal Shear Load (13)');
    disp('LoadID_p5');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p5');
    %
elseif strcmp(stepID, 'p6')
    disp('Analysis with a Positive Transverse Shear Load');
    disp('LoadID_p6');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p6');
    %
elseif strcmp(stepID, 'n1')
    disp('Analysis with a Negative Longitudinal Normal Load');
    disp('LoadID_n1');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n1');
    %
elseif strcmp(stepID, 'n2')
    disp('Analysis with a Negative Transverse Normal Load (22)');
    disp('LoadID_n2');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n2');
    %
elseif strcmp(stepID, 'n3')
    disp('Analysis with a Negative Transverse Normal Load (33)');
    disp('LoadID_n3');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n3');
    %
elseif strcmp(stepID, 'n4')
    disp('Analysis with a Negative Longitudinal Shear Load (12)');
    disp('LoadID_n4');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n4');
    %
elseif strcmp(stepID, 'n5')
    disp('Analysis with a Negative Longitudinal Shear Load (13)');
    disp('LoadID_n5');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n5');
    %
elseif strcmp(stepID, 'n6')
    disp('Analysis with a Negative Transverse Shear Load');
    disp('LoadID_n6');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n6');
    %
elseif strcmp(stepID, 'n3p6')
    disp('Analysis with a Transverse Shear + Compression Load');
    disp('LoadID_n3p6');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_n3p6');
    %
elseif strcmp(stepID, 'p2p3')
    disp('Analysis with Biaxial Stress State S22 + S33');
    disp('LoadID_p2p3');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p2p3');
    %
elseif strcmp(stepID, 'p2p4')
    disp('Analysis with Transverse Tension + Longitudinal Shear');
    disp('LoadID_p2p4');
    disp(' ');
    %
    dir_fea_name = strcat(dir_name,'/LoadID_p2p4');
    %
end
%
%% Write lines that are common to every load case:
%
if exist(dir_fea_name,'dir') == 0, mkdir(dir_fea_name); end
%
if explicit_option == 1
    file_name = strcat(dir_fea_name,'/run_explicit.inp');
else
    file_name = strcat(dir_fea_name,'/run_implicit.inp');
end
fid = fopen(file_name,'wt');
%
% Heading
%
fprintf(fid,'** ***************************************************************************\n');
fprintf(fid,'**   Created by A.R. Melro and M.A. Bessa on %s\n',datestr(now));
fprintf(fid,'**\n');
if explicit_option == 1
    fprintf(fid,'** INPUT file for EXPLICIT analysis\n');
else
    fprintf(fid,'** INPUT file for IMPLICIT analysis\n');
end
fprintf(fid,'**\n');
fprintf(fid,'** ***************************************************************************\n');
fprintf(fid,'**\n');
%
% Including other files
%
fprintf(fid,'** Include file with RVE mesh:\n');
fprintf(fid,'*INCLUDE, INPUT=../mesh.inp\n');
fprintf(fid,'** Amplitude for the load curve:\n');
fprintf(fid,'*Amplitude, name=AMP-1, DEFINITION=SMOOTH STEP\n');
fprintf(fid,'0.,   0.,     1.0,      1.0\n');
fprintf(fid,'** Include file with ALL material properties:\n');
if plastic_option == 1 % Running plasticity analysis
    fprintf(fid,'*INCLUDE, INPUT=../matprops_p.inp\n');
else % Running elastic analysis
    fprintf(fid,'*INCLUDE, INPUT=../matprops_e.inp\n');
end
fprintf(fid,'** Include file with node sets:\n');
fprintf(fid,'*INCLUDE, INPUT=../spnsets.inp\n');
if pbcs_option == 1
    fprintf(fid,'** Include file defining constraints for Periodic Boundary Conditions:\n');
    fprintf(fid,'*INCLUDE, INPUT=../pbcs.inp\n');
else % NO Periodic Boundary Conditions
    fprintf(fid,'**** COMMENTED OUT the Periodic Boundary Conditions:\n');
    fprintf(fid,'***INCLUDE, INPUT=../pbcs.inp\n');
end
%
if pbcs_option == 0 % Using PSEUDO-periodic boundary conditions
    %
    if strcmp(stepID, 'p1')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p2')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'***Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'**Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p3')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'***Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'**Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p4')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p5')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p6')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'***Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'**Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n1')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n2')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'***Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'**Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n3')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'***Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'**Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n4')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n5')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n6')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'***Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'**Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'n3p6')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'*Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'***Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'**Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p2p3')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'*Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'***Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'**Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'***Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'**Entire_FACE06, Entire_FACE05\n');
        %
    elseif strcmp(stepID, 'p2p4')
        %
        fprintf(fid,'** USE PSEUDO-periodic boundary conditions in faces that are not loaded:\n');
        fprintf(fid,'***Tie, name=Tie_Face01Face03, adjust=no, position tolerance=%6.6f\n',2*c);
        fprintf(fid,'**Entire_FACE01, Entire_FACE03\n');
        fprintf(fid,'***Tie, name=Tie_Face02Face04, adjust=no, position tolerance=%6.6f\n',2*a);
        fprintf(fid,'**Entire_FACE02, Entire_FACE04\n');
        fprintf(fid,'*Tie, name=Tie_Face06Face05, adjust=no, position tolerance=%6.6f\n',2*b);
        fprintf(fid,'Entire_FACE06, Entire_FACE05\n');
        %
    end
    %
end
fprintf(fid,'**\n');
%
% STEP
%
if explicit_option == 1 % Running explicit analysis
    fprintf(fid,'** EXPLICIT Load Step 1 --------------------------------------\n');
    fprintf(fid,'*Step, name=Step-1, nlgeom=NO\n');
    fprintf(fid,'*Dynamic, Explicit\n');
    fprintf(fid,', 1.,, 0.0001\n');
    fprintf(fid,'*Bulk Viscosity\n');
    fprintf(fid,'0.06, 1.2\n');
    fprintf(fid,'** Mass Scaling: Semi-Automatic\n');
    fprintf(fid,'**               Whole Model\n');
    fprintf(fid,'*Variable Mass Scaling, dt=1e-05, type=below min, frequency=1\n');
else % Running implicit analysis
    fprintf(fid,'** IMPLICIT Load Step 1 --------------------------------------\n');
    fprintf(fid,'*STEP, INC=1000000, UNSYMM=YES\n');
    fprintf(fid,'*DYNAMIC, HAFTOL=1E-3\n');
    fprintf(fid,' 0.0001,1,1E-15,0.0001\n');
end
%
fprintf(fid,'** FOR PERIODIC BOUNDARY CONDITIONS\n');
fprintf(fid,'** Note: Node %9d is "dummy1", i.e. Epsilon_1i\n',dummy1);
fprintf(fid,'** Note: Node %9d is "dummy2", i.e. Epsilon_2j\n',dummy2);
fprintf(fid,'** Note: Node %9d is "dummy3", i.e. Epsilon_3k\n',dummy3);
%
%% Set the Periodic Boundary Conditions
%
if strcmp(stepID, 'p1')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_11 = 0.03\n');
        fprintf(fid,' %9d, 1, 1, 0.03\n',dummy1);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_11 = 0.03\n');
        fprintf(fid,'** %9d, 1, 1, 0.03\n',dummy1);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE01,1,1,%6,4f\n',0.030*c);
        fprintf(fid,'ALLnodes_FACE03,1,1\n');
    end
    %
elseif strcmp(stepID, 'p2')
    %    
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_22 = 0.03\n');
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy2);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_22 = 0.03\n');
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy2);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,2,2,%6.4f\n',0.030*a);
        fprintf(fid,'ALLnodes_FACE04,2,2\n');
    end
    %
elseif strcmp(stepID, 'p3')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_33 = 0.03\n');
        fprintf(fid,' %9d, 3, 3, 0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_33 = 0.03\n');
        fprintf(fid,'** %9d, 3, 3, 0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE06,3,3,%6.4f\n',0.030*b);
        fprintf(fid,'ALLnodes_FACE05,3,3\n');
    end
    %
elseif strcmp(stepID, 'p4')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_12 = Epsilon_21 = 0.03\n');
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy1);
        fprintf(fid,' %9d, 1, 1, 0.03\n',dummy2);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_12 = Epsilon_21 = 0.03\n');
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy1);
        fprintf(fid,'** %9d, 1, 1, 0.03\n',dummy2);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,1,1,%6.4f\n',0.030*a);
        fprintf(fid,'ALLnodes_FACE04,1,3\n');
    end
    %
elseif strcmp(stepID, 'p5')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_13 = Epsilon_31 = 0.03\n');
        fprintf(fid,' %9d, 3, 3, 0.03\n',dummy1);
        fprintf(fid,' %9d, 1, 1, 0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_13 = Epsilon_31 = 0.03\n');
        fprintf(fid,'** %9d, 3, 3, 0.03\n',dummy1);
        fprintf(fid,'** %9d, 1, 1, 0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE06,1,1,%6.4f\n',0.030*b);
        fprintf(fid,'ALLnodes_FACE05,1,3\n');
    end
    %
elseif strcmp(stepID, 'p6')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_23 = Epsilon_32 = 0.03\n');
        fprintf(fid,' %9d, 3, 3, 0.03\n',dummy2);
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_23 = Epsilon_32 = 0.03\n');
        fprintf(fid,'** %9d, 3, 3, 0.03\n',dummy2);
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy3);
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_13 = Epsilon_31 = 0.03\n');
        fprintf(fid,'** %9d, 3, 3, 0.03\n',dummy1);
        fprintf(fid,'** %9d, 1, 1, 0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,3,3,%6.4f\n',0.030*a);
        fprintf(fid,'ALLnodes_FACE04,1,3\n');
    end
    %
elseif strcmp(stepID, 'n1')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_11 = -0.03\n');
        fprintf(fid,' %9d, 1, 1, -0.03\n',dummy1);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_11 = -0.03\n');
        fprintf(fid,'** %9d, 1, 1, -0.03\n',dummy1);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE01,1,1,%6.4f\n',-0.030*c);
        fprintf(fid,'ALLnodes_FACE03,1,1\n');
    end
    %
elseif strcmp(stepID, 'n2')
    %    
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_22 = -0.03\n');
        fprintf(fid,' %9d, 2, 2, -0.03\n',dummy2);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_22 = -0.03\n');
        fprintf(fid,'** %9d, 2, 2, -0.03\n',dummy2);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,2,2,%6.4f\n',-0.030*a);
        fprintf(fid,'ALLnodes_FACE04,2,2\n');
    end
    %
elseif strcmp(stepID, 'n3')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_33 = -0.03\n');
        fprintf(fid,' %9d, 3, 3, -0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_33 = -0.03\n');
        fprintf(fid,'** %9d, 3, 3, -0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE06,3,3,%6.4f\n',-0.030*b);
        fprintf(fid,'ALLnodes_FACE05,3,3\n');
    end
    %
elseif strcmp(stepID, 'n4')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_12 = Epsilon_21 = -0.03\n');
        fprintf(fid,' %9d, 2, 2, -0.03\n',dummy1);
        fprintf(fid,' %9d, 1, 1, -0.03\n',dummy2);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_12 = Epsilon_21 = -0.03\n');
        fprintf(fid,'** %9d, 2, 2, -0.03\n',dummy1);
        fprintf(fid,'** %9d, 1, 1, -0.03\n',dummy2);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,1,1,%6.4f\n',-0.030*a);
        fprintf(fid,'ALLnodes_FACE04,1,3\n');
    end
    %
elseif strcmp(stepID, 'n5')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_13 = Epsilon_31 = -0.03\n');
        fprintf(fid,' %9d, 3, 3, -0.03\n',dummy1);
        fprintf(fid,' %9d, 1, 1, -0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_13 = Epsilon_31 = -0.03\n');
        fprintf(fid,'** %9d, 3, 3, -0.03\n',dummy1);
        fprintf(fid,'** %9d, 1, 1, -0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE06,1,1,%6.4f\n',-0.030*b);
        fprintf(fid,'ALLnodes_FACE05,1,3\n');
    end
    %
elseif strcmp(stepID, 'n6')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_23 = Epsilon_32 = -0.03\n');
        fprintf(fid,' %9d, 3, 3, -0.03\n',dummy2);
        fprintf(fid,' %9d, 2, 2, -0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_23 = Epsilon_32 = -0.03\n');
        fprintf(fid,'** %9d, 3, 3, -0.03\n',dummy2);
        fprintf(fid,'** %9d, 2, 2, -0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,3,3,%6.4f\n',-0.030*a);
        fprintf(fid,'ALLnodes_FACE04,1,3\n');
    end
    %
elseif strcmp(stepID, 'n3p6')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_33 = -0.03\n');
        fprintf(fid,' %9d, 3, 3, -0.03\n',dummy3);
        fprintf(fid,'** Note: Epsilon_23 = Epsilon_32 = 0.03\n');
        fprintf(fid,' %9d, 3, 3, 0.03\n',dummy2);
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_33 = -0.03\n');
        fprintf(fid,'** %9d, 3, 3, -0.03\n',dummy3);
        fprintf(fid,'**** Note: Epsilon_23 = Epsilon_32 = 0.03\n');
        fprintf(fid,'** %9d, 3, 3, 0.03\n',dummy2);
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE06,3,3,%6.4f\n',-0.030*b);
        fprintf(fid,'ALLnodes_FACE06,2,2,%6.4f\n',0.030*b);
        fprintf(fid,'ALLnodes_FACE05,1,3\n');
    end
    %
elseif strcmp(stepID, 'p2p3')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_22 = 0.03\n');
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy2);
        fprintf(fid,'** Note: Epsilon_33 = 0.03\n');
        fprintf(fid,' %9d, 3, 3, 0.03\n',dummy3);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_22 = 0.03\n');
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy2);
        fprintf(fid,'**** Note: Epsilon_33 = 0.03\n');
        fprintf(fid,'** %9d, 3, 3, 0.03\n',dummy3);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,2,2,%6.4f\n',0.030*a);
        fprintf(fid,'ALLnodes_FACE04,2,2\n');
        fprintf(fid,'ALLnodes_FACE06,3,3,%6.4f\n',0.030*b);
        fprintf(fid,'ALLnodes_FACE05,3,3\n');
    end
    %
elseif strcmp(stepID, 'p2p4')
    %
    if pbcs_option == 1
        fprintf(fid,'*BOUNDARY,OP=NEW\n');
        fprintf(fid,'** Note: Epsilon_22 = 0.03\n');
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy2);
        fprintf(fid,'** Note: Epsilon_12 = Epsilon_21 = 0.03\n');
        fprintf(fid,' %9d, 2, 2, 0.03\n',dummy1);
        fprintf(fid,' %9d, 1, 1, 0.03\n',dummy2);
    else % If using PSEUDO-periodic boundary conditions
        fprintf(fid,'***BOUNDARY,OP=NEW\n');
        fprintf(fid,'**** Note: Epsilon_22 = 0.03\n');
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy2);
        fprintf(fid,'**** Note: Epsilon_12 = Epsilon_21 = 0.03\n');
        fprintf(fid,'** %9d, 2, 2, 0.03\n',dummy1);
        fprintf(fid,'** %9d, 1, 1, 0.03\n',dummy2);
        fprintf(fid,'** Since PSEUDO-PERIODIC BOUNDARY CONDITIONS ARE BEING APPLIED:\n');
        fprintf(fid,'*BOUNDARY,OP=NEW, AMPLITUDE=Amp-1\n');
        fprintf(fid,'ALLnodes_FACE02,1,2,%6.4f\n',0.030*a);
        fprintf(fid,'ALLnodes_FACE04,1,2\n');
    end
    %
end
%
%% Write remaining lines that are common to every load case
fprintf(fid,'*OUTPUT, FIELD, TIME INTERVAL=0.01\n');
fprintf(fid,'*ELEMENT OUTPUT\n');
fprintf(fid,' S, E, IVOL, EVOL\n');
fprintf(fid,'*ELEMENT OUTPUT, ELSET=MATRIX\n');
fprintf(fid,' SDV\n');
fprintf(fid,'*ELEMENT OUTPUT, ELSET=FIBRE\n');
fprintf(fid,' SDV\n');
if cohesive_choice == 1
    fprintf(fid,'*ELEMENT OUTPUT, ELSET=COHESI\n');
end
fprintf(fid,' SDEG\n');
fprintf(fid,'*NODE OUTPUT\n');
fprintf(fid,' U\n');
fprintf(fid,'*OUTPUT, HISTORY, TIME INTERVAL=0.01\n');
fprintf(fid,'*ENERGY OUTPUT\n');
fprintf(fid,' ALLIE, ALLAE, ALLCD, ALLWK, ALLKE, ALLPD, ALLSE, ETOTAL\n');
if explicit_option == 0 % Running implicit analysis
    fprintf(fid,'*CONTROLS,PARAMETERS=FIELD\n');
    fprintf(fid,' 0.01, , , , 0.05\n');
    fprintf(fid,'*CONTROLS,PARAMETERS=TIME INCREMENTATION\n');
    fprintf(fid,' 10, 10, 4, 500, 10, 4, 30, 50, 10, 3\n');
    fprintf(fid,' .25, .5, .75, .85\n');
end
fprintf(fid,'*CONTROLS,PARAMETERS=FIELD,FIELD=DISPLACEMENT\n');
fprintf(fid,' 0.10, 1., , , 0.10\n');
fprintf(fid,'*END STEP\n');
%
fclose(fid);
    %
%     tmp=pwd;
%     cd(dir_fea_name);
%     copyfile('../../../dir_umat_final_thick.f','dir_umat.f');
%     copyfile('../../../yield_c.inp');
%     copyfile('../../../yield_t.inp');
%     copyfile('../../../run.pbs');
%     copyfile('../../../abaqus_v6_to_copy.env','abaqus_v6.env');
%
% unix(strcat(abaqus_path,' j=run ask_delete=OFF int cpus=4'));    cd(tmp);
%
%
%
disp(' ');
disp('run.inp file created');
disp('Elapsed Time [min]: ');
disp(toc/60);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%