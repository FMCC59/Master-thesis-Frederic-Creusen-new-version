%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%            Matlab Routine to perform finite element analyses            %
%               on generated random distributions of fibres               %
%                                                                         %
%                                   v2.0                                  %
%                August 31, 2015                                          %
%                                                                         %
%                Code initially developed by:                             %
%                Antonio Rui Melro - antonio.melro@fe.up.pt               %
%                from University of Porto                                 %
%                                                                         %
%                and further extended and debugged by:                    %
%                Miguel A. Bessa - mbessa@u.northwestern.edu              %
%                                                                         %
%                Key publications from A.R. Melro:                        %
%   1-> To understand function "RAND_PER_uSTRU_GEN_3D":                   %
%       A.R. Melro, P.P. Camanho and  S.T. Pinho., 2008. Generation of    %
%       random distribution of fibres in long-fibre reinforced composites.%
%       Composites Science and Technology, 68, 2092-2102.                 %
%                                                                         %
%   2-> To understand the periodic boundary conditions implemented in     %
%       function "f_mesh_quad_per_3D":                                    %
%       Melro, A. R., Camanho, P. P., Pires, F. A., & Pinho, S. T. (2013).%
%       Micromechanical analysis of polymer composites reinforced by      %
%       unidirectional fibres: Part II–Micromechanical analyses.          %
%       International Journal of Solids and Structures, 50(11), 1906-1915.%
%                                                                         %
%                Key publications from the Northwestern team:             %
%   1-> To understand new plasticity law:                                 %
%       X. Bai, M.A. Bessa, et al. (2015). High-fidelity micro-scale      %
%       modeling of the thermo-visco-plastic behavior of carbon fiber     %
%       polymer matrix composites. Composite Structures, 134, 132-141.    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   APPENDED FUNCTIONS                                                    %
%                                                                         %
%         RAND_PER_uSTRU_GEN_3D - Generates a random distribution of      %
%                                 circles that represent the real fibre   %
%                                 distribution                            %
%         f_image_per - Generates images of fibre distribution in .BMP    %
%         f_characterise_per - Quantitative characterisation of the       %
%                              randomness of the fibre distribution       %
%         f_output - Generates output files in both .MAT and .PY formats  %
%         f_mesh_quad_per_3D_FINAL - Generates HEX + Wedge mesh with      %
%                              Periodic Boundary Conditions               %
%         f_fea_3D_FINAL - Generates global Input file for Abaqus,        %
%                             depending on the loading applied            %
%         f_post_proc - Performs post-processing analysis with the        %
%                       purpose to obtain global stresses and strains.    %
%                       Applies only to elastic analyses.                 %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clc; close all; clear all; tic; fclose('all');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global R delta Vol_fibre_req DISTMIN Max_fibres a b c N_guesses_max ...
    N_cycles_max N_change Square_size Square_inc inter_option ...
    quant_option image_option finim_option vorim_option mesh_option ...
    FEA_option mat_name dir_name s11 s22 s33 s12 s13 s23 S_base stepID;
global Fibre_pos dir_fea_name cohesive_choice geometry_option ...
    abaqus_path explicit_option plastic_option pbcs_option ...
    hardening_curves N R_A R_B linewidth thickness flowmultiplier;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Parameters                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
mat_name = 'NU_Daniel_Material';      % Name of the material input file (without ".mat"): mat_name.mat
%
hardening_curves = 'NU_Daniel_Material_hardening_curves'; % Name of file with hardening curves for this material (without ".txt")
%
abaqus_path = '/var/DassaultSystemes/SIMULIA/Commands/abaqus';  % Specify path to ABAQUS command
%/var/DassaultSystemes/SIMULIA/Commands
file_name = strcat(cd,'/Input/',mat_name,'.mat');
fid = fopen(file_name,'r');
%
for i=1:1:13, fgetl(fid); end
R = str2num(fgetl(fid))/2; fgetl(fid);  % fibre radius
Vol_fibre_req = str2num(fgetl(fid));    % fibre volume we wish to obtain
%
fclose(fid);
%
delta = 6;             % ratio a/R
DISTMIN = 2*R;         % Minimum distance between fibre centres %%%%% 2.07
Max_fibres = (delta/2)^2;  % Maximum number of fibres 2*R
N_guesses_max = 50E3;  % number of guesses in the random generation of fibre position
N_cycles_max = 15;     % maximum number of cycles that the routine runs
N_change = 3;          % number of iterations before changing criteria on First Heuristic; MinValue = 3
Square_size = 3*R;     % initial size of square in Second Heuristic
Square_inc = (8.5-10*Vol_fibre_req)*R; % increment to be given to the square size in Second Heuristic
S_base = 2*pi/60*R;    % Average of element size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option Parameters                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Essential options:
%
initdist = 1;           % Number to be given to the first distribution of fibers. This
                        % number will appear in the name of the DIST folder. For
                        % example, initdist=2 => folder with name "DIST02".
                        %
distcount = 1;          % Number of distributions to run. Each fiber distribution will create
                        % a folder inside the "mat_name" folder. For example, if
                        % initdist=2 and distcount=3 then the following folders
                        % will be created: DIST02, DIST03, DIST04
                        %
geometry_option = 0;    % Generate the random fiber distribution (RVE geometry); 0-NO, 1-YES
                        % if selected 0-NO then you must have a file named
                        % 'ALLvariables_RAND_FiberGeneration.mat' in each
                        % fiber distribution folder DIST01, DIST02, ...
                        %
mesh_option = 1;        % Generate mesh and auxiliary input files; 0-NO, 1-YES
                        %
FEA_option = 1;         % Generates Finite Element main input (run.inp); 0-NO, 1-YES
                        %
explicit_option = 10;    % Run explicit analysis? 0-NO (runs implicit), 1-YES 
                        %
plastic_option = 0;     % Run plasticity analysis? 0-NO, 1-YES
                        % Note: IMPLICIT analysis DOES NOT have a UMAT,
                        % we only wrote a VUMAT (for EXPLICIT)
                        %
pbcs_option = 1;        % Use Periodic Boundary Conditions? 0-NO, 1-YES
                        % Important Note: in EXPLICIT using periodic
                        % boundary conditions is EXTREMELY slow
                        % (maybe impossible?)
                        %
cohesive_choice = 0;    % Include cohesive elements on the mesh; 0-NO, 1-YES
%
%
% Other options:
%
% inter_option = 0;       % (NOT IMPLEMENTED) model the interlaminar area; 0-NO, 1-YES
quant_option = 0;       % quantitave analysis of distribution; 0-NO, 1-YES
image_option = 1;       % print an image of fibre distribution in the end of each step; 0-NO, 1-YES
finim_option = 1;       % plot the final image with fibre distribution; 0-NO, 1-YES
vorim_option = 0;       % draw in the final output image the Voronoi polygons; 0-NO, 1-YES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of analysis to run                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
step_type = {'p2' 'n2' 'p6' 'p2p3'};
% Step Type Id:
%               p1 --> Positive Longitudinal Normal Load (+EPS11);
%               p2 --> Positive Transverse Normal Load (+EPS22);
%               p3 --> Positive Transverse Normal Load (+EPS33);
%               p4 --> Positive Longitudinal Shear Load (+EPS12 & +EPS21);
%               p5 --> Positive Longitudinal Shear Load (+EPS13 & +EPS31);
%               p6 --> Positive Transverse Shear Load (+EPS23 & +EPS32);
%               n1 --> Negative Longitudinal Normal Load (-EPS11);
%               n2 --> Negative Transverse Normal Load (-EPS22);
%               n3 --> Negative Transverse Normal Load (-EPS33);
%               n4 --> Negative Longitudinal Shear Load (-EPS12 & -EPS21);
%               n5 --> Negative Longitudinal Shear Load (-EPS13 & -EPS31);
%               n6 --> Negative Transverse Shear Load (-EPS23 & -EPS32);
%             p2p3 --> Biaxial Stress State (+EPS22 & +EPS33);
%             p2p4 --> Transv. Tensile + Long. Shear (+EPS22 & -EPS12 & -EPS21);
%             n3p6 --> Biaxial Stress State (-EPS33 & +EPS23 & +EPS32);
%
% NOTE: We didn't code other combinations... It is trivial to change the
% input file (run.inp) to apply any desired load.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determination of RVE Dimensions                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = delta*R;                  % RVE width (yy direction)
b = delta*R;                  % RVE height (zz direction, aka, thickness direction)
c = round(1e5*2*S_base)/1e5;  % RVE thickness (xx direction, aka, fibre direction)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate Fibre Distribution                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for j=initdist:1:(initdist-1+distcount)
    status_mesh = 0;
    % Create a folder for this distribution, if it doesn't already exist:
    dir_name = strcat(cd,'/',mat_name,'/DIST0',int2str(j));
    if exist(dir_name,'dir') == 0, mkdir(dir_name); end
    %
    while status_mesh == 0
        status = 0;
        if geometry_option == 1
            while status == 0
                status = RAND_PER_uSTRU_GEN_3D; % Run random fiber generation
                if status == 1 % We converged to a valid fiber distribution
%                     % Save all variables from workspace to file
%                     file_name = strcat(dir_name,...
%                         '/AllVariables_RAND_FiberGeneration.mat');
%                     save(file_name)
                    break;
                end
            end
        else % if not asked to generate geometry, then load variables
            %previously generated by RAND_PER_uSTRU_GEN_3D
            file_name = strcat(dir_name,...
                '/AllVariables_RAND_FiberGeneration.mat');
            % Load all variables:
            load(file_name);
        end
        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternative Fibre distribution for FDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        N=1 %number of roads in x and y        
        thickness = 0.2
        linewidth = 0.4
        flowmultiplier = 1.1
        R_A = linewidth/2*flowmultiplier
        R_B = thickness/2*flowmultiplier
        A_road = pi*R_A*R_B
        a = N^0.5*linewidth
        b = N^0.5*thickness
        c = a

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Postition of fibres in a N*N matrix %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        Fibre_pos = zeros(N,4)
        j=2 %x coordinates
        for i =0:(N^0.5)-1
            for k =1:(N^0.5)
                Fibre_pos(k+(i*N^0.5),j)= linewidth*(0.5+(k-1))
            end
        end
        j=3 %y coordinates
        for i =0:(N^0.5)-1
            for k =1:(N^0.5)
                Fibre_pos(k+(i*N^0.5),j)= thickness*(0.5+(i))
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Prints the Image                                                        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if finim_option == 1
            index1 = 0;
            f_image_per(index1);
        end
        drawnow;
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Quantitative Characterization of Randomness of Fibre Distribution       %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if quant_option == 1
            f_characterise_per;
        end
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Generates Files with Mesh, PBCs and Material Properties                 %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        if mesh_option == 1
%             dir_name = strcat(cd,'/',mat_name,'/DIST0',int2str(j));
%             if exist(dir_name,'dir') == 0, mkdir(dir_name); end
            %
            status_mesh = f_mesh_quad_per_3D_FINAL;
%             if cohesive_choice == 1
%                 status_mesh = f_mesh_quad_per_3D; % Mesh with cohesive elements
%             else % Mesh without cohesive elements
%                 status_mesh = f_mesh_quad_per_3D_noCohesive;
%             end
            %
            if status_mesh == 0, continue, end
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Runs the Finite Element Analysis                                        %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            if FEA_option == 1
                for i=1:1:length(step_type)
                    stepID = step_type(i);
                    f_fea_3D_FINAL(step_type(i));
%                     dir_fea_name = strcat(dir_name,'/ID0',num2str(i));
%                     if stepID <= 6 && stepID >= 1
%                         f_post_proc(step_type(i));
%                         Dijkl(1,stepID) = s11; Dijkl(2,stepID) = s22;
%                         Dijkl(3,stepID) = s33; Dijkl(4,stepID) = s12;
%                         Dijkl(5,stepID) = s13; Dijkl(6,stepID) = s23;
%                     end
                end
                %
                % Determine Engineering Properties
                %
%                 if stepID <= 6 && stepID >= 1
%                     D11 = Dijkl(1,1);
%                     D22 = (Dijkl(2,2)+Dijkl(3,3))/2;
%                     D23 = (Dijkl(2,3)+Dijkl(3,2))/2;
%                     D12 = (Dijkl(1,2)+Dijkl(2,1)+Dijkl(1,3)+Dijkl(3,1))/4;
%                     %
%                     E1 = D11-2*D12*D12/(D22+D23);
%                     v12 = D12/(D22+D23);
%                     E2 = (D11*(D22+D23)-2*D12*D12)*(D22-D23)/(D11*D22-D12*D12);
%                     v23 = (D11*D23-D12*D12)/(D11*D22-D12*D12);
%                     G12 = Dijkl(4,4);
%                     G13 = Dijkl(5,5);
%                     G23 = (D22-D23)/2;
%                 end
%                 save(strcat(dir_name,'/',mat_name,'.mat'));
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Routine                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
