%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Check for Incompatibility of the Fibre Position       %
%                                                                         %
%  This file is part of the RAND_uSTRU_GEN.m file                         %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v2.0                                   %
%                                                                         %
%                Ant�nio Rui Melro - antonio.melro@fe.up.pt               %
%                                June 2006                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [fibre_overlap,min] = f_overlap(X_TMP,Y_TMP,fibre2test)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global N_fibre Fibre_pos DISTMIN a
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Verify Distribution for Fibre Overlaps                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
fibre_overlap = 0;
min = a*2;
xlim0 = X_TMP - 4*DISTMIN; xlim1 = X_TMP + 4*DISTMIN;
ylim0 = Y_TMP - 4*DISTMIN; ylim1 = Y_TMP + 4*DISTMIN;
for k=1:N_fibre
    if k ~= fibre2test
        xx = Fibre_pos(k,2); yy = Fibre_pos(k,3);
        if xx > xlim0 && xx < xlim1 && yy > ylim0 && yy < ylim1
            new_dist = sqrt((xx-X_TMP)^2 + (yy-Y_TMP)^2);
            if new_dist < DISTMIN
                fibre_overlap = 1;
            else
                if new_dist < min
                    min = new_dist;
                end
            end
        end
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

