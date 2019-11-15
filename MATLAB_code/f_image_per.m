%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Generate an image of the RVE in .eps format           %
%                                                                         %
%  This file is part of the XXXX_PER_uSTRU_GEN.m sequence of files        %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v2.0                                   %
%                                                                         %
%                António Rui Melro - antonio.melro@fe.up.pt               %
%                             October  2008                               %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  Legend for variable "index1":                                          %
%       0 - image in the end of routine RAND_PER_uSTRU_GEN.m              %
%       1 - image in the end of step 1 of routine RAND_PER_uSTRU_GEN.m    %
%       2 - image in the end of step 2 of routine RAND_PER_uSTRU_GEN.m    %
%       3 - image in the end of step 3 of routine RAND_PER_uSTRU_GEN.m    %
%       4 - image for SPER_uSTRU_GEN.m                                    %
%       5 - image for HPER_uSTRU_GEN.m                                    %
%       6 - image for WONG_uSTRU_GEN.m                                    %
%       9 - image for animation sequence of routine RAND_PER_uSTRU_GEN.m  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function f_image_per(index1)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global a b R N_cycles Fibre_pos N_fibre vorim_option Vol_fibre A_1_fibre ...
    R_A R_B linewidth thickness N flowmultiplier;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Opens Figure and Adjusts Properties                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure(1);
hold on;
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'TickDir','out');
rectangle('Position',[0,0,a,b],'EdgeColor','w');
axis([-R a+R -R b+R]);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates the Title                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if index1 == 1 || index1 == 2 || index1 == 3
    title1 = num2str(N_cycles);
    title2 = num2str(index1);
    title3 = num2str(Vol_fibre*100);
    title(['Iteration: ', title1,';   Step: ', title2,';   v_f=', title3,'%']);
elseif index1 == 0
    title1 = num2str(N_cycles-1);
    title3 = num2str(Vol_fibre*100);
    title(['Iteration: ', title1,';   Step: End;   v_f=', title3,'%']);
elseif index1 == 4 || index1 == 5
    title3 = num2str(Vol_fibre*100);
    title(['Fibre Volume = ', title3,'%']);
elseif index1 == 6
    title1 = num2str(N_cycles);
    title2 = num2str(Vol_fibre*100);
    title(['Iteration: ', title1,';   Fibre Volume = ', title2,'%']);
elseif index1 == 9
    title3 = num2str(Vol_fibre*100);
    title(['v_f=', title3,'%']);
end
whitebg('black');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots the Position of all Fibres                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fibre_scale = 25.4^2 * 72^2 * (210-3*a/R)/40;
% for i=1:N_fibre
   % scatter(Fibre_pos(i,2),Fibre_pos(i,3),A_1_fibre*Fibre_scale,'w','filled');
%     %      text(Fibre_pos(i,2),Fibre_pos(i,3),num2str(i),'Color','r','FontWeight',...
%     %          'bold','HorizontalAlignment','center');
% end
% if vorim_option == 1
%     voronoi(Fibre_pos(:,2),Fibre_pos(:,3),'w-');
% end
% hold off;

axis auto %[0 linewidth*3* 0 thickness*3*])
pbaspect ([1 1 1])
for i=1:N
xp = Fibre_pos(i,2)
yp = Fibre_pos(i,3)
 t = linspace(0,2*pi);              %linearly spaced vector (2pi)
 X = R_A*cos(t);                      %Definition of X axis of circumpherence
 Y = R_B*sin(t);                      %Definition of Y axis of circumpherence
 x = xp + X; %Definition of X axis of circumpherence
 y = yp + Y; %Definition of Y axis of circumpherence
 plot(x,y,'r-')
 hold on
end
%scatter(Fibre_pos(:,2),Fibre_pos(:,3),200,"r");
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saves the Image in .eps Format                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if index1 == 1
    name = num2str(3*N_cycles-2,'%03d');
    saveas(1,name,'eps');
%     print -dps2 -append anim;
    close(1);
end
if index1 == 2
    name = num2str(3*N_cycles-1,'%03d');
    saveas(1,name,'eps');
%     print -dps2 -append anim;
    close(1);
end
if index1 == 3
    name = num2str(3*N_cycles,'%03d');
    saveas(1,name,'eps');
%     print -dps2 -append anim;
    close(1);
end
%N.m sequence of files        %
%  and should not be used separately from it.
if index1 == 4
    name = 'SPER_uSTRU';
    saveas(1,name,'eps');
    close(1);
end
if index1 == 5
    name = 'HPER_uSTRU';
    saveas(1,name,'eps');
    close(1);
end
if index1 == 6
    name = 'WONG_uSTRU';
    saveas(1,name,'eps');
    close(1);
end
if index1 == 9
    print -dps2 -append animation;
    close(1);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%