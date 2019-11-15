%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       Function to Quantitatively Characterise Randomness of             %
%       Fibre Distribution                                                %
%                                                                         %
%  This file is part of the XXXX_uSTRU_GEN.m sequence of files            %
%  and should not be used separately from it.                             %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                                  v5.0                                   %
%                                                                         %
%                António Rui Melro - antonio.melro@fe.up.pt               %
%                               March 2007                                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%   MODIFICATION LOG:                                                     %
%                                                                         %
%       2007/03/23 - Added First Neighbour Orientation                    %
%       2006/07/18 - Added Pair Distribution Function                     %
%       2006/07/17 - Added Second Order Intensity Function                %
%       2006/07/14 - Added Neighbour Distances Calculation                %
%       2006/07/13 - First Release                                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function f_characterise_per
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of Global Variables                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
global a b R Fibre_pos N_fibre DISTMIN L_h h_R g_h neigh neigh_o;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1 - Voronoi Polygon Areas                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for i=1:N_fibre
    X(i,:) = [Fibre_pos(i,2) Fibre_pos(i,3)];
end
[V,C] = voronoin(X);
pArea = zeros(N_fibre,1);
k_poly = 0;
for i=1:N_fibre
    Top_lim = length(C{i});
    GO = 0;
    for j=1:Top_lim
        X_TMP = V(C{i}(1,j),1); Y_TMP = V(C{i}(1,j),2);
        if X_TMP < -R || X_TMP > a+R || Y_TMP < -R || Y_TMP > b+R
            GO = 1;
        end
    end
    if GO == 0
        XX = zeros(Top_lim,1); YY = zeros(Top_lim,1);
        for j=1:Top_lim
            XX(j) = V(C{i}(1,j),1);
            YY(j) = V(C{i}(1,j),2);
        end
        k_poly = k_poly + 1;
        pArea(k_poly) = polyarea(XX,YY);
    end
end
for i=1:k_poly
    pArea1(i)=pArea(i);
end
coeff_var(1) = std(pArea1) / mean(pArea1);
disp('Coefficient of Variation of Voronoi polygon areas: ');
disp(coeff_var(1)); disp (' ');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2 - Voronoi Polygon Distances                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
k_dist = 0;
for i=1:1:N_fibre-1
    Top_lim = length(C{i});
    xlim0 = X(i,1) - 4*DISTMIN; xlim1 = X(i,1) + 4*DISTMIN;
    ylim0 = X(i,2) - 4*DISTMIN; ylim1 = X(i,2) + 4*DISTMIN;
    for j=i+1:N_fibre
        if X(j,1) > xlim0 && X(j,1) < xlim1 && X(j,2) > ylim0 && X(j,2) < ylim1
            Top_lim1 = length(C{j});
            GO = 0;
            for k=1:Top_lim
                vert_test = C{i}(1,k);
                if V(C{i}(1,k),1) > -R && V(C{i}(1,k),1) < a+R && V(C{i}(1,k),2) > -R && V(C{i}(1,k),2) < b+R
                    for l=1:Top_lim1
                        if vert_test == C{j}(1,l) && V(C{j}(1,l),1) > -R && V(C{j}(1,l),1) < a+R && V(C{j}(1,l),2) > -R && V(C{j}(1,l),2) < b+R
                            GO = GO + 1;
                        end
                    end
                end
            end
        end
        if GO == 2
            k_dist = k_dist + 1;
            fdist(k_dist) = sqrt((X(i,1)-X(j,1))^2 + (X(i,2)-X(j,2))^2);
        end
        GO = 0;
    end
end
coeff_var(2) = std(fdist)/mean(fdist);
disp('Coefficient of Variation of Voronoi polygon distances: ');
disp(coeff_var(2)); disp (' ');
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 - Second Order Intensity Function                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
h = [DISTMIN:0.5*R:0.3*a];
for i=1:N_fibre
    X(i,1) = Fibre_pos(i,2);
    X(i,2) = Fibre_pos(i,3);
end
bound = [-2*R -2*R;a+2*R -2*R;a+2*R b+2*R;-2*R b+2*R];
K_h = cskhat(X, bound, h);
% %
% % K_h for Hard-core point set
% %
% for i=1:1:length(h)
%     K_h_HC(i) = pi*(h(i)-DISTMIN)^2;
% end
% %
% % K_h for Poisson point set
% %
% for i=1:1:length(h)
%     K_h_Poisson(i) = pi*(h(i))^2;
% end
% %
% figure(2);
% hold on;
% box on;
% title(['Second Order Intensity Function']);
% whitebg('white');
% plot(h,K_h,'ko-');
% plot(h,K_h_HC,'k--');
% plot(h,K_h_Poisson,'k:');
% axis([0 h(length(h))+R 0 max(K_h)*1.1]);
% legend({'Present Model','Hard Core','Poisson'},...
%     'EdgeColor',[0 0 0],'TextColor',[0 0 0],'Location','NW');
% xlabel('Distance [mm]');
% ylabel('K(r)');
% hold off;
%
for i=1:length(h)
    L_h(i) = sqrt(K_h(i)/pi)-h(i);
end
h_R = h/R;
%
% figure(3);
% hold on;
% box on;
% % title(['Deviation of K(r) from Poisson Process']);
% whitebg('white');
% plot(h_R,L_h,'ko-');
% xlabel('h/R');
% ylabel('L(h)');
% hold off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 - Pair Distribution Function                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% p = polyfit(h,K_h,2);
% dp(1) = p(1)*2; dp(2)=p(2);
p = polyfit(h,K_h,10);
dp(1) = p(1)*10; dp(2)=p(2)*9; dp(3)=p(3)*8; dp(4)=p(4)*7; dp(5)=p(5)*6; dp(6)=p(6)*5;
dp(7)=p(7)*4; dp(8)=p(8)*3; dp(9)=p(9)*2; dp(10)=p(10);
for i=1:length(h)
    g_h(i) = polyval(dp,h(i))/(2*pi*h(i));
end
% %
% % figure(2); hold on; plot(h,polyval(p,h),'b'); hold off;
% %
% figure(4);
% hold on;
% box on;
% title(['Pair Distribution Function']);
% whitebg('white');
% plot(h,g_h,'ko-');
% plot(h,1,'k:');
% xlabel('Distance [mm]');
% ylabel('g(r)');
% hold off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5 - Neighbour Distances                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
Neighbour = zeros(N_fibre,5);
for i=1:N_fibre
    Neighbour(i,1) = i;
    Neighbour(i,2) = 2*a;
    Neighbour(i,3) = 3*a;
    Neighbour(i,4) = 4*a;
    XC = Fibre_pos(i,2); YC = Fibre_pos(i,3);
    xlim0 = XC - 4*DISTMIN; xlim1 = XC + 4*DISTMIN;
    ylim0 = YC - 4*DISTMIN; ylim1 = YC + 4*DISTMIN;
    for j=1:N_fibre
        if i ~= j
            X_TMP = Fibre_pos(j,2); Y_TMP = Fibre_pos(j,3);
            if X_TMP > xlim0 && X_TMP < xlim1 && Y_TMP > ylim0 && Y_TMP < ylim1
                dist = sqrt((X_TMP-XC)^2 + (Y_TMP-YC)^2);
                if dist < Neighbour(i,4)
                    if dist < Neighbour(i,3)
                        if dist < Neighbour(i,2)
                            Neighbour(i,2) = dist;
                            if Y_TMP >= YC
                                Neighbour(i,5) = acos((X_TMP-XC)/dist);
                            else
                                Neighbour(i,5) = acos((X_TMP-XC)/dist)+pi;
                            end
                        else
                            Neighbour(i,3) = dist;
                        end
                    else
                        Neighbour(i,4) = dist;
                    end
                end
            end
        end
    end
end
%
% figure(5);
% hold on;
% box on;
% title(['Neighbour Distances']);
% whitebg('white');
for i=2:4
    min_d = DISTMIN; max_d = [0 3*R 4.5*R 8*R]; n_step_d = 12;
%     min_d = min(Neighbour(:,i)); max_d = max(Neighbour(:,i)); n_step_d =
%     12;
    step_d = (max_d(i)-min_d)/n_step_d; l = 0;
    for j=min_d:step_d:max_d(i)-step_d
        l = l + 1;
        neigh(l,2*i-3) = j + step_d/2; neigh (l,2*i-2) = 0;
        for k=1:N_fibre
            if Neighbour(k,i) >= j && Neighbour(k,i) < j+step_d
                neigh(l,2*i-2) = neigh(l,2*i-2) + 1;
            end
        end
    end
%     neigh(l,2*i-2) = neigh(l,2*i-2) + 1;
    neigh(:,2*i-2) = neigh(:,2*i-2)/N_fibre;
%     %
%     if i == 2
%         plot(neigh(:,2*i-3),neigh(:,2*i-2),'ko-');
%     elseif i == 3
%         plot(neigh(:,2*i-3),neigh(:,2*i-2),'ks-');
%     else
%         plot(neigh(:,2*i-3),neigh(:,2*i-2),'k^-');
%     end
end
% % h=[0.005:0.0005:0.015];
% % a1 = (a-2*R)^2; a2 = 4*(2*R)*(a-2*R); a3 = 4*(2*R)^2; at = a1+a2+a3; af = pi*R*R; A=a*b;
% % lambda_poisson = Vol_fibre*(a1*a1+a2*a2+a3*a3)/(af*A*at);
% % D_h = lambda_poisson.*2.*pi.*h.*exp(-1*pi.*lambda_poisson.*h.*h);
% % plot(h,D_h,'k--');
% axis([min(neigh(:,1))-step_d max(neigh(:,5))+step_d...
%     0 0.1+max([max(neigh(:,2)) max(neigh(:,4)) max(neigh(:,6))])]);
% legend({'1st Neighbour Distance','2nd Neighbour Distance','3rd Neighbour Distance'},...%'Poisson Distribution'},...
%     'EdgeColor',[0 0 0],'TextColor',[0 0 0],'Position',[0.53 0.75 0.35 0.14]);
% xlabel('Distance [mm]');
% ylabel('Probability Density Function');
% hold off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6 - Neighbour Orientations                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% figure(6);
% hold on;
% box on;
% title(['First Neighbour Orientations']);
% whitebg('white');
min_d = 0; max_d = 2*pi; n_step_d = 24;
step_d = (max_d-min_d)/n_step_d;
neigh_o(1,1) = 0; neigh_o(1,2) = 0;
l = 1;
for j=min_d:step_d:max_d-step_d
    l = l + 1;
    neigh_o(l,1) = j+step_d; neigh_o (l,2) = 0;
    for k=1:N_fibre
        if Neighbour(k,5) >= j && Neighbour(k,5) < j+step_d
            neigh_o(l,2) = neigh_o(l,2) + 1;
        end
    end
end
neigh_o(:,2) = neigh_o(:,2)/N_fibre;
for i=2:l
    neigh_o(i,2) = neigh_o(i,2) + neigh_o(i-1,2);
end
% %
% plot([0 2*pi],[0 1],'k:');
% plot(neigh_o(:,1),neigh_o(:,2),'ko-');
% axis([0 2*pi+0.2 0 0.05+max(neigh_o(:,2))]);
% xlabel('Orientation [rad]');
% ylabel('Cumulative Distribution Function');
% hold off;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of Function                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%