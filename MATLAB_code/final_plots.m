clc; clear all; close all;
% 
%
% Figure for ID11
%

load(strcat(cd,'/Silenka_Fiedler/DIST06/ID11/stress-strain.mat'),'E','S');
E_01_11 = E;
S_01_11 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST07/ID11/stress-strain.mat'),'E','S');
E_02_11 = E;
S_02_11 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST08/ID11/stress-strain.mat'),'E','S');
E_03_11 = E;
S_03_11 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST09/ID11/stress-strain.mat'),'E','S');
E_04_11 = E;
S_04_11 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST010/ID11/stress-strain.mat'),'E','S');
E_05_11 = E;
S_05_11 = S;
load(strcat(cd,'/Silenka_Fiedler_20/DIST02/ID11/stress-strain.mat'),'E','S');
E_01_20_11 = E;
S_01_20_11 = S;

clear E S;

max01 = 0; imax01 = 1; max02 = 0; imax02 = 1; max03 = 0; imax03 = 1;
max04 = 0; imax04 = 1; max05 = 0; imax05 = 1; max06 = 0; imax06 = 1;
for i=2:1:length(S_01_11), if S_01_11(i) > max01, imax01=i; max01=S_01_11(i); end, end
for i=2:1:length(S_02_11), if S_02_11(i) > max02, imax02=i; max02=S_02_11(i); end, end
for i=2:1:length(S_03_11), if S_03_11(i) > max03, imax03=i; max03=S_03_11(i); end, end
for i=2:1:length(S_04_11), if S_04_11(i) > max04, imax04=i; max04=S_04_11(i); end, end
for i=2:1:length(S_05_11), if S_05_11(i) > max05, imax05=i; max05=S_05_11(i); end, end
for i=2:1:length(S_01_20_11), if S_01_20_11(i) > max06, imax06=i; max06=S_01_20_11(i); end, end

imax01 = imax01 + 3; imax02 = imax02 + 3; imax03 = imax03 + 3;
imax04 = imax04 + 3; imax05 = imax05 + 3; imax06 = imax06 + 3;

for i=length(S_01_11):-1:imax01, S_01_11(i)=[];E_01_11(i)=[]; end
for i=length(S_02_11):-1:imax02, S_02_11(i)=[];E_02_11(i)=[]; end
for i=length(S_03_11):-1:imax03, S_03_11(i)=[];E_03_11(i)=[]; end
for i=length(S_04_11):-1:imax04, S_04_11(i)=[];E_04_11(i)=[]; end
for i=length(S_05_11):-1:imax05, S_05_11(i)=[];E_05_11(i)=[]; end
for i=length(S_01_20_11):-1:imax06, S_01_20_11(i)=[];E_01_20_11(i)=[]; end

figure(1);
hold on;
set(gca,'FontSize',13);
xlabel('Far-field Strain \epsilon_{11}^o (%)','Position',[0.035/2 -130]);
ylabel('\sigma_{11}^o (MPa)');
rectangle('Position',[0.0001,min([max01,max02,max03,max04,max06]),0.05,...
    max([max01,max02,max03,max04,max06])-...
    min([max01,max02,max03,max04,max06])],'Curvature',[0,0],...
    'Linestyle','none','FaceColor',[0.9 0.9 0.9]);
plot(E_01_11,S_01_11,'k-');
plot(E_01_20_11,S_01_20_11,'ro-');
plot(E_05_11,S_05_11,'b-x');
plot(E_02_11,S_02_11,'k-');
plot(E_04_11,S_04_11,'k-');
plot(E_03_11,S_03_11,'k-');
axis([0 0.035 0 1400]);
legend('CASES 1,2,3,4','CASE 6','CASE 5','Location','SouthEast');
legend boxoff;
hold off;
 
%
% Figure for ID13
%

load(strcat(cd,'/Silenka_Fiedler/DIST01/ID13/stress-strain.mat'),'E','S');
E_01_13 = E;
S_01_13 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST02/ID13/stress-strain.mat'),'E','S');
E_02_13 = E;
S_02_13 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST03/ID13/stress-strain.mat'),'E','S');
E_03_13 = E;
S_03_13 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST04/ID13/stress-strain.mat'),'E','S');
E_04_13 = E;
S_04_13 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST05/ID13/stress-strain.mat'),'E','S');
E_05_13 = E;
S_05_13 = S;
load(strcat(cd,'/Silenka_Fiedler_20/DIST01/ID13/stress-strain.mat'),'E','S');
E_01_20_13 = E;
S_01_20_13 = S;

clear E S;

max01 = 0; imax01 = 1; max02 = 0; imax02 = 1; max03 = 0; imax03 = 1;
max04 = 0; imax04 = 1; max05 = 0; imax05 = 1; max06 = 0; imax06 = 1;
for i=2:1:length(S_01_13)
    if S_01_13(i) > max01, imax01=i; max01=S_01_13(i); end
    if S_02_13(i) > max02, imax02=i; max02=S_02_13(i); end
    if S_03_13(i) > max03, imax03=i; max03=S_03_13(i); end
    if S_04_13(i) > max04, imax04=i; max04=S_04_13(i); end
    if S_05_13(i) > max05, imax05=i; max05=S_05_13(i); end
    if S_01_20_13(i) > max06, imax06=i; max06=S_01_20_13(i); end
end
imax01 = imax01 + 3; imax02 = imax02 + 3; imax03 = imax03 + 3;
imax04 = imax04 + 3; imax05 = imax05 + 3; imax06 = imax06 + 3;

for i=length(S_01_13):-1:imax01, S_01_13(i)=[];E_01_13(i)=[]; end
for i=length(S_02_13):-1:imax02, S_02_13(i)=[];E_02_13(i)=[]; end
for i=length(S_03_13):-1:imax03, S_03_13(i)=[];E_03_13(i)=[]; end
for i=length(S_04_13):-1:imax04, S_04_13(i)=[];E_04_13(i)=[]; end
for i=length(S_05_13):-1:imax05, S_05_13(i)=[];E_05_13(i)=[]; end
for i=length(S_01_20_13):-1:imax06, S_01_20_13(i)=[];E_01_20_13(i)=[]; end

figure(2);
hold on;
set(gca,'FontSize',13);
xlabel('Far-field Strain \epsilon_{22}^o (%)','Position',[0.015/2 -7.4]);
ylabel('\sigma_{22}^o (MPa)');
rectangle('Position',[0.0001,min([max01,max02,max03,max04,max06]),0.035,...
    max([max01,max02,max03,max04,max06])-...
    min([max01,max02,max03,max04,max06])],'Curvature',[0,0],...
    'Linestyle','none','FaceColor',[0.9 0.9 0.9]);
plot(E_01_13,S_01_13,'k-');
plot(E_01_20_13,S_01_20_13,'ro-');
plot(E_05_13,S_05_13,'b-x');
plot(E_02_13,S_02_13,'k-');
plot(E_03_13,S_03_13,'k-');
plot(E_04_13,S_04_13,'k-');
axis([0 0.015 0 80]);
legend('CASES 1,2,3,4','CASE 6','CASE 5','Location','SouthEast');
legend boxoff;
hold off;

% 
% Figure for ID14
% 

load(strcat(cd,'/Silenka_Fiedler/DIST01/ID14/stress-strain.mat'),'E','S');
E_01_14 = E;
S_01_14 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST02/ID14/stress-strain.mat'),'E','S');
E_02_14 = E;
S_02_14 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST03/ID14/stress-strain.mat'),'E','S');
E_03_14 = E;
S_03_14 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST04/ID14/stress-strain.mat'),'E','S');
E_04_14 = E;
S_04_14 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST05/ID14/stress-strain.mat'),'E','S');
E_05_14 = E;
S_05_14 = S;
load(strcat(cd,'/Silenka_Fiedler_20/DIST01/ID14/stress-strain.mat'),'E','S');
E_01_20_14 = E;
S_01_20_14 = S;

clear E S;

max01 = 0; imax01 = 1; max02 = 0; imax02 = 1; max03 = 0; imax03 = 1;
max04 = 0; imax04 = 1; max05 = 0; imax05 = 1; max06 = 0; imax06 = 1;
for i=2:1:length(S_01_14)
    if S_01_14(i) > max01, imax01=i; max01=S_01_14(i); end
    if S_02_14(i) > max02, imax02=i; max02=S_02_14(i); end
    if S_03_14(i) > max03, imax03=i; max03=S_03_14(i); end
    if S_04_14(i) > max04, imax04=i; max04=S_04_14(i); end
    if S_05_14(i) > max05, imax05=i; max05=S_05_14(i); end
    if S_01_20_14(i) > max06, imax06=i; max06=S_01_20_14(i); end
end
imax01 = imax01 + 3; imax02 = imax02 + 3; imax03 = imax03 + 3;
imax04 = imax04 + 3; imax05 = imax05 + 3; imax06 = imax06 + 2;

for i=length(S_01_14):-1:imax01, S_01_14(i)=[];E_01_14(i)=[]; end
for i=length(S_02_14):-1:imax02, S_02_14(i)=[];E_02_14(i)=[]; end
for i=length(S_03_14):-1:imax03, S_03_14(i)=[];E_03_14(i)=[]; end
for i=length(S_04_14):-1:imax04, S_04_14(i)=[];E_04_14(i)=[]; end
for i=length(S_05_14):-1:imax05, S_05_14(i)=[];E_05_14(i)=[]; end
for i=length(S_01_20_14):-1:imax06, S_01_20_14(i)=[];E_01_20_14(i)=[]; end

figure(3);
hold on;
set(gca,'FontSize',13);
xlabel('Far-field Strain \epsilon_{12}^o (%)','Position',[0.015/2 -4.6]);
ylabel('\sigma_{12}^o (MPa)');
rectangle('Position',[0.0001,min([max01,max02,max03,max04,max06]),0.015,...
    max([max01,max02,max03,max04,max06])-...
    min([max01,max02,max03,max04,max06])],'Curvature',[0,0],...
    'Linestyle','none','FaceColor',[0.9 0.9 0.9]);
plot(E_01_14,S_01_14,'k-');
plot(E_01_20_14,S_01_20_14,'ro-');
plot(E_05_14,S_05_14,'b-x');
plot(E_02_14,S_02_14,'k-');
plot(E_03_14,S_03_14,'k-');
plot(E_04_14,S_04_14,'k-');
axis([0 0.015 0 50]);
legend('CASES 1,2,3,4','CASE 6','CASE 5','Location','SouthEast');
legend boxoff;
hold off;


%
% Figure for ID16
%

load(strcat(cd,'/Silenka_Fiedler/DIST01/ID16/stress-strain.mat'),'E','S');
E_01_16 = E;
S_01_16 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST02/ID16/stress-strain.mat'),'E','S');
E_02_16 = E;
S_02_16 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST03/ID16/stress-strain.mat'),'E','S');
E_03_16 = E;
S_03_16 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST04/ID16/stress-strain.mat'),'E','S');
E_04_16 = E;
S_04_16 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST05/ID16/stress-strain.mat'),'E','S');
E_05_16 = E;
S_05_16 = S;
load(strcat(cd,'/Silenka_Fiedler_20/DIST01/ID16/stress-strain.mat'),'E','S');
E_01_20_16 = E;
S_01_20_16 = S;

clear E S;

max01 = 0; imax01 = 1; max02 = 0; imax02 = 1; max03 = 0; imax03 = 1;
max04 = 0; imax04 = 1; max05 = 0; imax05 = 1; max06 = 0; imax06 = 1;
for i=2:1:length(S_01_16)
    if S_01_16(i) > max01, imax01=i; max01=S_01_16(i); end
    if S_02_16(i) > max02, imax02=i; max02=S_02_16(i); end
    if S_03_16(i) > max03, imax03=i; max03=S_03_16(i); end
    if S_04_16(i) > max04, imax04=i; max04=S_04_16(i); end
    if S_05_16(i) > max05, imax05=i; max05=S_05_16(i); end
    if S_01_20_16(i) > max06, imax06=i; max06=S_01_20_16(i); end
end
imax01 = imax01 + 5; imax02 = imax02 + 2; imax03 = imax03 + 3;
imax04 = imax04 + 3; imax05 = imax05 + 3; imax06 = imax06 + 3;

for i=length(S_01_16):-1:imax01, S_01_16(i)=[];E_01_16(i)=[]; end
for i=length(S_02_16):-1:imax02, S_02_16(i)=[];E_02_16(i)=[]; end
for i=length(S_03_16):-1:imax03, S_03_16(i)=[];E_03_16(i)=[]; end
for i=length(S_04_16):-1:imax04, S_04_16(i)=[];E_04_16(i)=[]; end
for i=length(S_05_16):-1:imax05, S_05_16(i)=[];E_05_16(i)=[]; end
for i=length(S_01_20_16):-1:imax06, S_01_20_16(i)=[];E_01_20_16(i)=[]; end

figure(4);
hold on;
set(gca,'FontSize',13);
xlabel('Far-field Strain \epsilon_{23}^o (%)','Position',[0.015/2 -4.6]);
ylabel('\sigma_{23}^o (MPa)');
rectangle('Position',[0.0001,min([max01,max02,max03,max04,max06]),0.015,...
    max([max01,max02,max03,max04,max06])-...
    min([max01,max02,max03,max04,max06])],'Curvature',[0,0],...
    'Linestyle','none','FaceColor',[0.9 0.9 0.9]);
plot(E_01_16,S_01_16,'k-');
plot(E_01_20_16,S_01_20_16,'ro-');
plot(E_05_16,S_05_16,'b-x');
plot(E_02_16,S_02_16,'k-');
plot(E_03_16,S_03_16,'k-');
plot(E_04_16,S_04_16,'k-');
axis([0 0.015 0 50]);
legend('CASES 1,2,3,4','CASE 6','CASE 5','Location','SouthEast');
legend boxoff;
hold off;


% 
% Figure for ID23
% 

load(strcat(cd,'/Silenka_Fiedler/DIST01/ID23/stress-strain.mat'),'E','S');
E_01_23 = E;
S_01_23 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST02/ID23/stress-strain.mat'),'E','S');
E_02_23 = E;
S_02_23 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST03/ID23/stress-strain.mat'),'E','S');
E_03_23 = E;
S_03_23 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST04/ID23/stress-strain.mat'),'E','S');
E_04_23 = E;
S_04_23 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST05/ID23/stress-strain.mat'),'E','S');
E_05_23 = E;
S_05_23 = S;
load(strcat(cd,'/Silenka_Fiedler_20/DIST01/ID23/stress-strain.mat'),'E','S');
E_01_20_23 = E;
S_01_20_23 = S;

clear E S;

max01 = 0; imax01 = 1; max02 = 0; imax02 = 1; max03 = 0; imax03 = 1;
max04 = 0; imax04 = 1; max05 = 0; imax05 = 1; max06 = 0; imax06 = 1;
for i=2:1:length(S_01_23)
    if S_01_23(i) < max01, imax01=i; max01=S_01_23(i); end
    if S_02_23(i) < max02, imax02=i; max02=S_02_23(i); end
    if S_03_23(i) < max03, imax03=i; max03=S_03_23(i); end
    if S_04_23(i) < max04, imax04=i; max04=S_04_23(i); end
    if S_05_23(i) < max05, imax05=i; max05=S_05_23(i); end
    if S_01_20_23(i) < max06, imax06=i; max06=S_01_20_23(i); end
end
imax01 = imax01 + 4; imax02 = imax02 + 4; imax03 = imax03 + 4;
imax04 = imax04 + 3; imax05 = imax05 + 6; imax06 = imax06 + 4;

for i=length(S_01_23):-1:imax01, S_01_23(i)=[];E_01_23(i)=[]; end
for i=length(S_02_23):-1:imax02, S_02_23(i)=[];E_02_23(i)=[]; end
for i=length(S_03_23):-1:imax03, S_03_23(i)=[];E_03_23(i)=[]; end
for i=length(S_04_23):-1:imax04, S_04_23(i)=[];E_04_23(i)=[]; end
for i=length(S_05_23):-1:imax05, S_05_23(i)=[];E_05_23(i)=[]; end
for i=length(S_01_20_23):-1:imax06, S_01_20_23(i)=[];E_01_20_23(i)=[]; end

figure(5);
hold on;
set(gca,'FontSize',13);
xlabel('Far-field Strain |\epsilon_{22}^o| (%)','Position',[0.02/2 -12]);
ylabel('|\sigma _{22}^o| (MPa)');
rectangle('Position',[0.0001,min(abs([max01,max02,max03,max04,max06])),0.02,...
    max(abs([max01,max02,max03,max04,max06]))-...
    min(abs([max01,max02,max03,max04,max06]))],'Curvature',[0,0],...
    'Linestyle','none','FaceColor',[0.9 0.9 0.9]);
plot(-E_01_23,-S_01_23,'k-');
plot(-E_01_20_23,-S_01_20_23,'ro-');
plot(-E_05_23,-S_05_23,'b-x');
plot(-E_03_23,-S_03_23,'k-');
plot(-E_04_23,-S_04_23,'k-');
plot(-E_02_23,-S_02_23,'k-');
axis ([0 0.02 0 130]);
legend('CASES 1,2,3,4','CASE 6','CASE 5','Location','SouthEast');
legend boxoff;
hold off;

 
%
% Figure for ID236
%

load(strcat(cd,'/Silenka_Fiedler/DIST01/ID236/stress-strain.mat'),'E','S');
E_01_236 = E;
S_01_236 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST02/ID236/stress-strain.mat'),'E','S');
E_02_236 = E;
S_02_236 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST03/ID236/stress-strain.mat'),'E','S');
E_03_236 = E;
S_03_236 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST04/ID236/stress-strain.mat'),'E','S');
E_04_236 = E;
S_04_236 = S;
load(strcat(cd,'/Silenka_Fiedler/DIST05/ID236/stress-strain.mat'),'E','S');
E_05_236 = E;
S_05_236 = S;
load(strcat(cd,'/Silenka_Fiedler_20/DIST01/ID236/stress-strain.mat'),'E','S');
E_01_20_236 = E;
S_01_20_236 = S;

clear E S;

max01 = 0; imax01 = 1; max02 = 0; imax02 = 1; max03 = 0; imax03 = 1;
max04 = 0; imax04 = 1; max05 = 0; imax05 = 1; max06 = 0; imax06 = 1;
for i=2:1:length(S_01_236)
    if S_01_236(i) > max01, imax01=i; max01=S_01_236(i); end
    if S_02_236(i) > max02, imax02=i; max02=S_02_236(i); end
    if S_03_236(i) > max03, imax03=i; max03=S_03_236(i); end
    if S_04_236(i) > max04, imax04=i; max04=S_04_236(i); end
    if S_05_236(i) > max05, imax05=i; max05=S_05_236(i); end
    if S_01_20_236(i) > max06, imax06=i; max06=S_01_20_236(i); end
end
imax01 = imax01 + 3; imax02 = imax02 + 5; imax03 = imax03 + 6;
imax04 = imax04 + 6; imax05 = imax05 + 6; imax06 = imax06 + 2;

for i=length(S_01_236):-1:imax01, S_01_236(i)=[];E_01_236(i)=[]; end
for i=length(S_02_236):-1:imax02, S_02_236(i)=[];E_02_236(i)=[]; end
for i=length(S_03_236):-1:imax03, S_03_236(i)=[];E_03_236(i)=[]; end
for i=length(S_04_236):-1:imax04, S_04_236(i)=[];E_04_236(i)=[]; end
for i=length(S_05_236):-1:imax05, S_05_236(i)=[];E_05_236(i)=[]; end
for i=length(S_01_20_236):-1:imax06, S_01_20_236(i)=[];E_01_20_236(i)=[]; end

figure(6);
hold on;
set(gca,'FontSize',13);
xlabel('Far-field Strain |\epsilon_{22}^o| (%)','Position',[0.01/2 -4.6]);
ylabel('|\sigma_{22}^o| (MPa)');
rectangle('Position',[0.0001,min([max01,max02,max03,max04,max06]),0.015,...
    max([max01,max02,max03,max04,max06])-...
    min([max01,max02,max03,max04,max06])],'Curvature',[0,0],...
    'Linestyle','none','FaceColor',[0.9 0.9 0.9]);
plot(E_01_236,S_01_236,'k-');
plot(E_01_20_236,S_01_20_236,'ro-');
plot(E_05_236,S_05_236,'b-x');
plot(E_02_236,S_02_236,'k-');
plot(E_03_236,S_03_236,'k-');
plot(E_04_236,S_04_236,'k-');
axis([0 0.01 0 50]);
legend('CASES 1,2,3,4','CASE 6','CASE 5','Location','SouthEast');
legend boxoff;
hold off;


 

saveas(1,'RVE_11_se.eps','psc2');
saveas(2,'RVE_13_se.eps','psc2');
saveas(3,'RVE_14_se.eps','psc2');
saveas(4,'RVE_16_se.eps','psc2');
saveas(5,'RVE_23_se.eps','psc2');
saveas(6,'RVE_236_se.eps','psc2');


