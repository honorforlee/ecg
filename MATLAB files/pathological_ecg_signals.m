%% Pathological ECG signals
clear all; close all;

ecg_AF = load('ecg_AF.mat'); 
ecg_VF = load('ecg_VF.mat');
ecg_SSS = load('ecg_SSS.mat');
ecg_PVC = load('ecg_PVC.mat');

ecg_abnorm_AF = ecg_AF.ecg;
ecg_abnorm_VF = ecg_VF.ecg;
ecg_abnorm_SSS = ecg_SSS.ecg;
ecg_abnorm_PVC = ecg_PVC.ecg;

Fs_AF = ecg_AF.Fs;
Fs_VF = ecg_VF.Fs;
Fs_SSS = ecg_SSS.Fs;
Fs_PVC = ecg_PVC.Fs;

%% AF ECG

N_AF = length(ecg_abnorm_AF);
t = [0:N_AF-1]/Fs_AF;

figure
subplot 411
plot(t,ecg_abnorm_AF);
set(gca,'FontSize',16)
title('ECG AF','FontSize',16); 
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',14);
xlim([300 304]);

hold on;

plot(302.285,ecg_abnorm_AF(108014),'ro','markersize',12);
text(302.2, -250, 'Q', 'Color', 'r','fontsize',16,'fontweight','b');

plot(302.316,ecg_abnorm_AF(108085),'ro','markersize',12);
text(302.37, 300, 'R', 'Color', 'r','fontsize',16,'fontweight','b');

plot(302.332,ecg_abnorm_AF(108869),'ro','markersize',12);
text(302.37, -250, 'S', 'Color', 'r','fontsize',16,'fontweight','b');

%% VF ECG

N_VF = length(ecg_abnorm_VF);
t = [0:N_VF-1]/Fs_VF;

subplot 412
plot(t,ecg_abnorm_VF);
set(gca,'FontSize',16)
title('ECG VF','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',14);
xlim([216 230]);

%% SSS ECG

N_SSS = length(ecg_abnorm_SSS);
t = [0:N_SSS-1]/Fs_SSS;

subplot 413
plot(t,ecg_abnorm_SSS);
set(gca,'FontSize',16)
title('ECG SSS','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',14);
xlim([0 6]);

hold on;

plot(2.781,ecg_abnorm_SSS(1002),'ro','markersize',12);
text(2.88, 1135, 'R', 'Color', 'r','fontsize',16,'fontweight','b');

plot(2.803,ecg_abnorm_SSS(1010),'ro','markersize',12);
text(2.65, 908, 'S', 'Color', 'r','fontsize',16,'fontweight','b');

plot(3.078,ecg_abnorm_SSS(1109),'ro','markersize',12);
text(3.17, 1060, 'T', 'Color', 'r','fontsize',16,'fontweight','b');

%% PVC ECG

N_PVC = length(ecg_abnorm_PVC);
t = [0:N_PVC-1]/Fs_PVC;

subplot 414
plot(t,ecg_abnorm_PVC);
set(gca,'FontSize',16)
title('ECG PVC','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',14);
xlim([0 4]);

hold on;

plot(1.028,ecg_abnorm_PVC(371),'ro','markersize',12);
text(1.1, 260, 'R', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.069,ecg_abnorm_PVC(385),'ro','markersize',12);
text(1.12, -260, 'S', 'Color', 'r','fontsize',16,'fontweight','b');

hold off;