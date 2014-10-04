%% Normal ECG Signals
clear all; close all;

ecg_normal_1 = load('ecg_normal_1.mat'); % Loading of ecg file
ecg_normal_2 = load('ecg_normal_2.mat');
ecg_normal_3 = load('ecg_normal_3.mat');

ecg_norm_1 = ecg_normal_1.ecg; % Loading samples of ecg_normal_1
ecg_norm_2 = ecg_normal_2.ecg;
ecg_norm_3 = ecg_normal_3.ecg;

Fs_1 = ecg_normal_1.Fs; % Loading sample frequency of ecg_normal_1
Fs_2 = ecg_normal_2.Fs;
Fs_3 = ecg_normal_3.Fs;

%% Normal ECG 1 Signal

% Plot with ecg_display function

figure(1)
interval = 0:1/Fs_1:4; % Interval between 0 s and 4 s
ecg_display(interval,ecg_norm_1); % Please see ecg_display.m for ecg_display function
set(gca,'FontSize',16)
title('ECG 1','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);

%------------------------

N = length(ecg_norm_1); % Number of ecg samples
t = [0:N-1]/Fs_1; % Time interval set with samples of the ecg signal and its sample frequency

figure(2)
subplot 311
plot(t,ecg_norm_1); 
set(gca,'FontSize',16)
title('ECG 1','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]); % limit used to limit display 

hold on;

plot(0.2833,ecg_norm_1(103),'ro','markersize',12);
text(0.25, 95,'P', 'Color', 'r','fontsize',16,'fontweight','b');

plot(0.4167,ecg_norm_1(151),'ro','markersize',12);
text(0.3, -140, 'Q', 'Color', 'r','fontsize',16,'fontweight','b');

plot(0.45,ecg_norm_1(163),'ro','markersize',12);
text(0.5, 332, 'R', 'Color', 'r','fontsize',16,'fontweight','b');

plot(0.4667,ecg_norm_1(169),'ro','markersize',12);
text(0.52, -130, 'S', 'Color', 'r','fontsize',16,'fontweight','b');

plot(0.675,ecg_norm_1(244),'ro','markersize',12);
text(0.73, 70 , 'T', 'Color', 'r','fontsize',16,'fontweight','b');

hold off;


%% Normal ECG 2 Signal

N = length(ecg_norm_2);
t = [0:N-1]/Fs_2;

subplot 312
plot(t,ecg_norm_2);
set(gca,'FontSize',16)
%set(gca,'XTick',0:0.1:4)
title('ECG 2','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);

hold on;

plot(1.378,ecg_norm_2(497),'ro','markersize',12);
text(1.37, 210,'P', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.511,ecg_norm_2(545),'ro','markersize',12);
text(1.42, -230, 'Q', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.55,ecg_norm_2(559),'ro','markersize',12);
text(1.62, 760, 'R', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.625,ecg_norm_2(585),'ro','markersize',12);
text(1.61, -250, 'S', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.85,ecg_norm_2(667),'ro','markersize',12);
text(1.84, 290 , 'T', 'Color', 'r','fontsize',16,'fontweight','b');

hold off;

%% Normal ECG 3 Signal

N = length(ecg_norm_3);
t = [0:N-1]/Fs_3;

subplot 313
plot(t,ecg_norm_3);
set(gca,'FontSize',16)
title('ECG 3','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);

hold on;

plot(1.567,ecg_norm_3(565),'ro','markersize',12);
text(1.558, 200,'P', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.7,ecg_norm_3(612),'ro','markersize',12);
text(1.65, -290, 'Q', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.725,ecg_norm_3(622),'ro','markersize',12);
text(1.715, 750, 'R', 'Color', 'r','fontsize',16,'fontweight','b');

plot(1.747,ecg_norm_3(630),'ro','markersize',12);
text(1.78, -290, 'S', 'Color', 'r','fontsize',16,'fontweight','b');

plot(2.083,ecg_norm_3(751),'ro','markersize',12);
text(2.07, 200 , 'T', 'Color', 'r','fontsize',16,'fontweight','b');

hold off;