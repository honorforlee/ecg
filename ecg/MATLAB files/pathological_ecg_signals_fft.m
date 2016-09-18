%% Pathological ECG Spectrums
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

N = 1000;
N1 = N*10;
F_AF = ([0:N1-1]/N1-0.5)*Fs_AF;

fft_ecg_AF = abs(fft(ecg_abnorm_AF([300*Fs_AF:300*Fs_AF+N]),N1));
fft_ecg_AF = fftshift(fft_ecg_AF);

figure
subplot 221
plot(F_AF,fft_ecg_AF);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG AF spectrum','FontSize',16);
xlim([-50 50])

%% VF ECG

N = 1000;
N1 = N*10;
F_VF = ([0:N1-1]/N1-0.5)*Fs_VF;

fft_ecg_VF = abs(fft(ecg_abnorm_VF([216*Fs_VF:216*Fs_VF+N]),N1));
fft_ecg_VF = fftshift(fft_ecg_VF);

subplot 222
plot(F_VF,fft_ecg_VF);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG VF spectrum','FontSize',16);
xlim([-50 50])

%% SSS ECG

N = 1000;
N1 = N*10;
F_SSS = ([0:N1-1]/N1-0.5)*Fs_SSS;

fft_ecg_SSS = abs(fft(ecg_abnorm_SSS([1:N]),N1));
fft_ecg_SSS = fftshift(fft_ecg_SSS);

subplot 223
plot(F_SSS,fft_ecg_SSS);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG SSS spectrum','FontSize',16);
xlim([-50 50])

%% PVC ECG

N = 1000;
N1 = N*10;
F_PVC = ([0:N1-1]/N1-0.5)*Fs_PVC;

fft_ecg_PVC = abs(fft(ecg_abnorm_PVC([1:N]),N1));
fft_ecg_PVC = fftshift(fft_ecg_PVC);

subplot 224
plot(F_PVC,fft_ecg_PVC);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG PVC spectrum','FontSize',16);
xlim([-50 50])