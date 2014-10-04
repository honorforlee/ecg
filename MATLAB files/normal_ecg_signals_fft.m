%% Normal ECG Spectrums
clear all; close all;

ecg_normal_1 = load('ecg_normal_1.mat');
ecg_normal_2 = load('ecg_normal_2.mat');
ecg_normal_3 = load('ecg_normal_3.mat');

ecg_norm_1 = ecg_normal_1.ecg; 
ecg_norm_2 = ecg_normal_2.ecg;
ecg_norm_3 = ecg_normal_3.ecg;

Fs_1 = ecg_normal_1.Fs; 
Fs_2 = ecg_normal_2.Fs;
Fs_3 = ecg_normal_3.Fs;

%% Normal ECG 1 Spectrum

s = 15; % Use s seconds of samples
N1 = Fs_1*s;
F1 = ([0:N1-1]/N1-0.5)*Fs_1;
fft_ecg_1 = abs(fft(ecg_norm_1([1:N1]),N1)).^2;
fft_ecg_1 = fftshift(fft_ecg_1);

figure
subplot 311
plot(F1,fft_ecg_1);
set(gca,'FontSize',16)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG 1 spectrum','FontSize',16);
xlim([-30 30]);

%% Normal ECG 2 Spectrum

s = 15;
N2 = Fs_2*s;
F2 = ([0:N2-1]/N2-0.5)*Fs_2;
fft_ecg_2 = abs(fft(ecg_norm_2([1:N2]),N2)).^2;
fft_ecg_2 = fftshift(fft_ecg_2);

subplot 312
plot(F2,fft_ecg_2);
set(gca,'FontSize',16)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG 2 spectrum','FontSize',16);
xlim([-30 30]);

%% Normal ECG 3 Spectrum

s = 15;
N3 = Fs_3*s;
F3 = ([0:N3-1]/N3-0.5)*Fs_3;
fft_ecg_3 = abs(fft(ecg_norm_3([1:N3]),N3)).^2;
fft_ecg_3 = fftshift(fft_ecg_3);

subplot 313
plot(F3,fft_ecg_3);
set(gca,'FontSize',16)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',16);
title('ECG 3 spectrum','FontSize',16);
xlim([-30 30]);