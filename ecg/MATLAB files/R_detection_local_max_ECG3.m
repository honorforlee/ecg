%% ECG 3 R Peaks Detection / Method of local maxima
clear all; close all;

ecg_normal_3 = load('ecg_normal_3.mat');
ecg_norm_3 = ecg_normal_3.ecg;
Fs_3 = ecg_normal_3.Fs;

N = length(ecg_norm_3);
t = [0:N-1]/Fs_3;

figure(1)
subplot 211
plot(t,ecg_norm_3);
set(gca,'FontSize',14)
title('ECG 1 using findpeaks','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);

threshold = (max(ecg_norm_3([1:10000]))-mean(ecg_norm_3([1:10000])))/3 + mean(ecg_norm_3);

%% Using findpeaks

[R_peak_findpeaks,R_peak_loc_findpeaks] = findpeaks(ecg_norm_3,'MinPeakHeight',threshold); 
% peak is the amplitude and peak_loc is the location of the R waves

hold on;
plot(R_peak_loc_findpeaks/Fs_3,ecg_norm_3(R_peak_loc_findpeaks),'rv','MarkerFaceColor','r','MarkerSize',8);
plot([0 max(t)], threshold*[1 1],'r--')
legend('ECG signal','R waves','Threshold');

%% Using self-made function

R_peak_loc_self = find(ecg_norm_3 > threshold); % find ecg samples > threshold computed
S = length(R_peak_loc_self);

for i=1:S-1
    if R_peak_loc_self(i+1)-R_peak_loc_self(i)>50 % if the distance between two consecutive peaks is greater than 50 samples, R_peak_loc_self is the location of an R wave
        R_loc_detected(i) = R_peak_loc_self(i)-2;
    end
end

R_peak_loc_self_final = R_loc_detected(find(R_loc_detected > 0)); 

subplot 212
plot(t,ecg_norm_3);
set(gca,'FontSize',14)
title('ECG 1 using selfmade function','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);
hold on;
plot(R_peak_loc_self_final/Fs_3,ecg_norm_3(R_peak_loc_self_final),'rv','MarkerFaceColor','r','MarkerSize',8);
plot([0 max(t)], threshold*[1 1],'r--')
legend('ECG signal','R waves','Threshold');

%% Automatic identification of cardiac pathologies

% Tachycardia and Bradycardia detection as well as Heart Rate Variability
% detection are used here because of R waves detected needed to compute
% the heart rate.
 
%% Tachycardia/Bradycardia detection

for i = 1:length(R_peak_loc_findpeaks)-1
    d(i) = (R_peak_loc_findpeaks(i+1)-R_peak_loc_findpeaks(i))/Fs_3; % Time distance between two consecutive R waves
    d_mean = mean(d(i));
end

for i=1:length(d)
    
    bpm(i) = 60/d(i); % beats per minute
    bpm_mean = mean(bpm(i));
    
    if (bpm > 110)
        state = 'Tachycardia';
    elseif (bpm < 60);
        state = 'Bradycardia';
    else
        state = 'Normal';
    end
end

%% Heart Rate Variability detection

for i = 1:length(R_peak_loc_findpeaks)-2
    for j=R_peak_loc_findpeaks(i):R_peak_loc_findpeaks(i+1)
        v(j) = d(i) + ((d(i+1)-d(i))*(j-R_peak_loc_findpeaks(i)))/(R_peak_loc_findpeaks(i+1)-R_peak_loc_findpeaks(i)); % v(t) expression
    end
end

N1 = 1000000;
F=([0:N1-1]/N1-0.5)*Fs_3;
v_fft = abs(fft(v(1:20*Fs_3),N1)).^2;
v_fft = fftshift(v_fft);

v_ulf_fft = v_fft(1:500009); % Ultra low frequency part
v_vlf_fft = v_fft(500009:500112);  % Very low frequency part
v_lf_fft = v_fft(500112:500418);  % Low frequency part
v_hf_fft = v_fft(500418:501112);  % High frequency part

figure(2)
plot(F(1:500009),v_ulf_fft);
hold on;
plot(F(500009:500112),v_vlf_fft,'g');
plot(F(500112:500418),v_lf_fft,'r');
plot(F(500418:501112),v_hf_fft,'m');
set(gca,'FontSize',16)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('Spectrum amplitude','FontSize',16);
title('HRV spectrum for ECG 1','FontSize',16);
xlim([0 0.4]);
legend('0-0.003 Hz','0.003-0.04 Hz','0.04-0.15 Hz','0.15-0.4 Hz');