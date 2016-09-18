%% ECG 1 R Peaks Detection / Method of the derivative

clear all; close all;

ecg_normal_1 = load('ecg_normal_1.mat');
ecg_norm_1 = ecg_normal_1.ecg;
Fs_1 = ecg_normal_1.Fs;

N = length(ecg_norm_1);
t = [0:N-1]/Fs_1;

ecg_deriv = diff(ecg_norm_1);

threshold_deriv = (max(ecg_deriv([1:10000]))-mean(ecg_deriv([1:10000])))/2; % threshold for derived ecg
threshold = (max(ecg_norm_1([1:10000]))-mean(ecg_norm_1([1:10000])))/2; % threshold for ecg

t2 = [0:N-2]/Fs_1;
figure(1)
subplot 211
plot(t2, ecg_deriv,[0 N-2],[0 0],'g');
set(gca,'FontSize',14)
title('ECG 1 derived','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);

hold on;

plot([0 max(t)], threshold_deriv*[1 1],'r--')
plot([0 max(t)], -threshold_deriv*[1 1],'r--')
legend('ECG 1 derived','y = 0','Threshold');
xlim([0 4])

[peak_max,peak_loc_max] = findpeaks(ecg_deriv,'MinPeakHeight',threshold_deriv); % find positive peaks for derived ecg
ecg_deriv_inv = -ecg_deriv; 
[peak_min,peak_loc_min] = findpeaks(ecg_deriv_inv,'MinPeakHeight',threshold_deriv); % find negative peaks for derived ecg

hold on;
plot(peak_loc_max/Fs_1,ecg_deriv(peak_loc_max),'rv','MarkerFaceColor','r','MarkerSize',8);
plot(peak_loc_min/Fs_1,-ecg_deriv_inv(peak_loc_min),'rv','MarkerFaceColor','r','MarkerSize',8);

for i = 1:min(length(peak_loc_max),length(peak_loc_min))
    for j = peak_loc_max(i):peak_loc_min(i)
        R_loc(j+1) = 0.5*abs(sign(ecg_deriv(j+1))-sign(ecg_deriv(j)));
    end
end

R_peak_loc = find(R_loc == 1);
R_peak_amp = ecg_norm_1(R_peak_loc);

for i = 1:length(R_peak_amp) % if peaks are below threshold, they are deleted
    if R_peak_amp(i) < threshold
       R_peak_loc(i) = 0;
    end
end

R_peak_loc = R_peak_loc(find(R_peak_loc > 0));

subplot 212
plot(t,ecg_norm_1);
set(gca,'FontSize',14)
title('ECG 1','FontSize',16);
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
xlim([0 4]);

hold on;
plot(R_peak_loc/Fs_1,ecg_norm_1(R_peak_loc),'rv','MarkerFaceColor','r','MarkerSize',8);
legend('ECG signal','R waves');