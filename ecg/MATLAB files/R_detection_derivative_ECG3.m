%% ECG 3 R Peaks Detection / Method of the derivative

clear all; close all;

ecg_normal_3 = load('ecg_normal_3.mat');
ecg_norm_3 = ecg_normal_3.ecg;
Fs_3 = ecg_normal_3.Fs;

N = length(ecg_norm_3);
t = [0:N-1]/Fs_3;

ecg_deriv = diff(ecg_norm_3);

threshold_deriv = (max(ecg_deriv([1:10000]))-mean(ecg_deriv([1:10000])))/3;
threshold = (max(ecg_norm_3([1:10000]))-mean(ecg_norm_3([1:10000])))/2;

t2 = [0:N-2]/Fs_3;
figure(1)
subplot 211
plot(t2, ecg_deriv,[0 N-2],[0 0],'g');
set(gca,'FontSize',12)
title('ECG 2 derived','FontSize',14);
xlabel('Time (s)','FontSize',14);
ylabel('ECG signal amplitude','FontSize',14);

hold on;

plot([0 max(t)], threshold_deriv*[1 1],'r--')
plot([0 max(t)], -threshold_deriv*[1 1],'r--')
legend('ECG 1 derived','y = 0','Threshold');
xlim([0 10])

[peak_max,peak_loc_max] = findpeaks(ecg_deriv,'MinPeakHeight',threshold_deriv);
ecg_deriv_inv = -ecg_deriv; 
[peak_min,peak_loc_min] = findpeaks(ecg_deriv_inv,'MinPeakHeight',threshold_deriv);

hold on;
plot(peak_loc_max/Fs_3,ecg_deriv(peak_loc_max),'rv','MarkerFaceColor','r');
plot(peak_loc_min/Fs_3,-ecg_deriv_inv(peak_loc_min),'rv','MarkerFaceColor','r');

for i = 1:min(length(peak_loc_max),length(peak_loc_min))
    for j=peak_loc_max(i):peak_loc_min(i)
        R_loc(j+1) = 0.5*abs(sign(ecg_deriv(j+1))-sign(ecg_deriv(j)));
    end
end

R_peak_loc = find(R_loc == 1);
R_peak_amp = ecg_norm_3(R_peak_loc);

for i=1:length(R_peak_amp) 
    if R_peak_amp(i) < threshold
       R_peak_loc(i) = 0;
    end
end

R_peak_loc = R_peak_loc(find(R_peak_loc > 0));

subplot 212
plot(t,ecg_norm_3);
set(gca,'FontSize',12)
title('ECG 2','FontSize',14);
xlabel('Time (s)','FontSize',14);
ylabel('ECG signal amplitude','FontSize',14);
xlim([0 10]);

hold on;
plot(R_peak_loc/Fs_3,ecg_norm_3(R_peak_loc),'rv','MarkerFaceColor','r');
legend('ECG signal','R waves');