clear all; close all;

%% Pan Tompkins algorithm 

ecg_normal_1 = load('ecg_normal_1');
ecg_norm_1 = ecg_normal_1.ecg;
Fs = ecg_normal_1.Fs;

N = length(ecg_norm_1);
t = [0:N-1]/Fs;

%% Filtering

%% Low-pass filter
l_b = [1 zeros(1,5) -2 zeros(1,5) 1];
l_a = [1,-2,1];

% Low-pass filter frequency response
[H_l,F] = freqz(l_b,l_a,1024,Fs);

figure(1)
subplot 411
plot(F,db(H_l)); % Amplitude
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('H_low filter amplitude |H_{low}(f)|','FontSize',9);
title('Low-pass filter frequency response','FontSize',16);

subplot 412
plot(F,unwrap(angle(H_l)));
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('H_low filter phase arg(H_{low}(f))','FontSize',9); % Phase

ecg_low_filtered = filter(l_b,l_a,ecg_norm_1); % ECG signal after low-pass filtering

subplot 413
plot(t,ecg_low_filtered);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16); 
ylabel('ECG signal amplitude','FontSize',9);
title('Low-pass filtered ECG 1','FontSize',16);
xlim([0 4]);

s = 15; 
N1 = Fs * s;
F=([0:N1-1]/N1-0.5)*Fs;
fft_low_ecg = abs(fft(ecg_low_filtered([1:N1]),N1)).^2;
fft_low_ecg = fftshift(fft_low_ecg);

subplot 414
plot(F,fft_low_ecg);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',9);
title('Low-pass filtered ECG 1 spectrum','FontSize',16);
xlim([-30 30]);

%% High-pass filter
h_b = [-1 zeros(1,15) 32,-32 zeros(1,14) 1];
h_a = [1,1];

% High-pass filter frequency response
[H_h,F] = freqz(h_b,h_a,1024,Fs);

figure(2)
subplot 411
plot(F,db(H_h)); % Amplitude
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('H_high filter amplitude |H_{high}(f)|','FontSize',9);
title('High-pass filter frequency response','FontSize',16);

subplot 412
plot(F,unwrap(angle(H_h))); 
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('H_high filter phase arg(H_{high}(f))','FontSize',9); % Phase

ecg_high_filtered = filter(h_b,h_a,ecg_low_filtered); % ECG signal after high-pass filtering

subplot 413
plot(t,ecg_high_filtered);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16); 
ylabel('ECG signal amplitude','FontSize',9);
title('High-pass filtered ECG 1','FontSize',16);
xlim([0 4]);

F=([0:N1-1]/N1-0.5)*Fs;
fft_high_ecg = abs(fft(ecg_high_filtered([1:N1]),N1)).^2;
fft_high_ecg = fftshift(fft_high_ecg);

subplot 414
plot(F,fft_high_ecg);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',9);
title('High-pass filtered ECG 1 spectrum','FontSize',16);
xlim([-30 30]);

ecg_filtered = ecg_high_filtered; % ECG signal after band-pass filtering

%% Deriving

% Differentiation filter
b_diff = [Fs/8,2*Fs/8,0,-2*Fs/8,-Fs/8];
a_diff = [1];

% Frequency response
[H_diff,F] = freqz(b_diff,a_diff,1024,Fs,'whole');

figure(3)
subplot 411
plot(F,db(H_diff)); % Amplitude
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('H_diff filter amplitude |H_{diff}(f)|','FontSize',9);
title('Differentiation filter frequency response','FontSize',16);

subplot 412
plot(F,unwrap(angle(H_diff))); 
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('H_diff filter phase arg(H_{diff}(f))','FontSize',9);

ecg_derived = filter(b_diff,a_diff,ecg_filtered); % ECG signal derived

subplot 413
plot(t,ecg_derived);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16); 
ylabel('ECG signal amplitude','FontSize',9);
title('Derived ECG 1','FontSize',16);
xlim([0 4]);

F=([0:N1-1]/N1-0.5)*Fs;
fft_derived_ecg = abs(fft(ecg_derived([1:N1]),N1)).^2;
fft_derived_ecg = fftshift(fft_derived_ecg);

subplot 414
plot(F,fft_derived_ecg);
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('ECG spectrum amplitude','FontSize',9);
title('Derived ECG 1 spectrum','FontSize',16);
xlim([-30 30]);

%% Squaring

ecg_squared = abs(ecg_derived).^2; % ECG signal squared 

figure(4)
plot(t,ecg_squared);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
title('Squared ECG 1','FontSize',16);
xlim([0 4])

%% Moving-window integration

Win_width = 40;

% Moving window integration filter
b_MWI = ones(1,Win_width);
a_MWI = [Win_width];

ecg_MWI = filter(b_MWI,a_MWI,ecg_squared); % ECG after the Moving Window Integration

figure(5)
plot(t,ecg_MWI);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
title('ECG 1 after the Moving Window Integration','FontSize',16);
xlim([0 4])

%% Thresholding

ecg_MWI_deriv = diff(ecg_MWI);
t2 = [0:N-2]/Fs;

figure(6)
subplot 211
plot(t2, ecg_MWI_deriv);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
title('Derivate of ECG 1 after Moving-Window Integration','FontSize',16);
xlim([0 4])

max_vect = [max(abs(ecg_MWI_deriv([1*Fs:10*Fs]))) max(abs(ecg_MWI_deriv([40*Fs:50*Fs]))) max(abs(ecg_MWI_deriv([90*Fs:100*Fs]))) ]; % 3 intervals used to obtain an estimate of the maximum of the signal
threshold = mean(max_vect)/1.8;

[peak_max,peak_loc_max] = findpeaks(ecg_MWI_deriv,'MinPeakHeight',threshold);
ecg_MWI_deriv_inv = -ecg_MWI_deriv;
[peak_min,peak_loc_min] = findpeaks(ecg_MWI_deriv_inv,'MinPeakHeight',threshold);

ecg_MWI_slope_begin = peak_loc_max;
ecg_MWI_slope_end = peak_loc_min;

for i = 1:min(length(ecg_MWI_slope_begin),length(ecg_MWI_slope_end))
    for j = ecg_MWI_slope_begin(i):ecg_MWI_slope_end(i)
        if (ecg_MWI(j+1) > ecg_MWI(j))
            R_peak_loc_range_mean(j) = (ecg_MWI_slope_begin(i)+ecg_MWI_slope_end(i))/2; % Each peak range is obtained by computing the mean of the beginning and the end for each slope
        end
    end
end

R_peak_loc_range = R_peak_loc_range_mean(find(R_peak_loc_range_mean > 0));

for i = 1:length(R_peak_loc_range)-1
    if R_peak_loc_range(i+1)-R_peak_loc_range(i)>50
       R_peak_loc(i) = R_peak_loc_range(i);
    end
end

R_peak_loc = round(R_peak_loc(find(R_peak_loc>0)));

subplot 212
plot(t,ecg_MWI);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
title('ECG 1 after Moving-Window Integration','FontSize',16);
xlim([0 4])

hold on;
plot(R_peak_loc/Fs,ecg_MWI(R_peak_loc),'rv','MarkerFaceColor','r');
legend('ECG signal','R waves');

R_peak_loc_delay = R_peak_loc - Win_width - 2; % Delay necessary because MWI introduces a time shift from the orginal signal 

for i=1:length(R_peak_loc_delay)
    if (ecg_norm_1(R_peak_loc_delay(i)) < 0)
       R_peak_loc_delay(i) = R_peak_loc_delay(i) - 8;
    end
end

%% Q and S waves detection

for i=1:length(R_peak_loc_delay)-1
    [S_peak(i) S_peak_loc(i)] = min(ecg_norm_1([R_peak_loc_delay(i):R_peak_loc_delay(i)+15]));
    S_peak_loc(i)=S_peak_loc(i)+R_peak_loc_delay(i);
end

for i=1:length(R_peak_loc_delay)-1
    [Q_peak(i) Q_peak_loc(i)] = min(ecg_norm_1([R_peak_loc_delay(i)-15:R_peak_loc_delay(i)]));
    Q_peak_loc(i)=R_peak_loc_delay(i)-15+Q_peak_loc(i);
end

%% P and T waves detection

for i=1:length(R_peak_loc_delay)-1
    [T_peak(i) T_peak_loc(i)] = max(ecg_norm_1([R_peak_loc_delay(i)+30:R_peak_loc_delay(i)+0.7*(R_peak_loc_delay(i+1)-R_peak_loc_delay(i))]));
    T_peak_loc(i)=T_peak_loc(i)+R_peak_loc_delay(i)+30;
end

for i=1:length(R_peak_loc_delay)-1
    Range(i) = round(R_peak_loc_delay(i)+0.7*(R_peak_loc_delay(i+1)-R_peak_loc_delay(i)));
end

for i=1:length(R_peak_loc_delay)-1
    [P_peak(i) P_peak_loc(i)] = max(ecg_norm_1([Range(i):R_peak_loc_delay(i+1)-20]));
    P_peak_loc(i) = P_peak_loc(i)+Range(i);
end

%% G1 differiator filter
G1_b = [1 zeros(1,5) -1];
G1_a = [1];

% G1 differiator filter frequency response
[H_G1,F] = freqz(G1_b,G1_a,1024,Fs);

figure(7)
subplot 211
plot(F,db(H_G1)); % Amplitude
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('G1 filter amplitude |H_{G1}(f)|','FontSize',16);
title('G1 diferiator filter frequency response','FontSize',16);

subplot 212
plot(F,unwrap(angle(H_G1))); 
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('G1 filter phase arg(H_{G1}(f))','FontSize',16); % Phase

%% G2 low-pass filter
G2_b = [1 zeros(1,7) -1];
G2_a = [1,-1];

% G1 differiator filter frequency response
[H_G2,F] = freqz(G2_b,G2_a,1024,Fs);

figure(8)
subplot 211
plot(F,db(H_G2)); % Amplitude
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('G2 filter amplitude |H_{G2}(f)|','FontSize',16);
title('G2 Low-pass filter frequency response','FontSize',16);

subplot 212
plot(F,unwrap(angle(H_G2))); 
set(gca,'FontSize',14)
xlabel('Frequency (Hz)','FontSize',16);
ylabel('G2 filter phase arg(H_{G2}(f))','FontSize',16); % Phase

%% Final display

figure(9)
plot(t,ecg_norm_1);
set(gca,'FontSize',14)
xlabel('Time (s)','FontSize',16);
ylabel('ECG signal amplitude','FontSize',16);
title('ECG 1 with P,QRS and T waves display','FontSize',16);
xlim([0 4])

hold on;
plot(P_peak_loc/Fs,ecg_norm_1(P_peak_loc),'cv','MarkerFaceColor','c','MarkerSize',8);
plot(Q_peak_loc/Fs,ecg_norm_1(Q_peak_loc),'gv','MarkerFaceColor','g','MarkerSize',8);
plot(R_peak_loc_delay/Fs,ecg_norm_1(R_peak_loc_delay),'rv','MarkerFaceColor','r','MarkerSize',8);
plot(S_peak_loc/Fs,ecg_norm_1(S_peak_loc),'bv','MarkerFaceColor','b','MarkerSize',8);
plot(T_peak_loc/Fs,ecg_norm_1(T_peak_loc),'mv','MarkerFaceColor','m','MarkerSize',8);
legend('ECG signal','P waves','Q waves','R waves','S waves','T waves');

